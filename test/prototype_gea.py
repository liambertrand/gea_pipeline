#!/usr/bin/env python3
"""
numpy re-implementation of the pipeline's statistical core, used to sanity-check
the simulator (are planted signals recoverable?) and to pick sensible defaults.
It mirrors the R scripts: filtering -> ancestry-aware imputation -> PCA ->
env pruning (|r|, VIF) -> partial RDA + rdadapt -> LFMM2 (ridge) -> intersection
-> gene-level windows. It is NOT the production pipeline.
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.covariance import MinCovDet

sp_dir, sp = sys.argv[1], sys.argv[2]
gff = sys.argv[3]
K_STRUCT = 3

# ---------------- read VCF
chrom, pos, G = [], [], []
with open(f"{sp_dir}/{sp}.vcf") as fh:
    for line in fh:
        if line[0] == "#":
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
            continue
        f = line.rstrip().split("\t")
        if int(f[5]) < 30:
            continue
        chrom.append(f[0]); pos.append(int(f[1]))
        G.append([np.nan if g == "./." else g.count("1") for g in f[9:]])
G = np.array(G, dtype=float).T  # n x L
chrom = np.array(chrom); pos = np.array(pos)
maf = np.nanmean(G, 0) / 2; maf = np.minimum(maf, 1 - maf)
miss = np.isnan(G).mean(0)
keep = (maf >= 0.05) & (miss <= 0.1)
G, chrom, pos = G[:, keep], chrom[keep], pos[keep]
ids = np.char.add(np.char.add(chrom, ":"), pos.astype(str))
print(f"{sp}: {G.shape[0]} ind x {G.shape[1]} SNPs after filtering")

# ---------------- structure: PCA on thinned set (every 10th SNP as LD-prune stand-in)
thin = np.arange(0, G.shape[1], 10)
Gt = G[:, thin].copy()
Gt = np.where(np.isnan(Gt), np.nanmean(Gt, 0), Gt)
Gs = (Gt - Gt.mean(0)) / (Gt.std(0) + 1e-9)
U, S, Vt = np.linalg.svd(Gs, full_matrices=False)
PCs = U[:, :K_STRUCT] * S[:K_STRUCT]
clusters = KMeans(K_STRUCT, n_init=10, random_state=1).fit_predict(PCs)

# impute: mode genotype within ancestry cluster
for k in range(K_STRUCT):
    idx = clusters == k
    sub = G[idx]
    for j in np.where(np.isnan(sub).any(0))[0]:
        vals = sub[~np.isnan(sub[:, j]), j]
        mode = stats.mode(vals, keepdims=False).mode if len(vals) else np.nanmedian(G[:, j])
        col = sub[:, j]; col[np.isnan(col)] = mode; sub[:, j] = col
    G[idx] = sub

# ---------------- environment
env = pd.read_csv(f"{sp_dir}/{sp}_samples.csv").set_index("sample_id").loc[samples]
collapse = {21: "developed", 22: "developed", 41: "forest", 42: "forest", 43: "forest",
            52: "shrub", 71: "grass", 81: "agri", 82: "agri", 31: "barren"}
lc = env["nlcd_class"].map(collapse)
dummies = pd.get_dummies(lc, prefix="lc").astype(float)
dummies = dummies.loc[:, dummies.mean() >= 0.05]
cand = env[[f"bio{i}" for i in range(1, 20)]].join(dummies)
priority = ["bio1", "bio12", "bio4", "bio15", "bio2", "bio7", "lc_forest", "lc_shrub", "lc_grass"]
order = [c for c in priority if c in cand] + [c for c in cand if c not in priority]
sel = []
R = cand.corr().abs()
for c in order:
    if all(R.loc[c, s] < 0.7 for s in sel):
        sel.append(c)
def vif(df):
    X = (df - df.mean()) / df.std()
    return pd.Series(np.diag(np.linalg.inv(np.corrcoef(X.values.T))), index=df.columns)
while True:
    v = vif(cand[sel])
    if v.max() < 5: break
    sel.remove(v.idxmax())
print("selected env:", sel)
X = cand[sel].values.astype(float); X = (X - X.mean(0)) / X.std(0)

# ---------------- partial RDA conditioned on PCs
Y = G - G.mean(0)
Z = np.c_[np.ones(len(X)), PCs]
def resid(A, Z): return A - Z @ np.linalg.lstsq(Z, A, rcond=None)[0]
Yr, Xr = resid(Y, Z), resid(X, Z)
Yhat = Xr @ np.linalg.lstsq(Xr, Yr, rcond=None)[0]
u, s, vt = np.linalg.svd(Yhat, full_matrices=False)
ev = s ** 2
print("RDA constrained eigen prop:", np.round(ev[:len(sel)] / ev[:len(sel)].sum(), 3))
K_rda = 2
load = vt[:K_rda].T
zs = (load - load.mean(0)) / load.std(0)
mcd = MinCovDet(random_state=1).fit(zs)
md = mcd.mahalanobis(zs)
lam = np.median(md) / stats.chi2.ppf(0.5, K_rda)
p_rda = stats.chi2.sf(md / lam, K_rda)

def bh(p):
    n = len(p); o = np.argsort(p); q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty(n); out[o] = np.minimum(q, 1); return out
q_rda = bh(p_rda)

# ---------------- LFMM2 ridge
def lfmm2(Y, X, K, lam=1e-5):
    n, d = X.shape
    Q, sv, _ = np.linalg.svd(X, full_matrices=True)
    dvec = np.ones(n); dvec[:d] = np.sqrt(lam / (lam + sv ** 2))
    Dm = np.diag(dvec); Dinv = np.diag(1 / dvec)
    uu, ss, vv = np.linalg.svd(Dm @ Q.T @ Y, full_matrices=False)
    W = Q @ Dinv @ (uu[:, :K] * ss[:K]) @ vv[:K]
    Ulat = Q @ Dinv @ (uu[:, :K] * ss[:K])
    return Ulat
Yc = G - G.mean(0)
Ulat = lfmm2(Yc, X, K_STRUCT)
C = np.c_[np.ones(len(X)), X, Ulat]
beta, *_ = np.linalg.lstsq(C, Yc, rcond=None)
res = Yc - C @ beta
dfree = len(X) - C.shape[1]
s2 = (res ** 2).sum(0) / dfree
XtXi = np.linalg.inv(C.T @ C)
se = np.sqrt(np.outer(np.diag(XtXi)[1:1 + X.shape[1]], s2))
z = beta[1:1 + X.shape[1]] / se  # d x L
gif = np.median(z ** 2, axis=1) / stats.chi2.ppf(0.5, 1)
p_lfmm = stats.chi2.sf(z ** 2 / gif[:, None], 1)
print("LFMM GIF per variable:", dict(zip(sel, np.round(gif, 2))))
q_lfmm = np.vstack([bh(p) for p in p_lfmm])
lfmm_hit = (q_lfmm < 0.05).any(0)
rda_hit = q_rda < 0.05
both = rda_hit & lfmm_hit

# ---------------- truth
truth = pd.read_csv(f"{sp_dir}/truth_snps.tsv", sep="\t")
causal = set(truth.loc[truth.role == "causal", "snp_id"])
tgenes = set(pd.read_csv(f"{sp_dir}/truth_genes.tsv", sep="\t").gene_id)
is_causal = np.isin(ids, list(causal))
is_true = np.isin(ids, list(truth.snp_id))
for name, hit in [("RDA", rda_hit), ("LFMM", lfmm_hit), ("both", both)]:
    print(f"{name:5s}: n={hit.sum():5d}  causal recovered={np.sum(hit & is_causal)}/{is_causal.sum()}"
          f"  FDR(vs causal+linked)={1 - np.sum(hit & is_true) / max(hit.sum(), 1):.2f}")

# ---------------- gene-level: candidates inside gene +/- 5 kb, binomial vs genome rate
g = pd.read_csv(gff, sep="\t", comment="#", header=None)
g = g[g[2] == "gene"]
gid = g[8].str.extract(r"ID=([^;]+)")[0].values
rate = both.mean()
rows = []
for (c, st, en, i) in zip(g[0], g[3], g[4], gid):
    m = (chrom == c) & (pos >= st - 5000) & (pos <= en + 5000)
    n, k = m.sum(), (both & m).sum()
    p = stats.binom.sf(k - 1, n, rate) if n and k else 1.0
    rows.append((i, n, k, p))
gt = pd.DataFrame(rows, columns=["gene", "n", "k", "p"])
gt["q"] = bh(gt.p.values)
hitg = set(gt.loc[(gt.q < 0.05) & (gt.k >= 2), "gene"])
print(f"genes: detected={len(hitg)}  true={len(tgenes)}  TP={len(hitg & tgenes)}  FDR={1 - len(hitg & tgenes) / max(len(hitg), 1):.2f}")

#!/usr/bin/env python3
"""
Simulate fake whole-genome GEA datasets for several species that share one
reference genome (same chromosomes and gene IDs), so the pipeline can be
tested against a known truth.

For each species the simulator produces:
  <out>/<species>/<species>.vcf           biallelic SNPs, GT only, some missing
  <out>/<species>/<species>_samples.csv   one row per individual: coords,
                                          elevation, bio1-bio19, NLCD class
  <out>/<species>/truth_snps.tsv          causal + linked SNPs (the answer key)
  <out>/<species>/truth_genes.tsv         causal genes and the driving variable
and once, shared by all species:
  <out>/reference/annotation.gff3         gene models on the shared reference

The generative model is intentionally simple but has the features that make
real GEA hard:
  * ancestry clusters arranged along latitude, so neutral structure is
    confounded with temperature (tests the structure correction),
  * isolation-by-distance noise among sampling sites,
  * 19 correlated bioclim variables (tests collinearity pruning),
  * a categorical NLCD land-cover class per point,
  * selected loci whose allele frequency tracks one environmental variable,
    with nearby linked SNPs so gene/window-level signal exists,
  * some causal genes shared across species and some species-specific
    (tests the cross-species comparison step).

Only numpy is required.
"""
import argparse
import json
import os

import numpy as np

# NLCD 2019 class codes (subset relevant to rattlesnake habitat)
NLCD_CLASSES = [21, 22, 41, 42, 43, 52, 71, 81, 82, 31]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True)
    p.add_argument("--species", nargs="+", default=["sppA", "sppB", "sppC", "sppD"])
    p.add_argument("--n-sites", type=int, default=20)
    p.add_argument("--ind-per-site", type=int, default=8)
    p.add_argument("--n-chrom", type=int, default=6)
    p.add_argument("--chrom-len", type=int, default=3_000_000)
    p.add_argument("--snp-density", type=float, default=1 / 600,
                   help="SNPs per bp (after filtering, roughly)")
    p.add_argument("--n-causal", type=int, default=24,
                   help="causal genes per species")
    p.add_argument("--n-shared-causal", type=int, default=8,
                   help="causal genes shared by all species")
    p.add_argument("--fst", type=float, default=0.08)
    p.add_argument("--sel-strength", type=float, default=2.5,
                   help="logit slope of causal allele freq on the standardized env variable")
    p.add_argument("--missing", type=float, default=0.02)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


# --------------------------------------------------------------------------
# Shared reference: chromosomes + genes
# --------------------------------------------------------------------------
def make_reference(rng, n_chrom, chrom_len):
    genes = []
    gid = 0
    for c in range(1, n_chrom + 1):
        pos = int(rng.integers(5_000, 30_000))
        while pos < chrom_len - 40_000:
            glen = int(rng.integers(4_000, 40_000))
            gid += 1
            genes.append(dict(chrom=f"chr{c}", start=pos, end=pos + glen,
                              id=f"CRO_g{gid:05d}", strand=rng.choice(["+", "-"])))
            pos += glen + int(rng.integers(5_000, 60_000))
    return genes


def write_gff(genes, path):
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for g in genes:
            fh.write(f"{g['chrom']}\tsim\tgene\t{g['start']}\t{g['end']}\t.\t{g['strand']}\t.\t"
                     f"ID={g['id']};Name={g['id']}\n")


# --------------------------------------------------------------------------
# Sites and environment
# --------------------------------------------------------------------------
def make_sites(rng, n_sites, species_idx):
    # each species occupies its own chunk of the SW US, with lat + elevation gradients
    lat0 = 29 + 2.0 * species_idx
    lon0 = -113 + 2.5 * species_idx
    lat = lat0 + rng.uniform(0, 8, n_sites)
    lon = lon0 + rng.uniform(0, 7, n_sites)
    # elevation only loosely related to latitude so the two can be separated
    elev = np.clip(rng.uniform(100, 2800, n_sites) + 60 * (lat - lat0), 0, 3400)
    return lat, lon, elev


def make_env(rng, lat, lon, elev):
    n = len(lat)
    z = lambda v: (v - v.mean()) / v.std()
    # latent climate drivers
    temp = 24 - 0.9 * (lat - 30) - 0.0062 * elev + rng.normal(0, 0.7, n)  # deg C
    precip = 300 + 55 * (lon - lon.min()) + 0.12 * elev + rng.normal(0, 40, n)  # mm
    seas = 5 + 0.35 * (lat - 30) + rng.normal(0, 0.4, n)  # temperature seasonality
    env = {}
    env["bio1"] = temp
    env["bio2"] = 14 + 0.8 * z(seas) + rng.normal(0, 0.6, n)
    env["bio3"] = 40 - 2 * z(seas) + rng.normal(0, 1.5, n)
    env["bio4"] = 100 * seas + rng.normal(0, 20, n)
    env["bio5"] = temp + 13 + 0.5 * z(seas) + rng.normal(0, 0.5, n)
    env["bio6"] = temp - 13 - 1.2 * z(seas) + rng.normal(0, 0.5, n)
    env["bio7"] = env["bio5"] - env["bio6"]
    env["bio8"] = temp + rng.normal(4, 2.5, n)
    env["bio9"] = temp + rng.normal(-3, 2.5, n)
    env["bio10"] = temp + 8 + 0.3 * z(seas) + rng.normal(0, 0.4, n)
    env["bio11"] = temp - 8 - 0.6 * z(seas) + rng.normal(0, 0.4, n)
    env["bio12"] = precip
    env["bio13"] = precip * 0.17 + rng.normal(0, 6, n)
    env["bio14"] = precip * 0.02 + rng.normal(0, 2, n)
    env["bio15"] = 60 - 0.04 * precip + rng.normal(0, 6, n)
    env["bio16"] = precip * 0.42 + rng.normal(0, 12, n)
    env["bio17"] = precip * 0.08 + rng.normal(0, 5, n)
    env["bio18"] = precip * 0.30 + rng.normal(0, 25, n)
    env["bio19"] = precip * 0.20 + rng.normal(0, 25, n)
    return env


def assign_nlcd(rng, precip_ind, elev_ind):
    """Point NLCD class, more forest in wet/high places, shrub in dry places."""
    zp = (precip_ind - precip_ind.mean()) / precip_ind.std()
    ze = (elev_ind - elev_ind.mean()) / elev_ind.std()
    out = []
    for p, e in zip(zp, ze):
        w = np.array([
            0.05,                     # 21 developed open
            0.02,                     # 22 developed low
            np.exp(0.9 * p + 0.8 * e) * 0.15,  # 41 deciduous forest
            np.exp(0.7 * p + 1.2 * e) * 0.20,  # 42 evergreen forest
            np.exp(0.8 * p + 0.5 * e) * 0.05,  # 43 mixed forest
            np.exp(-1.2 * p) * 0.45,  # 52 shrub/scrub
            np.exp(-0.2 * p) * 0.25,  # 71 grassland
            0.04,                     # 81 pasture
            0.03,                     # 82 crops
            np.exp(-1.0 * p) * 0.03,  # 31 barren
        ])
        out.append(int(rng.choice(NLCD_CLASSES, p=w / w.sum())))
    return np.array(out)


# --------------------------------------------------------------------------
# Genotypes
# --------------------------------------------------------------------------
def logistic(x):
    return 1 / (1 + np.exp(-x))


def simulate_species(rng, sp, sp_idx, args, genes, causal_shared):
    n_sites, n_per = args.n_sites, args.ind_per_site
    lat, lon, elev = make_sites(rng, n_sites, sp_idx)
    site_env = make_env(rng, lat, lon, elev)

    # individuals: jitter coordinates / elevation slightly around site
    site_of = np.repeat(np.arange(n_sites), n_per)
    n_ind = len(site_of)
    ilat = lat[site_of] + rng.normal(0, 0.02, n_ind)
    ilon = lon[site_of] + rng.normal(0, 0.02, n_ind)
    ielev = np.clip(elev[site_of] + rng.normal(0, 40, n_ind), 0, None)
    ienv = {k: v[site_of] + rng.normal(0, 0.03 * v.std(), n_ind) for k, v in site_env.items()}
    nlcd = assign_nlcd(rng, ienv["bio12"], ielev)
    forest = np.isin(nlcd, [41, 42, 43]).astype(float)

    # ---- ancestry: 3 clusters arranged along latitude (confounded with bio1)
    K = 3
    centers = np.quantile(lat, [0.15, 0.5, 0.85])
    d = -((lat[:, None] - centers[None, :]) ** 2) / 2.0
    q_site = np.exp(d) / np.exp(d).sum(1, keepdims=True)

    # ---- SNP positions
    snps = []
    for c in range(1, args.n_chrom + 1):
        n = rng.poisson(args.chrom_len * args.snp_density)
        pos = np.sort(rng.choice(np.arange(1_000, args.chrom_len - 1_000), n, replace=False))
        snps += [(f"chr{c}", int(p)) for p in pos]
    L = len(snps)
    chrom_arr = np.array([s[0] for s in snps])
    pos_arr = np.array([s[1] for s in snps])

    # ---- neutral allele frequencies (Balding-Nichols)
    p_anc = rng.beta(0.7, 0.7, L) * 0.9 + 0.05
    a = p_anc * (1 - args.fst) / args.fst
    b = (1 - p_anc) * (1 - args.fst) / args.fst
    p_clu = rng.beta(a[None, :].repeat(K, 0), b[None, :].repeat(K, 0))  # K x L
    p_site = q_site @ p_clu  # sites x L
    # isolation-by-distance style noise, spatially smooth via random site effects
    coords = np.c_[lat, lon]
    D = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(-1))
    C = np.exp(-D / 2.0)
    Lc = np.linalg.cholesky(C + 1e-6 * np.eye(n_sites))
    ibd = (Lc @ rng.normal(0, 0.35, (n_sites, L)))
    p_site = logistic(np.log(p_site / (1 - p_site)) + ibd)

    # ---- causal genes: shared + species-specific
    n_spec = args.n_causal - len(causal_shared)
    remaining = [g for g in range(len(genes)) if g not in causal_shared]
    spec = list(rng.choice(remaining, n_spec, replace=False))
    causal_gene_idx = list(causal_shared) + spec
    drivers_pool = ["bio1"] * 10 + ["bio12"] * 7 + ["bio4"] * 3 + ["forest"] * 4
    # deterministic driver for shared genes (same across species) → convergence
    drivers = {}
    for i, g in enumerate(causal_gene_idx):
        if g in causal_shared:
            drivers[g] = ["bio1", "bio12", "bio1", "bio4", "bio12", "bio1", "forest", "bio12"][
                causal_shared.index(g) % 8]
        else:
            drivers[g] = drivers_pool[int(rng.integers(len(drivers_pool)))]

    site_driver = {
        "bio1": site_env["bio1"], "bio12": site_env["bio12"], "bio4": site_env["bio4"],
    }
    # forest proportion by site (from individual classes)
    site_driver["forest"] = np.array([forest[site_of == s].mean() for s in range(n_sites)])

    truth_rows, gene_rows = [], []
    causal_flag = np.zeros(L, dtype=bool)
    linked_to = {}  # snp idx -> causal snp idx
    for g in causal_gene_idx:
        G = genes[g]
        in_gene = np.where((chrom_arr == G["chrom"]) & (pos_arr >= G["start"]) & (pos_arr <= G["end"]))[0]
        if len(in_gene) == 0:
            # force one SNP into the gene
            j = int(np.argmin(np.abs(np.where(chrom_arr == G["chrom"], pos_arr, 10**12) - (G["start"] + G["end"]) // 2)))
            in_gene = np.array([j])
        cidx = int(rng.choice(in_gene))
        env_v = site_driver[drivers[g]]
        zv = (env_v - env_v.mean()) / (env_v.std() + 1e-9)
        sign = rng.choice([-1, 1])
        base = rng.normal(0, 0.4)
        p_site[:, cidx] = np.clip(logistic(base + sign * args.sel_strength * zv), 0.01, 0.99)
        causal_flag[cidx] = True
        gene_rows.append((G["id"], G["chrom"], G["start"], G["end"], drivers[g],
                          "shared" if g in causal_shared else "specific"))
        truth_rows.append((f"{chrom_arr[cidx]}:{pos_arr[cidx]}", G["id"], drivers[g], "causal"))
        # linked SNPs: within 10 kb
        near = np.where((chrom_arr == G["chrom"]) & (np.abs(pos_arr - pos_arr[cidx]) <= 10_000))[0]
        for j in near:
            if j != cidx and not causal_flag[j]:
                linked_to[j] = cidx

    # ---- draw genotypes
    geno = np.empty((n_ind, L), dtype=np.int8)
    for s in range(n_sites):
        idx = np.where(site_of == s)[0]
        geno[idx] = rng.binomial(2, p_site[s][None, :].repeat(len(idx), 0))
    # linkage: linked SNPs copy the causal genotype with prob r2-ish decaying with distance
    for j, c in linked_to.items():
        dist = abs(pos_arr[j] - pos_arr[c])
        r = 0.9 * np.exp(-dist / 6_000)
        copy = rng.random(n_ind) < r
        geno[copy, j] = geno[copy, c]
        truth_rows.append((f"{chrom_arr[j]}:{pos_arr[j]}", None, None, "linked"))
    # fill gene/driver for linked
    causal_info = {r[0]: r for r in truth_rows if r[3] == "causal"}
    truth_final = []
    for r in truth_rows:
        if r[3] == "linked":
            j = [k for k in linked_to if f"{chrom_arr[k]}:{pos_arr[k]}" == r[0]][0]
            ci = causal_info[f"{chrom_arr[linked_to[j]]}:{pos_arr[linked_to[j]]}"]
            truth_final.append((r[0], ci[1], ci[2], "linked"))
        else:
            truth_final.append(r)

    # ---- missing data
    miss = rng.random(geno.shape) < args.missing

    # ---- write
    odir = os.path.join(args.out, sp)
    os.makedirs(odir, exist_ok=True)
    samples = [f"{sp}_{s+1:02d}_{i+1:02d}" for s in range(n_sites) for i in range(n_per)]
    gt_str = {0: "0/0", 1: "0/1", 2: "1/1"}
    with open(os.path.join(odir, f"{sp}.vcf"), "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for c in range(1, args.n_chrom + 1):
            fh.write(f"##contig=<ID=chr{c},length={args.chrom_len}>\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        bases = np.array(list("ACGT"))
        for j in range(L):
            ref, alt = rng.choice(bases, 2, replace=False)
            qual = int(rng.integers(15, 999))
            col = geno[:, j]
            gts = ["./." if miss[i, j] else gt_str[int(col[i])] for i in range(n_ind)]
            fh.write(f"{chrom_arr[j]}\t{pos_arr[j]}\t.\t{ref}\t{alt}\t{qual}\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")

    with open(os.path.join(odir, f"{sp}_samples.csv"), "w") as fh:
        cols = ["sample_id", "species", "site", "lat", "lon", "elevation"] + [f"bio{i}" for i in range(1, 20)] + ["nlcd_class"]
        fh.write(",".join(cols) + "\n")
        for i in range(n_ind):
            row = [samples[i], sp, f"{sp}_site{site_of[i]+1:02d}", f"{ilat[i]:.5f}", f"{ilon[i]:.5f}", f"{ielev[i]:.1f}"]
            row += [f"{ienv[f'bio{k}'][i]:.3f}" for k in range(1, 20)]
            row += [str(nlcd[i])]
            fh.write(",".join(row) + "\n")

    with open(os.path.join(odir, "truth_snps.tsv"), "w") as fh:
        fh.write("snp_id\tgene_id\tdriver\trole\n")
        for r in truth_final:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(odir, "truth_genes.tsv"), "w") as fh:
        fh.write("gene_id\tchrom\tstart\tend\tdriver\tclass\n")
        for r in gene_rows:
            fh.write("\t".join(map(str, r)) + "\n")

    return dict(species=sp, n_ind=n_ind, n_snps=L, n_causal=len(causal_gene_idx),
                n_linked=len(linked_to))


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    os.makedirs(os.path.join(args.out, "reference"), exist_ok=True)
    genes = make_reference(rng, args.n_chrom, args.chrom_len)
    write_gff(genes, os.path.join(args.out, "reference", "annotation.gff3"))
    causal_shared = list(rng.choice(len(genes), args.n_shared_causal, replace=False))
    summary = []
    for i, sp in enumerate(args.species):
        summary.append(simulate_species(np.random.default_rng(args.seed + 100 * (i + 1)),
                                        sp, i, args, genes, [int(g) for g in causal_shared]))
    summary = dict(n_genes=len(genes), species=summary,
                   shared_causal_genes=[genes[g]["id"] for g in causal_shared],
                   params={k: v for k, v in vars(args).items()})
    with open(os.path.join(args.out, "simulation_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()

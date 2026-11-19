# Environment-by-SNP (GEA) pipeline

Maps environmentally associated SNPs, genomic windows and genes **within each
species** from whole-genome resequencing, then compares the resulting gene sets
across species. Built on Capblancq & Forester (2021) RDA plus genome-wide LFMM2,
with results intersected.

```
VCF ─► filter (bcftools/plink2) ─┬─► LD-pruned ─► PCA + sNMF (K) ─┐
                                 └─► all SNPs ─► ancestry-aware imputation
samples.csv (bioclim, NLCD, coords) ─► NLCD dummies ─► |r| + VIF pruning
                                                                   │
            ┌──────────────── partial RDA + rdadapt ◄──────────────┤
            │                 LFMM2 ridge (genome-wide) ◄──────────┘
            ▼
   tier1 = RDA ∩ LFMM  ─► sliding windows + genes (binomial enrichment)
            ─► per-species report ─► cross-species overlap (hypergeometric)
```

## Repository contents

```
gea-pipeline/
├── README.md                     this file
├── run_viridis.sh                single-species wrapper: input checks, dry run, real run (C. viridis, 8 threads)
├── run_crovir3_test.sh           earlier wrapper for the CroVir3 test VCF (multi-species Snakefile)
├── config/
│   ├── config.yaml               multi-species settings (4 CVOS species, cross-species step)
│   └── config_viridis.yaml       single-species settings used for the C. viridis run
├── envs/
│   ├── snakemake.yaml            conda env to run Snakemake (create once)
│   └── gea.yaml                  conda env for all analysis steps (bcftools, plink2, R + vegan, LEA, qvalue, robust)
├── profiles/
│   ├── local/config.yaml         Snakemake profile for a workstation
│   └── slurm/config.yaml         Snakemake profile for a SLURM cluster
├── workflow/
│   ├── Snakefile                 multi-species workflow (per-species GEA + cross-species comparison)
│   ├── Snakefile_viridis         single-species workflow (no cross-species steps)
│   └── scripts/
│       ├── common.R              shared helpers (argument parsing, config, q-values, rdadapt)
│       ├── env_prep.R            NLCD habitat dummies, |r| + VIF predictor pruning, scaling
│       ├── structure.R           PCA + sNMF on LD-pruned SNPs, choice of K, ancestry-aware imputation
│       ├── rda.R                 variance partitioning, partial RDA, rdadapt outliers
│       ├── lfmm.R                genome-wide ridge LFMM2 (LEA), GIF calibration, q-values
│       ├── candidates_genes.R    tier-1/tier-2 SNPs, window and gene enrichment (binomial + BH)
│       ├── plots.R               allele frequency vs environment, regional zooms, gene x predictor heatmap
│       ├── report.R              per-species report.md
│       ├── cross_species.R       gene overlap across species and figures
│       ├── evaluate.R            simulation only: recall and FDR against the known truth
│       ├── replot_with_legends.R standalone: redraws the main figures with full legends (PNG + PDF)
│       └── snp_annotation_classes.R  standalone: CDS / UTR / intron / flank / intergenic classes + gene homolog table
├── sim/
│   └── simulate_gea.py           simulator of four species with known causal genes
└── test/
    ├── config_test.yaml          settings for the simulation test
    ├── run_test.sh               end-to-end test on simulated data
    └── prototype_gea.py          numpy re-implementation of the core statistics (benchmark)
```

Not tracked: input data (VCFs, environment CSVs, reference annotation) and `results*/`
folders, which are large; see "Run on real data" for where they go.

## Design decisions (and why)

| Step | Choice | Reason |
|---|---|---|
| Scope | Each species analysed alone; same settings | Within-species GEA; the cross-species step only compares outputs |
| Structure | sNMF K (min cross-entropy) + PCs on LD-pruned SNPs | Neutral structure along the latitude gradient mimics climate selection |
| Lat / lon / elevation | Used as **spatial conditioning**, not as predictors | They proxy climate *and* geography/ancestry; bioclim carries the environment |
| NLCD | Point class → collapsed habitat dummies (or buffer % cover) | Categorical codes can't enter RDA/LFMM as numbers |
| Collinearity | Priority-ordered \|r\| < 0.7 then VIF < 5, per species | Collinearity differs between species' ranges |
| RDA | pRDA `G ~ env + Condition(structure + geo)`, rdadapt q-values | Capblancq & Forester 2021; the q-values give an FDR, unlike ±3 SD |
| LFMM2 | Run **genome-wide**, independent of RDA | Its GIF calibration/FDR assumes the full p-value distribution; running it on RDA hits only would be circular |
| Candidates | tier1 = both methods; tier2 = one method | Consensus sharply lowers false positives (see test results) |
| SNP → gene | gene ± 5 kb; binomial test for excess tier1 SNPs | WGS SNPs are in LD; counting single SNPs overstates signal and favours long genes |
| Cross-species | Shared reference gene IDs, or an orthology table | Pairwise hypergeometric overlap within genes testable in both species |

## Install

```bash
conda env create -f envs/snakemake.yaml   # snakemake 8 + slurm plugin
conda activate snakemake
# analysis tools (bcftools, plink2, R, vegan, LEA, qvalue, robust) are
# installed automatically by snakemake from envs/gea.yaml
```

## Status (2026-11-18)
- **Simulation benchmark** (`test/prototype_gea.py`, numpy re-implementation of the core steps, one simulated species: 160 individuals, ~27k SNPs, 24 causal genes):
  - RDA alone: 21/24 causal SNPs, FDR 0.23
  - LFMM2 alone: 24/24 causal SNPs, FDR 0.06
  - tier1 (RDA ∩ LFMM2): 21/24 causal SNPs, FDR 0.02
  - gene level: 18/24 causal genes, 0 false genes
- **Real-data test run completed (16 Nov 2026)** with `workflow/Snakefile_viridis` on the lab workstation: 100 *C. viridis* (40 localities), 5,325,812 SNPs on CroVir3, K = 7, 6 predictors, 504 tier-1 SNPs, 18 genes and 61 windows called. RDA permutations were switched off (`rda.permutations: 0`) because they exceeded 125 GB of memory.
- Still to do: run the full Snakemake simulation test (`test/run_test.sh`); robustness re-runs (kinship, K, Z chromosome, leave-one-locality-out); RDA permutations on a SNP subset; the four-species run on CroVir4.

## Software versions
Minimum versions are pinned in `envs/*.yaml` (Snakemake ≥ 8, Python ≥ 3.11, bcftools ≥ 1.17,
R ≥ 4.3, vegan ≥ 2.6, LEA ≥ 3.0). Record the exact versions used in a run with
`conda list -p .snakemake/conda/<env>` and `snakemake --version`.

## Test on simulated data

```bash
SMALL=1 bash test/run_test.sh     # ~2 species, a few minutes
bash test/run_test.sh             # 4 species
```

`sim/simulate_gea.py` builds four fake species on one shared reference
(6 chromosomes, ~320 genes, ~30k SNPs, 160 individuals across 20 sites). Each
species has:
- ancestry clusters arranged along latitude, so structure is confounded with bio1
- isolation by distance
- 19 correlated bioclim variables
- NLCD point classes
- 24 causal genes, 8 of them shared across species, each driven by bio1, bio12, bio4 or forest cover
- linked SNPs around each causal site

`evaluation.tsv` reports recall and FDR at the SNP and gene level, and whether
the right driver was attributed.

## Run on real data

1. Put one samples CSV per species in `data/`. Columns: `sample_id, lat, lon,
   bio1…bio19, nlcd_class` (the QGIS point-sampling output). IDs must match the
   VCF sample names. Extra columns are ignored unless you list them in
   `env.candidate_vars`.
2. Edit `config/config.yaml`: species blocks, VCF paths, the GFF, and
   `env.priority`. The priority list is the variables you consider biologically
   most meaningful; it decides which member of a correlated pair is kept.
3. Dry run: `snakemake -s workflow/Snakefile -n`
4. Bridges-2: set `slurm_account` in `profiles/slurm/config.yaml`, then
   `snakemake -s workflow/Snakefile --profile profiles/slurm`

### Single-species run on a workstation (as used for *C. viridis*)
1. Edit the paths at the top of `config/config_viridis.yaml` (VCF, environment CSV, GTF) and
   `csv_columns` (the sample-ID, latitude and longitude column names in your CSV).
2. `bash run_viridis.sh` checks the files, sample IDs and contig names, then does a dry run.
3. `bash run_viridis.sh --run` runs the pipeline (start it inside `screen` or `tmux`); results go to
   `results_viridis/Cviridis/`, starting with `report.md`.
4. Optional extras, using the pipeline's R environment (`ENV=$(ls -d .snakemake/conda/*_ | head -1)`):
   - `$ENV/bin/Rscript workflow/scripts/replot_with_legends.R results_viridis/Cviridis`
   - `$ENV/bin/Rscript workflow/scripts/snp_annotation_classes.R results_viridis/Cviridis <annotation.gtf>`

### NLCD as habitat proportions (recommended)
A single point class is noisy for a snake's home range. In QGIS, compute the %
cover of each NLCD class within a buffer (e.g. 500 m–1 km), add columns such as
`lc_forest, lc_shrub, lc_grass`, and set:
```yaml
env:
  nlcd:
    mode: proportion
    prop_cols: [lc_forest, lc_shrub, lc_grass, lc_developed]
```

### Scaling to WGS
RDA and LFMM2 hold the full genotype matrix in memory: individuals × SNPs × 8 bytes.
For example, 150 × 5 M SNPs ≈ 6 GB, plus overhead. The SLURM profile gives these
jobs 96–128 GB. For more SNPs, pre-filter the VCF (depth, mappability, repeats),
or split by chromosome and merge. That extension is on the to-do list.

## Outputs (`results/<species>/`)
| File | Content |
|---|---|
| `genes.tsv` | Every gene: SNPs, tier1 SNPs, enrichment p/q, `called`, associated predictors |
| `candidate_snps.tsv` | SNPs significant in RDA and/or LFMM2, with tier, predictors, gene / nearest gene |
| `windows.tsv.gz` | Sliding-window enrichment |
| `report.md` | Summary: K, selected predictors, variance partitioning, GIFs, top genes |
| `structure/`, `env/`, `rda/`, `lfmm/` | Intermediate tables and PDFs |
| `plots/allele_freq_vs_env.pdf` | For each top gene, the allele frequency of its lead SNP plotted against the gene's main predictor. Points are populations (`plots.pop_col`, sized by n) or individuals, with a logistic fit |
| `plots/regional_zoom.pdf` | One page per top gene: RDA and LFMM2 −log10 p across the gene ± 100 kb, coloured by tier, with a gene-model track |
| `plots/gene_env_heatmap.pdf` | Called genes × selected predictors, filled by the number of tier1 SNPs (tables: `gene_env_matrix.tsv`, `lead_snps.tsv`) |
| `results/cross_species/` | Gene × species overlap and pairwise tests, plus figures: `upset.pdf` (shared-gene intersections), `pairwise_heatmap.pdf` (fold enrichment vs. chance, shared counts, * = q < 0.05), `shared_gene_drivers.pdf` (main predictor of each shared gene in each species), `variance_partitioning.pdf` (environment / geography / structure / confounded, one bar per species) |

## Next extensions
- R.SamBada as an optional third method
- GO enrichment of called genes (gene-length-aware null)
- Genomic offset (gradient forest / RDA-based) under future climate
- Per-chromosome chunking for very large SNP sets

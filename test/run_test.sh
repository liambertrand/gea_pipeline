#!/usr/bin/env bash
# End-to-end test on simulated data (run from anywhere).
#   bash test/run_test.sh              # 4 species (~30k SNPs, 160 ind each)
#   SMALL=1 bash test/run_test.sh      # quick: 2 species, 3 chromosomes
#   CORES=8 bash test/run_test.sh      # more cores
# Extra arguments are passed to snakemake (e.g. -n for a dry run).
set -euo pipefail
cd "$(dirname "$0")/.."

RUN=()
if [[ "${SMALL:-0}" == "1" ]]; then
  [[ -f test/data/sppA/sppA.vcf ]] || python3 sim/simulate_gea.py --out test/data --species sppA sppB \
      --n-chrom 3 --n-causal 16 --n-shared-causal 6 --seed 7
  RUN=(--config "run_species=[sppA,sppB]")
else
  [[ -f test/data/sppD/sppD.vcf ]] || python3 sim/simulate_gea.py --out test/data
fi

# Snakemake >= 8 uses --sdm conda; Snakemake 7 uses --use-conda
MAJOR=$(snakemake --version | cut -d. -f1)
if (( MAJOR >= 8 )); then SDM=(--sdm conda); else SDM=(--use-conda); fi
[[ "${NO_CONDA:-0}" == "1" ]] && SDM=()   # tools already on PATH

snakemake -s workflow/Snakefile --configfile test/config_test.yaml --cores "${CORES:-4}" \
  "${SDM[@]}" "${RUN[@]}" -p "$@"

echo; echo "==== evaluation against simulated truth ===="
for f in results_test/*/evaluation.tsv; do echo "## $f"; column -t -s $'\t' "$f"; echo; done
if [[ -f results_test/cross_species/pairwise_tests.tsv ]]; then
  echo "## cross-species overlap"; column -t -s $'\t' results_test/cross_species/pairwise_tests.tsv
fi

#!/usr/bin/env bash
# =====================================================================
# run_viridis.sh -- single-species GEA run (C. viridis, CroVir3), max 8 threads
#
#   1. Fill in the EDIT ME block in config/config_viridis.yaml
#   2. bash run_viridis.sh          # checks your files + dry run (runs nothing)
#   3. bash run_viridis.sh --run    # real run (start it inside tmux)
# =====================================================================
set -euo pipefail
cd "$(dirname "$0")"
CONFIG="config/config_viridis.yaml"
CORES=8                                   # total CPU threads the pipeline may use

eval "$(conda shell.bash hook)"
conda activate snakemake || { echo "ERROR: conda env 'snakemake' missing (see guide Step 1)"; exit 1; }
command -v bcftools >/dev/null || { echo "ERROR: bcftools missing: conda install -n snakemake -c conda-forge -c bioconda bcftools htslib"; exit 1; }

# Read settings from the config
read_cfg() { python3 -c "import yaml,sys; d=yaml.safe_load(open('$CONFIG')); v=d
for k in sys.argv[1].split('.'): v=v[k]
print(v)" "$1"; }
VCF=$(read_cfg vcf); CSV=$(read_cfg samples_csv); ANNOT=$(read_cfg annotation)
ID_COL=$(read_cfg csv_columns.sample_id); EXPECT=$(read_cfg expected_samples)

echo "== Input files"
fail=0
for f in "$VCF" "$CSV" "$ANNOT"; do
  if [[ -s "$f" ]]; then echo "  ok       $f"; else echo "  MISSING  $f"; fail=1; fi
done
(( fail == 0 )) || { echo "Fix the paths in $CONFIG"; exit 1; }

echo "== Sample IDs (CSV column '$ID_COL' vs VCF)"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
bcftools query -l "$VCF" | sort > "$TMP/vcf"
python3 - "$CSV" "$ID_COL" > "$TMP/csv" <<'PY'
import csv, sys
with open(sys.argv[1], newline="", encoding="utf-8-sig") as fh:
    r = csv.DictReader(fh)
    if sys.argv[2] not in r.fieldnames:
        sys.exit(f"ERROR: column '{sys.argv[2]}' not in CSV. Columns: {r.fieldnames}")
    for row in r:
        if row[sys.argv[2]].strip(): print(row[sys.argv[2]].strip())
PY
sort -o "$TMP/csv" "$TMP/csv"
n_csv=$(wc -l < "$TMP/csv"); n_vcf=$(wc -l < "$TMP/vcf"); n_both=$(comm -12 "$TMP/csv" "$TMP/vcf" | wc -l)
echo "  CSV samples: $n_csv | VCF samples: $n_vcf | in both (will be analysed): $n_both | expected: $EXPECT"
if (( n_both < n_csv )); then
  echo "  In CSV but NOT in the VCF (these will be skipped):"; comm -23 "$TMP/csv" "$TMP/vcf" | head -20 | sed 's/^/    /'
fi
(( n_both > 0 )) || { echo "ERROR: no CSV IDs match VCF sample names"; exit 1; }
(( n_both == EXPECT )) || echo "  WARNING: $n_both samples will be analysed, not $EXPECT"

echo "== Chromosome names (VCF vs annotation)"
bcftools view -h "$VCF" | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' | sort -u > "$TMP/vchr"
gzip -cdf "$ANNOT" | grep -v '^#' | cut -f1 | sort -u > "$TMP/achr"
n_sh=$(comm -12 "$TMP/vchr" "$TMP/achr" | wc -l)
echo "  shared contig names: $n_sh"
(( n_sh > 0 )) || { echo "ERROR: VCF and annotation share no chromosome names"; exit 1; }

# Snakemake >= 8 uses --sdm conda; Snakemake 7 uses --use-conda
if (( $(snakemake --version | cut -d. -f1) >= 8 )); then SDM=(--sdm conda); else SDM=(--use-conda); fi
SMK=(snakemake -s workflow/Snakefile_viridis --configfile "$CONFIG" --cores "$CORES"
     "${SDM[@]}" --printshellcmds --rerun-incomplete)

if [[ "${1:-}" == "--run" ]]; then
  mkdir -p logs; LOG="logs/viridis_$(date +%F_%H%M).log"
  echo "== Running (max $CORES threads); log: $LOG"
  "${SMK[@]}" 2>&1 | tee "$LOG"
  echo "== Done. Start with: $(read_cfg outdir)/$(read_cfg species_name)/report.md"
else
  echo "== Dry run (nothing executed)"
  "${SMK[@]}" -n
  echo "== Checks passed. Run for real with:  bash run_viridis.sh --run"
fi

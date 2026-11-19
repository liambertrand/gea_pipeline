#!/usr/bin/env bash
# =====================================================================
# run_crovir3_test.sh
# Runs the GEA pipeline on the CroVir3 C. viridis test VCF.
#
#   bash run_crovir3_test.sh          # checks your files + dry run (runs nothing)
#   bash run_crovir3_test.sh --run    # checks your files + real run
# =====================================================================
set -euo pipefail

# ######################## EDIT ME ################################
# Use FULL paths (run `echo ~` to see your home folder; "~" won't work here).

# >>> INSERT YOUR VCF FILE HERE (multi-sample .vcf, .vcf.gz or .bcf)
VCF="/home/YOUR_USERNAME/env_by_snps/VCF/YOUR_VIRIDIS_FILE.vcf.gz"

# >>> INSERT YOUR QGIS CSV FILE HERE (the file exported in Step 6 of the guide)
CSV="/home/YOUR_USERNAME/env_by_snps/Env_files/YOUR_QGIS_EXPORT.csv"

# Column in the QGIS CSV holding the sequencing IDs (must match VCF sample names exactly)
ID_COL="sample_id"

# Columns in the QGIS CSV holding decimal-degree coordinates
LAT_COL="lat"
LON_COL="lon"

# >>> INSERT THE CroVir3 GENE ANNOTATION HERE (GTF or GFF3, .gz is fine; same release as the reference)
ANNOT="/home/YOUR_USERNAME/env_by_snps/Ref_files/CroVir3_annotation.gtf"

# Land cover: "proportion" if the CSV has lc_forest, lc_shrub, ... columns (1 km buffer);
#             "class" if it only has the point code column nlcd_class
NLCD_MODE="proportion"
NLCD_PROP_COLS="[lc_forest, lc_shrub, lc_grass, lc_developed, lc_agri]"

# Column grouping individuals into localities (used for allele-frequency plots); "" if none
SITE_COL="site"

# Z / W scaffold names to exclude, comma-separated (snakes are ZW). "" to keep everything.
SEX_SCAFFOLDS=""          # e.g. "scaffold-Z1,scaffold-W"

# CPU cores to use (about 3/4 of `nproc`)
CORES=8
# ##################################################################

SPECIES="Cviridis_crovir3"
OUTDIR="results_crovir3"
PIPE_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$PIPE_DIR"

die() { echo "ERROR: $*" >&2; exit 1; }
say() { echo; echo "== $*"; }

# ---------------------------------------------------------------- conda
eval "$(conda shell.bash hook)"
conda activate snakemake || die "conda env 'snakemake' not found - run the one-time setup first"
command -v bcftools >/dev/null || die "bcftools missing: conda install -n snakemake -c conda-forge -c bioconda bcftools htslib"
[[ -f workflow/Snakefile ]] || die "save this script inside the gea-pipeline folder"

# ---------------------------------------------------------------- files exist
say "Checking input files"
for f in "$VCF" "$CSV" "$ANNOT"; do [[ -s "$f" ]] || die "file not found or empty: $f"; done
echo "VCF: $VCF"; echo "CSV: $CSV"; echo "Annotation: $ANNOT"

# ---------------------------------------------------------------- CSV columns
say "Checking CSV columns"
HEADER=$(head -1 "$CSV" | tr -d '\r"')
missing=""
need_cols=("$ID_COL" "$LAT_COL" "$LON_COL")
for i in $(seq 1 19); do need_cols+=("bio$i"); done
[[ "$NLCD_MODE" == "class" ]] && need_cols+=("nlcd_class")
if [[ "$NLCD_MODE" == "proportion" ]]; then
  for c in $(echo "$NLCD_PROP_COLS" | tr -d '[] ' | tr ',' ' '); do need_cols+=("$c"); done
fi
[[ -n "$SITE_COL" ]] && need_cols+=("$SITE_COL")
for c in "${need_cols[@]}"; do [[ ",$HEADER," == *",$c,"* ]] || missing+=" $c"; done
[[ -z "$missing" ]] || die "CSV is missing column(s):$missing
   CSV header is: $HEADER"
echo "all required columns present"

# ---------------------------------------------------------------- sample IDs
say "Matching sample IDs (VCF vs CSV)"
TMP=$(mktemp -d)
bcftools query -l "$VCF" | sort > "$TMP/vcf_ids"
awk -F',' -v col="$ID_COL" '
  NR==1 { gsub(/\r|"/, ""); for (i=1; i<=NF; i++) if ($i==col) c=i; next }
        { gsub(/\r|"/, ""); print $c }' "$CSV" | sort > "$TMP/csv_ids"
n_vcf=$(wc -l < "$TMP/vcf_ids"); n_csv=$(wc -l < "$TMP/csv_ids")
n_both=$(comm -12 "$TMP/vcf_ids" "$TMP/csv_ids" | wc -l)
dups=$(uniq -d "$TMP/csv_ids" | head -5 | tr '\n' ' ')
echo "VCF samples: $n_vcf | CSV rows: $n_csv | in both: $n_both"
[[ -z "$dups" ]] || die "duplicate IDs in CSV: $dups"
if (( n_both < n_vcf )); then
  echo "In VCF but NOT in CSV (pipeline will stop on these):"
  comm -23 "$TMP/vcf_ids" "$TMP/csv_ids" | head -20
  die "every VCF sample needs a CSV row. Fix the names, or drop samples from the VCF with:
   bcftools view -S keep_ids.txt -Oz -o subset.vcf.gz $VCF"
fi
(( n_csv > n_vcf )) && echo "note: $((n_csv - n_vcf)) CSV rows have no VCF sample; they are ignored"
(( n_vcf < 20 )) && echo "WARNING: only $n_vcf individuals - fine for a test, too few for biology"

# ---------------------------------------------------------------- blanks in CSV
blank=$(awk -F',' 'NR>1 { gsub(/\r/, ""); for (i=1; i<=NF; i++) if ($i=="" || $i=="NA") { print NR; break } }' "$CSV" | head -5 | tr '\n' ' ')
[[ -z "$blank" ]] || die "empty/NA cells in CSV on line(s): $blank (point on NoData? fix or remove)"

# ---------------------------------------------------------------- chromosome names
say "Checking annotation and chromosome names (VCF vs annotation)"
if gzip -cdf "$ANNOT" | grep -v '^#' | head -2000 | grep -q 'gene_id "'; then
  ANNOT_FMT="GTF"; GENE_ID_ATTR="gene_id"
else
  ANNOT_FMT="GFF3"; GENE_ID_ATTR="ID"
fi
n_gene_lines=$(gzip -cdf "$ANNOT" | awk -F'\t' '$3=="gene"' | wc -l)
echo "annotation format: $ANNOT_FMT | 'gene' lines: $n_gene_lines"
(( n_gene_lines > 0 )) || echo "note: no 'gene' lines; gene spans will be built from transcript/exon lines"
bcftools view -h "$VCF" | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' | sort -u > "$TMP/vcf_chr"
gzip -cdf "$ANNOT" | grep -v '^#' | cut -f1 | sort -u > "$TMP/gff_chr"
n_shared=$(comm -12 "$TMP/vcf_chr" "$TMP/gff_chr" | wc -l)
echo "contigs in VCF: $(wc -l < "$TMP/vcf_chr") | in annotation: $(wc -l < "$TMP/gff_chr") | shared: $n_shared"
(( n_shared > 0 )) || die "no shared chromosome names. VCF: $(head -3 "$TMP/vcf_chr" | tr '\n' ' ') annotation: $(head -3 "$TMP/gff_chr" | tr '\n' ' ')
   Rename with: bcftools annotate --rename-chrs map.txt (guide Step 5)"

# ---------------------------------------------------------------- config for this run
say "Writing config_crovir3.yaml"
EXTRA=""
[[ -n "$SEX_SCAFFOLDS" ]] && EXTRA="-t ^${SEX_SCAFFOLDS}"
cat > config_crovir3.yaml <<YAML
# Generated by run_crovir3_test.sh - overrides config/config.yaml for this test
outdir: $OUTDIR
run_species: [$SPECIES]
species:
  $SPECIES:
    vcf: "$VCF"
    samples: "$CSV"
gff: "$ANNOT"            # GTF or GFF3; format detected automatically
gene_id_attr: $GENE_ID_ATTR
filter:
  extra_bcftools: "$EXTRA"
  max_snp_missing: 0.20     # a bit looser than default for a test file
structure:
  K_max: 6
env:
  sample_id_col: "$ID_COL"
  geo_vars: [$LAT_COL, $LON_COL]
  nlcd:
    column: nlcd_class
    mode: $NLCD_MODE
    prop_cols: $NLCD_PROP_COLS
rda:
  permutations: 0           # fast test; use 999 for real results
plots:
  pop_col: "$SITE_COL"
YAML
cat config_crovir3.yaml

SMK=(snakemake -s workflow/Snakefile --configfile config/config.yaml config_crovir3.yaml
     --profile profiles/local --cores "$CORES")

# ---------------------------------------------------------------- run
if [[ "${1:-}" == "--run" ]]; then
  say "Running pipeline (log: logs/crovir3_$(date +%F_%H%M).log)"
  mkdir -p logs
  "${SMK[@]}" 2>&1 | tee "logs/crovir3_$(date +%F_%H%M).log"
  say "Done. Start with $OUTDIR/$SPECIES/report.md"
else
  say "Dry run (nothing is executed)"
  "${SMK[@]}" -n
  say "Checks passed. To run for real:  bash run_crovir3_test.sh --run"
fi
rm -rf "$TMP"

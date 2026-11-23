#!/usr/bin/env bash
# Stream LDlink's annotation collections from the S3 mongodump (and the 2025 GWAS Catalog backup) to TSV, keeping
# dump order (column 1 = ord) and only the fields LDlink reads. Needs the LDlink account credentials in ../../.env.
# usage (from experiments/): ./fetch_annotations.sh [collection ...]   → out/annot_src/<collection>.tsv + .log
#   no arguments: every collection; otherwise only the named ones (e.g. ./fetch_annotations.sh gwas_catalog)
set -euo pipefail
set -a; source ../../.env; set +a
S3=s3://nci-cbiit-dceg-data-nonprod/ldlink
OUT=out/annot_src; mkdir -p "$OUT"
POS="chromosome_grch37 position_grch37 chromosome_grch38 position_grch38"
NUMQ="position_grch37,position_grch38,txStart,txEnd,begin_grch37,end_grch37,begin_grch38,end_grch38"
ONLY=("$@")
want() { [ ${#ONLY[@]} -eq 0 ] && return 0; local x; for x in "${ONLY[@]}"; do [ "$x" = "$1" ] && return 0; done; return 1; }

declare -A F=(
  [refGene_grch37]="bin name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds score name2 cdsStartStat cdsEndStat exonFrames"
  [refGene_grch38]="bin name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds score name2 cdsStartStat cdsEndStat exonFrames"
  [recomb]="rate $POS"
  [regulome]="score $POS"
  [forge_score]="snp_id score"
  [genes_name_coords]="name chromosome_grch37 begin_grch37 end_grch37 chromosome_grch38 begin_grch38 end_grch38"
  [gtex_tissues]="tissueSiteDetailId tissueSiteDetail tissueSiteDetailAbbr tissueSite colorHex colorRgb datasetId hasEGenes eGeneCount rnaSeqAndGenotypeSampleCount uberonId"
  [gtex_genes]="gene_id_v26 gene_name_v26 ncbi_id"
  [platforms]="platform code"
  [snp_col]="data $POS"
  [gtex_tissue_eqtl]="variant_id gene_id tss_distance ma_samples ma_count maf pval_nominal slope slope_se pval_nominal_threshold min_pval_nominal pval_beta tissueSiteDetailId $POS"
)
for c in refGene_grch37 refGene_grch38 recomb regulome forge_score genes_name_coords gtex_tissues gtex_genes platforms snp_col gtex_tissue_eqtl; do
  want "$c" || continue
  echo "== $c $(date -u +%T)"
  # LDlink queries these numerically ($gte/$lte): numeric-looking strings are written "s:..." so they stay un-keyed (B5)
  aws s3 cp "$S3/deploy/mongodb/$c.bson.gz" - | NUMERIC_QUERY="$NUMQ" node bson_stream.mjs ${F[$c]} > "$OUT/$c.tsv" 2> "$OUT/$c.log"
  tail -c 400 "$OUT/$c.log"; echo
done
# GWAS Catalog: the 2025-10-24 backup of the live collection (newer than the 2023 dump); mongoexport JSON array
if want gwas_catalog; then
  GWAS="DATE ADDED TO CATALOG|PUBMEDID|FIRST AUTHOR|DATE|JOURNAL|STUDY|DISEASE/TRAIT|MAPPED_TRAIT|MAPPED_GENE|STRONGEST SNP-RISK ALLELE|SNPS|MERGED|SNP_ID_CURRENT|CONTEXT|RISK ALLELE FREQUENCY|P-VALUE|PVALUE_MLOG|P-VALUE (TEXT)|OR or BETA|95% CI (TEXT)|STUDY ACCESSION|chromosome|position_grch37|position_grch38"
  IFS='|' read -r -a GF <<< "$GWAS"
  echo "== gwas_catalog (2025 backup) $(date -u +%T)"
  aws s3 cp "$S3/backups/LDLink.gwas_catalog20251024.json" - | NUMERIC_QUERY="$NUMQ" node json_array_stream.mjs "${GF[@]}" > "$OUT/gwas_catalog.tsv" 2> "$OUT/gwas_catalog.log"
  tail -c 600 "$OUT/gwas_catalog.log"; echo
  printf '%s\n' "${GF[@]}" > "$OUT/gwas_catalog.fields"
fi
echo "done $(date -u +%T)"

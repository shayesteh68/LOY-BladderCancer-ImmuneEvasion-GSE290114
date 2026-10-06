#!/usr/bin/env bash
# ==============================================================================
#  restore_loy_inputs.sh      version 1.0.0
#  Project : LOY / GSE290114   (repo: github.com/shayesteh68)
#
#  WHAT IT DOES (one run, one command):
#    1. finds the LOY repository on this machine and verifies its git remote,
#    2. takes an md5 snapshot of the tracked results/ directory  (READ-ONLY),
#    3. downloads the missing raw input  data/raw/GSE290114_gene_count.txt.gz
#       from NCBI GEO FTP (file listing on the server verified 2026-10-06),
#    4. rebuilds the missing       metadata/samplesheet.csv
#       from the official GEO series matrix (sample ids are taken from the
#       downloaded count matrix header, conditions from the GEO annotations),
#    5. rebuilds         data/reference/mouse_gene2entrez.tsv  if it is absent
#       from the official NCBI Mus musculus gene_info table,
#    6. validates EVERY input the scripts in scripts/ require, cross-checks the
#       new files against the existing results/, and
#       re-checks the results/ md5 snapshot (must be unchanged),
#    7. commits ONLY the restored input files, and pushes if --push is given.
#
#  It never overwrites or recomputes anything inside results/.
#  It stops with a FATAL message instead of guessing when a mapping is unclear.
#
#  USAGE
#    bash restore_loy_inputs.sh            # restore + validate + commit locally
#    bash restore_loy_inputs.sh --push     # same, and push to origin afterwards
#    bash restore_loy_inputs.sh --repo /full/path/to/repo   # only if auto-detect fails
#
#  REQUIREMENTS: bash 4+, curl, gzip, awk, sort, md5sum, git.  No R needed here.
# ==============================================================================
set -Eeuo pipefail

VERSION="1.0.0"
TS="$(date +%Y%m%d_%H%M%S)"
LOG="${HOME}/loy_restore_inputs_${TS}.log"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/loy_restore.XXXXXX")"
trap 'rm -rf "${TMP}"' EXIT

# --- verified source URLs (NCBI FTP listing checked 2026-10-06) ---------------
COUNTS_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE290nnn/GSE290114/suppl/GSE290114_gene_count.txt.gz"
MATRIX_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE290nnn/GSE290114/matrix/GSE290114_series_matrix.txt.gz"
GENEINFO_URL="https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz"

# --- repository candidates (existence-checked; nothing is assumed) ------------
REPO_GUESS_1="digital_home/LOY-BladderCancer-ImmuneEvasion-GSE290114"
REPO_GUESS_2="${HOME}/LOY-BladderCancer-ImmuneEvasion-GSE290114"

REQUIRED_ANNOT_COLS="gene_id gene_name gene_chr gene_biotype gene_description"
EXPECTED_CONDITIONS="Y_Scr Y_KO"

DO_PUSH=0
REPO_ARG=""

# --- output helper ------------------------------------------------------------
mkdir -p "$(dirname "${LOG}")"
exec > >(tee -a "${LOG}") 2>&1

log()  { printf '[%s] %s\n'          "$(date +%H:%M:%S)" "$*"; }
note() { printf '           %s\n'    "$*"; }
hdr()  { printf '\n%s\n  %s\n%s\n' "==============================================================================" "$*" "=============================================================================="; }
die()  { printf '\n[%s] FATAL: %s\n' "$(date +%H:%M:%S)" "$*"; printf 'Log kept at: %s\n' "${LOG}"; exit 1; }

usage() {
  sed -n '2,40p' "$0"
}

while [ $# -gt 0 ]; do
  case "$1" in
    --push) DO_PUSH=1 ;;
    --repo) shift; REPO_ARG="${1:-}" ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'Unknown argument: %s\n' "$1" >&2; exit 2 ;;
  esac
  shift
done

hdr "LOY / GSE290114  -  INPUT RESTORATION  (script v${VERSION})"
log "host=$HOSTNAME  shell=$BASH_VERSION  date=$(date -Iseconds)"
log "log file: ${LOG}"

# ==============================================================================
# STEP 0 - locate the repository
# ==============================================================================
hdr "STEP 0  locating the LOY repository"

resolve_repo() {
  local c hit g r
  if [ -n "${REPO_ARG}" ]; then
    [ -d "${REPO_ARG}/.git" ] || die "--repo does not point to a git checkout: ${REPO_ARG}"
    printf '%s\n' "${REPO_ARG}"; return 0
  fi
  for c in "${REPO_GUESS_1}" "${REPO_GUESS_2}"; do
    [ -n "${c}" ] || continue
    if [ -d "${c}/.git" ]; then printf '%s\n' "${c}"; return 0; fi
  done
  hit="$(find "${HOME}" -maxdepth 4 -type d -name .git 2>/dev/null | while IFS= read -r g; do
           r="$(dirname -- "${g}")"
           if git -C "${r}" remote -v 2>/dev/null | grep -q 'shayesteh68'; then printf '%s\n' "${r}"; break; fi
         done | head -n 1)"
  [ -n "${hit}" ] && printf '%s\n' "${hit}"
}

REPO="$(resolve_repo)"
[ -n "${REPO}" ] || die "Could not find the LOY repository under ${HOME}. Re-run: bash restore_loy_inputs.sh --repo /full/path/to/repo"
REPO="$(cd "${REPO}" && pwd)"
cd "${REPO}"
log "repository found: ${REPO}"
note "HEAD   : $(git rev-parse --short HEAD 2>/dev/null || echo '?')  ($(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo '?'))"
note "remote : $(git remote -v 2>/dev/null | awk 'NR==1{print $2}')"
git remote -v 2>/dev/null | grep -q 'shayesteh68' || log "WARNING: remote does not contain 'shayesteh68' - check that this is the right clone."

for d in data/raw data/reference metadata scripts results; do
  if [ -d "${d}" ]; then note "present: ${d}/"; else note "MISSING: ${d}/"; fi
done

COUNTS="${REPO}/data/raw/GSE290114_gene_count.txt.gz"
SHEET="${REPO}/metadata/samplesheet.csv"
REFS="${REPO}/data/reference/mouse_gene2entrez.tsv"

mkdir -p "${REPO}/data/raw" "${REPO}/data/reference" "${REPO}/metadata"

# ==============================================================================
# STEP 1 - freeze the existing results (md5 snapshot, read-only contract)
# ==============================================================================
hdr "STEP 1  snapshotting results/ (must not change during this run)"

snap() {
  if [ -d "${REPO}/results" ]; then
    find "${REPO}/results" -type f -print0 2>/dev/null | sort -z | xargs -0 -r md5sum 2>/dev/null | sort -k2
  fi
  return 0
}
snap > "${TMP}/results_before.md5"
if [ -s "${TMP}/results_before.md5" ]; then
  log "results/ files snapshotted: $(wc -l < "${TMP}/results_before.md5")"
  while read -r sum path; do note "${sum}  ${path#${REPO}/}"; done < "${TMP}/results_before.md5"
else
  log "WARNING: results/ is empty or absent - nothing to protect."
fi

# ==============================================================================
# STEP 2 - download the raw count matrix
# ==============================================================================
hdr "STEP 2  raw count matrix (missing input #1)"

fetch() {
  local url="$1" dest="$2"
  if [ -s "${dest}" ] && gzip -t "${dest}" 2>/dev/null; then
    log "already present and valid, download skipped: ${dest#${REPO}/}"
    return 0
  fi
  log "downloading $(basename "${dest}")"
  note "${url}"
  curl -fL --retry 3 --retry-delay 3 --connect-timeout 20 --progress-bar -o "${dest}.part" "${url}" \
    || die "download failed: ${url}"
  mv -f "${dest}.part" "${dest}"
  gzip -t "${dest}" 2>/dev/null || die "downloaded file is not a valid gzip: ${dest}"
  log "saved: ${dest#${REPO}/}   size=$(du -h "${dest}" | cut -f1)"
}

fetch "${COUNTS_URL}" "${COUNTS}"
log "md5(${COUNTS#${REPO}/}) = $(md5sum "${COUNTS}" | cut -d' ' -f1)"

# --- header / schema of the count matrix --------------------------------------
gzip -dc "${COUNTS}" | head -n 1 > "${TMP}/counts_header.txt"
HEADER="$(cat "${TMP}/counts_header.txt")"
NCOL="$(awk -F'\t' '{print NF; exit}' "${TMP}/counts_header.txt")"
log "count matrix: ${NCOL} columns in the header"

idx_of_col() {
  awk -F'\t' -v want="$1" '{for(i=1;i<=NF;i++) if($i==want){print i; exit}}' "${TMP}/counts_header.txt"
}

GENE_ID_COL=""
for c in ${REQUIRED_ANNOT_COLS}; do
  i="$(idx_of_col "${c}")"
  [ -n "${i}" ] || die "count matrix header has no column named '${c}' - scripts/01_deseq2_analysis.R requires it."
  note "annotation column ok: ${c}  (position ${i})"
  [ "${c}" = "gene_id" ] && GENE_ID_COL="${i}"
done

: > "${TMP}/sample_cols.txt"
awk -F'\t' -v keep="$(printf '%s\t' ${REQUIRED_ANNOT_COLS})" '
  {for(i=1;i<=NF;i++){ if(index(keep,"\t" $i "\t")==0 && (i==1 || $i!=prev_keep)) print $i }
   exit}' "${TMP}/counts_header.txt" > "${TMP}/raw_cols.txt"
# robust selection: print every header name that is not one of the required annotation names
awk -F'\t' -v keep="$(printf '%s\t' ${REQUIRED_ANNOT_COLS})" '{
  for(i=1;i<=NF;i++){ k=1; n=split(keep,a,"\t"); for(j=1;j<=n;j++){ if($i==a[j]) k=0 } if(k==1) print $i }
}' "${TMP}/counts_header.txt" > "${TMP}/sample_cols.txt"

NSAMPLE_COLS="$(wc -l < "${TMP}/sample_cols.txt")"
[ "${NSAMPLE_COLS}" -ge 1 ] || die "no sample columns detected in the count matrix header."
log "sample columns detected in the count matrix: ${NSAMPLE_COLS}"
sed 's/^/           - /' "${TMP}/sample_cols.txt"

DUPES="$(sort "${TMP}/sample_cols.txt" | uniq -d | paste -sd, - || true)"
[ -z "${DUPES}" ] || die "duplicate sample column names in the count matrix header: ${DUPES}"

NROW="$(gzip -dc "${COUNTS}" | wc -l)"
NENS="$(gzip -dc "${COUNTS}" | awk -F'\t' -v c="${GENE_ID_COL}" '$c ~ /^ENSMUSG/ {n++} END{print n+0}')"
log "data rows (incl. header): ${NROW}   rows whose gene_id starts with ENSMUSG: ${NENS}"
[ "${NENS}" -gt 0 ] || log "WARNING: no ENSMUSG* values found in the gene_id column - inspect the file before running the pipeline."

# ==============================================================================
# STEP 3 - samplesheet (missing input #2), built from the GEO series matrix
# ==============================================================================
hdr "STEP 3  metadata/samplesheet.csv (missing input #2)"

log "downloading the official series matrix"
note "${MATRIX_URL}"
curl -fL --retry 3 --retry-delay 3 --connect-timeout 20 -o "${TMP}/series_matrix.txt.gz" "${MATRIX_URL}" \
  || die "download failed: ${MATRIX_URL}"
gzip -t "${TMP}/series_matrix.txt.gz" || die "series matrix is not a valid gzip file"
gzip -dc "${TMP}/series_matrix.txt.gz" > "${TMP}/series_matrix.txt"
[ -s "${TMP}/series_matrix.txt" ] || die "series matrix is empty"

# >>> SERIES_AWK
awk -F'\t' '
  function clean(v){ gsub(/^"|"$/,"",v); return v }
  $1=="!Sample_geo_accession"        { n=0; for(i=2;i<=NF;i++){ gsm[i-1]=clean($i); n++ } next }
  $1=="!Sample_title"                { for(i=2;i<=NF;i++) ttl[i-1]=clean($i); next }
  $1=="!Sample_source_name_ch1"      { for(i=2;i<=NF;i++) src[i-1]=clean($i); next }
  $1=="!Sample_characteristics_ch1"  { for(i=2;i<=NF;i++){ v=clean($i); ch[i-1] = (ch[i-1]=="" ? v : ch[i-1] " | " v) } next }
  END{
    for(i=1;i<=n;i++) printf "%s\t%s\t%s\t%s\n", gsm[i], ttl[i], src[i], ch[i]
  }
' "${TMP}/series_matrix.txt" > "${TMP}/geo.tsv"
# <<< SERIES_AWK

NGEO="$(wc -l < "${TMP}/geo.tsv")"
[ "${NGEO}" -ge 1 ] || die "the series matrix contained no !Sample_geo_accession rows."
log "GEO samples in the series matrix: ${NGEO}"
while IFS=$'\t' read -r g t s c; do note "${g} | title=${t} | source=${s} | ${c}"; done < "${TMP}/geo.tsv"

# >>> MAP_COND
map_condition() {
  local hay
  hay="$(printf '%s' "$1" | tr '[:upper:]' '[:lower:]')"
  if printf '%s' "${hay}" | grep -Eq '(^|[^a-z])(ko|knockout|knock-out|knock out|crispr|mutant|mut)([^a-z]|$)'; then
    printf 'Y_KO\n'; return 0
  fi
  if printf '%s' "${hay}" | grep -Eq '(scr|scrambl|control|ctrl|(^|[^a-z])(wt|wildtype|wild)[^a-z]|non-?target)'; then
    printf 'Y_Scr\n'; return 0
  fi
  printf 'UNKNOWN\n'
}
# <<< MAP_COND

: > "${TMP}/map.tsv"
if [ -s "${SHEET}" ]; then
  log "samplesheet already exists - it will be VALIDATED, not overwritten"
  awk -F',' 'NR>1{ v=$1; gsub(/^"|"$/,"",v); print v }' "${SHEET}" | sed '/^$/d' > "${TMP}/sheet_ids.txt"
else
  log "samplesheet absent - rebuilding it from the GEO annotations"
  while IFS= read -r col; do
    [ -n "${col}" ] || continue
    lc="$(printf '%s' "${col}" | tr '[:upper:]' '[:lower:]')"
    line="$(awk -F'\t' -v k="${lc}" 'tolower($1)==k || tolower($2)==k {print; exit}' "${TMP}/geo.tsv")"
    how="exact"
    if [ -z "${line}" ]; then
      line="$(awk -F'\t' -v k="${lc}" 'index(tolower($0),k)>0 {print; exit}' "${TMP}/geo.tsv")"
      how="substring"
    fi
    [ -n "${line}" ] || die "cannot map count-matrix column '${col}' to any GEO sample. GSMs available: $(cut -f1 "${TMP}/geo.tsv" | paste -sd, -)"
    cond="$(map_condition "${line}")"
    [ "${cond}" != "UNKNOWN" ] || die "column '${col}' maps to GEO row '$(printf '%s' "${line}" | cut -f1,2)' but its condition is not identifiable (no KO / Scr signal). Review ${TMP}/geo.tsv and the GEO page, then decide manually."
    printf '%s\t%s\t%s\t%s\t%s\n' "${col}" "${cond}" "$(printf '%s' "${line}" | cut -f1)" "$(printf '%s' "${line}" | cut -f2-)" "${how}" >> "${TMP}/map.tsv"
  done < "${TMP}/sample_cols.txt"

  csv_field() { case "$1" in *[,\"]*) printf '"%s"' "$(printf '%s' "$1" | sed 's/"/""/g')" ;; *) printf '%s' "$1" ;; esac }
  {
    printf 'sample_id,condition\n'
    while IFS=$'\t' read -r sid cond rest; do
      printf '%s,%s\n' "$(csv_field "${sid}")" "$(csv_field "${cond}")"
    done < "${TMP}/map.tsv"
  } > "${SHEET}"
  {
    printf 'sample_id\tcondition\tgsm\ttitle\tsource_name\tcharacteristics\tmatch\n'
    while IFS=$'\t' read -r sid cond gsm rest how; do
      printf '%s\t%s\t%s\t%s\t%s\n' "${sid}" "${cond}" "${gsm}" "${rest}" "${how}"
    done < "${TMP}/map.tsv"
  } > "${REPO}/metadata/GEO_annotation_source.tsv"
  log "written: metadata/samplesheet.csv"
  log "written: metadata/GEO_annotation_source.tsv (provenance, not read by the pipeline)"
  log "condition mapping produced:"
  while IFS=$'\t' read -r sid cond gsm rest how; do note "${sid}  ->  ${cond}   (${gsm}, ${how})"; done < "${TMP}/map.tsv"
  cut -d',' -f1 "${SHEET}" | tail -n +2 | sed '/^$/d' | sed 's/^"//; s/"$//' > "${TMP}/sheet_ids.txt"
fi

# ==============================================================================
# STEP 4 - Entrez reference map (only if absent)
# ==============================================================================
hdr "STEP 4  data/reference/mouse_gene2entrez.tsv"

if [ -s "${REFS}" ]; then
  log "already present - left untouched ($(wc -l < "${REFS}") lines): ${REFS#${REPO}/}"
else
  log "absent - rebuilding from the official NCBI table (columns: Symbol, GeneID)"
  note "${GENEINFO_URL}"
  curl -fL --retry 3 --retry-delay 3 --connect-timeout 20 -o "${TMP}/gene_info.gz" "${GENEINFO_URL}" \
    || die "download failed: ${GENEINFO_URL}"
  {
    printf 'Symbol\tGeneID\n'
    gzip -dc "${TMP}/gene_info.gz" | awk -F'\t' '$1 !~ /^#/ && $2!="" && $3!="" {print $3 "\t" $2}' | sort -u
  } > "${REFS}"
  log "written: ${REFS#${REPO}/}  ($(wc -l < "${REFS}") lines)"
fi

# ==============================================================================
# STEP 5 - full input validation against what scripts/ actually require
# ==============================================================================
hdr "STEP 5  validation"

FAIL=0
ok()   { printf '  [ OK ]   %s\n' "$*"; }
warn() { printf '  [WARN]   %s\n' "$*"; }
bad()  { printf '  [FAIL]   %s\n' "$*"; FAIL=$((FAIL+1)); }

# 5.1 columns of the count matrix
for c in ${REQUIRED_ANNOT_COLS}; do
  [ -n "$(idx_of_col "${c}")" ] && ok "count matrix column present: ${c}" || bad "count matrix column missing: ${c}"
done

# 5.2 samplesheet schema and cross-check with the count matrix
if [ -s "${SHEET}" ]; then
  ok "samplesheet exists: metadata/samplesheet.csv ($(($(wc -l < "${SHEET}")-1)) data rows)"
  shdr="$(head -n 1 "${SHEET}" | tr -d '\r' | tr -d '"')"
  case "${shdr}" in
    sample_id,condition*) ok "samplesheet header starts with: ${shdr}" ;;
    *) bad "samplesheet header is '${shdr}' - scripts/01_deseq2_analysis.R expects the first columns to be sample_id,condition" ;;
  esac
  cut -d',' -f2 "${SHEET}" | tail -n +2 | sed 's/^"//; s/"$//' | sed '/^$/d' | sort -u > "${TMP}/conds.txt"
  while IFS= read -r cc; do
    case " ${EXPECTED_CONDITIONS} " in *" ${cc} "*) ok "condition value used: ${cc}" ;; *) bad "unexpected condition value: '${cc}' (expected ${EXPECTED_CONDITIONS})" ;; esac
    n=$(grep -c ",\"?${cc}\"?$" "${SHEET}" || true); note "n(${cc}) = ${n}"
    [ "${n}" -ge 2 ] || bad "condition '${cc}' has only ${n} sample(s); DESeq2 needs at least 2 per condition."
  done < "${TMP}/conds.txt"
  sort "${TMP}/sheet_ids.txt" > "${TMP}/a.txt"
  sort "${TMP}/sample_cols.txt" > "${TMP}/b.txt"
  if diff -q "${TMP}/a.txt" "${TMP}/b.txt" >/dev/null; then
    ok "sample ids in the samplesheet are an exact set match with the count matrix columns (${NSAMPLE_COLS} samples)"
  else
    bad "samplesheet sample ids do not match the count matrix columns exactly:"
    diff "${TMP}/a.txt" "${TMP}/b.txt" | sed 's/^/           /' || true
  fi
else
  bad "metadata/samplesheet.csv is missing"
fi

# 5.3 the four pipeline scripts and the reference they consume
for s in 01_deseq2_analysis.R 02_tf_enrichment.R 03_immune_deconvolution.R 04_pathway_analysis.R; do
  if [ -s "${REPO}/scripts/${s}" ]; then ok "script present: scripts/${s}"; else warn "script not found (name may differ): scripts/${s}"; fi
done
[ -s "${REFS}" ] && ok "entrez reference present: data/reference/mouse_gene2entrez.tsv" || bad "entrez reference missing"
if [ -s "${REFS}" ]; then
  head -n 1 "${REFS}" | grep -q 'Symbol' && head -n 1 "${REFS}" | grep -q 'GeneID' \
    && ok "entrez reference has the Symbol / GeneID columns script 02 merges on" \
    || bad "entrez reference header must contain the columns Symbol and GeneID"
fi

# 5.4 new inputs vs. the results that were produced earlier
DEG="${REPO}/results/deseq2_deg_significant.csv"
if [ -s "${DEG}" ] && [ "${NENS}" -gt 0 ]; then
  grep -oE 'ENSMUSG[0-9]+' "${DEG}" | sort -u | head -n 200 > "${TMP}/deg_ids.txt"
  nprobe="$(wc -l < "${TMP}/deg_ids.txt")"
  if [ "${nprobe}" -gt 0 ]; then
    gzip -dc "${COUNTS}" | cut -f"${GENE_ID_COL}" | grep -Fx -f "${TMP}/deg_ids.txt" | sort -u > "${TMP}/deg_found.txt"
    nfound="$(wc -l < "${TMP}/deg_found.txt")"
    if [ "${nfound}" -eq "${nprobe}" ]; then
      ok "independent cross-check: all ${nprobe} gene ids sampled from the existing results/deseq2_deg_significant.csv are present in the downloaded count matrix"
    elif [ "${nfound}" -gt 0 ]; then
      warn "cross-check: ${nfound}/${nprobe} gene ids from the existing results were found in the downloaded count matrix"
    else
      warn "cross-check: none of the sampled gene ids from the existing results were found - verify this is the correct GEO file"
    fi
  fi
else
  warn "cross-check skipped (results file or gene ids unavailable)"
fi

# 5.5 results directory unchanged
snap > "${TMP}/results_after.md5"
if diff -q "${TMP}/results_before.md5" "${TMP}/results_after.md5" >/dev/null; then
  ok "results/ is byte-identical to the snapshot taken at the start of this run"
else
  bad "results/ changed during this run - investigate immediately:"
  diff "${TMP}/results_before.md5" "${TMP}/results_after.md5" | sed 's/^/           /' || true
fi
[ "${FAIL}" -eq 0 ] && ok "ALL INPUT CHECKS PASSED" || die "${FAIL} validation check(s) failed - see the [FAIL] lines above. Nothing was committed."

# ==============================================================================
# STEP 6 - commit the restored inputs (results/ untouched)
# ==============================================================================
hdr "STEP 6  commit"

git add -f -- "data/raw/$(basename "${COUNTS}")" "metadata/samplesheet.csv" 2>/dev/null || true
if [ -s "${REPO}/metadata/GEO_annotation_source.tsv" ]; then git add -f -- "metadata/GEO_annotation_source.tsv" 2>/dev/null || true; fi
if [ -s "${REFS}" ]; then git add -f -- "data/reference/mouse_gene2entrez.tsv" 2>/dev/null || true; fi

log "staged files:"
git diff --cached --name-only | sed 's/^/           /' || true

if [ -n "$(git diff --cached --name-only)" ]; then
  git -c user.name="$(git config user.name || echo restore-tool)" \
      -c user.email="$(git config user.email || echo restore-tool@localhost)" \
      commit -m "inputs: restore GSE290114 raw counts matrix, samplesheet and Entrez map (NCBI-verified sources, all input checks passed)" \
    || die "git commit failed"
  log "commit created: $(git rev-parse --short HEAD)"
else
  log "nothing new to commit (files already tracked and unchanged)"
fi

if [ "${DO_PUSH}" -eq 1 ]; then
  log "pushing to origin"
  git push || die "git push failed - resolve it and re-run with --push"
  ok "pushed: $(git rev-parse --abbrev-ref HEAD) -> origin"
else
  note "push not requested. To publish: cd \"${REPO}\" && git push"
fi

# ==============================================================================
# STEP 7 - R environment status (informational; nothing here needs R)
# ==============================================================================
hdr "STEP 7  R environment status for the pipeline scripts"
if command -v Rscript >/dev/null 2>&1; then
  log "Rscript on PATH: $(Rscript --version 2>&1 | head -n 1)"
else
  warn "Rscript is NOT on PATH in this shell."
fi
if command -v conda >/dev/null 2>&1; then
  log "conda detected; environments whose name matches /r_/ :"
  conda env list 2>/dev/null | awk '/(^r_|deseq)/{print "           " $0}' || true
  note "if the pipeline environment exists, run:  conda activate r_deseq_env"
else
  note "conda not on PATH in this shell."
fi
for p in DESeq2 clusterProfiler pheatmap ggplot2 org.Mm.eg.db; do
  if command -v Rscript >/dev/null 2>&1; then
    v="$(Rscript -e "cat(if(requireNamespace('${p}',quietly=TRUE)) as.character(packageVersion('${p}')) else 'NOT-INSTALLED')" 2>/dev/null || echo 'check-failed')"
    printf '           %-16s %s\n' "${p}" "${v}"
  fi
done

# ==============================================================================
# SUMMARY
# ==============================================================================
hdr "SUMMARY"
log "new/verified inputs:"
note "data/raw/$(basename "${COUNTS}")   md5=$(md5sum "${COUNTS}" | cut -d' ' -f1)   size=$(du -h "${COUNTS}" | cut -f1)"
[ -s "${SHEET}" ] && note "metadata/samplesheet.csv            md5=$(md5sum "${SHEET}" | cut -d' ' -f1)   rows=$(($(wc -l < "${SHEET}")-1))"
[ -s "${REFS}" ]  && note "data/reference/mouse_gene2entrez.tsv md5=$(md5sum "${REFS}" | cut -d' ' -f1)"
log "results/ : untouched and md5-verified against the pre-run snapshot."
log "full transcript: ${LOG}"
log "next (only when you decide to re-run the pipeline):"
note "cd \"${REPO}\" && conda activate r_deseq_env && Rscript scripts/01_deseq2_analysis.R"
log "DONE"

#!/usr/bin/env bash
# ==============================================================================
#  restore_loy_inputs.sh      version 1.0.5
#  Project : LOY / GSE290114   (repo: github.com/shayesteh68)
#
#  CHANGELOG
  #    1.0.5 - root cause of the same STEP 3 abort.
  #            * the annotation filter that 1.0.4 added just before the locked map is
  #              now applied to the header scan itself, so the count-matrix sample
  #              list is correct the moment it is built (STEP 2 now reports 6 sample
  #              columns, not 7).
  #            * header field names are sanitised (invisible bytes such as CR and
  #              whitespace, plus letter case, removed) before they are compared both
  #              with ANNOT_COLS_ALL and with the locked map. A trailing CR or space
  #              on the last header field - which defeats an exact string comparison -
  #              can no longer misclassify the annotation column 'tf_family' as a
  #              sample. Verified on synthetic headers: tf_family+CR, tf_family+space,
  #              TF_family and tf_family are all excluded, while an unknown extra
  #              column is still kept and therefore still stops the run at the locked
  #              map guard.
  #            * the header line is read with CR stripped, and the sample count used
  #              by the STEP 5 set-match report is recomputed after the filter so that
  #              every consumer sees the same list.
  #    1.0.4 - samplesheet-rebuild fix for the observed STEP 3 abort.
  #            * the count-matrix header scan let the trailing annotation column
  #              'tf_family' through as a sample column, so STEP 3 tried to map it
  #              through the locked table and stopped with "no entry in the locked
  #              GEO map". Any detected column whose name (case- and space-
  #              insensitive) occurs in the authoritative annotation list
  #              ANNOT_COLS_ALL is now excluded before the locked map is applied.
  #              A column that is neither annotation nor part of the locked map
  #              still stops the run, so nothing is ever guessed.
#    1.0.3 - unblocking release.
#            * the library-name cross-check no longer aborts the run. It still runs
#              and its result is recorded in metadata/GEO_annotation_source.tsv, but
#              when the GEO description field does not literally contain the matrix
#              column name (for example the title is "CRISPR-YScr Replicate 1"
#              rather than "RC_1") the script WARNS and proceeds with the locked
#              map. The sample group itself is still cross-verified from the GEO
#              genotype field, so a wrong group can never be written silently.
#    1.0.2 - GEO column-mapping release.
#            * the annotation-column list now covers all 10 non-sample columns of
#              GSE290114_gene_count.txt (gene_start/end/strand/length and tf_family
#              were missing, so 11 columns were mistaken for samples instead of 6).
#            * STEP 3 maps matrix columns through an authoritative, verified table
#              (Library name -> GSM -> group) and cross-checks it against the series
#              matrix; the word "crispr" is no longer treated as a KO signal.
#    1.0.1 - terminal-safety release, fixes the observed exit code 141.
#            * the logging line "exec > >(tee -a LOG) 2>&1" is gone. A process
#              substitution on stdout plus "set -o pipefail" is what allowed a
#              downstream SIGPIPE to end the whole terminal shell with 141.
#            * every "writer | head" pipe is gone, because such a pipe makes the
#              writer die of SIGPIPE (128+13 = 141):
#                 gzip -dc counts.gz | head -n 1
#                 grep -oE ... | sort -u | head -n 200
#                 find ~ | while ... | head -n 1
#              All three now read their input to the end, which cannot raise
#              SIGPIPE. This was the line that killed the run.
#            * logging appends to the log file directly, with no pipe at all.
#            * curl no longer paints the \r progress bar, which overwrote the
#              real message on the terminal line.
#            * an ERR trap now prints the exact failing command and line.
#            * refuses to run when source'd / dot-sourced.
#    1.0.0 - first release.
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
# USAGE_BEGIN
#    bash restore_loy_inputs.sh            # restore + validate + commit locally
#    bash restore_loy_inputs.sh --push     # same, and push to origin afterwards
#    bash restore_loy_inputs.sh --repo /full/path/to/repo   # only if auto-detect fails
#
#    Run it with "bash", never with "source" or "." :
#      a source'd script runs inside YOUR terminal shell, so a single failure
#      in it terminates the terminal itself. That is the exact mechanism behind
#      "The terminal process /bin/bash terminated with exit code: 141".
# USAGE_END
#
#  REQUIREMENTS: bash 4+, curl, gzip, awk, sort, md5sum, git.  No R needed here.
# ==============================================================================

# --- refuse to run inside the interactive shell (source / . script.sh) --------
if [ "${BASH_SOURCE[0]}" != "$0" ]; then
  printf 'FATAL: this script must be run as a command, not source'"'"'d.\n'
  printf 'Use:  bash %s [--push] [--repo /path/to/repo]\n' "${BASH_SOURCE[0]}"
  return 1 2>/dev/null || exit 1
fi

set -Eeuo pipefail

VERSION="1.0.6"
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
# v1.0.2: EVERY non-sample column of GSE290114_gene_count.txt must be listed here.
ANNOT_COLS_ALL="gene_id gene_name gene_chr gene_start gene_end gene_strand gene_length gene_biotype gene_description tf_family"
EXPECTED_CONDITIONS="Y_Scr Y_KO"
# v1.0.2: authoritative column -> GSM -> group map, locked from the GEO series matrix
# (Library name: RC_* = CRISPR-YScr, RP_* = CRISPR Y-KO).
LOCKED_MAP="RC_1 GSM8806400 Y_Scr
RC_2 GSM8806401 Y_Scr
RC_3 GSM8806402 Y_Scr
RP_1 GSM8806403 Y_KO
RP_2 GSM8806404 Y_KO
RP_3 GSM8806405 Y_KO"

DO_PUSH=0
REPO_ARG=""
DYING=0

# --- logging: append to the log file directly, never through a pipe -----------
: > "${LOG}" 2>/dev/null || true
_emit() { printf '%s\n' "$1"; printf '%s\n' "$1" >> "${LOG}" 2>/dev/null || true; }
log()  { _emit "[$(date +%H:%M:%S)] $*"; }
note() { _emit "           $*"; }
out()  { local __l; while IFS= read -r __l; do _emit "${__l}"; done; }
hdr()  { _emit ""; _emit "=============================================================================="; _emit "  $*"; _emit "=============================================================================="; }
die()  { DYING=1; _emit ""; _emit "[$(date +%H:%M:%S)] FATAL: $*"; _emit "Log kept at: ${LOG}"; exit 1; }

on_err() {
  local rc=$?
  [ "${DYING}" -eq 1 ] && return 0
  _emit ""
  _emit "[$(date +%H:%M:%S)] ERROR: a command exited with status ${rc} (line ${LINENO}): ${BASH_COMMAND}"
  _emit "           the script stops here; nothing beyond this point was committed."
  _emit "Log kept at: ${LOG}"
  return 0
}
trap on_err ERR

usage() { awk '/^# USAGE_BEGIN/{f=1;next} /^# USAGE_END/{f=0} f' "$0"; }

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
log "host=${HOSTNAME}  shell=${BASH_VERSION}  date=$(date -Iseconds)"
log "log file: ${LOG}"

# ==============================================================================
# STEP 0 - locate the repository
# ==============================================================================
hdr "STEP 0  locating the LOY repository"

resolve_repo() {
  local c hit="" g r
  if [ -n "${REPO_ARG}" ]; then
    [ -d "${REPO_ARG}/.git" ] || die "--repo does not point to a git checkout: ${REPO_ARG}"
    printf '%s\n' "${REPO_ARG}"; return 0
  fi
  for c in "${REPO_GUESS_1}" "${REPO_GUESS_2}"; do
    [ -n "${c}" ] || continue
    if [ -d "${c}/.git" ]; then printf '%s\n' "${c}"; return 0; fi
  done
  find "${HOME}" -maxdepth 4 -type d -name .git 2>/dev/null > "${TMP}/gitdirs.txt" || true
  while IFS= read -r g; do
    r="$(dirname -- "${g}")"
    if git -C "${r}" remote -v 2>/dev/null | grep -q 'shayesteh68'; then hit="${r}"; break; fi
  done < "${TMP}/gitdirs.txt"
  [ -n "${hit}" ] && printf '%s\n' "${hit}"
  return 0
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
  curl -fL --retry 3 --retry-delay 3 --connect-timeout 20 --no-progress-meter -o "${dest}.part" "${url}" \
    || die "download failed: ${url}"
  mv -f "${dest}.part" "${dest}"
  gzip -t "${dest}" 2>/dev/null || die "downloaded file is not a valid gzip: ${dest}"
  log "saved: ${dest#${REPO}/}   size=$(du -h "${dest}" | cut -f1)"
}

fetch "${COUNTS_URL}" "${COUNTS}"
log "md5(${COUNTS#${REPO}/}) = $(md5sum "${COUNTS}" | cut -d' ' -f1)"

# --- header / schema of the count matrix (no early-exit reader: SIGPIPE-safe) --
gzip -dc "${COUNTS}" | sed -n '1p' | tr -d '\r' > "${TMP}/counts_header.txt"   # v1.0.5: CR-safe header
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

awk -F'\t' -v keep="$(printf '%s\t' ${ANNOT_COLS_ALL})" '
  # v1.0.5: sanitise every header field before comparing it with the annotation
  #          list, so invisible bytes (CR/whitespace) or letter case can no longer
  #          let an annotation column through as a sample column.
  function sanitise(v,   i,c,out){ out=""; for(i=1;i<=length(v);i++){ c=substr(v,i,1); if(c ~ /[A-Za-z0-9_.-]/) out=out c } return out }
  function lower(v){ return tolower(v) }
  BEGIN{ nk=split(keep,a,"\t"); for(j=1;j<=nk;j++) kk[j]=lower(sanitise(a[j])) }
  {
    for(i=1;i<=NF;i++){
      name=sanitise($i)
      if(name=="") continue
      key=lower(name); isanno=0
      for(j=1;j<=nk;j++){ if(kk[j]!="" && key==kk[j]){ isanno=1; break } }
      if(isanno==0) print name
    }
  }' "${TMP}/counts_header.txt" > "${TMP}/sample_cols.txt"

NSAMPLE_COLS="$(wc -l < "${TMP}/sample_cols.txt")"
[ "${NSAMPLE_COLS}" -ge 1 ] || die "no sample columns detected in the count matrix header."
log "sample columns detected in the count matrix: ${NSAMPLE_COLS}"
sed 's/^/           - /' "${TMP}/sample_cols.txt" | out

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
curl -fL --retry 3 --retry-delay 3 --connect-timeout 20 --no-progress-meter -o "${TMP}/series_matrix.txt.gz" "${MATRIX_URL}" \
  || die "download failed: ${MATRIX_URL}"
gzip -t "${TMP}/series_matrix.txt.gz" || die "series matrix is not a valid gzip file"
gzip -dc "${TMP}/series_matrix.txt.gz" > "${TMP}/series_matrix.txt"
[ -s "${TMP}/series_matrix.txt" ] || die "series matrix is empty"

awk -F'\t' '
  function clean(v){ gsub(/^"|"$/,"",v); return v }
  $1=="!Sample_geo_accession"        { n=0; for(i=2;i<=NF;i++){ gsm[i-1]=clean($i); n++ } next }
  $1=="!Sample_title"                { for(i=2;i<=NF;i++) ttl[i-1]=clean($i); next }
  $1=="!Sample_source_name_ch1"      { for(i=2;i<=NF;i++) src[i-1]=clean($i); next }
  $1=="!Sample_characteristics_ch1"  { for(i=2;i<=NF;i++){ v=clean($i); ch[i-1] = (ch[i-1]=="" ? v : ch[i-1] " | " v) } next }
  $1=="!Sample_description"          { for(i=2;i<=NF;i++){ v=clean($i); lib[i-1] = (lib[i-1]=="" ? v : lib[i-1] " | " v) } next }
  END{
    for(i=1;i<=n;i++) printf "%s\t%s\t%s\t%s\t%s\n", gsm[i], ttl[i], src[i], ch[i], lib[i]
  }
' "${TMP}/series_matrix.txt" > "${TMP}/geo.tsv"

NGEO="$(wc -l < "${TMP}/geo.tsv")"
[ "${NGEO}" -ge 1 ] || die "the series matrix contained no !Sample_geo_accession rows."
log "GEO samples in the series matrix: ${NGEO}"
while IFS=$'\t' read -r g t s c l; do note "${g} | title=${t} | source=${s} | ${c} | ${l}"; done < "${TMP}/geo.tsv"

map_condition() {
  local hay
  hay="$(printf '%s' "$1" | tr '[:upper:]' '[:lower:]')"
  # v1.0.2: "crispr" is NOT a knockout signal - both arms of GSE290114 are CRISPR
  # screens, so it used to mislabel the YScr controls as Y_KO.
  if printf '%s' "${hay}" | grep -Eq '(y[-_ ]?ko|[^a-z]ko[^a-z]|knockout|knock-?out|knock out|(^|[^a-z])mut(ant)?([^a-z]|$))'; then
    printf 'Y_KO\n'; return 0
  fi
  if printf '%s' "${hay}" | grep -Eq '(y[-_ ]?scr|yscr|[^a-z]scr|scrambl|control|ctrl|non-?target)'; then
    printf 'Y_Scr\n'; return 0
  fi
  printf 'UNKNOWN\n'
}

: > "${TMP}/map.tsv"
if [ -s "${SHEET}" ]; then
  log "samplesheet already exists - it will be VALIDATED, not overwritten"
  awk -F',' 'NR>1{ v=$1; gsub(/^"|"$/,"",v); print v }' "${SHEET}" | sed '/^$/d' > "${TMP}/sheet_ids.txt"
else
  log "samplesheet absent - rebuilding it from the GEO annotations"
  # v1.0.4: the count matrix has 10 annotation columns + 6 sample columns. The
  #          header scan can let an annotation column through as a "sample"
  #          (observed: the trailing 'tf_family' column of
  #          GSE290114_gene_count.txt); the locked map has no entry for it, so the
  #          whole run used to abort. Any detected column whose name matches the
  #          authoritative annotation list is dropped here. A column that is
  #          neither annotation nor locked still stops at the guard below.
  : > "${TMP}/sample_cols.clean.txt"
  while IFS= read -r col_chk; do
    [ -n "${col_chk}" ] || continue
    chk_key="$(printf '%s' "${col_chk}" | tr -d '[:space:]' | tr '[:upper:]' '[:lower:]')"
    anno_key="$(printf '%s' "${ANNOT_COLS_ALL}" | tr '[:upper:]' '[:lower:]')"
    case " ${anno_key} " in
      *" ${chk_key} "*) log "excluded annotation column '${col_chk}' from the samplesheet"; continue ;;
    esac
    printf '%s\n' "${col_chk}" >> "${TMP}/sample_cols.clean.txt"
  done < "${TMP}/sample_cols.txt"
  mv -f "${TMP}/sample_cols.clean.txt" "${TMP}/sample_cols.txt"
  NSAMPLE_COLS="$(wc -l < "${TMP}/sample_cols.txt")"
  while IFS= read -r col; do
    [ -n "${col}" ] || continue
    lock=""
    while IFS=' ' read -r m_col m_gsm m_cond; do
      if [ "${m_col}" = "${col}" ]; then lock="${m_col} ${m_gsm} ${m_cond}"; break; fi
    done <<< "${LOCKED_MAP}"
    [ -n "${lock}" ] || die "count-matrix column '${col}' has no entry in the locked GEO map (RC_1..RC_3 = CRISPR-YScr, RP_1..RP_3 = CRISPR Y-KO). Refusing to guess a group."
    gsm="${lock#* }"; cond="${gsm##* }"; gsm="${gsm%% *}"
    geo_row="$(awk -F'\t' -v g="${gsm}" '$1==g{print; exit}' "${TMP}/geo.tsv")"
    [ -n "${geo_row}" ] || die "locked GSM ${gsm} for column '${col}' is absent from the GEO series matrix."
    lib="$(printf '%s\n' "${geo_row}" | cut -f5 | tr '[:upper:]' '[:lower:]')"
    liball="$(printf '%s\n' "${geo_row}" | tr '[:upper:]' '[:lower:]')"
    lib_confirmed="no"
    case "${lib}" in
      *"${col,,}"*) lib_confirmed="yes" ;;
    esac
    if [ "${lib_confirmed}" = "no" ]; then
      case "${liball}" in
        *"${col,,}"*) lib_confirmed="yes" ;;
      esac
    fi
    if [ "${lib_confirmed}" = "no" ]; then
      log "WARNING: the GEO series matrix does not literally name column '${col}' for ${gsm}; proceeding on the locked map. The group is still cross-checked against the GEO genotype below. GEO description: '$(printf '%s\n' "${geo_row}" | cut -f5)'"
    fi
    geo_cond="$(map_condition "${geo_row}")"
    [ "${geo_cond}" = "${cond}" ] || die "group disagreement for column '${col}': locked=${cond}, GEO=${geo_cond} (${gsm}, '$(printf '%s\n' "${geo_row}" | cut -f2)'). Refusing to guess."
    prov="$(printf '%s\n' "${geo_row}" | cut -f2) | $(printf '%s\n' "${geo_row}" | cut -f3) | $(printf '%s\n' "${geo_row}" | cut -f5)"
    printf '%s\t%s\t%s\t%s\t%s\n' "${col}" "${cond}" "${gsm}" "${prov}" "locked+geo-verified(lib=${lib_confirmed})" >> "${TMP}/map.tsv"
  done < "${TMP}/sample_cols.txt"

  csv_field() {
    local v="$1"
    case "${v}" in
      *","*|*'"'*)
        local w="${v//\"/\"\"}"
        printf '"%s"' "${w}"
        ;;
      *)
        printf '%s' "${v}"
        ;;
    esac
  }
  {
    printf 'sample_id,condition\n'
    while IFS="$(printf '\t')" read -r sid cond gsm rest how; do
      printf '%s,%s\n' "$(csv_field "${sid}")" "$(csv_field "${cond}")"
    done < "${TMP}/map.tsv"
  } > "${SHEET}"
  {
    printf 'sample_id\tcondition\tgsm\tprovenance\tmatch\n'
    while IFS="$(printf '\t')" read -r sid cond gsm rest how; do
      printf '%s\t%s\t%s\t%s\t%s\n' "${sid}" "${cond}" "${gsm}" "${rest}" "${how}"
    done < "${TMP}/map.tsv"
  } > "${REPO}/metadata/GEO_annotation_source.tsv"
  log "written: metadata/samplesheet.csv"
  log "written: metadata/GEO_annotation_source.tsv (provenance, not read by the pipeline)"
  log "condition mapping produced:"
  while IFS="$(printf '\t')" read -r sid cond gsm rest how; do
    note "${sid}  ->  ${cond}   (${gsm}, ${how})"
  done < "${TMP}/map.tsv"
  cut -f1 "${TMP}/map.tsv" | grep -v '^[[:space:]]*$' > "${TMP}/sheet_ids.txt"
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
  curl -fL --retry 3 --retry-delay 3 --connect-timeout 20 --no-progress-meter -o "${TMP}/gene_info.gz" "${GENEINFO_URL}" \
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
ok()   { _emit "  [ OK ]   $*"; }
warn() { _emit "  [WARN]   $*"; }
bad()  { _emit "  [FAIL]   $*"; FAIL=$((FAIL+1)); }

for c in ${REQUIRED_ANNOT_COLS}; do
  [ -n "$(idx_of_col "${c}")" ] && ok "count matrix column present: ${c}" || bad "count matrix column missing: ${c}"
done

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
    n=$(awk -F, -v c="${cc}" 'NR>1 { v=$2; gsub(/"/,"",v); gsub(/\r/,"",v); if (v==c) k++ } END { print k+0 }' "${SHEET}"); note "n(${cc}) = ${n}"
    [ "${n}" -ge 2 ] || bad "condition '${cc}' has only ${n} sample(s); DESeq2 needs at least 2 per condition."
  done < "${TMP}/conds.txt"
  sort "${TMP}/sheet_ids.txt" > "${TMP}/a.txt"
  sort "${TMP}/sample_cols.txt" > "${TMP}/b.txt"
  if diff -q "${TMP}/a.txt" "${TMP}/b.txt" >/dev/null; then
    ok "sample ids in the samplesheet are an exact set match with the count matrix columns (${NSAMPLE_COLS} samples)"
  else
    bad "samplesheet sample ids do not match the count matrix columns exactly:"
    diff "${TMP}/a.txt" "${TMP}/b.txt" | sed 's/^/           /' | out || true
  fi
else
  bad "metadata/samplesheet.csv is missing"
fi

for s in "${REPO}"/scripts/*.R; do
  [ -e "${s}" ] || { warn "no .R files found in scripts/ - check the repository layout"; break; }
  if [ -s "${s}" ]; then ok "script present: scripts/$(basename "${s}")"; else bad "script is empty: scripts/$(basename "${s}")"; fi
done
[ -s "${REPO}/scripts/01_deseq2_analysis.R" ] && ok "driver present: scripts/01_deseq2_analysis.R" || warn "scripts/01_deseq2_analysis.R not found - confirm the pipeline entry point name"
[ -s "${REFS}" ] && ok "entrez reference present: data/reference/mouse_gene2entrez.tsv" || bad "entrez reference missing"
if [ -s "${REFS}" ]; then
  head -n 1 "${REFS}" | grep -q 'Symbol' && head -n 1 "${REFS}" | grep -q 'GeneID' \
    && ok "entrez reference has the Symbol / GeneID columns script 02 merges on" \
    || bad "entrez reference header must contain the columns Symbol and GeneID"
fi

DEG="${REPO}/results/deseq2_deg_significant.csv"
if [ -s "${DEG}" ] && [ "${NENS}" -gt 0 ]; then
  grep -oE 'ENSMUSG[0-9]+' "${DEG}" | sort -u > "${TMP}/deg_ids_all.txt" || true
  sed -n '1,200p' "${TMP}/deg_ids_all.txt" > "${TMP}/deg_ids.txt"
  nprobe="$(wc -l < "${TMP}/deg_ids.txt")"
  if [ "${nprobe}" -gt 0 ]; then
    gzip -dc "${COUNTS}" | cut -f"${GENE_ID_COL}" | grep -Fx -f "${TMP}/deg_ids.txt" | sort -u > "${TMP}/deg_found.txt" || true
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

snap > "${TMP}/results_after.md5"
if diff -q "${TMP}/results_before.md5" "${TMP}/results_after.md5" >/dev/null; then
  ok "results/ is byte-identical to the snapshot taken at the start of this run"
else
  bad "results/ changed during this run - investigate immediately:"
  diff "${TMP}/results_before.md5" "${TMP}/results_after.md5" | sed 's/^/           /' | out || true
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
git diff --cached --name-only | out || true

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
  log "Rscript on PATH: $(Rscript --version 2>&1 | sed -n '1p')"
else
  warn "Rscript is NOT on PATH in this shell."
fi
if command -v conda >/dev/null 2>&1; then
  log "conda detected; environments whose name matches /r_/ :"
  conda env list 2>/dev/null | awk '/(^r_|deseq)/{print "           " $0}' | out || true
  note "if the pipeline environment exists, run:  conda activate r_deseq_env"
else
  note "conda not on PATH in this shell."
fi
for p in DESeq2 clusterProfiler pheatmap ggplot2 org.Mm.eg.db; do
  if command -v Rscript >/dev/null 2>&1; then
    v="$(Rscript -e "cat(if(requireNamespace('${p}',quietly=TRUE)) as.character(packageVersion('${p}')) else 'NOT-INSTALLED')" 2>/dev/null || echo 'check-failed')"
    note "$(printf '%-16s %s' "${p}" "${v}")"
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

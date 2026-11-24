#!/usr/bin/env bash
set -euo pipefail

TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR=${PIPELINE_DIR:-$(cd "$TEST_DIR/.." && pwd)}
RUN="$PIPELINE_DIR/run_pipeline.sh"
SHARE=${CURRENT_SHARE_ROOT:-$PIPELINE_DIR/share_package}

fail() { echo "[FAIL] $*" >&2; exit 1; }
TMP=$(mktemp -d "$TEST_DIR/.tmp_download_routing.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

CURRENT_SHARE_ROOT="$SHARE" conda run -n jj_mfm bash "$RUN" current-validate-package > "$TMP/validate.log" 2>&1
rg -F 'active PacBio HiFi ready donors: 294/294' "$TMP/validate.log" >/dev/null || fail "current package validation output mismatch"

CURRENT_SHARE_ROOT="$SHARE" OUT_ROOT="$TMP/illumina" PLATFORM_FILTER=Illumina_30x INCLUDE_RNA=0 \
  conda run -n jj_mfm bash "$RUN" current-download-list 01 > "$TMP/illumina.log" 2>&1
[[ -s "$TMP/illumina/phase_01/all.urls" ]] || fail "current Illumina URL list missing"
[[ "$(wc -l < "$TMP/illumina/phase_01/all.urls")" == 1168 ]] || fail "current Illumina URL count mismatch"

CURRENT_SHARE_ROOT="$SHARE" OUT_ROOT="$TMP/mixed" PLATFORM_FILTER=all INCLUDE_RNA=0 \
  conda run -n jj_mfm bash "$RUN" current-download-list 01 > "$TMP/mixed.log" 2>&1
[[ -s "$TMP/mixed/phase_01/all.urls" ]] || fail "current mixed URL list missing"
[[ "$(wc -l < "$TMP/mixed/phase_01/all.urls")" == 1453 ]] || fail "current mixed URL count mismatch"
rg -F $'HG03125\tIllumina_30x\t1KGP_ILLUMINA' "$TMP/mixed/phase_01/samples.tsv" >/dev/null || fail "HG03125 Illumina routing missing"
rg -F $'HG03125\twgs_cram\thttps://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989126/HG03125.final.cram' "$TMP/mixed/phase_01/urls.with_kind.tsv" >/dev/null || fail "HG03125 CRAM URL missing"
if rg -F $'HG03125\twgs_pacbio' "$TMP/mixed/phase_01/urls.with_kind.tsv" >/dev/null; then fail "HG03125 CLR was routed in active package"; fi
rg -F $'__REFERENCE__\treference_fasta_gz' "$TMP/mixed/phase_01/urls.with_kind.tsv" >/dev/null || fail "exact reference routing missing"

LEGACY="$TMP/legacy"
mkdir -p "$LEGACY/scripts"
printf '%s\n' '#!/usr/bin/env bash' 'printf "%s\\n" "$@" > "${LEGACY_CAPTURE:?}"' > "$LEGACY/scripts/00_make_phase_assignments.sh"
chmod +x "$LEGACY/scripts/00_make_phase_assignments.sh"
LEGACY_CAPTURE="$TMP/legacy.args" DOWNLOAD_1000G_ROOT="$LEGACY" \
  conda run -n jj_mfm bash "$RUN" download-plan marker > "$TMP/legacy.log" 2>&1
rg -F 'historical 1KGP 3,202 downloader' "$TMP/legacy.log" >/dev/null || fail "legacy compatibility warning missing"
rg -Fx 'marker' "$TMP/legacy.args" >/dev/null || fail "legacy alias did not forward arguments"

for mode in current-validate-package current-download-list current-download current-verify-download current-download-vcf current-preprocess-phase; do
  rg -F "$mode)" "$RUN" >/dev/null || fail "missing current mode: $mode"
done

echo "[PASS] current/legacy download routing and frozen PacBio resource contract"

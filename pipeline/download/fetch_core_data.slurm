#!/bin/bash
#SBATCH --job-name=fetch-core-data
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=END,FAIL
#
# Pull a sequencing core's data release (HTTPS directory listing behind
# basic auth, e.g. UMN Genomics Center) onto HiPerGator, then verify it.
# Account/QoS/log path are passed at submit time so this file stays generic:
#
#   sbatch --account=<GROUP> --qos=<GROUP> \
#          --output=<DEST>/logs/fetch-%j.out \
#          pipeline/download/fetch_core_data.slurm <URL> <DEST> [CUT_DIRS]
#
# Credentials are never on the command line or in this repo. Before
# submitting, create a private wget config holding them (you type the
# password; it is not echoed):
#
#   umask 077
#   read -rp "core username: " U; read -rsp "core password: " P; echo
#   printf 'user = %s\npassword = %s\n' "$U" "$P" > ~/.core_download.wgetrc
#   unset U P; chmod 600 ~/.core_download.wgetrc
#
# Delete it once the download is verified: rm ~/.core_download.wgetrc
#
# Re-running is safe: -N skips files already complete and re-fetches any
# that are partial or changed.

set -euo pipefail

URL="${1:?usage: fetch_core_data.slurm <URL> <DEST> [CUT_DIRS]}"
DEST="${2:?usage: fetch_core_data.slurm <URL> <DEST> [CUT_DIRS]}"
CUT_DIRS="${3:-2}"   # leading URL path components to drop (core's own wget line says 2)
CREDS="${CORE_WGETRC:-${HOME}/.core_download.wgetrc}"

if [[ ! -f "${CREDS}" ]]; then
  echo "ERROR: ${CREDS} not found -- create it first (see header)." >&2
  exit 1
fi
if [[ "$(stat -c %a "${CREDS}")" != "600" ]]; then
  echo "ERROR: ${CREDS} must be chmod 600." >&2
  exit 1
fi

mkdir -p "${DEST}/logs"
cd "${DEST}"

echo "[$(date)] host=$(hostname) url=${URL} dest=${DEST}"

# --no-check-certificate matches the core's own documented command.
# -R index.html* drops the auto-generated directory listing pages.
wget --config="${CREDS}" \
  -nH -np -N -r --cut-dirs "${CUT_DIRS}" \
  --no-check-certificate \
  -R "index.html*" \
  --tries=20 --waitretry=30 --read-timeout=120 \
  --progress=dot:giga \
  "${URL}"

echo "[$(date)] download finished; verifying"

# Find the top-level folder wget created (last component of the URL path).
TOP="$(basename "${URL%/}")"
cd "${DEST}/${TOP}"

# 1) Manifest of everything received, with sizes.
find . -type f ! -name 'MANIFEST*' -printf '%s\t%p\n' | sort -k2 > MANIFEST.sizes.tsv
echo "files: $(wc -l < MANIFEST.sizes.tsv)  total: $(du -sh . | cut -f1)"

# 2) If the core shipped checksum files, verify against them.
mapfile -t SUMS < <(find . -type f \( -iname '*md5*' -o -iname '*.sha256' \) ! -name 'MANIFEST*')
if (( ${#SUMS[@]} )); then
  for s in "${SUMS[@]}"; do
    echo "checking ${s}"
    tool=md5sum; [[ "${s,,}" == *sha256* ]] && tool=sha256sum
    ( cd "$(dirname "${s}")" && "${tool}" -c --quiet "$(basename "${s}")" )
  done
else
  echo "no core checksum files found -- skipping md5 verification"
fi

# 3) gzip integrity of every FASTQ (catches truncated transfers).
find . -type f -name '*.fastq.gz' -print0 \
  | xargs -0 -n1 -P "${SLURM_CPUS_PER_TASK:-4}" sh -c 'gzip -t "$0" || echo "CORRUPT: $0"' \
  | tee gzip_check.txt
if grep -q CORRUPT gzip_check.txt; then
  echo "ERROR: corrupt FASTQ(s) listed above -- re-run this job to re-fetch." >&2
  exit 1
fi

# 4) Our own checksums, for later integrity checks after moves/copies.
find . -type f -name '*.fastq.gz' -print0 | sort -z \
  | xargs -0 -n1 -P "${SLURM_CPUS_PER_TASK:-4}" md5sum > MANIFEST.md5

echo "[$(date)] all checks passed"

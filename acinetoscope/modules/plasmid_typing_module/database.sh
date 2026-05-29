#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
TEMP_DIR=$(mktemp -d)

echo "Creating data directory if not exists: ${DATA_DIR}"
mkdir -p "${DATA_DIR}"

echo "Cloning AcinetobacterPlasmidTyping repository into temporary directory..."
git clone --depth 1 https://github.com/MehradHamidian/AcinetobacterPlasmidTyping.git "${TEMP_DIR}"

echo "Listing files in cloned repository for debugging:"
ls -la "${TEMP_DIR}"

# Try to find the DNA rep sequences file using a more flexible pattern
DB_FILE=$(find "${TEMP_DIR}" -maxdepth 1 -type f -name "*rep_DNA-seqs*.fasta" | head -1)
if [ -z "${DB_FILE}" ]; then
    echo "Error: Could not find rep DNA sequences file."
    echo "Available files:"
    ls "${TEMP_DIR}"
    rm -rf "${TEMP_DIR}"
    exit 1
fi

echo "Copying ${DB_FILE} to ${DATA_DIR}/apt_reps_v3.fasta"
cp "${DB_FILE}" "${DATA_DIR}/apt_reps_v3.fasta"

echo "Cleaning up temporary directory..."
rm -rf "${TEMP_DIR}"

echo "Database downloaded successfully: ${DATA_DIR}/apt_reps_v3.fasta"

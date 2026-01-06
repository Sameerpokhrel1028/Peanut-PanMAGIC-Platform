#!/usr/bin/env bash
#SBATCH --job-name=hifi_assembly
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=1000G
#SBATCH --time=190:00:00
#SBATCH --output=logs/hifi_assembly.%j.out
#SBATCH --error=logs/hifi_assembly.%j.err


#Minimal script for a single file
#hifiadapterfilt.sh -t 10 -p filename filename.hifi.fastq.gz
#hifiasm -t 10 -l0 -o filename.hifiasm filename_HIFI.filt.fastq.gz # lo flag represents homozygosity

set -euo pipefail

# -----------------------
# Requirements (examples)
# -----------------------
# Use module load if your cluster supports modules, otherwise ensure these are on PATH.
# module load HiFiAdapterFilt ( https://github.com/sheinasim-USDA/HiFiAdapterFilt)
# module load hifiasm (https://github.com/chhylp123/hifiasm)

# -----------------------
# User settings
# -----------------------
THREADS="${SLURM_CPUS_PER_TASK:-30}"

# Update these paths for your system (or load from a config file)
INDIR="${INDIR:-/path/to/hifi_fastq}"
WORKDIR="${WORKDIR:-/path/to/workdir/Genome_Assembly}"
ADAPTER_OUT="${ADAPTER_OUT:-/path/to/workdir/Adapter_Filter}"

mkdir -p "$WORKDIR" "$ADAPTER_OUT" logs

# -----------------------
# Choosing in batches due to higher number
# -----------------------
# Batch 1
FILES=(
  "IAC322_HIFI.fastq.gz"
  "GP-NC-WS17_HIFI.fastq.gz"
  "GeorgiaGreen_HIFI.fastq.gz"
  "GA12Y_HIFI.fastq.gz"
  "CC812_HIFI.fastq.gz"
  "C76_16_HIFI.fastq.gz"
)

# Batch 2
# FILES=(
#   "TifNV_HIFI.fastq.gz"
#   "ICG1471_HIFI.fastq.gz"
#   "York_HIFI.fastq.gz"
#   "Bailey_HIFI.fastq.gz"
#   "C99R_HIFI.fastq.gz"
#   "Georganic_HIFI.fastq.gz"
# )

# Batch 3
# FILES=(
#   "Lariat_HIFI.fastq.gz"
#   "Marc1_HIFI.fastq.gz"
#   "C431_HIFI.fastq.gz"
#   "Florida07_HIFI.fastq.gz"
#   "NC94022.ccs.fastq.gz"
# )

echo "Input dir:   $INDIR"
echo "Work dir:    $WORKDIR"
echo "Adapter out: $ADAPTER_OUT"
echo "Threads:     $THREADS"
echo

# -----------------------
# Main loop
# -----------------------
cd "$INDIR"

for FILE in "${FILES[@]}"; do
  if [[ ! -f "$FILE" ]]; then
    echo "ERROR: Missing input file: $INDIR/$FILE" >&2
    exit 1
  fi

  # Sample naming: supports *.fastq.gz and *.ccs.fastq.gz
  BASENAME="$(basename "$FILE")"
  SAMPLE="${BASENAME%.fastq.gz}"
  SAMPLE="${SAMPLE%.ccs}"

  echo "=== ${SAMPLE} ==="

  # Copy to workdir (keeps original inputs untouched)
  cp -f "$FILE" "$WORKDIR/"
  cd "$WORKDIR"

  echo "[1/3] Adapter filtering: ${SAMPLE}"
  # HiFiAdapterFilt expects files in the current directory for the given prefix
  # If your naming differs, adjust the -p prefix and/or input file naming convention.
  bash hifiadapterfilt.sh -t "$THREADS" -p "$SAMPLE" -o "$ADAPTER_OUT"

  # Remove the copied raw fastq from workdir to save space (original remains in $INDIR)
  rm -f "$WORKDIR/$BASENAME"

  echo "[2/3] Assembly (hifiasm): ${SAMPLE}"
  hifiasm -o "${SAMPLE}.hifiasm" -t "$THREADS" -l0 \
    "$ADAPTER_OUT/${SAMPLE}_HIFI.filt.fastq.gz"

  echo "[3/3] Convert primary contigs GFA -> FASTA: ${SAMPLE}"
  awk '/^S/{print ">"$2"\n"$3}' "${SAMPLE}.hifiasm.bp.p_ctg.gfa" > "${SAMPLE}.hifiasm.bp.p_ctg.fa"

  echo "Done: ${SAMPLE}"
  echo

  cd "$INDIR"
done

echo "All assemblies finished."

#!/bin/bash
#SBATCH --job-name=rdeval_array          # Job name
#SBATCH --output=rdeval_%A_%a.out       # One log per array task
#SBATCH --error=rdeval_%A_%a.err
#SBATCH --partition=vgl_a               # Partition/queue to use
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --acctg-freq=task=30s

# Usage: sbatch --array=0-(N-1) rdeval_array.sh filelist.txt

if [ $# -ne 1 ]; then
    echo "Usage: $0 <filelist.txt>"
    exit 1
fi

FILELIST="$1"

# Get the file corresponding to this array index
FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$FILELIST")

if [ -z "$FILE" ]; then
    echo "[ERROR] No file found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Script for rdeval_v.0.0.8 (array mode)"
echo "[DEBUG] starting at $(date)"
echo "[DEBUG] SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "[DEBUG] processing file: $FILE"


# Output dir
OUTDIR="results_v0.0.8"
mkdir -p "$OUTDIR"

# Summary TSV (shared by all array tasks)
SUMMARY_TSV="${OUTDIR}/rdeval_summary.tsv"

# Only task 0 creates/overwrites the header
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    echo -e "file\tTotal_read_length\tMaxRSS_kB\tElapsed_wallclock" > "$SUMMARY_TSV"
else
    # Wait until header file exists (in case task 0 is a bit slower to start)
    while [ ! -f "$SUMMARY_TSV" ]; do
        sleep 2
    done
fi

start_time=$(date +%s)

time_log="${OUTDIR}/$(basename "${FILE}").time.log"
result_file="${OUTDIR}/$(basename "${FILE}").results.rdeval"

# Run rdeval with /usr/bin/time -v
/usr/bin/time -v -o "$time_log" \
    ./Software/rdeval_0.0.8/rdeval_static_binary \
    "$FILE" --cmd -j 4 > "$result_file"

end_time=$(date +%s)
elapsed_seconds=$((end_time - start_time))

echo "[DEBUG] $FILE finished in ${elapsed_seconds}s"

# Extract Total read length from rdeval result
total_len=$(awk -F': ' '/^Total read length/ {print $2}' "$result_file")

# Extract MaxRSS and wall clock time from /usr/bin/time -v
max_rss=$(awk -F': +' '/Maximum resident set size \(kbytes\)/ {print $2}' "$time_log")
wallclock=$(awk -F': ' '/Elapsed \(wall clock\) time/ {print $2}' "$time_log")

# Append a row to the summary TSV
echo -e "${FILE}\t${total_len}\t${max_rss}\t${wallclock}" >> "$SUMMARY_TSV"

echo "[DEBUG] ending at $(date)"
echo "[DEBUG] Summary file is: $SUMMARY_TSV"

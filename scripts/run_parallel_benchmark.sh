#!/bin/bash
set -euo pipefail

DATASETS=("berlin52" "d198" "pr439" "pr1002")
PROCS=(2 4 8)
GEN=200
SEED=42

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BIN="${ROOT_DIR}/tsp_ga"
DATA_DIR="${ROOT_DIR}/data"
OUTDIR="${ROOT_DIR}/outputs/parallel_runs"
LOGFILE="${OUTDIR}/runtime_log.csv"

if [[ ! -x "${BIN}" ]]; then
  echo "❌ Cannot find executable at ${BIN}. Build the project first (make clean && make)." >&2
  exit 1
fi

if ! command -v mpirun >/dev/null 2>&1; then
  echo "❌ mpirun not found on PATH. Install OpenMPI/MPICH and ensure mpirun is available." >&2
  exit 1
fi

mkdir -p "${OUTDIR}"
echo "timestamp,dataset,np,best_length,real_time_sec" > "${LOGFILE}"

for dataset in "${DATASETS[@]}"; do
  dataset_path="${DATA_DIR}/${dataset}.tsp"
  if [[ ! -f "${dataset_path}" ]]; then
    echo "⚠️ Skipping ${dataset}: missing file ${dataset_path}" >&2
    continue
  fi

  for np in "${PROCS[@]}"; do
    outfile="${OUTDIR}/${dataset}_p${np}.txt"
    echo "🚀 Running: ${dataset} with -np ${np}"

    timestamp="$(date +"%Y-%m-%d %H:%M:%S")"

    {
      echo "----- RUN START: ${timestamp} -----"
      # --oversubscribe avoids failures on laptops with fewer hardware slots than requested
      /usr/bin/time -p mpirun --oversubscribe -np "${np}" "${BIN}" "${dataset_path}" --generations "${GEN}" --seed "${SEED}"
      echo "----- RUN END: $(date +"%Y-%m-%d %H:%M:%S") -----"
    } &> "${outfile}"

    best_length="$(awk -F': ' '/Best tour length/ {val=$2} END {gsub(/^[[:space:]]+|[[:space:]]+$/, "", val); if (val=="") print "NA"; else print val}' "${outfile}")"
    real_time="$(awk '/^real / {val=$2} END {if (val=="") print "NA"; else print val}' "${outfile}")"

    echo "${timestamp},${dataset},${np},${best_length},${real_time}" >> "${LOGFILE}"
    echo "✔ Done ${dataset} -np ${np}: ${best_length} in ${real_time}s"
  done
done

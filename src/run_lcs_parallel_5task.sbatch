#!/bin/bash
#SBATCH -p cpu
#SBATCH --job-name=LCS_MPI_1ppn
#SBATCH --nodes=2
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=5
#SBATCH --mem=4G
#SBATCH --time=11:00:00
#SBATCH --output=log_mpi_%j.out

cd $SLURM_SUBMIT_DIR

echo "=== Compilando lcs_mpi ==="
make paralelo || { echo "Falha na compilação"; exit 1; }

INPUT_DIR="testes/inputs"
OUTPUT_DIR="testes/outputs_5tasks"
mkdir -p "${OUTPUT_DIR}"
rm -f "${OUTPUT_DIR}"/output_*.out

for fileA in "${INPUT_DIR}"/*A.in; do
  for fileB in "${INPUT_DIR}"/*B.in; do
    baseA=$(basename "${fileA}" .in)
    baseB=$(basename "${fileB}" .in)

    output_file="${OUTPUT_DIR}/output_${baseA}_${baseB}.out"

    echo "[$(date +'%H:%M:%S')] ${baseA} vs ${baseB}" | tee -a "$output_file"
    for i in $(seq 1 21); do
      echo "--- Execução $i ---" >> "$output_file"
      mpirun -np 10 --map-by ppr:5:node ./lcs_mpi "$fileA" "$fileB" \
        >> "$output_file" 2>&1
    done
  done
done

echo "=== Todos os outputs estão em ${OUTPUT_DIR} ==="


#!/bin/bash
#SBATCH -p cpu
#SBATCH --job-name=LCS_MPI_1ppn
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --mem=4G
#SBATCH --time=11:00:00
#SBATCH --output=log_mpi_%j.out

cd $SLURM_SUBMIT_DIR

echo "=== Compilando lcs_mpi ==="
make paralelo || { echo "Falha na compilação"; exit 1; }

INPUT_DIR="testes/inputs_resumidos"
OUTPUT_DIR="testes/outputs_4tasks"
mkdir -p "${OUTPUT_DIR}"
rm -f "${OUTPUT_DIR}"/output_*.out

for fileA in "${INPUT_DIR}"/*A.in; do
  for fileB in "${INPUT_DIR}"/*B.in; do
    baseA=$(basename "${fileA}" .in)
    baseB=$(basename "${fileB}" .in)

    # Pula se estivermos revertendo o par (B precede A lexicograficamente)
    if [[ "${baseB}" < "${baseA}" ]]; then
      continue
    fi

    output_file="${OUTPUT_DIR}/output_${baseA}_${baseB}.out"

    echo "[$(date +'%H:%M:%S')] ${baseA} vs ${baseB}" | tee -a "$output_file"
    for i in $(seq 1 21); do
      echo "--- Execução $i ---" >> "$output_file"
      mpirun -np 8 --map-by ppr:4:node ./lcs_mpi "$fileA" "$fileB" \
        >> "$output_file" 2>&1
    done
  done
done

echo "=== Todos os outputs estão em ${OUTPUT_DIR} ==="


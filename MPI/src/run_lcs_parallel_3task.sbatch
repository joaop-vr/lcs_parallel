#!/bin/bash
#SBATCH -p cpu
#SBATCH --job-name=LCS_MPI
#SBATCH --nodes=2
#SBATCH --ntasks=6
#SBATCH --ntasks-per-node=3
#SBATCH --mem=4G
#SBATCH --cpu-freq=3.0-3.0:performance
#SBATCH --time=12:00:00
#SBATCH --output=log_mpi_%j.out

cd $SLURM_SUBMIT_DIR

# 1) Compilar
echo "=== Compilando lcs_mpi ==="
make paralelo || { echo "Falha na compilação"; exit 1; }

# 2) Diretórios
INPUT_DIR="testes/inputs_resumidos"
OUTPUT_DIR="testes/outputs_3tasks"
mkdir -p "${OUTPUT_DIR}"
rm -f "${OUTPUT_DIR}"/output_*.out

# 3) Loop sobre pares A vs B
for fileA in "${INPUT_DIR}"/*A.in; do
  for fileB in "${INPUT_DIR}"/*B.in; do
    baseA=$(basename "${fileA}" .in)
    baseB=$(basename "${fileB}" .in)
    output_file="${OUTPUT_DIR}/output_${baseA}_${baseB}.out"

    echo "[$(date +'%H:%M:%S')] ${baseA} vs ${baseB}" | tee -a "$output_file"
    for i in {1..21}; do
      echo "--- Execução $i ---" >> "$output_file"
      # Corrigido: Usar informação do SLURM e remover mapeamento manual
      mpirun -np 6 --map-by ppr:3:node ./lcs_mpi "$fileA" "$fileB" >> "$output_file" 2>&1
    done
  done
done

echo "=== Todos os outputs estão em ${OUTPUT_DIR} ==="


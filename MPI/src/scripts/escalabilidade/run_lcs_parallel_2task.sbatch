#!/bin/bash
#SBATCH -p cpu
#SBATCH --job-name=LCS_MPI    
#SBATCH --nodes=2                  # 2 nós
#SBATCH --ntasks=4                 # 4 processos MPI no total
#SBATCH --ntasks-per-node=2
#SBATCH --mem=4G
#SBATCH --cpu-freq=3500000-3500000:performance
#SBATCH --time=07:00:00
#SBATCH --output=log_mpi_%j.out

cd $SLURM_SUBMIT_DIR

# 1) Compilar
echo "=== Compilando lcs_mpi ==="
make paralelo || { echo "Falha na compilação"; exit 1; }

# 2) Diretórios
INPUT_DIR="testes/inputs_resumidos"
OUTPUT_DIR="testes/outputs_2tasks"
mkdir -p "${OUTPUT_DIR}"
rm -f "${OUTPUT_DIR}"/output_*.out

# 3) Loop sobre pares A vs B (sem filtrar aqui, roda todos)
for fileA in "${INPUT_DIR}"/*A.in; do
  for fileB in "${INPUT_DIR}"/*B.in; do
    baseA=$(basename "${fileA}" .in)
    baseB=$(basename "${fileB}" .in)
    output_file="${OUTPUT_DIR}/output_${baseA}_${baseB}.out"

    echo "[$(date +'%H:%M:%S')] ${baseA} vs ${baseB}" | tee -a "$output_file"
    for i in $(seq 1 21); do
      echo "--- Execução $i ---" >> "$output_file"
      # 1 por nó: --map-by ppr:1:node garante 1 rank em cada nó
      mpirun -np 4 --map-by ppr:2:node ./lcs_mpi "$fileA" "$fileB" \
        >> "$output_file" 2>&1
    done
  done
done

echo "=== Todos os outputs estão em ${OUTPUT_DIR} ==="


#!/bin/bash
#SBATCH -p cpu
#SBATCH --job-name=LCS_MPI_1ppn
#SBATCH --nodes=2
#SBATCH --ntasks=6
#SBATCH --ntasks-per-node=3
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=log_mpi_%j.out

cd $SLURM_SUBMIT_DIR

# 1) Compilar
echo "=== Compilando lcs_mpi ==="
make paralelo || { echo "Falha na compilação"; exit 1; }

# 2) Diretórios
INPUT_DIR="testes/inputs"
OUTPUT_DIR="testes/outputs_3tasks"
mkdir -p "${OUTPUT_DIR}"
rm -f "${OUTPUT_DIR}"/output_*.out

# 3) Lista restrita de arquivos A
for fileA in "${INPUT_DIR}/teste_40000A.in" "${INPUT_DIR}/teste_50000A.in" "${INPUT_DIR}/teste_60000A.in" "${INPUT_DIR}/teste_70000A.in"; do
  baseA=$(basename "${fileA}" .in)

  # Se for teste_40000A, apenas comparar com teste_70000B
  if [[ "$baseA" == "teste_40000A" ]]; then
    fileB="${INPUT_DIR}/teste_70000B.in"
    baseB=$(basename "${fileB}" .in)
    output_file="${OUTPUT_DIR}/output_${baseA}_${baseB}.out"

    echo "[$(date +'%H:%M:%S')] ${baseA} vs ${baseB}" | tee -a "$output_file"
    for i in {1..21}; do
      echo "--- Execução $i ---" >> "$output_file"
      mpirun -np 6 --map-by ppr:3:node ./lcs_mpi "$fileA" "$fileB" >> "$output_file" 2>&1
    done
  else
    # Caso contrário, comparar com todos os arquivos B
    for fileB in "${INPUT_DIR}"/*B.in; do
      baseB=$(basename "${fileB}" .in)
      output_file="${OUTPUT_DIR}/output_${baseA}_${baseB}.out"

      echo "[$(date +'%H:%M:%S')] ${baseA} vs ${baseB}" | tee -a "$output_file"
      for i in {1..21}; do
        echo "--- Execução $i ---" >> "$output_file"
        mpirun -np 6 --map-by ppr:3:node ./lcs_mpi "$fileA" "$fileB" >> "$output_file" 2>&1
      done
    done
  fi
done

echo "=== Todos os outputs estão em ${OUTPUT_DIR} ==="


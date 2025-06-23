#!/bin/bash
#SBATCH -p cpu
#SBATCH --job-name=LCS_MPI
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --mem=4G
#SBATCH --cpu-freq=3500000-3500000:performance
#SBATCH --time=12:00:00
#SBATCH --output=log_mpi_%j.out

cd $SLURM_SUBMIT_DIR

# 1) Compilar
echo "=== Compilando lcs_mpi ==="
make paralelo || { echo "Falha na compilação"; exit 1; }

# 2) Diretórios
INPUT_DIR="testes/inputs_resumidos"
OUTPUT_DIR="testes/outputs_4tasks"
mkdir -p "${OUTPUT_DIR}"
rm -f "${OUTPUT_DIR}"/output_*.out

# 3) Loop apenas sobre pares com mesmo nome base (XA vs XB)
for fileA in "${INPUT_DIR}"/*A.in; do
  # Extrai a base do nome (ex: teste_1000)
  base=$(basename "${fileA}" A.in)
  fileB="${INPUT_DIR}/${base}B.in"

  # Verifica se o arquivo B correspondente existe
  if [[ -f "$fileB" ]]; then
    output_file="${OUTPUT_DIR}/output_${base}A_${base}B.out"

    echo "[$(date +'%H:%M:%S')] ${base}A vs ${base}B" | tee -a "$output_file"
    for i in {1..21}; do
      echo "--- Execução $i ---" >> "$output_file"
      # Corrigido: Usar informação do SLURM e remover mapeamento manual
      mpirun -np 8 --map-by ppr:4:node ./lcs_mpi "$fileA" "$fileB" >> "$output_file" 2>&1
    done
  else
    echo "Arquivo correspondente não encontrado para ${fileA}"
  fi
done

echo "=== Todos os outputs estão em ${OUTPUT_DIR} ==="


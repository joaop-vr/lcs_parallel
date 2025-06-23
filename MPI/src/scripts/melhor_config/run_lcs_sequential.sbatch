#!/bin/bash
#SBATCH -p cpu                 # partição
#SBATCH --job-name=LCS_SEQ     # nome do job
#SBATCH --nodes=1              # 1 nó
#SBATCH --cpus-per-task=1      # 1 CPU por tarefa
#SBATCH --mem=4G               # memória
#SBATCH --time=07:00:00        # tempo máximo
#SBATCH --output=log_seq_%j.out  # arquivo de log (%j = job ID)

cd $SLURM_SUBMIT_DIR

# 1) Compila o executável sequencial
echo "=== Compilando lcs_sequencial ==="
make sequencial || { echo "Falha na compilação"; exit 1; }

# 2) Prepara diretório de outputs
INPUT_DIR="testes/melhor_config/inputs_resumidos"
OUTPUT_DIR="testes/outputs"

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
    for i in $(seq 1 21); do
      echo "--- Execução $i ---" >> "$output_file"
      ./lcs_sequencial "$fileA" "$fileB" >> "$output_file"
    done
  else
    echo "Arquivo correspondente não encontrado para ${fileA}"
  fi
done


echo "=== Todos os outputs estão em ${OUTPUT_DIR} ==="

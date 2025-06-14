#!/bin/bash
#SBATCH -p cpu                 # partição
#SBATCH --job-name=LCS_SEQ     # nome do job
#SBATCH --nodes=1              # 1 nó
#SBATCH --ntasks=1             # 1 tarefa (processo)
#SBATCH --cpus-per-task=1      # 1 CPU por tarefa
#SBATCH --mem=4G               # memória
#SBATCH --time=03:00:00        # tempo máximo
#SBATCH --output=log_seq_%j.out  # arquivo de log (%j = job ID)

cd $SLURM_SUBMIT_DIR

# 1) Compila o executável sequencial
echo "=== Compilando lcs_sequencial ==="
make sequencial || { echo "Falha na compilação"; exit 1; }

# 2) Prepara diretório de outputs
INPUT_DIR="testes/inputs"
OUTPUT_DIR="testes/outputs"

mkdir -p "${OUTPUT_DIR}"
rm -f "${OUTPUT_DIR}"/output_*.out

# 3) Executa apenas os pares equivalentes: teste_XXXXXA.in vs teste_XXXXXB.in
for fileA in "${INPUT_DIR}"/*A.in; do
  for fileB in "${INPUT_DIR}"/*B.in; do
    baseA=$(basename "${fileA}" .in)
    baseB=$(basename "${fileB}" .in)

    # extrai apenas o número do nome (ex: 10000)
    numA="${baseA//[!0-9]/}"
    numB="${baseB//[!0-9]/}"

    # só executa se o número for o mesmo
    if [[ "$numA" != "$numB" ]]; then
      continue
    fi

    output_file="${OUTPUT_DIR}/output_${baseA}_${baseB}.out"
    echo "[$(date +'%H:%M:%S')] ${baseA} vs ${baseB}"
    for i in $(seq 1 21); do
      echo "--- Execução $i ---" >> "$output_file"
      ./lcs_sequencial "$fileA" "$fileB" >> "$output_file"
    done
  done
done

echo "=== Todos os outputs estão em ${OUTPUT_DIR} ==="


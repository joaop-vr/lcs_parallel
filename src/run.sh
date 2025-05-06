#!/bin/bash

# Diretório para armazenar os resultados (cria se não existir)
OUTPUT_DIR="resultados"
mkdir -p "$OUTPUT_DIR"

# Parâmetros a serem testados
INPUT_SIZES=(1000 5000 10000 50000 100000)  # Tamanhos das entradas
THREADS=(1 2 4 8)                           # Número de thre
schedule="static" 
CHUNKS=(64 128)    # Tamanhos de chunk

# Pergunta ao usuário se deseja continuar
read -p "Confirmar esse schedule: ${schedule}? (s/n) " confirmacao
if [[ "$confirmacao" != "s" && "$confirmacao" != "S" ]]; then
    echo "Execução cancelada pelo usuário."
    exit 1
fi

# Loop principal
for size in "${INPUT_SIZES[@]}"; do
    
    ./generate_input "$size" "$size" "teste_$size"
    
    for threads in "${THREADS[@]}"; do
        for chunk in "${CHUNKS[@]}"; do
            echo "[TESTE] size=$size | threads=$threads | schedule=${schedule} | chunk=$chunk" 
            
            # Nome do arquivo de saída
            OUTPUT_FILE="${OUTPUT_DIR}/size_${size}_threads_${threads}_sched_${schedule}_chunk_${chunk}.txt" 
            
            for run in {1..21}; do
                echo "--- Execução $run ---" >> "$OUTPUT_FILE"
                ./lcs "teste_${size}A.in" "teste_${size}B.in" "$threads" "$chunk" >> "$OUTPUT_FILE"
                echo "" >> "$OUTPUT_FILE"
            done
        done
    done
    
    #rm teste_${size}*.in
done

echo "Todos os testes concluídos! Resultados salvos em: $OUTPUT_DIR/"

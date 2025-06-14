#!/bin/bash

# Diretório para armazenar os resultados (cria se não existir)
OUTPUT_DIR="resultados"
mkdir -p "$OUTPUT_DIR"

# Parâmetros a serem testados
INPUT_SIZES=(1000 5000 10000 50000 100000)   # Tamanhos das entradas
THREADS=(1 2 4 8)                             # Número de threads
SCHEDULES=("static" "dynamic" "guided")      # Tipos de escalonamento
CHUNKS=(64 128)                              # Tamanhos de chunk

# Loop principal
for size in "${INPUT_SIZES[@]}"; do
    
    ./generate_input "$size" "$size" "teste_$size"
    
    for threads in "${THREADS[@]}"; do
        for chunk in "${CHUNKS[@]}"; do
            for schedule in "${SCHEDULES[@]}"; do
                echo "[TESTE] size=$size | threads=$threads | schedule=${schedule} | chunk=$chunk" 
                
                # Nome do arquivo de saída
                OUTPUT_FILE="${OUTPUT_DIR}/size_${size}_threads_${threads}_sched_${schedule}_chunk_${chunk}.txt" 
                
                for run in {1..21}; do
                    echo "--- Execução $run ---" >> "$OUTPUT_FILE"
                    ./lcs "teste_${size}A.in" "teste_${size}B.in" "$threads" "$chunk" "$schedule" >> "$OUTPUT_FILE"
                    echo "" >> "$OUTPUT_FILE"
                done
            done
        done
    done

    #rm teste_${size}*.in
done

echo "Todos os testes concluídos! Resultados salvos em: $OUTPUT_DIR/"

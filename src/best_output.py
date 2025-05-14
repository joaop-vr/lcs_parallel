import os
import re
import sys

def extrair_hiperparametros(nome_arquivo):
    # Remove prefixo e separa os hiperparâmetros
    nome_limpo = nome_arquivo.replace("results_", "")
    partes = nome_limpo.split("_")
    pares = [f"{partes[i]}={partes[i+1]}" for i in range(0, len(partes), 2)]
    return ", ".join(pares)

def extrair_tempos_paralelos(conteudo):
    match_init = re.search(r"Média Tempo init matriz par:\s+([\d.]+)", conteudo)
    match_lcs  = re.search(r"Média Tempo LCS par:\s+([\d.]+)", conteudo)
    if match_init and match_lcs:
        return float(match_init.group(1)), float(match_lcs.group(1))
    return None, None

def encontrar_melhor_configuracao(diretorio):
    melhor_tempo_total = float('inf')
    melhor_config = ""
    melhor_arquivo = ""

    for nome_arquivo in os.listdir(diretorio):
        if not nome_arquivo.startswith("results_"):
            continue

        caminho_arquivo = os.path.join(diretorio, nome_arquivo)
        if not os.path.isfile(caminho_arquivo):
            continue

        with open(caminho_arquivo, 'r') as f:
            conteudo = f.read()

        tempo_init, tempo_lcs = extrair_tempos_paralelos(conteudo)
        if tempo_init is None or tempo_lcs is None:
            continue

        tempo_total = tempo_init + tempo_lcs

        if tempo_total < melhor_tempo_total:
            melhor_tempo_total = tempo_total
            melhor_config = extrair_hiperparametros(nome_arquivo)
            melhor_arquivo = nome_arquivo

    if melhor_config:
        print("Melhor configuração encontrada:")
        print(f"Arquivo: {melhor_arquivo}")
        print(f"Hiperparâmetros: {melhor_config}")
        print(f"Tempo total paralelo (init + LCS): {melhor_tempo_total:.10f} segundos")
    else:
        print("Nenhum resultado paralelo válido encontrado.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python melhor_config.py <caminho_para_diretorio>")
    else:
        encontrar_melhor_configuracao(sys.argv[1])

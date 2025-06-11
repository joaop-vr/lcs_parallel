import os
import re
import sys
from collections import defaultdict
from statistics import mean, variance

def extrair_chave_config(nome_arquivo):
    """
    Extrai uma chave no formato: threads=8, sched=guided, chunk=128
    Ignora o 'size'
    """
    partes = nome_arquivo.replace("results_", "").split("_")
    config = {}
    for i in range(0, len(partes), 2):
        if partes[i] == "size":
            continue
        config[partes[i]] = partes[i+1]
    return f"threads={config.get('threads')}, sched={config.get('sched')}, chunk={config.get('chunk')}"

def extrair_tempos_paralelos(conteudo):
    match_init = re.search(r"Média Tempo init matriz par:\s+([\d.]+)", conteudo)
    match_lcs  = re.search(r"Média Tempo LCS par:\s+([\d.]+)", conteudo)
    if match_init and match_lcs:
        return float(match_init.group(1)) + float(match_lcs.group(1))
    return None

def calcular_estatisticas(diretorio):
    grupos = defaultdict(list)

    for nome_arquivo in os.listdir(diretorio):
        if not nome_arquivo.startswith("results_"):
            continue

        caminho_arquivo = os.path.join(diretorio, nome_arquivo)
        if not os.path.isfile(caminho_arquivo):
            continue

        with open(caminho_arquivo, 'r') as f:
            conteudo = f.read()

        tempo_total = extrair_tempos_paralelos(conteudo)
        if tempo_total is None:
            continue

        chave = extrair_chave_config(nome_arquivo)
        grupos[chave].append(tempo_total)

    if not grupos:
        print("Nenhum dado paralelo válido encontrado.")
        return

    print("Estatísticas por configuração (ignorando tamanho de entrada):\n")
    melhor_config = None
    melhor_media = float('inf')

    for config, tempos in grupos.items():
        media = mean(tempos)
        var = variance(tempos) if len(tempos) > 1 else 0.0
        print(f"Configuração: {config}")
        print(f"  Média dos tempos totais:     {media:.10f} segundos")
        print(f"  Variância dos tempos totais: {var:.10f}\n")

        if media < melhor_media:
            melhor_media = media
            melhor_config = config

    print("Melhor configuração média encontrada:")
    print(f"  Configuração: {melhor_config}")
    print(f"  Média total paralela: {melhor_media:.10f} segundos")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python agrupar_resultados.py <caminho_para_diretorio>")
    else:
        calcular_estatisticas(sys.argv[1])

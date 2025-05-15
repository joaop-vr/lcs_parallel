import re
import os
import sys

def calcular_media_variancia(valores):
    if not valores:
        return 0.0, 0.0
    media = sum(valores) / len(valores)
    variancia = sum((x - media) ** 2 for x in valores) / len(valores)
    return media, variancia

def calcular_tempos_medios_em_arquivos(diretorio):
    if not os.path.isdir(diretorio):
        print(f"Erro: '{diretorio}' não é um diretório válido.")
        return

    # Cria o diretório de saída, se não existir
    pasta_resultados = os.path.join(diretorio, "results_interpreted")
    os.makedirs(pasta_resultados, exist_ok=True)

    for nome_arquivo in os.listdir(diretorio):
        if nome_arquivo.startswith("results_"):
            continue

        caminho_arquivo = os.path.join(diretorio, nome_arquivo)
        if not os.path.isfile(caminho_arquivo):
            continue

        with open(caminho_arquivo, 'r') as f:
            conteudo = f.read()

        # Coleta dos tempos
        tempos_init_seq = [float(t) for t in re.findall(r"Tempo para init matriz seq: ([\d.]+)", conteudo)]
        tempos_init_par = [float(t) for t in re.findall(r"Tempo para init matriz par: ([\d.]+)", conteudo)]
        tempos_lcs_seq  = [float(t) for t in re.findall(r"Tempo para LCS seq: ([\d.]+)", conteudo)]
        tempos_lcs_par  = [float(t) for t in re.findall(r"Tempo para LCS par: ([\d.]+)", conteudo)]

        # Médias e variâncias
        media_init_seq, var_init_seq = calcular_media_variancia(tempos_init_seq)
        media_init_par, var_init_par = calcular_media_variancia(tempos_init_par)
        media_lcs_seq,  var_lcs_seq  = calcular_media_variancia(tempos_lcs_seq)
        media_lcs_par,  var_lcs_par  = calcular_media_variancia(tempos_lcs_par)

        # Caminho do novo arquivo
        nome_saida = "results_" + nome_arquivo
        caminho_saida = os.path.join(pasta_resultados, nome_saida)

        # Escrita
        with open(caminho_saida, 'w') as f:
            f.write(f"Média Tempo init matriz seq: {media_init_seq:.10f} segundos\n")
            f.write(f"Variância Tempo init matriz seq: {var_init_seq:.10f}\n\n")
            
            f.write(f"Média Tempo init matriz par: {media_init_par:.10f} segundos\n")
            f.write(f"Variância Tempo init matriz par: {var_init_par:.10f}\n\n")
            
            f.write(f"Média Tempo LCS seq:        {media_lcs_seq:.10f} segundos\n")
            f.write(f"Variância Tempo LCS seq:    {var_lcs_seq:.10f}\n\n")
            
            f.write(f"Média Tempo LCS par:        {media_lcs_par:.10f} segundos\n")
            f.write(f"Variância Tempo LCS par:    {var_lcs_par:.10f}\n")

        print(f"Processado: {nome_arquivo} → {caminho_saida}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python read_outputs.py <caminho_para_diretorio>")
    else:
        calcular_tempos_medios_em_arquivos(sys.argv[1])

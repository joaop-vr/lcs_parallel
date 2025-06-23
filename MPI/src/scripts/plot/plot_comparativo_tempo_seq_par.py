#!/usr/bin/env python3
import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# Regex para extrair N do nome do arquivo
SIZE_RE = re.compile(r'teste_(\d+)[A-Za-z]', re.IGNORECASE)
# Regex para capturar tempo total e desvio relativo
TOTAL_RE = re.compile(r'(?:tempo total|Tempo total)\s*[,: ]+\s*([0-9.]+)\s*[, ]+\s*([0-9.]+)', re.IGNORECASE)

def coletar(root_dir):
    """
    Procura por arquivos .out em root_dir ou seu primeiro subdiretório não-vazio,
    extrai (N, mean, std_abs) da linha 'tempo total, mean, rel_std'.
    Retorna três numpy.arrays: Ns, means, stds_abs.
    """
    # lista de .out no nível atual
    arquivos = [f for f in os.listdir(root_dir) if f.lower().endswith('.out')]
    if not arquivos:
        # busca em subdiretórios
        for sub in sorted(os.listdir(root_dir)):
            subp = os.path.join(root_dir, sub)
            if os.path.isdir(subp):
                out_sub = [f for f in os.listdir(subp) if f.lower().endswith('.out')]
                if out_sub:
                    root_dir = subp
                    arquivos = out_sub
                    print(f"> usando arquivos em {root_dir}")
                    break

    if not arquivos:
        raise RuntimeError(f'Nenhum arquivo .out encontrado em {root_dir}')

    dados = []
    for fn in sorted(arquivos):
        m = SIZE_RE.search(fn)
        if not m:
            print(f"⚠️ Não encontrou tamanho no nome: {fn}")
            continue
            
        N = int(m.group(1))
        file_path = os.path.join(root_dir, fn)
        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                txt = f.read()
        except Exception as e:
            print(f"⚠️ Erro ao ler {file_path}: {e}")
            continue
            
        mt = TOTAL_RE.search(txt)
        if not mt:
            print(f"⚠️ Não encontrou 'tempo total' em {fn}")
            continue
            
        try:
            mean = float(mt.group(1))
            rel_std = float(mt.group(2))
            std_abs = mean * rel_std
            dados.append((N, mean, std_abs))
            print(f"✅ {fn}: N={N}, tempo={mean:.4f} ± {std_abs:.4f}")
        except ValueError as e:
            print(f"⚠️ Erro de conversão em {fn}: {e}")

    if not dados:
        raise RuntimeError(f'Arquivos .out encontrados, mas nenhum com linha "tempo total" válida em {root_dir}')

    dados.sort(key=lambda x: x[0])
    Ns, means, stds = zip(*dados)
    return np.array(Ns), np.array(means), np.array(stds)


def plot_e_salvar(Ns, means, stds_abs, title, xlabel, ylabel, output_path):
    plt.figure(figsize=(10, 6))
    
    # Cria um array de posições igualmente espaçadas para os valores de N
    x_positions = np.arange(len(Ns))
    
    # Plot principal com barras de erro usando as posições definidas
    plt.errorbar(
        x_positions, means,  # Usa as posições em vez dos valores de N diretamente
        yerr=stds_abs, 
        fmt='o-',         
        linewidth=2,      
        markersize=8,     
        capsize=5,        
        capthick=2,       
        ecolor='#d62728', 
        elinewidth=2,     
        label='Média e Desvio Padrão'
    )
    
    # Preenchimento também usa as posições definidas
    plt.fill_between(
        x_positions, 
        means - stds_abs, 
        means + stds_abs,
        alpha=0.2,        
        color='#1f77b4',  
        label='Intervalo de Desvio Padrão'
    )
    
    plt.title(title, fontsize=16, pad=20)
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    
    # Configurações de grade e escala
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    
    # Define os ticks do eixo x nas posições criadas, com os rótulos sendo os valores de N
    plt.xticks(x_positions, [str(n) for n in Ns], fontsize=12, rotation=45)
    plt.yticks(fontsize=12)
    
    # Adiciona legenda
    plt.legend(fontsize=12, loc='best')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'✅ Gráfico salvo: {output_path}')

def main():
    parser = argparse.ArgumentParser(
        description='Plota tempo total ± desvio: sequencial e paralelo em arquivos separados.'
    )
    parser.add_argument('seq_dir', help='Diretório (ou pai) com resultados sequenciais')
    parser.add_argument('par_dir', help='Diretório (ou pai) com resultados paralelos')
    parser.add_argument('--seq-out',  default='tempo_seq.png',
                        help='PNG de saída para sequencial')
    parser.add_argument('--par-out',  default='tempo_par.png',
                        help='PNG de saída para paralelo')
    args = parser.parse_args()

    # Sequencial
    print("\nColetando dados sequenciais...")
    try:
        Ns_seq, mean_seq, std_seq = coletar(args.seq_dir)
        plot_e_salvar(
            Ns_seq, mean_seq, std_seq,
            title='Código Sequencial: Tempo Total e Desvio Padrão',
            xlabel='Tamanho de entrada $N$',
            ylabel='Tempo total (s)',
            output_path=args.seq_out
        )
    except Exception as e:
        print(f"❌ Erro no processamento sequencial: {e}")

    # Paralelo
    print("\nColetando dados paralelos...")
    try:
        Ns_par, mean_par, std_par = coletar(args.par_dir)
        plot_e_salvar(
            Ns_par, mean_par, std_par,
            title='Código Paralelo: Tempo Total e Desvio Padrão',
            xlabel='Tamanho de entrada $N$',
            ylabel='Tempo total (s)',
            output_path=args.par_out
        )
    except Exception as e:
        print(f"❌ Erro no processamento paralelo: {e}")

if __name__ == '__main__':
    main()
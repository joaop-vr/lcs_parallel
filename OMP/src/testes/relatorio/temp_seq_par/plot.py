import os
import re
import matplotlib.pyplot as plt
import sys

def parse_results(file_content):
    pattern = r"Média Tempo (.*?) seq:\s+([\d\.]+).*?\nVariância Tempo \1 seq:\s+([\d\.]+).*?\n\nMédia Tempo \1 par:\s+([\d\.]+).*?\nVariância Tempo \1 par:\s+([\d\.]+)"
    matches = re.findall(pattern, file_content, re.DOTALL)

    data = {}
    for label, mean_seq, var_seq, mean_par, var_par in matches:
        key = label.strip()
        data[key] = {
            'seq': {'mean': float(mean_seq), 'var': float(var_seq)},
            'par': {'mean': float(mean_par), 'var': float(var_par)}
        }

    return data

def extract_input_size(filename):
    match = re.search(r"size_(\d+)", filename)
    return int(match.group(1)) if match else None

def plot_total_times(entries, exec_type, output_path):
    entries.sort(key=lambda x: x['size'])
    x = [entry['size'] for entry in entries]
    y = [entry['total_mean'] for entry in entries]
    yerr = [entry['total_var'] for entry in entries]

    plt.figure(figsize=(10, 6))
    plt.errorbar(x, y, yerr=yerr, fmt='-o', capsize=5, ecolor='gray',
                 color='blue' if exec_type == 'seq' else 'green', label=exec_type.upper())
    plt.title(f"Tempo total {'sequencial' if exec_type == 'seq' else 'paralelo'}")
    plt.xlabel("Tamanho da entrada")
    plt.ylabel("Tempo total (s)")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_comparative(seq_entries, par_entries, output_path):
    # Ordenar e garantir que ambos tenham os mesmos tamanhos
    seq_entries.sort(key=lambda x: x['size'])
    par_entries.sort(key=lambda x: x['size'])

    x_seq = [entry['size'] for entry in seq_entries]
    y_seq = [entry['total_mean'] for entry in seq_entries]
    yerr_seq = [entry['total_var'] for entry in seq_entries]

    x_par = [entry['size'] for entry in par_entries]
    y_par = [entry['total_mean'] for entry in par_entries]
    yerr_par = [entry['total_var'] for entry in par_entries]

    plt.figure(figsize=(10, 6))
    plt.errorbar(x_seq, y_seq, yerr=yerr_seq, fmt='-o', capsize=5, label='Sequencial', color='blue')
    plt.errorbar(x_par, y_par, yerr=yerr_par, fmt='-o', capsize=5, label='Paralelo', color='green')
    plt.title("Comparativo: Tempo total Sequencial vs Paralelo")
    plt.xlabel("Tamanho da entrada")
    plt.ylabel("Tempo total (s)")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def main(directory):
    seq_entries = []
    par_entries = []

    for filename in os.listdir(directory):
        if filename.startswith("results_") and filename.endswith(".txt"):
            size = extract_input_size(filename)
            if size is None:
                continue

            filepath = os.path.join(directory, filename)
            with open(filepath, 'r') as f:
                content = f.read()
                results = parse_results(content)

                if 'init matriz' in results and 'LCS' in results:
                    seq_total = results['init matriz']['seq']['mean'] + results['LCS']['seq']['mean']
                    seq_var = results['init matriz']['seq']['var'] + results['LCS']['seq']['var']

                    par_total = results['init matriz']['par']['mean'] + results['LCS']['par']['mean']
                    par_var = results['init matriz']['par']['var'] + results['LCS']['par']['var']

                    seq_entries.append({'size': size, 'total_mean': seq_total, 'total_var': seq_var})
                    par_entries.append({'size': size, 'total_mean': par_total, 'total_var': par_var})

    # Gera os gráficos individuais
    plot_total_times(seq_entries, 'seq', os.path.join(directory, "tempo_total_sequencial.png"))
    plot_total_times(par_entries, 'par', os.path.join(directory, "tempo_total_paralelo.png"))

    # Gera o gráfico comparativo
    plot_comparative(seq_entries, par_entries, os.path.join(directory, "comparativo_seq_vs_par.png"))

    print("Gráficos gerados com sucesso!")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python script.py <diretório>")
    else:
        main(sys.argv[1])

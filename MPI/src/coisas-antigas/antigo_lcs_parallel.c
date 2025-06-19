#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif

typedef unsigned short mtype;

/* Mapeia A,C,G,T para 0..3; retorna -1 se outro char */
static inline int char_idx(char c) {
    switch(c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return -1;
    }
}

/* 
 * Função para ler arquivos com as sequências (cada linha, retirando '\n' e '\r')
 * Apenas o processo rank 0 lê do disco; depois a string será broadcasted.
 */
char* read_seq(const char *fname, int *out_len) {
    FILE *fseq = fopen(fname, "rt");
    if (!fseq) {
        fprintf(stderr, "Error reading file %s\n", fname);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    fseek(fseq, 0L, SEEK_END);
    long size = ftell(fseq);
    rewind(fseq);

    char *seq = (char*) calloc(size + 1, sizeof(char));
    if (!seq) {
        fprintf(stderr, "Erro allocating memory for sequence %s.\n", fname);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int i = 0;
    while (!feof(fseq)) {
        int c = fgetc(fseq);
        if (c == EOF) break;
        if (c == '\n' || c == '\r') continue;
        seq[i++] = (char)c;
    }
    seq[i] = '\0';
    *out_len = i;
    fclose(fseq);
    return seq;
}

int main(int argc, char *argv[]) {
    int rank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc < 3) {
        if (rank == 0) {
            fprintf(stderr, "Uso: %s <arquivoA> <arquivoB>\n", argv[0]);
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    char *seq_A = NULL, *seq_B = NULL;
    int size_A = 0, size_B = 0;

    /* Processo 0 lê as duas sequências e obtém seus tamanhos */
    if (rank == 0) {
        seq_A = read_seq(argv[1], &size_A);
        seq_B = read_seq(argv[2], &size_B);
    }

    /* Broadcast do tamanho das sequências para todos os ranks */
    MPI_Bcast(&size_A, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&size_B, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* Agora todos os ranks sabem size_A e size_B; alocam espaço para armazenar as strings */
    if (rank != 0) {
        seq_A = (char*) malloc((size_A + 1) * sizeof(char));
        seq_B = (char*) malloc((size_B + 1) * sizeof(char));
        if (!seq_A || !seq_B) {
            fprintf(stderr, "[rank %d] Falha ao alocar espaço para sequências\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    /* Broadcast das próprias sequências (incluindo o '\0' no final) */
    MPI_Bcast(seq_A, size_A + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(seq_B, size_B + 1, MPI_CHAR, 0, MPI_COMM_WORLD);

    // -----------------------------
    // 1) Construção da tabela P4
    // -----------------------------
    // Medição de tempo apenas em rank 0:
    double tP_start = 0.0, tP_end = 0.0;
    if (rank == 0) {
        tP_start = MPI_Wtime();
    }

    int *P4 = NULL;
    if (rank == 0) {
        // Aloca P4[4][size_A+1]
        P4 = (int*) malloc(4 * (size_A + 1) * sizeof(int));
        if (!P4) {
            fprintf(stderr, "[rank %d] Falha ao alocar P4 table\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        // Inicializa P4[idx][0] = 0 para idx=0..3
        for (int idx = 0; idx < 4; ++idx) {
            P4[idx * (size_A + 1) + 0] = 0;
        }
        // (Opcional: construir idx_A para acelerar char_idx de seq_A)
        int *idx_A = (int*) malloc((size_A + 1) * sizeof(int));
        if (!idx_A) {
            fprintf(stderr, "[rank %d] Falha ao alocar idx_A\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        idx_A[0] = -1;
        for (int j = 1; j <= size_A; ++j) {
            idx_A[j] = char_idx(seq_A[j-1]);
        }
        // Preenche para j = 1..size_A
        for (int j = 1; j <= size_A; ++j) {
            int idx_here = idx_A[j]; // 0..3 ou -1
            for (int idx = 0; idx < 4; ++idx) {
                if (idx == idx_here) {
                    P4[idx * (size_A + 1) + j] = j;
                } else {
                    P4[idx * (size_A + 1) + j] = P4[idx * (size_A + 1) + (j-1)];
                }
            }
        }
        free(idx_A);
    } else {
        // Aloca espaço para receber P4 via broadcast
        P4 = (int*) malloc(4 * (size_A + 1) * sizeof(int));
        if (!P4) {
            fprintf(stderr, "[rank %d] Falha ao alocar P4 table (receptor)\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    // Broadcast de P4 para todos
    MPI_Bcast(P4, 4 * (size_A + 1), MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        tP_end = MPI_Wtime();
    }
    // -----------------------------
    // 2) Inicialização de outras estruturas
    // -----------------------------
    // Medição de tempo em rank 0:
    double tInit_start = 0.0, tInit_end = 0.0;
    MPI_Barrier(MPI_COMM_WORLD); // sincroniza antes de medir inicialização
    if (rank == 0) {
        tInit_start = MPI_Wtime();
    }

    // Pré-computar idx_A globalmente para DP direto
    int *idx_A_glob = (int*) malloc((size_A + 1) * sizeof(int));
    if (!idx_A_glob) {
        fprintf(stderr, "[rank %d] Falha ao alocar idx_A_glob\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    idx_A_glob[0] = -1;
    for (int j = 1; j <= size_A; ++j) {
        idx_A_glob[j] = char_idx(seq_A[j-1]); // 0..3 ou -1
    }

    // Buffers para DP: prev_row e curr_row
    mtype *rowA = (mtype*) malloc((size_A + 1) * sizeof(mtype));
    mtype *rowB = (mtype*) malloc((size_A + 1) * sizeof(mtype));
    if (!rowA || !rowB) {
        fprintf(stderr, "[rank %d] Falha ao alocar rowA/rowB\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    mtype *prev_row = rowA;
    mtype *curr_row = rowB;
    // Inicializa prev_row[j] = 0
    for (int j = 0; j <= size_A; ++j) {
        prev_row[j] = 0;
    }

    // Preparar arrays para Allgatherv: sendcounts e displs
    int *sendcounts = (int*) malloc(nprocs * sizeof(int));
    int *displs     = (int*) malloc(nprocs * sizeof(int));
    if (!sendcounts || !displs) {
        fprintf(stderr, "[rank %d] Falha ao alocar sendcounts/displs\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int base_cols = size_A / nprocs;
    int rem_cols  = size_A % nprocs;
    int offset = 0;
    for (int p = 0; p < nprocs; ++p) {
        int cnt = (p < rem_cols) ? (base_cols + 1) : base_cols;
        sendcounts[p] = cnt;
        displs[p]     = offset;
        offset += cnt;
    }
    // full_row não mais necessário se usarmos swap de rowA/rowB e Allgatherv diretamente em curr_row.

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        tInit_end = MPI_Wtime();
    }
    // -----------------------------
    // 3) Loop principal DP: medir tempo em rank 0
    // -----------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    double tDP_start = 0.0, tDP_end = 0.0;
    if (rank == 0) {
        tDP_start = MPI_Wtime();
    }

    // Loop i = 1..size_B
    for (int i = 1; i <= size_B; ++i) {
        curr_row[0] = 0;

        // partição local de colunas
        int local_count = sendcounts[rank];
        int j_offset    = displs[rank];
        int local_j_start = j_offset + 1;
        int local_j_end   = j_offset + local_count;

        int idx_c = char_idx(seq_B[i-1]); // 0..3 ou -1

        // Computa curr_row[j] localmente
        if (idx_c >= 0) {
            // Se caractere válido, usamos P4 e idx_A_glob
            for (int j = local_j_start; j <= local_j_end; ++j) {
                if (idx_c == idx_A_glob[j]) {
                    curr_row[j] = (mtype)(prev_row[j-1] + 1);
                } else {
                    int pos = P4[idx_c * (size_A + 1) + j];
                    if (pos == 0) {
                        curr_row[j] = prev_row[j];
                    } else {
                        mtype v = (mtype)(prev_row[pos-1] + 1);
                        curr_row[j] = (mtype) max(prev_row[j], v);
                    }
                }
            }
        } else {
            // Caractere inesperado: tratamos como sem match => curr_row[j] = prev_row[j]
            for (int j = local_j_start; j <= local_j_end; ++j) {
                curr_row[j] = prev_row[j];
            }
        }

        // Allgatherv para montar a linha completa em curr_row[1..size_A]
        MPI_Allgatherv(
            /*sendbuf*/    &curr_row[local_j_start],
            /*sendcount*/  local_count,
            /*sendtype*/   MPI_UNSIGNED_SHORT,
            /*recvbuf*/    &curr_row[1],
            /*recvcounts*/ sendcounts,
            /*displs*/     displs,
            /*recvtype*/   MPI_UNSIGNED_SHORT,
            MPI_COMM_WORLD
        );
        // curr_row[1..size_A] agora contém a linha i completa; curr_row[0]=0

        // Swap prev_row <-> curr_row
        mtype *tmp = prev_row;
        prev_row = curr_row;
        curr_row = tmp;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        tDP_end = MPI_Wtime();
    }

    int lcs_score = prev_row[size_A];
    if (rank == 0) {
        double time_P    = tP_end - tP_start;
        double time_Init = tInit_end - tInit_start;
        double time_DP   = tDP_end - tDP_start;
        printf("LCS = %d\n", lcs_score);
        printf("Tempo de construção de P:         %.6f segundos\n", time_P);
        printf("Tempo de inicialização (dados):   %.6f segundos\n", time_Init);
        printf("Tempo de cálculo DP (linha-a-linha): %.6f segundos\n", time_DP);
        printf("Tempo total (aprox., P + Init + DP): %.6f segundos\n", time_P + time_Init + time_DP);
    }

    // Libera recursos
    free(seq_A);
    free(seq_B);
    free(P4);
    free(idx_A_glob);
    free(rowA);
    free(rowB);
    free(sendcounts);
    free(displs);

    MPI_Finalize();
    return EXIT_SUCCESS;
}

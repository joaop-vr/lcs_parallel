#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif

typedef unsigned short mtype;

/* 
 * Função para ler arquivos com as sequências (cada linha, retirando '\n')
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
        if (c == '\n') continue;
        seq[i++] = (char)c;
    }
    seq[i] = '\0';
    *out_len = i;
    fclose(fseq);
    return seq;
}

/* Aloca matriz (size_B+1)×(size_A+1) de mtype */
mtype** allocate_matrix(int size_A, int size_B) {
    mtype **M = (mtype**) malloc((size_B + 1) * sizeof(mtype*));
    if (!M) {
        fprintf(stderr, "Falha ao alocar ponteiros das linhas.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (int i = 0; i <= size_B; ++i) {
        M[i] = (mtype*) malloc((size_A + 1) * sizeof(mtype));
        if (!M[i]) {
            fprintf(stderr, "Falha ao alocar linha %d da matriz.\n", i);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    return M;
}

/* Inicializa primeira linha e primeira coluna com zero */
void init_matrix(mtype **M, int size_A, int size_B) {
    for (int j = 0; j <= size_A; ++j) {
        M[0][j] = 0;
    }
    for (int i = 1; i <= size_B; ++i) {
        M[i][0] = 0;
    }
}

/* Libera matriz */
void free_matrix(mtype **M, int size_B) {
    for (int i = 0; i <= size_B; ++i) {
        free(M[i]);
    }
    free(M);
}

/*
 * Versão MPI do cálculo paralelo de LCS (sem cache blocking).
 * Cada processo mantém a matriz inteira em memória.
 * Para cada diagonal “diag” (i+j = diag), cada rank:
 *   1) calcula os seus próprios elementos da diagonal,
 *   2) preenche um vetor local diag_vals[k] (k = 0..len-1) com os valores que calculou,
 *      deixando 0 em todas as outras posições,
 *   3) faz MPI_Allreduce(diag_vals, recv_diag_vals, len, MPI_INT, MPI_MAX) para obter
 *      o valor de cada posição k da diagonal correta,
 *   4) grava recv_diag_vals[k] de volta em score_matrix[i][j], de modo que todos os processos
 *      terminem a iteração da diagonal com a mesma informação na matriz.
 * No final, o valor LCS = score_matrix[size_B][size_A].
 */
int LCS_mpi(mtype **score_matrix,
            int size_A, int size_B,
            const char *seq_A, const char *seq_B,
            int rank, int nprocs)
{
    int max_len = (size_A < size_B ? size_A : size_B);
    // Aloca buffers para cada diagonal (valores locais e valores reduzidos globalmente)
    int *diag_vals      = (int*) malloc(max_len * sizeof(int));
    int *recv_diag_vals = (int*) malloc(max_len * sizeof(int));
    if (!diag_vals || !recv_diag_vals) {
        fprintf(stderr, "[rank %d] Falha ao alocar buffers de diagonal\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Percorre todas as diagonais (a soma i+j = diag vai de 2 até size_A+size_B)
    for (int diag = 2; diag <= size_A + size_B; ++diag) {
        // limites de i na diagonal corrente
        int start = (diag > size_A ? (diag - size_A) : 1);
        int end   = (diag > size_B ? size_B : (diag - 1));
        int len   = end - start + 1;  // número de elementos nessa diagonal

        // Inicializa diag_vals[0..len-1] = 0
        for (int k = 0; k < len; ++k) {
            diag_vals[k] = 0;
        }

        // Divide a faixa de k = 0..len-1 entre os nprocs (quociente “base” e resto “rem”)
        int base = len / nprocs;
        int rem  = len % nprocs;
        int local_k_start, local_count;
        if (rank < rem) {
            local_count    = base + 1;
            local_k_start  = rank * (base + 1);
        } else {
            local_count   = base;
            local_k_start = rem * (base + 1) + (rank - rem) * base;
        }
        int local_k_end = local_k_start + local_count - 1;

        // Se não houver nenhum elemento para este rank nessa diagonal, local_count pode ser 0
        if (local_count > 0) {
            for (int k = local_k_start; k <= local_k_end; ++k) {
                int i = start + k;
                int j = diag - i;
                // calcula score_matrix[i][j] 
                if (seq_A[j - 1] == seq_B[i - 1]) {
                    score_matrix[i][j] = score_matrix[i - 1][j - 1] + 1;
                } else {
                    mtype cima      = score_matrix[i - 1][j];
                    mtype esquerda  = score_matrix[i][j - 1];
                    score_matrix[i][j] = (cima > esquerda ? cima : esquerda);
                }
                // grava no vetor local de redução
                diag_vals[k] = (int) score_matrix[i][j];
            }
        }
        // Processos que não calcularam k alguma nessa diagonal mantêm diag_vals[k]=0,
        // o que não causa problema porque o valor correto será >= 0.

        // Roda Allreduce para obter, em recv_diag_vals[k], o valor máximo de diag_vals[k] entre todos os ranks.
        // Como somente quem de fato calculou aquela célula colocou o valor “real” em diag_vals[k],
        // e os demais deixaram 0, MPI_MAX preserva o valor correto. Em situações em que o valor LCS seja 0,
        // todos os processos também terão colocado 0, logo MAX(0,0,...,0) é 0 — perfeito.
        MPI_Allreduce(
            /*sendbuf*/ diag_vals,
            /*recvbuf*/ recv_diag_vals,
            /*count*/   len,
            /*datatype*/MPI_INT,
            /*op*/      MPI_MAX,
            MPI_COMM_WORLD
        );

        // Agora atualiza score_matrix[i][j] para TODOS os processos
        for (int k = 0; k < len; ++k) {
            int i = start + k;
            int j = diag - i;
            score_matrix[i][j] = (mtype) recv_diag_vals[k];
        }
        // Fim de iteração sobre a diagonal “diag”. Próxima iteração utiliza apenas dados de diagonais anteriores,
        // que já estão sincronizados em cada processo.
    }

    free(diag_vals);
    free(recv_diag_vals);
    // Retorna o valor em score_matrix[size_B][size_A] (todos os ranks têm o mesmo valor no final)
    return (int) score_matrix[size_B][size_A];
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

    /* Cada rank aloca e inicializa a matriz (size_B+1)×(size_A+1) */
    mtype **score_matrix = allocate_matrix(size_A, size_B);
    init_matrix(score_matrix, size_A, size_B);

    /* Mede tempo (MPI_Wtime) apenas no rank 0 e broadcast do tempo final, se quiser */
    double t0 = MPI_Wtime();
    int lcs_score = LCS_mpi(score_matrix, size_A, size_B, seq_A, seq_B, rank, nprocs);
    double t1 = MPI_Wtime();

    /* Apenas rank 0 imprime o resultado e o tempo decorrido */
    if (rank == 0) {
        printf("LCS = %d\n", lcs_score);
        printf("Tempo total (MPI): %.6f segundos\n", t1 - t0);
    }

    /* Libera recursos e finaliza MPI */
    free_matrix(score_matrix, size_B);
    free(seq_A);
    free(seq_B);
    MPI_Finalize();
    return EXIT_SUCCESS;
}

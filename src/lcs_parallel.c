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

    // *** Construção da tabela P sobre seq_A:
    //    P[c][j] = última posição p ≤ j em seq_A onde seq_A[p-1] == c, ou 0 se não houver.
    //    Usamos dimensão 256 para todos os bytes. Cada entrada é int.
    int *P = NULL;
    if (rank == 0) {
        P = (int*) malloc(256 * (size_A + 1) * sizeof(int));
        if (!P) {
            fprintf(stderr, "[rank %d] Falha ao alocar P table\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        // Inicializa P[c][0] = 0 para todo c
        for (int c = 0; c < 256; ++c) {
            P[c * (size_A + 1) + 0] = 0;
        }
        // Preenche para j = 1..size_A
        for (int j = 1; j <= size_A; ++j) {
            unsigned char ch = (unsigned char) seq_A[j-1];
            for (int c = 0; c < 256; ++c) {
                if ((unsigned char)c == ch) {
                    P[c * (size_A + 1) + j] = j;
                } else {
                    P[c * (size_A + 1) + j] = P[c * (size_A + 1) + (j-1)];
                }
            }
        }
    } else {
        // Aloca espaço para receber P via broadcast
        P = (int*) malloc(256 * (size_A + 1) * sizeof(int));
        if (!P) {
            fprintf(stderr, "[rank %d] Falha ao alocar P table (receptor)\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    // Broadcast de P para todos
    MPI_Bcast(P, 256 * (size_A + 1), MPI_INT, 0, MPI_COMM_WORLD);

    // *** Preparação para DP linha-a-linha
    // prev_row e curr_row de tamanho size_A+1
    mtype *prev_row = (mtype*) malloc((size_A + 1) * sizeof(mtype));
    mtype *curr_row = (mtype*) malloc((size_A + 1) * sizeof(mtype));
    if (!prev_row || !curr_row) {
        fprintf(stderr, "[rank %d] Falha ao alocar prev_row/curr_row\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    // Inicializa prev_row[j] = 0
    for (int j = 0; j <= size_A; ++j) {
        prev_row[j] = 0;
    }

    // *** Preparar arrays para Allgatherv: sendcounts e displacements
    int *sendcounts = (int*) malloc(nprocs * sizeof(int));
    int *displs     = (int*) malloc(nprocs * sizeof(int));
    if (!sendcounts || !displs) {
        fprintf(stderr, "[rank %d] Falha ao alocar sendcounts/displs\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    // Calcula partição estática de colunas 1..size_A entre processos
    // Cada processo terá sendcounts[p] colunas (p.ex., local_count), e displs[p] é deslocamento em elementos
    int base = size_A / nprocs;
    int rem  = size_A % nprocs;
    int offset = 0;
    for (int p = 0; p < nprocs; ++p) {
        int cnt = (p < rem) ? (base + 1) : base;
        sendcounts[p] = cnt;
        displs[p]     = offset; // offset em unidades de mtype, mas Allgatherv exige desloc em elemento
        offset += cnt;
    }
    // Note: sum sendcounts[p] == size_A

    // Medição de tempo: apenas rank 0
    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    // Loop i = 1..size_B: iterações sobre prefixos de seq_B
    for (int i = 1; i <= size_B; ++i) {
        // prev_row[0] deve ser 0
        curr_row[0] = 0;

        // Determina local_j_start e local_count usando a mesma partição
        int local_count = sendcounts[rank];
        int j_offset    = displs[rank]; // offset de 0..(size_A-1) para colunas 1..size_A
        int local_j_start = j_offset + 1;         // índice j global inicial
        int local_j_end   = j_offset + local_count; // índice j global final

        unsigned char c = (unsigned char) seq_B[i-1]; // caractere atual de seq_B

        // Computa curr_row[j] para j = local_j_start..local_j_end
        for (int j = local_j_start; j <= local_j_end; ++j) {
            if ((unsigned char)seq_A[j-1] == c) {
                // match direto
                curr_row[j] = (mtype)(prev_row[j-1] + 1);
            } else {
                int pos = P[c * (size_A + 1) + j]; // última ocorrência de c em seq_A até j
                if (pos == 0) {
                    curr_row[j] = prev_row[j];
                } else {
                    // prev_row[pos-1] + 1 em mtype
                    mtype v = (mtype)(prev_row[pos-1] + 1);
                    curr_row[j] = (mtype) max(prev_row[j], v);
                }
            }
        }

        // Agora reunimos todas as partes curr_row em full_curr_row compartilhada em cada processo.
        // Podemos usar MPI_Allgatherv: cada processo envia curr_row[ local_j_start .. local_j_end ]
        // Destino: buffer full_curr_row[1..size_A], index 0 deixamos =0.
        // Para convenção, podemos usar mesmo curr_row como buffer de envio/recepção:
        //   - Precisamos um buffer temporário full_row de mtype[size_A+1].
        mtype *full_row = (mtype*) malloc((size_A + 1) * sizeof(mtype));
        if (!full_row) {
            fprintf(stderr, "[rank %d] Falha ao alocar full_row\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        // full_row[0] = 0
        full_row[0] = 0;
        // MPI_Allgatherv espera ponteiro para o início do buffer local: &curr_row[local_j_start]
        // Mas curr_row[0..size_A], offset local_j_start: index local_j_start.
        MPI_Allgatherv(
            /*sendbuf*/    &curr_row[local_j_start],
            /*sendcount*/  local_count,
            /*sendtype*/   MPI_UNSIGNED_SHORT,
            /*recvbuf*/    &full_row[1],
            /*recvcounts*/ sendcounts,
            /*displs*/     displs,
            /*recvtype*/   MPI_UNSIGNED_SHORT,
            MPI_COMM_WORLD
        );
        // Note: full_row[1 + displs[p] .. 1 + displs[p] + sendcounts[p] - 1] receberá dados do processo p.
        // full_row[0] já definido = 0.

        // Troca prev_row <-> full_row para próxima iteração
        // Copiamos full_row em prev_row (poderíamos trocar ponteiros, mas prev_row foi alocado e usado, melhor copiar).
        for (int j = 0; j <= size_A; ++j) {
            prev_row[j] = full_row[j];
        }
        free(full_row);
        // Passa para próxima i
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();

    int lcs_score = prev_row[size_A];

    /* Apenas rank 0 imprime o resultado e o tempo decorrido */
    if (rank == 0) {
        printf("LCS = %d\n", lcs_score);
        printf("Tempo total (MPI, P-based) : %.6f segundos\n", t1 - t0);
    }

    // Libera recursos
    free(prev_row);
    free(curr_row);
    free(P);
    free(sendcounts);
    free(displs);
    free(seq_A);
    free(seq_B);

    MPI_Finalize();
    return EXIT_SUCCESS;
}


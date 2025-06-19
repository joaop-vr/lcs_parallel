#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif

typedef unsigned short mtype;

static inline int char_idx(char c) {
    switch(c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return -1;
    }
}

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

    if (rank == 0) {
        seq_A = read_seq(argv[1], &size_A);
        seq_B = read_seq(argv[2], &size_B);
    }

    MPI_Bcast(&size_A, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&size_B, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        seq_A = (char*) malloc((size_A + 1) * sizeof(char));
        seq_B = (char*) malloc((size_B + 1) * sizeof(char));
        if (!seq_A || !seq_B) {
            fprintf(stderr, "[rank %d] Falha ao alocar espaço para sequências\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    MPI_Bcast(seq_A, size_A + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(seq_B, size_B + 1, MPI_CHAR, 0, MPI_COMM_WORLD);

    // 1. Construção da tabela P4
    double tP_start = 0.0, tP_end = 0.0;
    if (rank == 0) tP_start = MPI_Wtime();

    int *P4 = (int*) malloc(4 * (size_A + 1) * sizeof(int));
    if (!P4) {
        fprintf(stderr, "[rank %d] Falha ao alocar P4 table\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (rank == 0) {
        for (int idx = 0; idx < 4; ++idx) {
            P4[idx * (size_A + 1) + 0] = 0;
        }
        int *idx_A = (int*) malloc((size_A + 1) * sizeof(int));
        if (!idx_A) {
            fprintf(stderr, "[rank %d] Falha ao alocar idx_A\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        idx_A[0] = -1;
        for (int j = 1; j <= size_A; ++j) {
            idx_A[j] = char_idx(seq_A[j-1]);
        }
        for (int j = 1; j <= size_A; ++j) {
            int idx_here = idx_A[j];
            for (int idx = 0; idx < 4; ++idx) {
                if (idx == idx_here) {
                    P4[idx * (size_A + 1) + j] = j;
                } else {
                    P4[idx * (size_A + 1) + j] = P4[idx * (size_A + 1) + (j-1)];
                }
            }
        }
        free(idx_A);
    }

    MPI_Bcast(P4, 4 * (size_A + 1), MPI_INT, 0, MPI_COMM_WORLD);
    if (rank == 0) tP_end = MPI_Wtime();

    // 2. Configuração do grid 2D
    double tInit_start = 0.0, tInit_end = 0.0;
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) tInit_start = MPI_Wtime();

    int dims[2] = {0, 0};
    MPI_Dims_create(nprocs, 2, dims);
    int periods[2] = {0, 0};
    MPI_Comm grid_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &grid_comm);
    
    int grid_rank;
    MPI_Comm_rank(grid_comm, &grid_rank);
    int coords[2];
    MPI_Cart_coords(grid_comm, grid_rank, 2, coords);
    int pi = coords[0], pj = coords[1];
    
    // Divisão das sequências
    int block_size_i = (size_B + dims[0] - 1) / dims[0];
    int block_size_j = (size_A + dims[1] - 1) / dims[1];
    
    int start_i = pi * block_size_i;
    int end_i = (pi == dims[0]-1) ? size_B : (pi+1) * block_size_i;
    int size_i = end_i - start_i;
    
    int start_j = pj * block_size_j;
    int end_j = (pj == dims[1]-1) ? size_A : (pj+1) * block_size_j;
    int size_j = end_j - start_j;
    
    // Alocação do bloco local
    mtype **local_block = (mtype**) malloc((size_i+1) * sizeof(mtype*));
    for (int i = 0; i <= size_i; i++) {
        local_block[i] = (mtype*) calloc(size_j+1, sizeof(mtype));
    }

    // Buffers para comunicação
    mtype *top_border = (mtype*) calloc(size_j+1, sizeof(mtype));
    mtype *left_border = (mtype*) calloc(size_i+1, sizeof(mtype));
    mtype top_left_corner = 0;
    
    // Tags para comunicação
    #define TOP_TAG 0
    #define LEFT_TAG 1
    #define CORNER_TAG 2

    // Wavefront calculation
    MPI_Barrier(MPI_COMM_WORLD);
    double tDP_start = 0.0, tDP_end = 0.0;
    if (rank == 0) tDP_start = MPI_Wtime();

    for (int wave = 0; wave < dims[0] + dims[1] - 1; wave++) {
        for (int i = 0; i < dims[0]; i++) {
            for (int j = 0; j < dims[1]; j++) {
                if (pi + pj == wave && pi == i && pj == j) {
                    // Receber bordas
                    if (pi > 0) {
                        int top_coords[2] = {pi-1, pj};
                        int top_rank;
                        MPI_Cart_rank(grid_comm, top_coords, &top_rank);
                        MPI_Recv(top_border, size_j+1, MPI_UNSIGNED_SHORT, 
                                 top_rank, TOP_TAG, grid_comm, MPI_STATUS_IGNORE);
                    }
                    
                    if (pj > 0) {
                        int left_coords[2] = {pi, pj-1};
                        int left_rank;
                        MPI_Cart_rank(grid_comm, left_coords, &left_rank);
                        MPI_Recv(left_border, size_i+1, MPI_UNSIGNED_SHORT, 
                                 left_rank, LEFT_TAG, grid_comm, MPI_STATUS_IGNORE);
                    }
                    
                    if (pi > 0 && pj > 0) {
                        int corner_coords[2] = {pi-1, pj-1};
                        int corner_rank;
                        MPI_Cart_rank(grid_comm, corner_coords, &corner_rank);
                        MPI_Recv(&top_left_corner, 1, MPI_UNSIGNED_SHORT, 
                                 corner_rank, CORNER_TAG, grid_comm, MPI_STATUS_IGNORE);
                    }

                    // Inicializar bordas do bloco
                    if (pi == 0 && pj == 0) {
                        for (int jj = 0; jj <= size_j; jj++) local_block[0][jj] = 0;
                        for (int ii = 0; ii <= size_i; ii++) local_block[ii][0] = 0;
                    } else {
                        if (pi > 0) {
                            for (int jj = 0; jj <= size_j; jj++) local_block[0][jj] = top_border[jj];
                        } else {
                            for (int jj = 0; jj <= size_j; jj++) local_block[0][jj] = 0;
                        }
                        if (pj > 0) {
                            for (int ii = 0; ii <= size_i; ii++) local_block[ii][0] = left_border[ii];
                        } else {
                            for (int ii = 0; ii <= size_i; ii++) local_block[ii][0] = 0;
                        }
                        if (pi > 0 && pj > 0) {
                            local_block[0][0] = top_left_corner;
                        }
                    }

                    // Calcular bloco
                    for (int ii = 1; ii <= size_i; ii++) {
                        int global_i = start_i + ii;
                        int idx_c = char_idx(seq_B[global_i-1]);
                        
                        for (int jj = 1; jj <= size_j; jj++) {
                            int global_j = start_j + jj;
                            
                            if (seq_A[global_j-1] == seq_B[global_i-1]) {
                                local_block[ii][jj] = local_block[ii-1][jj-1] + 1;
                            } else {
                                if (idx_c >= 0) {
                                    int pos = P4[idx_c * (size_A + 1) + global_j];
                                    mtype v = 0;
                                    if (pos > 0 && pos-1 >= start_j && pos-1 < end_j) {
                                        int local_col = pos - 1 - start_j;
                                        v = local_block[ii-1][local_col] + 1;
                                    }
                                    local_block[ii][jj] = max(local_block[ii-1][jj], 
                                                            max(local_block[ii][jj-1], v));
                                } else {
                                    local_block[ii][jj] = max(local_block[ii-1][jj], local_block[ii][jj-1]);
                                }
                            }
                        }
                    }

                    // Enviar bordas
                    if (pi < dims[0]-1) {
                        int bottom_coords[2] = {pi+1, pj};
                        int bottom_rank;
                        MPI_Cart_rank(grid_comm, bottom_coords, &bottom_rank);
                        MPI_Send(local_block[size_i], size_j+1, MPI_UNSIGNED_SHORT, 
                                bottom_rank, TOP_TAG, grid_comm);
                    }
                    
                    if (pj < dims[1]-1) {
                        int right_coords[2] = {pi, pj+1};
                        int right_rank;
                        MPI_Cart_rank(grid_comm, right_coords, &right_rank);
                        mtype *right_border = (mtype*) malloc((size_i+1) * sizeof(mtype));
                        for (int ii = 0; ii <= size_i; ii++) {
                            right_border[ii] = local_block[ii][size_j];
                        }
                        MPI_Send(right_border, size_i+1, MPI_UNSIGNED_SHORT, 
                                right_rank, LEFT_TAG, grid_comm);
                        free(right_border);
                    }
                    
                    if (pi < dims[0]-1 && pj < dims[1]-1) {
                        int bottom_right_coords[2] = {pi+1, pj+1};
                        int bottom_right_rank;
                        MPI_Cart_rank(grid_comm, bottom_right_coords, &bottom_right_rank);
                        mtype corner = local_block[size_i][size_j];
                        MPI_Send(&corner, 1, MPI_UNSIGNED_SHORT, 
                                bottom_right_rank, CORNER_TAG, grid_comm);
                    }
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) tDP_end = MPI_Wtime();

    // Coleta do resultado final
    int global_lcs = 0;
    if (pi == dims[0]-1 && pj == dims[1]-1) {
        global_lcs = local_block[size_i][size_j];
    }

    int bottom_right_coords[2] = {dims[0]-1, dims[1]-1};
    int bottom_right_rank;
    MPI_Cart_rank(grid_comm, bottom_right_coords, &bottom_right_rank);
    MPI_Bcast(&global_lcs, 1, MPI_INT, bottom_right_rank, grid_comm);

    if (rank == 0) {
        tInit_end = MPI_Wtime();
        double time_P = tP_end - tP_start;
        double time_Init = tInit_end - tInit_start;
        double time_DP = tDP_end - tDP_start;
        printf("LCS = %d\n", global_lcs);
        printf("Tempo de construção de P:         %.6f segundos\n", time_P);
        printf("Tempo de inicialização (dados):   %.6f segundos\n", time_Init);
        printf("Tempo de cálculo DP (wavefront):  %.6f segundos\n", time_DP);
        printf("Tempo total (aprox.):             %.6f segundos\n", time_P + time_Init + time_DP);
    }

    // Liberação de recursos
    free(seq_A);
    free(seq_B);
    free(P4);
    for (int i = 0; i <= size_i; i++) {
        free(local_block[i]);
    }
    free(local_block);
    free(top_border);
    free(left_border);
    MPI_Comm_free(&grid_comm);
    MPI_Finalize();
    return EXIT_SUCCESS;
}

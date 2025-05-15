#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

//#define DEBUGMATRIX

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

typedef unsigned short mtype;


/******************************************
 *            FUNÇÕES AUXILIARES          *
 ******************************************/

/* Função para medir o tempo */
void print_time(double start, double end, const char* stage) {
    printf("Tempo para %s: %.6f segundos\n", stage, end - start);
}

/* Função para ler arquivos com as sequencias */
char* read_seq(char *fname) {
    FILE *fseq = NULL;
    long size = 0;
    char *seq = NULL;
    int i = 0;

    fseq = fopen(fname, "rt");
    if (fseq == NULL) {
        printf("Error reading file %s\n", fname);
        exit(1);
    }

    fseek(fseq, 0L, SEEK_END);
    size = ftell(fseq);
    rewind(fseq);

    seq = (char *) calloc(size + 1, sizeof(char));
    if (seq == NULL) {
        printf("Erro allocating memory for sequence %s.\n", fname);
        exit(1);
    }

    while (!feof(fseq)) {
        seq[i] = fgetc(fseq);
        if ((seq[i] != '\n') && (seq[i] != EOF))
            i++;
    }
    seq[i] = '\0';

    fclose(fseq);

    return seq;
}

void print_matrix(char * seq_A, char * seq_B, mtype ** score_matrix, int size_A, int size_B) {
    int i, j;
    printf("Score Matrix:\n");
    printf("========================================\n");

    printf("    ");
    printf("%5c   ", ' ');

    for (j = 0; j < size_A; j++)
        printf("%5c   ", seq_A[j]);
    printf("\n");
    for (i = 0; i < size_B + 1; i++) {
        if (i == 0)
            printf("    ");
        else
            printf("%c   ", seq_B[i - 1]);
        for (j = 0; j < size_A + 1; j++) {
            printf("%5d   ", score_matrix[i][j]);
        }
        printf("\n");
    }
    printf("========================================\n");
}

void free_matrix(mtype **score_matrix, int size_B) {
    int i;
    for (i = 0; i < (size_B + 1); i++)
        free(score_matrix[i]);
    free(score_matrix);
}

mtype ** allocate_matrix(int size_A, int size_B) {
    int i;
    mtype ** score_matrix = (mtype **) malloc((size_B + 1) * sizeof(mtype *));
    for (i = 0; i < (size_B + 1); i++)
        score_matrix[i] = (mtype *) malloc((size_A + 1) * sizeof(mtype));
    return score_matrix;
}

void init_matrix(mtype ** score_matrix, int size_A, int size_B) {
    int i, j;
    for (j = 0; j < (size_A + 1); j++)
        score_matrix[0][j] = 0;

    for (i = 1; i < (size_B + 1); i++)
        score_matrix[i][0] = 0;
}


/******************************************
 *           FUNÇÃO PARALELA              *
 ******************************************/

// LCS com paralelismo otimizado
int LCS_parallel(mtype **score_matrix, int size_A, int size_B, char *seq_A, char *seq_B) {
    #pragma omp parallel num_threads(8)
    {
        for (int diag = 2; diag <= size_A + size_B; diag++) {
            int start = (diag > size_A) ? (diag - size_A) : 1;
            int end = (diag > size_B) ? size_B : (diag - 1);
            #pragma omp for schedule(guided, 128)
            for (int i = start; i <= end; i++) {
                int j = diag - i;
                if (seq_A[j - 1] == seq_B[i - 1]) {
                    score_matrix[i][j] = score_matrix[i - 1][j - 1] + 1;
                } else {
                    score_matrix[i][j] = max(score_matrix[i - 1][j], score_matrix[i][j - 1]);
                }
            }
        }
    }

    return score_matrix[size_B][size_A];
}



int main(int argc, char ** argv) {

    char *seq_A, *seq_B;
    int size_A, size_B;

    if (argc < 3 ) {
        printf("Uso: %s <arquivoA> <arquivoB>\n", argv[0]);
        return EXIT_FAILURE;
    }

    seq_A = read_seq(argv[1]);
    seq_B = read_seq(argv[2]);
    size_A = strlen(seq_A);
    size_B = strlen(seq_B);

    // ########## PARALELO ##########
    mtype ** score_matrix = allocate_matrix(size_A, size_B);

    init_matrix(score_matrix, size_A, size_B);
    mtype score_par = LCS_parallel(score_matrix, size_A, size_B, seq_A, seq_B);

    #ifdef DEBUGMATRIX
        print_matrix(seq_A, seq_B, score_matrix, size_A, size_B);
    #endif

    printf("%d\n\n", score_par);

    free_matrix(score_matrix, size_B);

    return EXIT_SUCCESS;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

//#define DEBUGMATRIX

#define NUM_THREADS 4
#define CACHE_LINE_SIZE 64

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

typedef unsigned short mtype;

/* Função para medir o tempo */
void printTime(double start, const char* stage) {
    double end = omp_get_wtime();
    printf("Tempo para %s: %.6f segundos\n", stage, end - start);
}

/* Read sequence from a file to a char vector.
   Filename is passed as parameter */
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

// Alocação da matriz com alinhamento para evitar falso compartilhamento
mtype **allocateScoreMatrix(int sizeA, int sizeB) {
    mtype **scoreMatrix = (mtype **)malloc((sizeB + 1) * sizeof(mtype *));
    for (int i = 0; i <= sizeB; i++) {
        size_t row_size = (sizeA + 1) * sizeof(mtype);
        size_t padded_size = ((row_size + CACHE_LINE_SIZE - 1) / CACHE_LINE_SIZE) * CACHE_LINE_SIZE;
        scoreMatrix[i] = (mtype *)aligned_alloc(CACHE_LINE_SIZE, padded_size);
        if (scoreMatrix[i] == NULL) {
            fprintf(stderr, "Erro ao alocar linha %d\n", i);
            exit(EXIT_FAILURE);
        }
    }
    return scoreMatrix;
}

void initScoreMatrix(mtype ** scoreMatrix, int sizeA, int sizeB) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i <= sizeB; i++) {
        for (int j = 0; j <= sizeA; j++) {
            scoreMatrix[i][j] = 0;
        }
    }

}

// LCS com paralelismo otimizado
int LCS(mtype **scoreMatrix, int sizeA, int sizeB, char *seqA, char *seqB) {
    #pragma omp parallel
    {
        for (int d = 2; d <= sizeA + sizeB; d++) {
            int i_start = (d > sizeA) ? (d - sizeA) : 1;
            int i_end = (d > sizeB) ? sizeB : (d - 1);
            #pragma omp for schedule(dynamic, 64)
            for (int i = i_start; i <= i_end; i++) {
                int j = d - i;
                if (seqA[j - 1] == seqB[i - 1]) {
                    scoreMatrix[i][j] = scoreMatrix[i - 1][j - 1] + 1;
                } else {
                    scoreMatrix[i][j] = max(scoreMatrix[i - 1][j], scoreMatrix[i][j - 1]);
                }
            }
        }
    }

    return scoreMatrix[sizeB][sizeA];
}

void printMatrix(char * seqA, char * seqB, mtype ** scoreMatrix, int sizeA, int sizeB) {
    int i, j;
    printf("Score Matrix:\n");
    printf("========================================\n");

    printf("    ");
    printf("%5c   ", ' ');

    for (j = 0; j < sizeA; j++)
        printf("%5c   ", seqA[j]);
    printf("\n");
    for (i = 0; i < sizeB + 1; i++) {
        if (i == 0)
            printf("    ");
        else
            printf("%c   ", seqB[i - 1]);
        for (j = 0; j < sizeA + 1; j++) {
            printf("%5d   ", scoreMatrix[i][j]);
        }
        printf("\n");
    }
    printf("========================================\n");
}

void freeScoreMatrix(mtype **scoreMatrix, int sizeB) {
    int i;
    for (i = 0; i < (sizeB + 1); i++)
        free(scoreMatrix[i]);
    free(scoreMatrix);
}

int main(int argc, char ** argv) {

    double start_time;
    char *seqA, *seqB;
    int sizeA, sizeB;

    if (argc != 3 && argc != 4) {
        printf("Uso: %s <arquivoA> <arquivoB> [num_threads]\n", argv[0]);
        return EXIT_FAILURE;
    }

    // Recebe o número de threads como argumento
    int num_threads = 4; // valor padrão
    if (argc == 4) {
        num_threads = atoi(argv[3]);
        if (num_threads <= 0) {
            printf("Número de threads inválido. Usando o valor padrão (4).\n");
            num_threads = 4; 
        }
    }

    // Define o número de threads
    omp_set_num_threads(num_threads);

    seqA = read_seq(argv[1]);
    seqB = read_seq(argv[2]);
    sizeA = strlen(seqA);
    sizeB = strlen(seqB);

    mtype ** scoreMatrix = allocateScoreMatrix(sizeA, sizeB);

    start_time = omp_get_wtime();
    initScoreMatrix(scoreMatrix, sizeA, sizeB);
    printTime(start_time, "inicialização da matriz de pontuação");

    start_time = omp_get_wtime();
    mtype score = LCS(scoreMatrix, sizeA, sizeB, seqA, seqB);
    printTime(start_time, "cálculo do LCS");

#ifdef DEBUGMATRIX
    printMatrix(seqA, seqB, scoreMatrix, sizeA, sizeB);
#endif

    printf("\nScore: %d\n", score);

    freeScoreMatrix(scoreMatrix, sizeB);
    return EXIT_SUCCESS;
}

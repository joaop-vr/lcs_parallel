#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>  // Adicionando a biblioteca para medir tempo

//#define DEBUGMATRIX

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

typedef unsigned short mtype;

/* Função para medir o tempo */
void printTime(clock_t start, const char* stage) {
    clock_t end = clock();
    double elapsed_time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Tempo para %s: %.6f segundos\n", stage, elapsed_time);
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

mtype ** allocateScoreMatrix(int sizeA, int sizeB) {
    int i;
    mtype ** scoreMatrix = (mtype **) malloc((sizeB + 1) * sizeof(mtype *));
    for (i = 0; i < (sizeB + 1); i++)
        scoreMatrix[i] = (mtype *) malloc((sizeA + 1) * sizeof(mtype));
    return scoreMatrix;
}

void initScoreMatrix(mtype ** scoreMatrix, int sizeA, int sizeB) {
    int i, j;
    for (j = 0; j < (sizeA + 1); j++)
        scoreMatrix[0][j] = 0;

    for (i = 1; i < (sizeB + 1); i++)
        scoreMatrix[i][0] = 0;
}

int LCS(mtype ** scoreMatrix, int sizeA, int sizeB, char * seqA, char *seqB) {
    int i, j;
    for (i = 1; i < sizeB + 1; i++) {
        for (j = 1; j < sizeA + 1; j++) {
            if (seqA[j - 1] == seqB[i - 1]) {
                scoreMatrix[i][j] = scoreMatrix[i - 1][j - 1] + 1;
            } else {
                scoreMatrix[i][j] = max(scoreMatrix[i-1][j], scoreMatrix[i][j-1]);
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
    clock_t start_time;

    char *seqA, *seqB;
    int sizeA, sizeB;

    // Leitura das sequências
    seqA = read_seq("fileA.in");
    seqB = read_seq("fileB.in");

    sizeA = strlen(seqA);
    sizeB = strlen(seqB);

    // Alocar matriz de pontuação
    mtype ** scoreMatrix = allocateScoreMatrix(sizeA, sizeB);

    // Medir o tempo de inicialização da matriz
    start_time = clock();
    initScoreMatrix(scoreMatrix, sizeA, sizeB);
    printTime(start_time, "inicialização da matriz de pontuação");

    // Medir o tempo do cálculo do LCS
    start_time = clock();
    mtype score = LCS(scoreMatrix, sizeA, sizeB, seqA, seqB);
    printTime(start_time, "cálculo do LCS");

#ifdef DEBUGMATRIX
    printMatrix(seqA, seqB, scoreMatrix, sizeA, sizeB);
#endif

    // Imprimir o resultado final
    printf("\nScore: %d\n", score);

    // Liberar a matriz de pontuação
    freeScoreMatrix(scoreMatrix, sizeB);

    return EXIT_SUCCESS;
}

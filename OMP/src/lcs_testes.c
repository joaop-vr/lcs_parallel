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


/******************************************
 *            FUNÇÕES AUXILIARES          *
 ******************************************/

/* Função para medir o tempo */
void printTime(double start, double end, const char* stage) {
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

/******************************************
 *          FUNÇÕES SEQUENCIAIS           *
 ******************************************/

int LCS_sequencial(mtype ** scoreMatrix, int sizeA, int sizeB, char * seqA, char *seqB) {
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


/******************************************
 *           FUNÇÕES PARALELAS            *
 ******************************************/

// LCS com paralelismo otimizado
int LCS_parallel(mtype **scoreMatrix, int sizeA, int sizeB, char *seqA, char *seqB, int n_threads, int chunk_size) {
    #pragma omp parallel num_threads(n_threads)
    {
        for (int diag = 2; diag <= sizeA + sizeB; diag++) {
            int start = (diag > sizeA) ? (diag - sizeA) : 1;
            int end = (diag > sizeB) ? sizeB : (diag - 1);
            #pragma omp for schedule(runtime)
            for (int i = start; i <= end; i++) {
                int j = diag - i;
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



int main(int argc, char ** argv) {

    double start_time, end_time;
    char *seqA, *seqB;
    int sizeA, sizeB;

    if (argc < 3 || argc > 6) {
        printf("Uso: %s <arquivoA> <arquivoB> [num_threads] [tam_chunk] [tipo_schedule]\n", argv[0]);
        printf("Exemplo: %s A.txt B.txt 4 128 dynamic\n", argv[0]);
        return EXIT_FAILURE;
    }

    // Config padrão
    int num_threads = 4;
    int chunk_size = 128;
    char *schedule_type = "static"; 

    // Leitura dos argumentos
    if (argc >= 4) num_threads = atoi(argv[3]);
    if (argc >= 5) chunk_size = atoi(argv[4]);
    if (argc == 6) schedule_type = argv[5];

    // Config do schedule
    omp_sched_t schedule_enum;
    if (strcmp(schedule_type, "static") == 0) {
        schedule_enum = omp_sched_static;
    } else if (strcmp(schedule_type, "dynamic") == 0) {
        schedule_enum = omp_sched_dynamic;
    } else if (strcmp(schedule_type, "guided") == 0) {
        schedule_enum = omp_sched_guided;
    } else {
        fprintf(stderr, "Tipo de escalonamento inválido: %s. Use static, dynamic, guided ou auto.\n", schedule_type);
        return EXIT_FAILURE;
    }
    omp_set_schedule(schedule_enum, chunk_size);

    seqA = read_seq(argv[1]);
    seqB = read_seq(argv[2]);
    sizeA = strlen(seqA);
    sizeB = strlen(seqB);

    printf("Tam A: %i\n", sizeA);
    printf("Tam B: %i\n", sizeB);


    // ########## SEQUENCIAL ##########

    mtype ** scoreMatrixSequencial = allocateScoreMatrix(sizeA, sizeB);

    start_time = omp_get_wtime();
    initScoreMatrix(scoreMatrixSequencial, sizeA, sizeB);
    end_time = omp_get_wtime();
    printTime(start_time, end_time, "init matriz seq");

    start_time = omp_get_wtime();
    mtype score_seq = LCS_sequencial(scoreMatrixSequencial, sizeA, sizeB, seqA, seqB);
    end_time = omp_get_wtime();
    printTime(start_time, end_time, "LCS seq");

    #ifdef DEBUGMATRIX
        printMatrix(seqA, seqB, scoreMatrixSequencial, sizeA, sizeB);
    #endif

    printf("Score seq: %d\n\n", score_seq);

    freeScoreMatrix(scoreMatrixSequencial, sizeB);

    // ########## PARALELO ##########

    mtype ** scoreMatrixParallel = allocateScoreMatrix(sizeA, sizeB);

    start_time = omp_get_wtime();
    initScoreMatrix(scoreMatrixParallel, sizeA, sizeB);
    end_time = omp_get_wtime();
    printTime(start_time, end_time, "init matriz par");

    start_time = omp_get_wtime();
    mtype score_par = LCS_parallel(scoreMatrixParallel, sizeA, sizeB, seqA, seqB, num_threads, chunk_size);
    end_time = omp_get_wtime();
    printTime(start_time, end_time, "LCS par");

    #ifdef DEBUGMATRIX
        printMatrix(seqA, seqB, scoreMatrixParallel, sizeA, sizeB);
    #endif

    printf("Score par: %d\n\n", score_par);

    freeScoreMatrix(scoreMatrixParallel, sizeB);

    return EXIT_SUCCESS;
}

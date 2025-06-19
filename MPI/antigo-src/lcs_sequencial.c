#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

typedef unsigned short mtype;

/* Read sequence from a file to a char vector.
 Filename is passed as parameter */

char* read_seq(char *fname) {
	//file pointer
	FILE *fseq = NULL;
	//sequence size
	long size = 0;
	//sequence pointer
	char *seq = NULL;
	//sequence index
	int i = 0;

	//open file
	fseq = fopen(fname, "rt");
	if (fseq == NULL ) {
		printf("Error reading file %s\n", fname);
		exit(1);
	}

	//find out sequence size to allocate memory afterwards
	fseek(fseq, 0L, SEEK_END);
	size = ftell(fseq);
	rewind(fseq);

	//allocate memory (sequence)
	seq = (char *) calloc(size + 1, sizeof(char));
	if (seq == NULL ) {
		printf("Erro allocating memory for sequence %s.\n", fname);
		exit(1);
	}

	//read sequence from file
	while (!feof(fseq)) {
		seq[i] = fgetc(fseq);
		if ((seq[i] != '\n') && (seq[i] != EOF))
			i++;
	}
	//insert string terminator
	seq[i] = '\0';

	//close file
	fclose(fseq);

	//return sequence pointer
	return seq;
}

mtype ** allocateScoreMatrix(int sizeA, int sizeB) {
	int i;
	//Allocate memory for LCS score matrix
	mtype ** scoreMatrix = (mtype **) malloc((sizeB + 1) * sizeof(mtype *));
	for (i = 0; i < (sizeB + 1); i++)
		scoreMatrix[i] = (mtype *) malloc((sizeA + 1) * sizeof(mtype));
	return scoreMatrix;
}

void initScoreMatrix(mtype ** scoreMatrix, int sizeA, int sizeB) {
	int i, j;
	//Fill first line of LCS score matrix with zeroes
	for (j = 0; j < (sizeA + 1); j++)
		scoreMatrix[0][j] = 0;

	//Do the same for the first collumn
	for (i = 1; i < (sizeB + 1); i++)
		scoreMatrix[i][0] = 0;
}

int LCS(mtype ** scoreMatrix, int sizeA, int sizeB, char * seqA, char *seqB) {
	int i, j;
	for (i = 1; i < sizeB + 1; i++) {
		for (j = 1; j < sizeA + 1; j++) {
			if (seqA[j - 1] == seqB[i - 1]) {
				/* if elements in both sequences match,
				 the corresponding score will be the score from
				 previous elements + 1*/
				scoreMatrix[i][j] = scoreMatrix[i - 1][j - 1] + 1;
			} else {
				/* else, pick the maximum value (score) from left and upper elements*/
				scoreMatrix[i][j] =
						max(scoreMatrix[i-1][j], scoreMatrix[i][j-1]);
			}
		}
	}
	return scoreMatrix[sizeB][sizeA];
}
void printMatrix(char * seqA, char * seqB, mtype ** scoreMatrix, int sizeA,
		int sizeB) {
	int i, j;

	//print header
	printf("Score Matrix:\n");
	printf("========================================\n");

	//print LCS score matrix allong with sequences

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


int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "Uso: %s <arquivoA.in> <arquivoB.in>\n", argv[0]);
        return EXIT_FAILURE;
    }

    // marcadores de tempo usando clock()
    clock_t t_start, t_after_init, t_after_lcs, t_end;
    double t_init, t_lcs, t_total;

    // Início do timer total
    t_start = clock();

    // leitura das sequências
    char *seqA = read_seq(argv[1]);
    if (!seqA) {
        fprintf(stderr, "Erro ao ler sequência A: %s\n", argv[1]);
        return EXIT_FAILURE;
    }
    char *seqB = read_seq(argv[2]);
    if (!seqB) {
        fprintf(stderr, "Erro ao ler sequência B: %s\n", argv[2]);
        free(seqA);
        return EXIT_FAILURE;
    }

    int sizeA = strlen(seqA);
    int sizeB = strlen(seqB);

    // alocação e inicialização da matriz
    mtype **scoreMatrix = allocateScoreMatrix(sizeA, sizeB);
    if (!scoreMatrix) {
        fprintf(stderr, "Falha ao alocar matriz de pontuações\n");
        free(seqA);
        free(seqB);
        return EXIT_FAILURE;
    }
    initScoreMatrix(scoreMatrix, sizeA, sizeB);

    // Marca o fim da inicialização
    t_after_init = clock();

    // cálculo do LCS
    mtype score = LCS(scoreMatrix, sizeA, sizeB, seqA, seqB);

    // Marca o fim do cálculo do LCS
    t_after_lcs = clock();

    // libera estruturas
    freeScoreMatrix(scoreMatrix, sizeB);
    free(seqA);
    free(seqB);

    // Marca o fim total
    t_end = clock();

    // converte em segundos
    t_init  = (double)(t_after_init - t_start)   / CLOCKS_PER_SEC;
    t_lcs   = (double)(t_after_lcs  - t_after_init) / CLOCKS_PER_SEC;
    t_total = (double)(t_end         - t_start)   / CLOCKS_PER_SEC;

    // imprime tempos e score
    printf("Tempos (segundos):\n");
    printf("  alocação+init: %.6f\n", t_init);
    printf("  LCS compute : %.6f\n", t_lcs);
    printf("  tempo total : %.6f\n", t_total);
    printf("Score: %d\n\n", score);

    return EXIT_SUCCESS;
}

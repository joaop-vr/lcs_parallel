#include <stdio.h>
#include <stdlib.h>
#include <time.h>

char random_base() {
    char bases[] = {'A', 'C', 'G', 'T'};
    return bases[rand() % 4];
}

void generate_sequence(FILE *file, int length) {
    for (int i = 0; i < length; i++) {
        fputc(random_base(), file);
    }
    fputc('\n', file); // newline no final
}

int main(int argc, char **argv) {
    if (argc != 4) {
        printf("Uso: %s <tamanho_seqA> <tamanho_seqB> <prefixo_saida>\n", argv[0]);
        return 1;
    }

    int lenA = atoi(argv[1]);
    int lenB = atoi(argv[2]);
    char *prefix = argv[3];

    if (lenA <= 0 || lenB <= 0) {
        printf("Tamanhos devem ser positivos.\n");
        return 1;
    }

    srand(time(NULL));

    char fileA[256], fileB[256];
    snprintf(fileA, sizeof(fileA), "%sA.in", prefix);
    snprintf(fileB, sizeof(fileB), "%sB.in", prefix);

    FILE *outA = fopen(fileA, "w");
    FILE *outB = fopen(fileB, "w");

    if (!outA || !outB) {
        perror("Erro ao abrir arquivos de saída");
        return 1;
    }

    generate_sequence(outA, lenA);
    generate_sequence(outB, lenB);

    fclose(outA);
    fclose(outB);

    printf("Arquivos gerados: %s e %s\n", fileA, fileB);
    return 0;
}

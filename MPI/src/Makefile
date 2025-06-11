# Variáveis para o compilador, flags e arquivos
CC = gcc
CFLAGS = -Wall -O3 -fopenmp   # Flags de compilação, incluindo OpenMP
LDFLAGS = -fopenmp            # Flags de linkedição, incluindo OpenMP
OBJ = lcs.o                   # Arquivo objeto
EXEC_LCS = lcs                # Executável principal
EXEC_GEN = generate_input     # Executável de geração de input

# Regra padrão (compila ambos)
all: $(EXEC_LCS) $(EXEC_GEN)

# Regra para gerar o lcs
$(EXEC_LCS): $(OBJ)
	$(CC) $(LDFLAGS) $(OBJ) -o $(EXEC_LCS)   # Linkando com OpenMP

# Regra para gerar o objeto
$(OBJ): lcs.c
	$(CC) $(CFLAGS) -c lcs.c

# Regra para generate_input
$(EXEC_GEN): generate_input.c
	$(CC) -Wall -o $(EXEC_GEN) generate_input.c

# Limpeza dos arquivos gerados
clean:
	rm -f $(OBJ) $(EXEC_LCS) $(EXEC_GEN)

# Compilação com debug
debug: CFLAGS += -DDEBUGMATRIX
debug: clean all

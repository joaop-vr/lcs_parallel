# Variáveis para o compilador, flags e arquivos
CC = gcc
CFLAGS = -Wall -O3 -fopenmp    # Flags de compilação, incluindo OpenMP
LDFLAGS = -fopenmp            # Flags de linkedição, incluindo OpenMP
OBJ = lcs.o                   # Arquivo objeto
EXEC = lcs                    # Nome do executável final

# Regra padrão (vai compilar o programa)
all: $(EXEC)

# Regra para gerar o executável
$(EXEC): $(OBJ)
	$(CC) $(LDFLAGS) $(OBJ) -o $(EXEC)   # Linkando com OpenMP

# Regra para gerar o arquivo objeto
$(OBJ): lcs.c
	$(CC) $(CFLAGS) -c lcs.c            # Compilando com OpenMP

# Limpeza dos arquivos gerados (executável e objetos)
clean:
	rm -f $(OBJ) $(EXEC)

# Regra de depuração (opcional)
debug: CFLAGS += -DDEBUGMATRIX
debug: clean $(EXEC)

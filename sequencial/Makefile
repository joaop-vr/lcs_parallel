# Variáveis
CC = gcc                # Compilador
CFLAGS = -Wall -g       # Flags de compilação (gerar informações de depuração)
SRC = lcs.c             # Arquivo fonte
OBJ = lcs.o             # Arquivo objeto
EXEC = lcs              # Executável gerado
INCLUDES =              # Incluir bibliotecas adicionais, se necessário

# Regra para compilar o programa
all: $(EXEC)

# Regra para gerar o executável
$(EXEC): $(OBJ)
	$(CC) $(OBJ) -o $(EXEC)

# Regra para gerar o arquivo objeto a partir do código fonte
$(OBJ): $(SRC)
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRC)

# Limpeza dos arquivos gerados (executável, objetos)
clean:
	rm -f $(OBJ) $(EXEC)

# Regra de depuração (opcional)
debug: CFLAGS += -DDEBUGMATRIX
debug: clean $(EXEC)

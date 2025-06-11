#############################
# Makefile para programa MPI
#############################

# Compilador e flags
CC      = mpicc
CFLAGS  = -std=c99 -Wall -O3

# Nome do executável
TARGET  = lcs_mpi

# Arquivos-fonte
SRCS    = lcs.c

# Regra padrão
all: $(TARGET)

# Como compilar o executável
$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -o $@ $^

# Para rodar rapidamente (exemplo):
# make run ARGS="fileA.txt fileB.txt"
run: $(TARGET)
	mpirun -np 4 ./$(TARGET) $(ARGS)

# Limpeza
.PHONY: clean
clean:
	rm -f $(TARGET) *.o

# Exibir ajuda
.PHONY: help
default: help
help:
	@echo "Uso: make [all|run|clean]"
	@echo "  all   : compila o executável ($(TARGET))"
	@echo "  run   : executa com mpirun (use 'make run ARGS=\"fileA fileB\"')"
	@echo "  clean : remove executáveis e objetos"

#############################
# Makefile para compilar em paralelo e sequencial
#############################

# Compiladores e flags
CC_MPI      = mpicc
CC_SEQ      = gcc
CFLAGS_SEQ      = -std=c99 -Wall
CFLAGS_MPI      = -std=c99 -Wall -O3

# Nomes dos executáveis
PAR_TARGET  = lcs_mpi
SEQ_TARGET  = lcs_sequencial

# Arquivos-fonte
PAR_SRCS    = lcs_parallel.c
SEQ_SRCS    = lcs_sequencial.c

# Letras especiais
.PHONY: all paralelo sequencial run clean help

default: help

# Regra para compilação paralela (MPI)
paralelo: $(PAR_TARGET)

$(PAR_TARGET): $(PAR_SRCS)
	$(CC_MPI) $(CFLAGS_MPI) -o $@ $^

# Regra para compilação sequencial
sequencial: $(SEQ_TARGET)

$(SEQ_TARGET): $(SEQ_SRCS)
	$(CC_SEQ) $(CFLAGS_SEQ) -o $@ $^

# Para rodar em paralelo (MPI)
run: paralelo
	mpirun -np 4 ./$(PAR_TARGET) $(ARGS)

# Limpeza de artefatos
clean:
	rm -f $(PAR_TARGET) $(SEQ_TARGET) *.o generate

# Ajuda
help:
	@echo "Uso: make [paralelo|sequencial|run|clean]"
	@echo "  paralelo    : compila com MPI ($(PAR_TARGET))"
	@echo "  sequencial  : compila apenas o lcs_sequencial ($(SEQ_TARGET))"
	@echo "  run         : executa $(PAR_TARGET) com mpirun (use 'make run ARGS=\"fileA fileB\"')"
	@echo "  clean       : remove executáveis e objetos"

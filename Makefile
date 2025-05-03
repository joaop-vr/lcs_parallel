CC = gcc
CFLAGS = -Wall -O2

all: generate_input

generate_input: generate_input.c
	$(CC) $(CFLAGS) -o generate_input generate_input.c

clean:
	rm -f generate_input *.in

CC=gcc-10
CFLAGS=-fopenmp

ballAlg: ballAlg.c
	$(CC) -o ballAlg ballAlg.c $(CFLAGS)
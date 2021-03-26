CC=gcc-10
CFLAGS=-fopenmp -O3

ballAlg: ballAlg.c
	$(CC) -o ballAlg ballAlg.c $(CFLAGS)

debug: ballAlg.c
	$(CC) -o ballAlg ballAlg.c $(CFLAGS) -DDEBUG
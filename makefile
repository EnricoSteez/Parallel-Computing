CC=gcc-10
CFLAGS=-fopenmp -O3

ballAlg: ballAlg.c
	$(CC) -o ballAlg ballAlg.c $(CFLAGS)

debug: ballAlg.c
	$(CC) -o ballAlg ballAlg.c $(CFLAGS) -DDEBUG

parallel: ballAlg-omp.c
	$(CC) -o ballAlg-omp ballAlg-omp.c $(CFLAGS)
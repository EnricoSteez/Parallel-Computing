CC=gcc-10
CFLAGS=-fopenmp -O3

sequential: ballAlg.c
	$(CC) -o ballAlg-seq ballAlg.c $(CFLAGS)

debug: ballAlg.c
	$(CC) -o ballAlg ballAlg.c $(CFLAGS) -DDEBUG

parallel: ballAlg-omp.c
	$(CC) -o ballAlg ballAlg-omp.c $(CFLAGS)
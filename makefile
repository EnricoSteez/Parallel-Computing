CC=gcc-10
CFLAGS=-fopenmp -O3

sequential: ballAlg.c
	$(CC) -o ballAlg-seq ballAlg.c $(CFLAGS)

parallel: ballAlg-omp.c
	$(CC) -o ballAlg ballAlg-omp.c $(CFLAGS)

parallelprofiller: ballAlg-omp.c
	kinst-ompp $(CC) -o ballAlg ballAlg-omp.c $(CFLAGS) 
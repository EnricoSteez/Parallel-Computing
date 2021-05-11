CC=gcc
CFLAGS=-fopenmp -O3 -lm

all: sequential parallel

sequential: ballAlg.c
	$(CC) -o ballAlg ballAlg.c $(CFLAGS)

parallel: ballAlg-omp.c
	$(CC) -o ballAlg-omp ballAlg-omp.c $(CFLAGS)

parallelprofiller: ballAlg-omp.c
	kinst-ompp $(CC) -o ballAlg ballAlg-omp.c $(CFLAGS)

ballquery: ballQuery.c
	$(CC) -o ballQuery ballQuery.c $(CFLAGS)

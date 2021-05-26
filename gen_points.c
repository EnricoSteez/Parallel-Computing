#include <stdio.h>
#include <stdlib.h>

#define RANGE 10

extern void print_point(double *, int);

double **create_array_pts(int n_dims, long np)
{
    double *_p_arr;
    double **p_arr;

    _p_arr = (double *) malloc(n_dims * np * sizeof(double));
    p_arr = (double **) malloc(np * sizeof(double *));
    if((_p_arr == NULL) || (p_arr == NULL)){
        printf("Error allocating array of points, exiting.\n");
        exit(4);
    }

    for(long i = 0; i < np; i++)
        p_arr[i] = &_p_arr[i * n_dims];

    return p_arr;
}


double **get_points(int argc, char *argv[], int *n_dims, long *np, int me, int nprocs)
{
    double **pt_arr;
    unsigned seed;
    long i, mynp;
    int j;

    if(me!=(nprocs-1)){
        mynp = atol(argv[2])/nprocs;
    } else {
        mynp = atol(argv[2]) - (atol(argv[2])/(nprocs))*(nprocs-1);
    }


    if(argc != 4){
        printf("Usage: %s <n_dims> <n_points> <seed>\n", argv[0]);
        exit(1);
    }

    *n_dims = atoi(argv[1]);
    if(*n_dims < 2){
        printf("Illegal number of dimensions (%d), must be above 1.\n", *n_dims);
        exit(2);
    }

    *np = atol(argv[2]);
    if(*np < 1){
        printf("Illegal number of points (%ld), must be above 0.\n", *np);
        exit(3);
    }

    seed = atoi(argv[3]);
    srandom(seed);

    pt_arr = (double **) create_array_pts(*n_dims, mynp);

    long myindex=0;
    double pointDim;

    for(i = 0; i < *np; i++)
        for(j = 0; j < *n_dims; j++){
            pointDim = RANGE * ((double) random()) / RAND_MAX;
            if(i>=(atol(argv[2])/nprocs)*me && i<((atol(argv[2])/nprocs)*me + mynp)))
                pt_arr[myindex++][j] = pointDim;
        }
                    


#ifdef DEBUG
    for(i = 0; i < *np; i++)
        print_point(pt_arr[i], *n_dims);
#endif

    return pt_arr;
}

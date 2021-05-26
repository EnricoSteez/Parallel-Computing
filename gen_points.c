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



void get_n_points(long array[], long indexes[], long npoints, int n, int first_proc, long first_point) {
    if(n == 1) {
        array[first_proc] = npoints;
        if(first_proc != 0) indexes[first_proc] = indexes[first_proc-1] + array[first_proc-1];
        else indexes[first_proc] = 0;
        return;
    }
    else {
        get_n_points(array, indexes, npoints/2, n/2, first_proc, first_point);
        get_n_points(array, indexes, npoints-npoints/2, n-n/2, first_proc+n/2, first_point+npoints/2);       
    }
}



double **get_points(int argc, char *argv[], int *n_dims, long *np, int me, int nprocs)
{
    double **pt_arr;
    unsigned seed;
    long i;
    int j;


    if(argc != 4){
        printf("Usage: %s <n_dims> <n_points> <seed>\n", argv[0]);
        exit(1);
    }

    *n_dims = atoi(argv[1]);
    if(*n_dims < 2){
        printf("Illegal number of dimensions (%d), must be above 1.\n", *n_dims);
        exit(2);
    }
    long number_of_points = atol(argv[2]);
    if(number_of_points < 1){
        printf("Illegal number of points (%ld), must be above 0.\n", *np);
        exit(3);
    }

    seed = atoi(argv[3]);
    srandom(seed);

    long array[nprocs];
    long indexes[nprocs];

    get_n_points(array, indexes, number_of_points, nprocs, 0, 0);

    pt_arr = (double **) create_array_pts(*n_dims, array[me]);

    *np = array[me];

    long myindex=0;
    double pointDim;

    for(i = 0; i < *np; i++)
        for(j = 0; j < *n_dims; j++){
            pointDim = RANGE * ((double) random()) / RAND_MAX;
            if(i >= indexes[me] && i < indexes[me] + array[me])
                pt_arr[myindex++][j] = pointDim;
        }
                    


#ifdef DEBUG
    for(i = 0; i < *np; i++)
        print_point(pt_arr[i], *n_dims);
#endif

    return pt_arr;
}

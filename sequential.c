#include <stdio.h>
#include <math.h>
#include "gen_points.c"

#define DIMENSIONS 3
#define NP 5000

typedef struct node {
    float coordinates[DIMENSIONS];
    float radius;
    struct node* left;
    struct node* right;
} node;

struct IndexCoord {
    long idx;
    double coord;
};



int find_extremes(node * a, node * b){
    
    //TODO

    return 0;
}

void print_point(double* point, int dim) {
    int j;
    printf("(");
    for(j = 0; j < dim; j++) {
        if(j != dim-1) printf("%lf, ", point[j]);

        else printf("%lf", point[j]);
    }
    printf(")\n");
}

static int compare (const void * a, const void * b)
{
    return ( (*(struct IndexCoord *) a).coord - (*(struct IndexCoord*) b).coord );
}

struct IndexCoord* project_on_dimension_and_sort(long n_points, int n_dim, int* indexes, double proj_table[n_points][n_dim]){
    struct IndexCoord* oneDim_projection = 
        (struct IndexCoord*) malloc(sizeof(struct IndexCoord)*n_points);

    int counter = 1;
    double coord = 0;

    
    for(int i = 0; i < n_dim; i++){
        for(long p = 0; p < n_points; p++){
            if((coord = proj_table[p][i]) == proj_table[p+1][i] && p < n_points - 1)
                    counter++;
            oneDim_projection[p].coord = coord;
            oneDim_projection[p].idx = p;
        }
        if(counter != n_points)
            break;
        counter = 1;
    }

    // #ifdef DEBUG 
    //     printf("\nPoints projected on one dimension BEFORE QSORT:\n");
    //     for (int n=0; n<n_points; n++){
    //     printf ("%f (%ld)\n",oneDim_projection[n].coord, oneDim_projection[n].idx);
    //     }
    //     printf("\n");
    // #endif

    qsort(oneDim_projection, n_points, sizeof(struct IndexCoord), compare); 

    #ifdef DEBUG 
        printf("\nPoints projected on one dimension AFTER QSORT:\n");
        for (int n=0; n<n_points; n++){
        printf ("%f (%ld)\n",oneDim_projection[n].coord, oneDim_projection[n].idx);
        }
        printf("\n");
    #endif

    return oneDim_projection;
}

void find_median(long n_points, int n_dim, struct IndexCoord oneDim_projection[n_points], 
    double proj_table[n_points][n_dim]){
    
    //return value for the median point
    double median_point[n_dim];
    
    if((n_points % 2) != 0){
        long idx_median = oneDim_projection[(n_points - 1) / 2].idx;
        for(int i = 0; i < n_dim; i++){
            median_point[i] = proj_table[idx_median][i];
        }
        // Divide the arrays
        long L[(n_points-1)/2];
        long R[(n_points+1)/2];

        for(int p = 0; p < ((n_points-1)/2); p++){
            L[p] = oneDim_projection[p].idx;
            R[p] = oneDim_projection[n_points-1-p].idx;
        }
        R[(n_points-1)/2] = oneDim_projection[(n_points-1)/2].idx;
    
    }
    else{
        long idx_median_1 = oneDim_projection[n_points / 2].idx;
        long idx_median_2 = oneDim_projection[(n_points / 2) - 1].idx;
        for(int i = 0; i < n_dim; i++){
            median_point[i] = (proj_table[idx_median_1][i] + proj_table[idx_median_2][i]) / 2;
        }
        long L[n_points/2];
        long R[n_points/2];
        for(int p = 0; p < (n_points / 2); p++){
            L[p] = oneDim_projection[p].idx;
            R[p] = oneDim_projection[n_points-1-p].idx;
        }
    }
    #ifdef DEBUG 
        printf("\nMEDIAN POINT (full coordinate):\n");

        //don't judge this code please
        for(int i = 0; i < n_dim; i++){
            if(i==0){
                printf("(%f", median_point[i]);
            } else if(i==--n_dim){
                printf(",%f)\n", median_point[i]);
            } else {
                printf(",%f", median_point[i]);
            }
        }
    #endif

    // We need to divide the points still and return as an input parameter
}

void orthogonal_projection(double **points_table, long n_points, int n_dim, int* indexes, double proj_table[n_points][n_dim]){
    double delta = 0, gamma = 0, phi = 0;
    double a, b, point;

    for(long p = 0; p < n_points; p++){
        for(int dim = 0; dim < n_dim; dim++){
            a = points_table[indexes[0]][dim];
            b = points_table[indexes[1]][dim];
            point = points_table[p][dim];
            delta += (point - a) * (b - a);
            gamma += pow((b - a), 2);
        }
        phi = delta / gamma;
        for(int dim = 0; dim < n_dim; dim++){
            a = points_table[indexes[0]][dim];
            b = points_table[indexes[1]][dim];
            proj_table[p][dim] = phi * (b - a) + a;
        }
        delta = 0;
        gamma = 0;
    }
}

double distance_between_points(double* point1, double* point2, int dim) {
    int i, sum = 0;
    for (i = 0; i < dim; i++) {
        sum += pow(point2[i] - point1[i], 2);
    }
    return sqrt(sum);
}

double find_radius(double *median_point, int* indexes, int n_points, int n_dim, double proj_table[n_points][n_dim]){
    
    double highest = 0, dist;
    for(int i=0; i<2; i++){
        dist = distance_between_points(median_point, proj_table[indexes[i]], n_dim);
        if(dist>highest)
            highest = dist;
    }
    return highest;
}

//returns the indices of the 2 furthest points
void furthest_points(double** points, int np, int dim, int furthest[]) {
    double* first_point = points[0];
    int i, max_index;
    double aux, max;

    for (i = 1; i < np; i++) {
        if((aux = distance_between_points(first_point, points[i], dim)) > max) {
            max = aux;
            max_index = i;
        }
    }

    int a = max_index;
    max = 0;
    for (i = 0; i < np; i++) {
        if((aux = distance_between_points(points[a], points[i], dim)) > max) {
            max = aux;
            max_index = i;
        }
    }
    furthest[0] = a;
    furthest[1] = max_index;
}

int main(int argc, char **argv){

    //1. get input sample points (use the function from the guide)
    int dim;
    long np;

    //double ** arr = create_array_pts(DIMENSIONS,NP);
    double ** points = get_points(argc, argv, &dim, &np);
    
    int j;

    //2. compute points a and b, furthest apart in the current set;
    
    int furthest[2];

    furthest_points(points, np, dim, furthest);
    int a = furthest[0];
    int b = furthest[1];

    printf("a ");

    print_point(points[a], dim);

    printf("b ");

    print_point(points[b], dim);

    //3. perform the orthogonal projection of all points onto line ab;

    double proj_table[np][dim]; 
    orthogonal_projection(points, np, dim, furthest, proj_table);
    for(int i = 0; i < np; i++){
        printf("Point %d ", (i));
        for(int k = 0; k < dim; k++){
            printf("Dim %d, value = %f ", k, points[i][k]);
        }
        printf("\n");
    }

    //proj_table are the points on line ab and this function projects them onto one single dimension
    struct IndexCoord* oneDim_projection = project_on_dimension_and_sort(np, dim, furthest, proj_table);
    
    //oneDim_projection are the points projected on one dimension, ready to find the median point
    find_median(np, dim, oneDim_projection, proj_table);



    //4. compute the center, defined as the median point over all projections;

free(oneDim_projection);
    //5. create two sets of points, L and R, deﬁned as the points whose projection lies to one side or the other of the center;

    return 0;
}

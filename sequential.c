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


static int compare (const void * a, const void * b)
{
    return ( (*(struct IndexCoord*)b).coord > (*(struct IndexCoord*)a).coord );
}

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

void find_median(long n_points, int n_dim, int* indexes, double proj_table[n_points][n_dim]){

    struct IndexCoord median_coords[n_points];
    int counter = 1;
    double coord = 0;
    double median_point[n_dim];

    for(int i = 0; i < n_dim; i++){
        for(long p = 0; p < n_points; p++){
            if((coord = proj_table[p][i]) == proj_table[p+1][i] && p < n_points - 1)
                    counter++;
            median_coords[p].coord = coord;
            median_coords[p].idx = p;
        }
        if(counter != n_points)
            break;
        counter = 1;
    }
    qsort(median_coords, n_points, sizeof(struct IndexCoord), compare); 
    
    //for (int n=0; n<n_points; n++){
    //    printf ("%f (%ld)\n",median_coords[n].coord, median_coords[n].idx);
    //}
    //printf("\n");
    
    if((n_points % 2) != 0){
        long idx_median = median_coords[(n_points - 1) / 2].idx;
        for(int i = 0; i < n_dim; i++){
            median_point[i] = proj_table[idx_median][i];
        }
    }
    else{
        long idx_median_1 = median_coords[n_points / 2].idx;
        long idx_median_2 = median_coords[(n_points / 2) - 1].idx;
        for(int i = 0; i < n_dim; i++){
            median_point[i] = (proj_table[idx_median_1][i] + proj_table[idx_median_2][i]) / 2;
        }
    }

    //printf("MEDIAN POINT\n");
    //for(int i = 0; i < n_dim; i++){
    //   printf("%f\n", median_point[i]);
    //}

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
            printf("Dim %d, value = %f ", k, proj_table[i][k]);
        }
        printf("\n");
    }

    find_median(np, dim, furthest, proj_table);
    //4. compute the center, defined as the median point over all projections;

    //5. create two sets of points, L and R, deﬁned as the points whose projection lies to one side or the other of the center;

    

    return 0;
}

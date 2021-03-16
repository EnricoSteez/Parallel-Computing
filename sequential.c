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

//Usually success = 0, failure = 1, we can change this, it's not relevant 

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


float orthogonal_projection(node n){
    float distance = 0;

    return distance;
}


double distance_between_points(double* point1, double* point2, int dim) {
    int i, sum = 0;
    for (i = 0; i < dim; i++) {
        sum += pow(point2[i] - point1[i], 2);
    }
    return sqrt(sum);
}


//returns the indeices of the 2 furthest points
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

    

    //4. compute the center, defined as the median point over all projections;

    //5. create two sets of points, L and R, deﬁned as the points whose projection lies to one side or the other of the center;

    

    return 0;
}
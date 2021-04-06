#include <stdio.h>
#include <math.h>
#include "gen_points.c"
#include <omp.h>
#include <string.h>

#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_RESET   "\x1b[0m"

typedef struct node {
    long id;
    double* coordinates;
    float radius;
    struct node* left;
    struct node* right;
} node;

struct ProjectedPoint {
    long idx;
    double *projectedCoords;
};

double ** points;
int dim;
long np;
long n_nodes;

void print_point(double* point, int dim) {
    int j;
    printf("(");
    for(j = 0; j < dim; j++) {
        if(j != dim-1) printf("%.1lf, ", point[j]);

        else printf("%.1lf", point[j]);
    }
    printf(")\n");
    fflush(stdin);
}

static int compare (const void * a, const void * b)
{
    double diff = ( (*(struct ProjectedPoint *) a).projectedCoords[0] - (*(struct ProjectedPoint*) b).projectedCoords[0] );
    if(diff < 0) {
        return -1;
    }
    else if(diff > 0) {
        return 1;
    }
    return 0;
}


double* find_center_and_rearrange_set(long current_set_size, long* current_set, struct ProjectedPoint* proj_table){

    long lindex=0;
    double* pointProjection;
    long rindex = current_set_size/2;
    double* center = (double *) malloc(dim*sizeof(double));

    if((current_set_size % 2) != 0){//ODD SET
        long idx_median = current_set_size / 2;
        memcpy(center, proj_table[idx_median].projectedCoords,dim*sizeof(double));
    }
    else{ //EVEN SET

        long idx_median_1 = current_set_size / 2 -1;
        long idx_median_2 = current_set_size / 2;

        //center to be returned
        // #pragma omp parallel for
        for(int i = 0; i < dim; i++){
            center[i] = (proj_table[idx_median_1].projectedCoords[i] + proj_table[idx_median_2].projectedCoords[i]) / 2;
        }
    }

    // #pragma omp parallel for
    for(long i=0; i < current_set_size; i++){
        *(current_set+i) = proj_table[i].idx;
    }
    
    return center;
}

//TODO PARALLELISE THIS MESS
void orthogonal_projection(long current_set_size, long* current_set, long* furthest_points, struct ProjectedPoint* proj_table){
    double delta = 0, gamma = 0, phi = 0;
    double a, b, point;
    int d;

    for(long p = 0; p < current_set_size; p++){
        for(d = 0; d < dim; d++){
            a = points[furthest_points[0]][d];
            b = points[furthest_points[1]][d];
            point = points[current_set[p]][d];
            delta += (point - a) * (b - a);
            gamma += pow((b - a), 2);
        }
        phi = delta / gamma;

        proj_table[p].projectedCoords = (double*) malloc(dim * sizeof(double));

        for(d = 0; d < dim; d++){
            a = points[furthest_points[0]][d];
            b = points[furthest_points[1]][d];
            proj_table[p].projectedCoords[d] = phi * (b - a) + a;
        }

        proj_table[p].idx = current_set[p];

        delta = 0;
        gamma = 0;
    }
}

double distance_between_points(double* point1, double* point2) {
    int i;
    double sum = 0;
    //#pragma omp parallel for reduction(+:sum)
    for (i = 0; i < dim; i++) {
        sum += pow(point2[i] - point1[i], 2);
    }
    return sqrt(sum*1.0);
}

double find_radius(double *center, long* current_set, long current_set_size){
    
    // //TODO ASK THE PROFESSOR
    // //Does this imply that all the computation is done after the barrier?
    // //Maybe using the critical region with nowait would do the comparison of the last thread that finishes after the loop
    // //since the other threads will do it as soon as they finish!
    // //so it would be probably more efficient
    double highest = 0, dist;
    double globalHighest=0;

    #pragma omp parallel 
    {
        #pragma omp for reduction(max:highest)
        for(long i=0; i<current_set_size; i++){
            dist = distance_between_points(center, points[current_set[i]]);
            if(dist>highest)
                highest = dist;
        }
        //fprintf(stderr,"Thread %d finished\n",omp_get_thread_num());
    }

    


    return highest;

    // double localMax=0, dist, globalMax=0;
    // #pragma omp parallel private(localMax, dist) shared(globalMax)
    // {
    //     #pragma omp for nowait
    //     for(long i=0; i<current_set_size; i++){
    //         dist = distance_between_points(center, points[current_set[i]]);
    //         if(dist>localMax)
    //         localMax = dist;
    //     }

    //     #pragma omp critical
    //     {
    //         globalMax = localMax > globalMax ? localMax : globalMax;
    //     }
    // }
    
    // return globalMax;


    // double highest = 0, dist;
    // double globalHighest=0;

    // for(long i=0; i<current_set_size; i++){
    //     dist = distance_between_points(center, points[current_set[i]]);
    //     if(dist>highest)
    //         highest = dist;
    // }

    // return highest;

}

long furthest_point_from_point(double* point, long* current_set, long current_set_size) {
    // long i, local_max_index=0, global_max_index=0;
    // double aux, local_max=0, global_max=0;
    
    // #pragma omp parallel private(aux, local_max, local_max_index) shared(global_max, global_max_index)
    // {
    //     #pragma omp for nowait
    //     for (i = 0; i < current_set_size; i++) {
    //         if((aux = distance_between_points(point, points[current_set[i]])) > local_max) {
    //             local_max = aux;
    //             local_max_index = current_set[i];
    //         }
    //     }

    //     #pragma omp critical
    //     {
    //         if(local_max > global_max){
    //             global_max = local_max;
    //             global_max_index = local_max_index;
    //         }
    //     }
    // }

    // return global_max_index;

    long i, local_max_index=0;
    double aux;
    double local_max=0;
    
    for (i = 0; i < current_set_size; i++) {

        if((aux = distance_between_points(point, points[current_set[i]])) > local_max) {
            local_max = aux;
            local_max_index = current_set[i];
        }
    }
    return local_max_index;
}

//sets furthest[] to the indices of the 2 furthest points
void furthest_points(long furthest[2], long* current_set, long current_set_size) {
    double* first_point = points[current_set[0]];
    long a = furthest_point_from_point(first_point, current_set, current_set_size);

    furthest[0] = a;
    furthest[1] =furthest_point_from_point(points[a], current_set, current_set_size);
}

struct node* build_tree(long node_index, long* current_set, long current_set_size) {
    
    if(current_set_size == 1) {
        //stop recursion
        struct node* res = (struct node*)malloc(sizeof(struct node));
        res->id = node_index;
        res->coordinates = points[current_set[0]];
        res->left = NULL;
        res->right = NULL;
        res->radius = 0;
        n_nodes++;

        return res;
    }    

    long a;
    long b;

    struct ProjectedPoint* proj_table;
    proj_table = (struct ProjectedPoint*) malloc (current_set_size* sizeof(struct ProjectedPoint));
    
    if(current_set_size > 2) {
        //compute points a and b, furthest apart in the current set;
        long furthest[2];
        // fprintf(stderr, "Calculate furthest points\n"); 

        furthest_points(furthest, current_set, current_set_size);
        a = furthest[0];
        b = furthest[1];
        //perform the orthogonal projection of all points onto line ab;
        // fprintf(stderr, "Make orth proj\n"); 

        orthogonal_projection(current_set_size, current_set, furthest, proj_table);
    }
    else {
        //if there are only 2 points, no need to make orthogonal projection
        a = current_set[0];
        b = current_set[1];
        proj_table[0].idx = a;
        proj_table[1].idx = b;

        proj_table[0].projectedCoords = (double*) malloc(dim * sizeof(double));
        proj_table[1].projectedCoords = (double*) malloc(dim * sizeof(double));

        for(int i = 0; i < dim; i++) {
            proj_table[0].projectedCoords[i] = points[a][i];
            proj_table[1].projectedCoords[i] = points[b][i];
        }
    }
    
    //proj_table are the points on line ab and this function projects them onto one single dimension
    //a and b are the indexes of the furthest points of this current_set


    qsort(proj_table, current_set_size, sizeof(struct ProjectedPoint), compare);

    //compute the center, defined as the median point over all projections;
    // fprintf(stderr, "find_center and update current_set"); 

    double* center = find_center_and_rearrange_set(current_set_size, current_set, proj_table);

    free(proj_table);

    long nextLeftSize = current_set_size/2;
    long nextRightSize = current_set_size%2==0 ? current_set_size/2 :current_set_size/2+1;

    struct node* res = (struct node*)malloc(sizeof(struct node));
    res->coordinates = center;
    res->radius = find_radius(center, current_set, current_set_size);
    res->id = node_index;
    res->left = build_tree(node_index + 1, current_set, nextLeftSize);
    res->right = build_tree(node_index +nextLeftSize* 2 , current_set+current_set_size/2, nextRightSize);
    n_nodes++;

    return res;
}

void dump_tree(struct node *node){
    printf("%ld ", node->id);

    if(node->left != NULL)
        printf("%ld %d ", node->left->id, node->right->id);
    else
        printf("-1 -1 ");

    printf("%lf ",node->radius);

    for(int i=0; i<dim; i++){
        printf("%lf ",node->coordinates[i]);
    }
    printf("\n"); fflush(stdin);

    if(node->left != NULL){
        dump_tree(node->left);
        dump_tree(node->right);
    }
}


int main(int argc, char **argv){
    double exec_time;

    // fprintf(stderr, "%d\n", omp_get_num_procs()); 
    exec_time = - omp_get_wtime();
    //get input sample points (use the function from the guide)

    points = get_points(argc, argv, &dim, &np);

    long* current_set = (long*) malloc(np * sizeof(long));

    // #pragma omp parallel for
    for(int j = 0; j < np; j++) {
        current_set[j] = j;
    }

    n_nodes = 0;

    struct node* tree = build_tree(0, current_set, np);

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1lf\n", exec_time); 

    printf("%d %d\n",dim,n_nodes);
    dump_tree(tree);

    return 0;
}

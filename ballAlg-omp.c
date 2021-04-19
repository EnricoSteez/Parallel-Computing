#include <stdio.h>
#include <math.h>
#include "gen_points.c"
#include <omp.h>
#include <string.h>

#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define NUM_THREADS 2

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
int nthreads;

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


double* find_center(long current_set_size, struct ProjectedPoint* proj_table){

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
        for(int i = 0; i < dim; i++){
            center[i] = (proj_table[idx_median_1].projectedCoords[i] + proj_table[idx_median_2].projectedCoords[i]) / 2;
        }
    }
    
    return center;
}

void rearrange_set(long current_set_size, long* current_set, struct ProjectedPoint* proj_table, long rec_level){

    //#pragma omp parallel for num_threads(nthreads/(pow(2,rec_lev)) if(rec_level<2)
    for(long i=0; i < current_set_size; i++){
        *(current_set+i) = proj_table[i].idx;
    }
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

    for (i = 0; i < dim; i++) {
        sum += pow(point2[i] - point1[i], 2);
    }
    return sqrt(sum*1.0);
}

double find_radius(double *center, long* current_set, long current_set_size, long rec_level, int n){

    double highest = 0, dist;
    double globalHighest=0;

    if(n < 0) {
        fprintf(stderr,"Find radius in set of %d with %d threads\n",current_set_size,n);
        #pragma omp parallel for reduction(max:highest) num_threads(n)
        for(long i=0; i<current_set_size; i++){
            dist = distance_between_points(center, points[current_set[i]]);
            if(dist>highest)
                highest = dist;
        }
    }
    else {
        for(long i=0; i<current_set_size; i++){
            dist = distance_between_points(center, points[current_set[i]]);
            if(dist>highest)
                highest = dist;
        }
    }

    return highest;

}

long furthest_point_from_point(double* point, long* current_set, long current_set_size,int n) {

    if(n<0){
        long i, local_max_index=0, global_max_index=0;
        double aux, local_max=0, global_max=0;
        
        #pragma omp parallel firstprivate(aux, local_max, local_max_index, i) shared(global_max, global_max_index) num_threads(n)
        {
            #pragma omp for nowait
            for (i = 0; i < current_set_size; i++) {
                if((aux = distance_between_points(point, points[current_set[i]])) > local_max) {
                    local_max = aux;
                    local_max_index = current_set[i];
                }
            }
            #pragma omp critical
            {
                if(local_max > global_max){
                    global_max = local_max;
                    global_max_index = local_max_index;
                    
                }
            }
        }

        return global_max_index;

    } else {
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
    
}

//sets furthest[] to the indices of the 2 furthest points
void furthest_points(long furthest[2], long* current_set, long current_set_size, int n) {
    double* first_point = points[current_set[0]];
    long a = furthest_point_from_point(first_point, current_set, current_set_size, n);

    furthest[0] = a;
    furthest[1] =furthest_point_from_point(points[a], current_set, current_set_size, n);
}

struct node* build_tree(long node_index, long* current_set, long current_set_size, long rec_level) {
    int n = nthreads/(pow(2,rec_level));

    double exec_findradius;

    if(rec_level < 1) {
        exec_findradius = -omp_get_wtime();
    }

    #pragma omp atomic
    n_nodes++;

    if(current_set_size == 1) {
        //stop recursion
        struct node* res = (struct node*)malloc(sizeof(struct node));
        res->id = node_index;
        res->coordinates = points[current_set[0]];
        res->left = NULL;
        res->right = NULL;
        res->radius = 0;

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

        furthest_points(furthest, current_set, current_set_size, n);
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

    qsort(proj_table, current_set_size, sizeof(struct ProjectedPoint), compare);

    double* center;

    center = find_center(current_set_size, proj_table);
    rearrange_set(current_set_size, current_set, proj_table, rec_level);

    free(proj_table);
    
    struct node* res = (struct node*)malloc(sizeof(struct node));
    res->coordinates = center;

    
    res->radius = find_radius(center, current_set, current_set_size,rec_level, n);


    res->id = node_index;

    long nextLeftSize = current_set_size/2;
    long nextRightSize = current_set_size%2==0 ? current_set_size/2 :current_set_size/2+1;


    // if(n > 1) {
    //     fprintf(stderr, "build_tree time in level %ld: %lf\n", rec_level, exec_findradius + omp_get_wtime());
    // }
    
    #pragma omp task if(n>1) 
    res->left = build_tree(node_index + 1, current_set, nextLeftSize, rec_level+1);
    #pragma omp task if(n>1) 
    res->right = build_tree(node_index +nextLeftSize* 2 , current_set+current_set_size/2, nextRightSize, rec_level+1);
    
    return res;
}

void dump_tree(struct node *node){
    printf("%ld ", node->id);

    if(node->left != NULL)
        printf("%ld %ld ", node->left->id, node->right->id);
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
    struct node* tree;

    if(argc == 5) {
        nthreads = atoi(argv[4]);
        argc = argc-1;
    }
    else {
        nthreads = NUM_THREADS;
    }
    omp_set_nested(1);
    omp_set_num_threads(nthreads);
    fprintf(stderr, "number of threads: %d\n", nthreads); 
    exec_time = - omp_get_wtime();
    //get input sample points (use the function from the guide)

    points = get_points(argc, argv, &dim, &np);

    long* current_set = (long*) malloc(np * sizeof(long));

    // #pragma omp parallel for
    for(int j = 0; j < np; j++) {
        current_set[j] = j;
    }

    n_nodes = 0;

    #pragma omp parallel
    #pragma omp single
    {
        tree = build_tree(0, current_set, np, 0);
    }

    

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1lf\n", exec_time); 

    printf("%d %ld\n",dim,n_nodes);
    dump_tree(tree);

    return 0;
}

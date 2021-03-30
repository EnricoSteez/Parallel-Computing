#include <stdio.h>
#include <math.h>
#include "gen_points.c"
#include <omp.h>

#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_RESET   "\x1b[0m"

typedef struct node {
    long id;
    double* coordinates;
    float radius;
    struct node* left;
    struct node* right;
} node;

struct IndexCoord {
    long idx;
    double coord;
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
}

static int compare (const void * a, const void * b)
{
    double diff = ( (*(struct IndexCoord *) a).coord - (*(struct IndexCoord*) b).coord );
    if(diff < 0) {
        return -1;
    }
    else if(diff > 0) {
        return 1;
    }
    return 0;
}

struct IndexCoord* project_on_dimension_and_sort(long current_set_size, long* current_set, double proj_table[][dim], long a, long b){
    struct IndexCoord* oneDim_projection = 
        (struct IndexCoord*) malloc(sizeof(struct IndexCoord)*current_set_size);

    double coord = 0;
    
    for(int i = 0; i < dim; i++){
        #pragma omp parallel for
        for(long p = 0; p < current_set_size; p++){
            oneDim_projection[p].coord = proj_table[p][i];
            oneDim_projection[p].idx = p;
        }

        if(current_set_size > 1 && proj_table[a][i] != proj_table[b][i])
            break;
    }

    qsort(oneDim_projection, current_set_size, sizeof(struct IndexCoord), compare); 

    return oneDim_projection;
}

double* find_median(long current_set_size, long* current_set, struct IndexCoord oneDim_projection[], double proj_table[][dim], long* R, long* L){
    
    //return value for the median point
    double* median_point = (double*)malloc(dim*sizeof(double));

    if((current_set_size % 2) != 0){//ODD SET
        long idx_median = oneDim_projection[(current_set_size - 1) / 2].idx;

            #pragma omp parallel for
            for(int i = 0; i < dim; i++){
                median_point[i] = proj_table[idx_median][i];
            }

            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    #pragma omp parallel for              
                    for(long p = 0; p < ((current_set_size-1)/2); p++){
                        *(L + p) = current_set[oneDim_projection[p].idx];
                    }
                }
                #pragma omp section
                {
                    #pragma omp parallel for
                    for(long p = 0; p < ((current_set_size-1)/2); p++){
                        *(R + p) = current_set[oneDim_projection[current_set_size-1-p].idx];
                    }
                }
            }
            //middle point belongs to Right set
            *(R + (current_set_size-1)/2) = current_set[oneDim_projection[(current_set_size-1)/2].idx];

    }
    else{ //EVEN SET

        long idx_median_1 = oneDim_projection[current_set_size / 2].idx;
        long idx_median_2 = oneDim_projection[(current_set_size / 2) - 1].idx;

        #pragma omp parallel for
        for(int i = 0; i < dim; i++){
            median_point[i] = (proj_table[idx_median_1][i] + proj_table[idx_median_2][i]) / 2;
        }

        #pragma omp parallel sections
        {
            #pragma omp section
            {
                #pragma omp parallel for
                for(long p = 0; p < (current_set_size / 2); p++){
                    *(L + p) = current_set[oneDim_projection[p].idx];
                }
            }
            #pragma omp section
            {
                #pragma omp parallel for
                for(long p = 0; p < (current_set_size / 2); p++){
                    *(R + p) = current_set[oneDim_projection[current_set_size-1-p].idx];
                }
            }
        }
    }

    return median_point;
    // We need to divide the points still and return as an input parameter
}

void orthogonal_projection(long current_set_size, long* current_set, long* furthest_points, double proj_table[][dim]){
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

        for(d = 0; d < dim; d++){
            a = points[furthest_points[0]][d];
            b = points[furthest_points[1]][d];
            proj_table[p][d] = phi * (b - a) + a;
            
        }

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

double find_radius(double *median_point, long* current_set, long current_set_size){
    
    double highest = 0, dist;
    for(long i=0; i<current_set_size; i++){
        dist = distance_between_points(median_point, points[current_set[i]]);
        if(dist>highest)
            highest = dist;
    }
    return highest;
}

long furthest_point_from_point(double* point, long* current_set, long current_set_size) {
    long i, max_index;
    double aux, max = 0;

    for (i = 0; i < current_set_size; i++) {
        if((aux = distance_between_points(point, points[current_set[i]])) > max) {
            max = aux;
            max_index = current_set[i];
        }
    }

    return max_index;
}

//returns the indices of the 2 furthest points
void furthest_points(long furthest[], long* current_set, long current_set_size) {
    double* first_point = points[current_set[0]];
    long a = furthest_point_from_point(first_point, current_set, current_set_size);

    furthest[0] = a;
    furthest[1] =furthest_point_from_point(points[a], current_set, current_set_size);
}

struct node* build_tree(long node_index, long* current_set, long current_set_size) {
    /**
     * 
     * NOTE:
     * current_set: array of size current_set_size containing the indices of the original array (points) that we are currently working on
     * proj_table: array of size current_set_size containing the corresponding projection to each point in point[current_set[i]]
     * oneDim_projection: array of size current_set_size containing the corresponding one dimension point to each point in proj_table
     *
     * */
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

    double proj_table[current_set_size][dim];
    long a;
    long b;
    if(current_set_size > 2) {
        //compute points a and b, furthest apart in the current set;
        long furthest[2];
        furthest_points(furthest, current_set, current_set_size);
        a = furthest[0];
        b = furthest[1];
        //perform the orthogonal projection of all points onto line ab;
        orthogonal_projection(current_set_size, current_set, furthest, proj_table);
    }
    else {
        //if there are only 2 points, no need to make orthogonal projection
        for(int i = 0; i < dim; i++) {
            proj_table[0][i] = points[current_set[0]][i];
            proj_table[1][i] = points[current_set[1]][i];
        }
        
        a = current_set[0];
        b = current_set[1];
    }

    // fprintf(stderr, "BEFORE PROJ ONE DIMENSION AND SORT\n"); 
    //proj_table are the points on line ab and this function projects them onto one single dimension
    //a and b are the indexes of the furthest points of this current_set
    struct IndexCoord* oneDim_projection = project_on_dimension_and_sort(current_set_size, current_set, proj_table, a, b);

    // fprintf(stderr, "AFTER PROJ ONE DIMENSION AND SORT\n");
    long* L;
    long* R;
    long L_size, R_size;
    if((current_set_size % 2) != 0){
        L_size = (current_set_size-1)/2;
        R_size = (current_set_size+1)/2;
    }
    else {
        L_size = current_set_size/2;
        R_size = current_set_size/2;
    }
    L = (long*)malloc(L_size * sizeof(long));
    R = (long*)malloc(R_size * sizeof(long));

    // fprintf(stderr, "AFTER CREATION OF L AND R\n");

    //compute the center, defined as the median point over all projections;
    //oneDim_projection are the points projected on one dimension, ready to find the median point
    double* median = find_median(current_set_size, current_set, oneDim_projection, proj_table, R, L);
    // fprintf(stderr, "AFTER FIND MEDIAN\n");

    free(oneDim_projection);
    //proj table should be dynamically allocated with malloc in order to free

    struct node* res = (struct node*)malloc(sizeof(struct node));
    res->coordinates = median;
    res->radius = find_radius(median, current_set, current_set_size);
    res->id = node_index;
    res->left = build_tree(node_index + 1, L, L_size);
    res->right = build_tree(node_index + 1 + (L_size*2-1), R, R_size);
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
    printf("\n");

    if(node->left != NULL){
        dump_tree(node->left);
        dump_tree(node->right);
    }
}


int main(int argc, char **argv){
    double exec_time;
    exec_time = - omp_get_wtime();
    //get input sample points (use the function from the guide)

    points = get_points(argc, argv, &dim, &np);

    long* current_set = (long*) malloc(np * sizeof(long));

    #pragma omp parallel for
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

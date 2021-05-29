#include <stdio.h>
#include <math.h>
#include "gen_points.c"
#include <omp.h>
#include <string.h>
#include <mpi.h>

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

void rearrange_set(long current_set_size, long* current_set, struct ProjectedPoint* proj_table, long rec_level, int n){

    if(n>1){
        #pragma omp parallel for num_threads(n)
        for(long i=0; i < current_set_size; i++){
            *(current_set+i) = proj_table[i].idx;
        }
    }
    else{
        for(long i=0; i < current_set_size; i++){
            *(current_set+i) = proj_table[i].idx;
        }
    }
}

double* orthogonal_projection(long current_set_size, long* current_set, long* furthest_points, double* proj_table, int n){
    double delta = 0, gamma = 0, phi = 0;
    double a, b, point;
    int d;

    double* center = (double*)malloc(dim*sizeof(double));
    double* center_aux = (double*)malloc(dim*sizeof(double));


    if(n > 1) {

        #pragma omp parallel for firstprivate(delta, gamma, phi) private(a, b, point, d)
        for(long p = 0; p < current_set_size; p++){
            for(d = 0; d < dim; d++){
                a = points[furthest_points[0]][d];
                b = points[furthest_points[1]][d];
                point = points[current_set[p]][d];
                delta += (point - a) * (b - a);
                gamma += pow((b - a), 2);
            }
            phi = delta / gamma;


            a = points[furthest_points[0]][0];
            b = points[furthest_points[1]][0];



            if(current_set_size%2==1 && p == current_set_size/2) {
                for(d = 0; d < dim; d++){
                    a = points[furthest_points[0]][d];
                    b = points[furthest_points[1]][d];
                    center[d] = phi * (b - a) + a;
                }
            }
            else if (current_set_size%2==0 && p == (current_set_size/2) - 1){
                for(d = 0; d < dim; d++){
                    a = points[furthest_points[0]][d];
                    b = points[furthest_points[1]][d];
                    center_aux[d] = phi * (b - a) + a;
                }
            }

            else if (current_set_size%2==0 && p == current_set_size/2){
                for(d = 0; d < dim; d++){
                    a = points[furthest_points[0]][d];
                    b = points[furthest_points[1]][d];
                    center[d] = (center_aux[d] + (phi * (b - a) + a))/2;
                }
            }

            proj_table[p] = phi * (b - a) + a;

            delta = 0;
            gamma = 0;
        }


    }

    else {
        for(long p = 0; p < current_set_size; p++){
            for(d = 0; d < dim; d++){
                a = points[furthest_points[0]][d];
                b = points[furthest_points[1]][d];
                point = points[current_set[p]][d];
                delta += (point - a) * (b - a);
                gamma += pow((b - a), 2);
            }
            phi = delta / gamma;

            a = points[furthest_points[0]][0];
            b = points[furthest_points[1]][0];
            proj_table[p] = phi * (b - a) + a;

            if(current_set_size%2==1 && p == current_set_size/2) {
                for(d = 0; d < dim; d++){
                    a = points[furthest_points[0]][d];
                    b = points[furthest_points[1]][d];
                    center[d] = phi * (b - a) + a;
                }
            }
            else if (current_set_size%2==0 && p == (current_set_size/2) - 1){
                for(d = 0; d < dim; d++){
                    a = points[furthest_points[0]][d];
                    b = points[furthest_points[1]][d];
                    center_aux[d] = phi * (b - a) + a;
                }
            }

            else if (current_set_size%2==0 && p == current_set_size/2){
                for(d = 0; d < dim; d++){
                    a = points[furthest_points[0]][d];
                    b = points[furthest_points[1]][d];
                    center[d] = (center_aux[d] + (phi * (b - a) + a))/2;
                }
            }

            //proj_table[p].idx = current_set[p];

            delta = 0;
            gamma = 0;
        }
    }

    return center;

    


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

    if(n>1) {
        #pragma omp parallel for reduction(max:highest) //num_threads(n)
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

    if(n>1){
        long i, local_max_index=0, global_max_index=0;
        double aux, local_max=0, global_max=0;
        
        #pragma omp parallel firstprivate(aux, local_max, local_max_index, i) shared(global_max, global_max_index) //num_threads(n)
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

void swap(double* a, double* b)
{
    double t = *a;
    *a = *b;
    *b = t;
}

void swap2(long* a, long* b)
{
    long t = *a;
    *a = *b;
    *b = t;
}
 
int partition (double* arr, long* other, int low, int high)
{
    int pivot = arr[high]; // pivot
    int i = (low - 1); // Index of smaller element and indicates the right position of pivot found so far
 
    for (int j = low; j <= high - 1; j++)
    {
        // If current element is smaller than the pivot
        if (arr[j] < pivot)
        {
            i++; // increment index of smaller element
            swap(&arr[i], &arr[j]);
            swap2(&other[i], &other[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    swap2(&other[i + 1], &other[high]);
    return (i + 1);
}

void quickSort(double* arr, long* other, int low, int high)
{
    if (low < high)
    {
        int pi = partition(arr, other, low, high);
        quickSort(arr, other, low, pi - 1);
        quickSort(arr, other, pi + 1, high);
    }
}

struct node* build_tree(long node_index, long* current_set, long current_set_size, long rec_level, int nprocs, int whichproc) {
    int p = nprocs/(pow(2,rec_level));
    int next_number_of_subtrees = (pow(2,rec_level+1));
    int n = nthreads/(pow(2,rec_level));

    double exec_findradius;

    if(n>1) {
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

    double* proj_table;
    proj_table = (double*) malloc (current_set_size* sizeof(double));

    double* center;
    if(current_set_size > 2) {
        //compute points a and b, furthest apart in the current set;
        long furthest[2];

        furthest_points(furthest, current_set, current_set_size, n);
        a = furthest[0];
        b = furthest[1];
        //perform the orthogonal projection of all points onto line ab;
        center = orthogonal_projection(current_set_size, current_set, furthest, proj_table, n);
    }
    else {

        //if there are only 2 points, no need to make orthogonal projection
        a = current_set[0];
        b = current_set[1];

        proj_table[0] = points[a][0];
        proj_table[1] = points[b][0];

        center = (double*)malloc(dim*sizeof(double));

        for(long d = 0; d < dim; d++){
            center[d] = (points[a][d] + points[b][d])/2;
        }
        
    }

    quickSort(proj_table, current_set, 0, current_set_size-1);



    //TODO
    // center = find_center(current_set_size, proj_table);

    //rearrange_set(current_set_size, current_set, proj_table, rec_level, n);

    free(proj_table);
    
    struct node* res = (struct node*)malloc(sizeof(struct node));
    res->coordinates = center;

    res->radius = find_radius(center, current_set, current_set_size,rec_level, n);

    res->id = node_index;

    long nextLeftSize = current_set_size/2;
    long nextRightSize = current_set_size%2==0 ? current_set_size/2 :current_set_size/2+1;

    //MPI
    res->left = NULL;
    res -> right = NULL;
    if(whichproc + pow(2, rec_level) < nprocs) {

        //I NEED THE POINTER
        MPI_Send( current_set , nextLeftSize , MPI_LONG ,  whichproc + pow(2, rec_level), 0 , MPI_COMM_WORLD);

        res->right = build_tree(node_index + nextLeftSize*2, current_set+current_set_size/2, nextRightSize, rec_level + 1, nprocs, whichproc); 
       
        //send left and right to two different processes
        //send 

        
    } else { //OPENMP
        #pragma omp task if(n>1) 
        res->left = build_tree(node_index + 1, current_set, nextLeftSize, rec_level + 1, nprocs, whichproc);
        #pragma omp task if(n>1) 
        res->right = build_tree(node_index +nextLeftSize* 2 , current_set+current_set_size/2, nextRightSize, rec_level + 1, nprocs, whichproc);
    }

    return res;
}

void dump_tree(struct node *node, int me){

    printf("%ld ", node->id);
//    fflush(stdout);

    if(node->left != NULL)
        printf("%ld %ld ", node->left->id, node->right->id);
    else{
        //left is null
        if(node->right == NULL) {
            //both null -> LEAF!
            printf("-1 -1 ");
        }
        else {
            //right is not null-> child is on the other processor!
            printf("%ld %ld ", node->id + 1, node->right->id);
        }
    }

    printf("%lf ",node->radius);

    for(int i=0; i<dim; i++){
        printf("%lf ",node->coordinates[i]);
    }
    printf("\n"); fflush(stdin);

    if(node->left != NULL){
        dump_tree(node->left, me);
    }
    if(node->right != NULL){
        dump_tree(node->right, me);
    }

}


int main(int argc, char **argv){
    double elapsed_time;
    int me, nprocs;
    struct node* tree;

    if(argc == 5) {
        nthreads = atoi(argv[4]);
        argc = argc-1;
        omp_set_num_threads(nthreads);
    }
    else {
        nthreads = omp_get_num_procs();
    }
    omp_set_nested(1);

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    //fprintf(stderr, "Hi from: [%d]\n", me); 
    
    //get input sample points (use the function from the guide)

    points = get_points(argc, argv, &dim, &np);
    
    //if(me==0){
    //    printf("%d %ld\n",dim,np*2-1);
//	fflush(stdout);
  //  }
    
    long* current_set = (long*) malloc(np * sizeof(long));

    // #pragma omp parallel for
    

    n_nodes = 0;

    int recv_size = np;

    // MAGIC THAT WORKS

    int aux;
    int level = -1;
    long id = 0;
    for (int i = 0; pow(2, i) <= me; i++) {
        level = i;
        aux = pow(2, i);
        if(recv_size % 2 == 0) {
            recv_size = recv_size/2;
            if((int)((me - aux) / aux) % 2 == 0) {
                id += 1;
            }
            else {
                id += (recv_size * 2);
            }
        }
        else {
            if((int)((me - aux) / aux) % 2 == 0) {
                recv_size = (recv_size - 1) /2;
                id += 1;
            }
            else {
                recv_size = (recv_size + 1) /2;
                id += ((recv_size-1) *2);
            }
        }
    }

    ///

    if(me==0)
        for(int j = 0; j < np; j++) current_set[j] = j;
    else {
        //fprintf(stderr, "[%d] is waiting for  processor %d\n",me, me - aux);
        MPI_Recv( current_set , recv_size , MPI_LONG , me - aux , 0 , MPI_COMM_WORLD , &status);
        
	//fprintf(stderr, "[%d] received: ",me);
	//for (int i = 0; i < recv_size; i++) {
        //    fprintf(stderr, "%ld, ", current_set[i]);
        //}
        //fprintf(stderr, "\n ");
	//fprintf(stderr, "[%d] ROOT NODE: %ld\n",me,id);
    }

    

    tree = build_tree(id, current_set, recv_size, level + 1, nprocs, me);
    // fprintf(stderr, "[%d] will dump tree %ld, pointer: %p\n",me, tree->id, tree);
    //fprintf(stderr, "process %d finished\n", me);
    fflush(stderr);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();
    if(me==0){
    	fprintf(stderr, "%.1lf\n", elapsed_time);
	fflush(stderr);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    //dump_tree(tree, me);
    //fprintf(stderr, "[%d] DUMP FINISHED!\n",me);
    //fflush(stderr);
    free(tree);
    //for(int i= 0; i < np; i++) {
    //	free(points[i]);
    //}
    //free(points);
    MPI_Finalize();

    return 0;
}

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

void orthogonal_projection(long current_set_size, long* current_set, long* furthest_points, struct ProjectedPoint* proj_table, int n){
    double delta = 0, gamma = 0, phi = 0;
    double a, b, point;
    int d;

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

struct node* build_tree(long node_index, long* current_set, long current_set_size, long rec_level, int nprocs, int whichproc) {
    //fprintf(stderr, "[%d] processor entered build_tree with set ", whichproc);
    // for (int i ; i < current_set_size; i++) {
    //     //fprintf(stderr, "%ld, ", current_set[i]);
    // }
    //fprintf(stderr, "\n");
    int p = nprocs/(pow(2,rec_level));
    int next_number_of_subtrees = (pow(2,rec_level+1));
    int n = nthreads/(pow(2,rec_level));

    // double exec_findradius;

    // if(n>1) {
    //     exec_findradius = -omp_get_wtime();
    // }

    #pragma omp atomic
    n_nodes++;

    if(current_set_size == 1) {
        //fprintf(stderr, "[%d] will stop recursion\n", whichproc);
        //stop recursion
        struct node* res = (struct node*)malloc(sizeof(struct node));
        res->id = node_index;
        res->coordinates = points[current_set[0]];
        res->left = NULL;
        res->right = NULL;
        res->radius = 0;

        //fprintf(stderr, "[%d] will return %ld\n", whichproc, res->id);
        return res;
    }

    //fprintf(stderr, "[%d] 1\n", whichproc);

    long a;
    long b;

    struct ProjectedPoint* proj_table;
    proj_table = (struct ProjectedPoint*) malloc (current_set_size* sizeof(struct ProjectedPoint));
    //fprintf(stderr, "[%d] 2\n", whichproc);
    if(current_set_size > 2) {
        //fprintf(stderr, "[%d] 2.5\n", whichproc);
        //compute points a and b, furthest apart in the current set;
        long furthest[2];

        furthest_points(furthest, current_set, current_set_size, n);
        a = furthest[0];
        b = furthest[1];
        //fprintf(stderr, "[%d] 3\n", whichproc);
        //perform the orthogonal projection of all points onto line ab;
        orthogonal_projection(current_set_size, current_set, furthest, proj_table, n);
        //fprintf(stderr, "[%d] 4\n", whichproc);
    }
    else {

        //fprintf(stderr, "[%d] 3.1\n", whichproc);
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

        //fprintf(stderr, "[%d] 4.1\n", whichproc);

    }

    qsort(proj_table, current_set_size, sizeof(struct ProjectedPoint), compare);
    //fprintf(stderr, "[%d] 5\n", whichproc);

    double* center;

    center = find_center(current_set_size, proj_table);
    //fprintf(stderr, "[%d] 5.center\n", whichproc);

    rearrange_set(current_set_size, current_set, proj_table, rec_level, n);
    //fprintf(stderr, "[%d] 5.rearrange\n", whichproc);

    free(proj_table);
    
    struct node* res = (struct node*)malloc(sizeof(struct node));
    res->coordinates = center;
    //fprintf(stderr, "[%d] 5.center = %lf\n", whichproc,center[0]);
    // fprintf(stderr, "[%d] Finding radius of set:\n", whichproc);
    // for(int i=0 ; i<current_set_size ; i++){
    //     fprintf(stderr, "[%d] %d:\n", whichproc,current_set[i]);
    // }

    
    res->radius = find_radius(center, current_set, current_set_size,rec_level, n);
    //fprintf(stderr, "[%d] 5.radius = %lf\n", whichproc, res->radius);


    res->id = node_index;

    long nextLeftSize = current_set_size/2;
    long nextRightSize = current_set_size%2==0 ? current_set_size/2 :current_set_size/2+1;

    //fprintf(stderr, "[%d] 6\n", whichproc);
    //MPI
    res->left = NULL;
    res -> right = NULL;
    if(whichproc + pow(2, rec_level) < nprocs) {

        //I NEED THE POINTER
        fprintf(stderr, "[%d] will send to processor %lf\n",whichproc, whichproc + pow(2, rec_level));
        MPI_Send( current_set , nextLeftSize , MPI_LONG ,  whichproc + pow(2, rec_level), 0 , MPI_COMM_WORLD);
        fprintf(stderr, "[%d] sent!\n", whichproc);

        res->right = build_tree(node_index + nextLeftSize*2, current_set+current_set_size/2, nextRightSize, rec_level + 1, nprocs, whichproc); 
       
        //send left and right to two different processes
        //send 

        
    } else { //OPENMP
        //fprintf(stderr, "[%d] will execute tasks\n",whichproc);
        #pragma omp task if(n>1) 
        res->left = build_tree(node_index + 1, current_set, nextLeftSize, rec_level + 1, nprocs, whichproc);
        #pragma omp task if(n>1) 
        res->right = build_tree(node_index +nextLeftSize* 2 , current_set+current_set_size/2, nextRightSize, rec_level + 1, nprocs, whichproc);
    }

    //fprintf(stderr, "[%d] will return %ld\n", whichproc, res->id);
    return res;
}

void dump_tree(struct node *node, int me){

    printf("%ld ", node->id);

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




struct node* build_tree_distributed(long node_index, long set_size, long rec_level, int nprocs, int whichproc){
    struct node* subtree;

// FURTHEST POINTS A and B
    double * A = calc_A_dist();
    double * B = calc_B_dist(A);

    struct ProjectedPoint* proj_table;
    proj_table = (struct ProjectedPoint*) malloc (current_set_size* sizeof(struct ProjectedPoint));
    //CALCULATE LOCAL ORTHOGONAL PROJECTION TABLE
    orthogonal_projection_v2(A, B, proj_table)
    
    //DISTRIBUTED SORTING ALGORITHM

    //FIND CENTER (MEDIAN POINT)
    //MEMORISE WHICH PROCESS HAS THE CENTER (ONE OF THE TWO CENTRAL PROCESSES(?))

    //FREE LOCAL ORTHOGONAL PROJECTION

    // if(me == "process that has the center"){
    //     SEND CENTER TO EVERYONE
    //     FIND LOCAL MAX DISTANT
    //     RECEIVE DISTANCES
    //     CALCULATE RADIUS = MAX
        
    // } else {
    //     RECEIVE CENTER 
    //     FIND LOCAL MAX DISTANT
    //     SENT TO "process that has the center"
    // }   


    return subtree;
}

double * local_furthest_from_point(double ** set, long set_size, double * point){
    double max=0;
    double distance;
    double * furthest;

    for(int i=0;i<set_size;i++){
        if((distance = distance_between_points(point,*(set+i))) > max){
            max = distance;
            furthest = *(set+i);
        }
    }
    return furthest;
}

double * calc_A_dist(){
    double * first = (double *) malloc(dim*sizeof(double));
    if(!me){
        first = points[0];
    }
    //broadcast point 0
    MPI_Bcast(first,dim,MPI_DOUBLE,0,MPI_COMM_WORLD);
    //now everyone has 'first' point
    //calculate local a, furthest from 'first' (new util function)
    double * a = local_furthest_from_point(points, np, first);
    double ** as;
    if(!me){
    //allocate gathering vector of points (to store as, only p0 needs to)
        as = (double **) malloc(nprocs * sizeof(double *));
        for(int i=0;i<nprocs;i++){
            *(as+i) = (double *) malloc(dim * sizeof(double));
        }
    }
    //gather everyone's 'a' at p0
    MPI_Gather(a, dim, MPI_DOUBLE, as, dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //calculate real A = furthest among 'a's
    double * A;
    if(!me)
        double * A = local_furthest_from_point(as, nprocs, first);

    //Bcast A
    MPI_Bcast(A, dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //now everyone has the real A of the set
    return A;
}

double * calc_B_dist(double * A){
    //calculate local b
    double * b = local_furthest_from_point(points, np, A);
    double ** bs;
    //allocate gathering vector of points (to store bs, only p0 needs to)
    if(!me){
        double ** bs = (double **) malloc(nprocs * sizeof(double *));
        for(int i=0;i<nprocs;i++){
            *(bs+i) = (double *) malloc(dim * sizeof(double));
        }
    }

    //gather everyone's 'b' at p0
    MPI_Gather(b, dim, MPI_DOUBLE, bs, dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //calculate real B = furthest among 'b's
    double * B;
    if(!me)
        double * B = local_furthest_from_point(bs, nprocs, A);
    //send B to everyone
    MPI_Bcast(B,dim,MPI_DOUBLE,0,MPI_COMM_WORLD);
    //now everyone has the real B of the set
    return B;
}

void orthogonal_projection_v2(double* A, double* B, struct ProjectedPoint* proj_table){
    double delta = 0, gamma = 0, phi = 0;
    double a, b, point;

    for(long p = 0; p < np; p++){
        for(int d = 0; d < dim; d++){
            a = A[d];
            b = B[d];
            point = points[p][d];
            delta += (point - a) * (b - a);
            gamma += pow((b - a), 2);
        }
        phi = delta / gamma;
        proj_table[p].projectedCoords = (double*) malloc(dim * sizeof(double));

        for(int d = 0; d < dim; d++){
            a = A[d];
            b = B[d];
            proj_table[p].projectedCoords[d] = phi * (b - a) + a;
        }
        proj_table[p].idx = p;
        delta = 0;
        gamma = 0;
    }
}

int main(int argc, char **argv){
    double elapsed_time;
    int me, nprocs;
    struct node* tree;
    struct node* subtree;


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
    
    fprintf(stderr, "Hi from: [%d]\n", me); 
    
    //get input sample points (use the function from the guide)

    points = get_points(argc, argv, &dim, &np, me, nprocs);
    
    if(me==0){
        printf("%d %ld\n",dim,np*2-1);
    }
    
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

    /// SUPER MAGIC



    // if(me==0)
        for(int j = 0; j < np; j++) current_set[j] = j;
    // else {
    //     fprintf(stderr, "[%d] is waiting for  processor %d\n",me, me - aux);
    //     MPI_Recv( current_set , recv_size , MPI_LONG , me - aux , 0 , MPI_COMM_WORLD , &status);
        
	// fprintf(stderr, "[%d] received: ",me);
	// for (int i = 0; i < recv_size; i++) {
    //         fprintf(stderr, "%ld, ", current_set[i]);
    //     }
    //     fprintf(stderr, "\n ");
	// fprintf(stderr, "[%d] ROOT NODE: %ld\n",me,id);
    // }

    
    subtree = build_tree_distributed(id, np, 0, nprocs, me);

    if(!me){
        //receive x points from right half of processes
    }
    else if (me==1){
        //receive y points from left half of processes
    } else if("I'm right sided"){
        //send to 0
    } else if ("I'm left sided"){
        //send to 1
    }

    tree = build_tree(id, current_set, recv_size, level + 1, nprocs, me);
    // fprintf(stderr, "[%d] will dump tree %ld, pointer: %p\n",me, tree->id, tree);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();
    if(me==0){
        printf("%d %ld\n",dim,np*2-1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    dump_tree(tree, me);
    fprintf(stderr, "[%d] DUMP FINISHED!\n",me);
    free(tree);
    MPI_Finalize();

    return 0;
}
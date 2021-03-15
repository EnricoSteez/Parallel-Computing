#include <stdio.h>

#define DIMENSIONS 4

typedef struct node{
    float coordinates[DIMENSIONS];
    float radius;
    node * left;
    node * right;
} node;

//Usually success = 0, failure = 1, we can change this, it's not relevant 

int find_extremes(node * a, node * b){
    

    //TODO

    return 0;
}

float orthogonal_projection(node n){
    float distance = 0;

    return distance;
}

int main(int argc, char **argv){

    //1. get input sample points (use the function from the guide)

    //2. compute points a and b, furthest apart in the current set;
    //passing pointers to nodes so that the function can set a and b to the extremes
    node a;
    node b;
    find_extremes(&a, &b);

    //3. perform the orthogonal projection of all points onto line ab;

    

    //4. compute the center, defined as the median point over all projections;

    //5. create two sets of points, L and R, deﬁned as the points whose projection lies to one side or the other of the center;

    

    return 0;
}
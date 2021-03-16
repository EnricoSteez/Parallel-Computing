//
//  main.cpp
//  Ball Tree Algorithm
//
//  Created by Miguel de Oliveira Guerreiro on 15/03/2021.
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "gen_points.c"



void calculate_distance(double **points_table, long n_points, int n_dim, int* indexes){
    double distance = 0, highest = 0;

    for(long pi = 0; pi < n_points - 1; pi++){
        for (long pj = pi + 1; pj < n_points; pj++){
            for(int dim = 0; dim < n_dim; dim++){
                distance += pow((points_table[pi][dim] - points_table[pj][dim]), 2);
            }
            if (distance > highest){
                highest = distance;
                indexes[0] = pi;
                indexes[1] = pj;
            }
            distance = 0;
        } 
    }
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

int main(int argc, char * argv[]) {
    int n_dim, random_seed;
    long n_points;
    
    if(argc != 4){
        printf("Usage: %s <n_dims> <n_points> <seed>\n", argv[0]);
        exit(1);
    }
    
    double ** points_table;
    points_table = get_points(argc, argv, &n_dim, &n_points);
    
    n_dim = atoi(argv[1]);
    n_points = atol(argv[2]);
    random_seed = atoi(argv[3]);
    
    for(int i = 0; i < n_points ; i++){
        printf("Point %d ", (i));
        for(int k = 0; k < n_dim; k++){
            printf("Dim %d, value = %f ", k, points_table[i][k]);
        }
        printf("\n");
    }

    int indexes[2];
    double proj_table[n_points][n_dim]; 
    calculate_distance(points_table, n_points, n_dim, indexes);
    orthogonal_projection(points_table, n_points, n_dim, indexes, proj_table);
    
    printf("%d, %d\n", indexes[0], indexes[1]);
    printf("PROJECTION TABLE\n");
    printf("%f, %f\n", proj_table[0][0], proj_table[1][1]);
    
    for(int i = 0; i < n_points ; i++){
        printf("Point %d ", (i));
        for(int k = 0; k < n_dim; k++){
            printf("Dim %d, value = %f ", k, proj_table[i][k]);
        }
        printf("\n");
    }
    return 0;
}

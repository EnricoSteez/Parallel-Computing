#include <stdio.h>
#include <math.h>
#include <string.h>

int main() {
    
    int np = 5;
    int me = 4;
    int recv_size = np;

    // MAGIC THAT WORKS
    int id = 0;
    int aux;
    int level = -1;
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

    printf("id: %d\n", id);
    return 0;
}
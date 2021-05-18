#include <stdio.h>
#include <math.h>

int main() {
    int aux;
    int level = -1;
    int me = 2;
    int recv_size = 5;
    for (int i = 0; pow(2, i) <= me; i++) {
        level = i;
        aux = pow(2, i);
        if(recv_size % 2 == 0) {
            recv_size = recv_size/2;
        }
        else {
            if((int)((me - aux) / aux) % 2 == 0) {
                recv_size = (recv_size - 1) /2;
            }
            else {
                recv_size = (recv_size + 1) /2;
            }
        }
    }
    printf("%d\n", recv_size);
}
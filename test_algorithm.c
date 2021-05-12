#include <math.h>
#include <stdio.h>

int main(int argc, char** argv) {
    int recv_size = 10;
    int me = 6;

    int aux;

    for (int i = 0; pow(2, i) <= me; i++) {
        aux = pow(2, i);
        if(recv_size % 2 == 0) {
            recv_size = recv_size/2;
        }
        else {
            if((int)((me - aux) / aux) % 2 == 0) {
                recv_size = (recv_size + 1) /2;
            }
            else {
                recv_size = (recv_size - 1) /2;
            }
        }
    }
    printf("recv_size : %d\n", recv_size);
    printf("parent : %d\n", me-aux);

    return 0;
}


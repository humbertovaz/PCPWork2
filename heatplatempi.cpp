nclude <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#define ITER 2
/*void iterationParallelSwapCpy(){
    int iter=0;
    while (iter < ITER){
        wtime = omp_get_wtime ();
        //#pragma omp parallel for
        for(int i=1; i<N-1; ++i){ 
            for(int j=1; j<M-1; ++j){
                G1[i][j] = 0.2*(
                    G2[i-1][j]+
                    G2[i+1][j]+
                    G2[i][j-1]+
                    G2[i][j+1]+
                    G2[i][j]);
            }
        }
            ++iter;
    }  
           
            //Copiar de volta para a Memória
            //#pragma omp parallel for
            for(int i=1; i<N-1;++i){
                    for(int j=1;j<M-1;++j)
                        G2[i][j]=G1[i][j];          
            }
}
*/
/*
void iterationSequentialCopSwap(){
    int iter=0;
    double** temp;
    while (iter < ITER){
        for(int i=1; i<N-1; i++){ 
            for(int j=1; j<M-1; j++){
                G1[i][j] = 0.2*(
                            G2[i-1][j]+
                            G2[i+1][j]+
                            G2[i][j-1]+
                            G2[i][j+1]+
                            G2[i][j]);
            }   
        }
    ++iter;
    temp = G2;
    G2 = G1;
    G1 = temp;
    }
              
}*/


int main( int argc, char *argv[]) {
int rank, msg;
MPI_Status status;
MPI::Init(&argc, &argv);
rank= MPI::MPI_COMM_WORLD.Get_rank();
/* Process 0 sends and Process 1 receives */
    if (rank == 0) {
        msg = 123456;
        //MPI::MPI_COMM_WORLD.Send( &msg, 1, MPI::INT, 1, 0);
    }/*
    else if (rank == 1) {
        MPI::MPI_COMM_WORLD.Recv(&msg, 1, MPI::INT, 0, 0)
        printf( "Received %d\n", msg);
    }*/
MPI::Finalize();
return 0;
}

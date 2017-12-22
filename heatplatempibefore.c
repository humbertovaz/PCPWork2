/****************************************************************************
 * FILE: mpi_heat2D.c
 * OTHER FILES: draw_heat.c  
 * DESCRIPTIONS:  
 *   HEAT2D Example - Parallelized C Version
 *   This example is based on a simplified two-dimensional heat 
 *   equation domain decomposition.  The initial temperature is computed to be 
 *   high in the middle of the domain and zero at the boundaries.  The 
 *   boundaries are held at zero throughout the simulation.  During the 
 *   time-stepping, an array containing two domains is used; these domains 
 *   alternate between old data and new data.
 *
 *   In this parallelized version, the grid is decomposed by the master
 *   process and then distributed by rows to the worker processes.  At each 
 *   time step, worker processes must exchange border data with neighbors, 
 *   because a grid point's current temperature depends upon it's previous
 *   time step value plus the values of the neighboring grid points.  Upon
 *   completion of all time steps, the worker processes return their results
 *   to the master process.
 *
 *   Two data files are produced: an initial data set and a final data set.
 *   An X graphic of these two states displays after all calculations have
 *   completed.
 * AUTHOR: Blaise Barney - adapted from D. Turner's serial C version. Converted
 *   to MPI: George L. Gusciora (1/95)
 * LAST REVISED: 06/12/13 Blaise Barney
 * 
 * Para executar:
 * mpiexec -np 4 ./ex2
 * 
 ****************************************************************************/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
extern void draw_heat(int nx, int ny);       /* X routine to create graph */

#define M      20                       /* x dimension of problem grid /Collumns*/
#define N      20                       /* y dimension of problem grid /Rows */
#define ITER       100                  /* number of time steps */
#define MAXWORKER   8                  /* maximum number of worker tasks */
#define MINWORKER   3                  /* minimum number of worker tasks */
#define BEGIN       1                  /* message tag */
#define MIDDLE      
#define END         99                 /*Final row */
#define MASTER      0                  /* taskid of first process */

double **G1;
double **G2;

struct Parms { 
  float cx;
  float cy;
} parms = {0.1, 0.1};

void fillMatrix(){
    int i,j;
        //CIMA
            for ( i = 0; i < M ; i++ ){
                G2[0][i] = 100;
                G1[0][i] = 100;
            }
        //BAIXO
            for ( i = 0; i < M ; i++ ){
                G2[N-1][i] = 100;
                G1[N-1][i] = 100;
            }
        //ESQUERDA
            for ( i = 1; i < N-1 ; i++ ){
                G2[i][0] = 100;
                G1[i][0] = 100;
            }
        //DIREITA    
            for ( i = 1; i < N -1 ; i++ ){
                G2[i][M-1] = 100;
                G1[i][M-1] = 100;
            }
        //RESTA    
            for ( i=1; i < N - 1 ; i++)
                for( j=1; j< M - 1 ; j++){
                    G2[i][j] = 0;
                }
    
}
void init(){
    int i;
    G1 = (double **) malloc(N*sizeof(double));
    G2 = (double **) malloc(N*sizeof(double));
    
    for(i = 0; i < N; i++){
        G1[i] = (double *) malloc(M*sizeof(double));
        G2[i] = (double *) malloc(M*sizeof(double));
    }
}

/**************************************************************************
 *  subroutine update
 ****************************************************************************/
void update(int start, int end, int ny, float *u1, float *u2)
{
   int ix, iy;
   for (ix = start; ix <= end; ix++) 
      for (iy = 1; iy <= ny-2; iy++) 
         *(u2+ix*ny+iy) = *(u1+ix*ny+iy)  + 
                          parms.cx * (*(u1+(ix+1)*ny+iy) +
                          *(u1+(ix-1)*ny+iy) - 
                          2.0 * *(u1+ix*ny+iy)) +
                          parms.cy * (*(u1+ix*ny+iy+1) +
                         *(u1+ix*ny+iy-1) - 
                          2.0 * *(u1+ix*ny+iy));
}

/*****************************************************************************
 *  subroutine inidat
 *****************************************************************************/
void inidat(int nx, int ny, float *u) {
int ix, iy;

for (ix = 0; ix <= nx-1; ix++) 
  for (iy = 0; iy <= ny-1; iy++)
     *(u+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
}

/**************************************************************************
 * subroutine prtdat
 **************************************************************************/
void prtdat(int nx, int ny, float *u1, char *fnam) {
int ix, iy;
FILE *fp;

fp = fopen(fnam, "w");
for (iy = ny-1; iy >= 0; iy--) {
  for (ix = 0; ix <= nx-1; ix++) {
    fprintf(fp, "%8.1f", *(u1+ix*ny+iy));
    if (ix != nx-1) 
      fprintf(fp, " ");
    else
      fprintf(fp, "\n");
    }
  }
fclose(fp);
}
void compute(double** linhas, double * linha){
            int i=1;
            int j;
            for(j=1; j<M-1; j++){
                linha[j] = 0.2*(
                            linhas[i-1][j]+
                            linhas[i+1][j]+
                            linhas[i][j-1]+
                            linhas[i][j+1]+
                            linhas[i][j]);
            }
}


int main (int argc, char *argv[]){
double* linha;
int	taskid,                     /* this task's unique id */
	numtasks,                   /* number of tasks */
	source,               /* to - from for message send-receive */
	index;              /* loop variables */
MPI_Status status;
int j,i;


/* First, find out my taskid and how many tasks are running */
   MPI_Init(NULL,NULL);
   MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
   MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
   printf("taskid:%d\n",taskid);
   if (taskid == MASTER) {
       printf("tid%d\n", taskid);
      /************************* master code *******************************/
      /* Initialize grid */
      init();
      fillMatrix();
      /* Distribute work to workers.  Must first figure out how many rows to */
      /* send and what to do with extra rows.  */
      for (i=1; i<=N-1; i++)
      {
         /*  Now send startup information to each worker  */
         index = i;
         MPI_Send(&i,1,MPI_INT,index,MASTER,MPI_COMM_WORLD); // Enviar indice para alocação
         MPI_Send(&G1[i-1], 3*M, MPI_DOUBLE,index,MASTER,MPI_COMM_WORLD); //Envia 3 linhas para o iésimo slave 
         printf("Sent to worker %d: rows=%d,%d,%d\n",index, index-1, index, index+1);
      }
      //Receber o que os workers mandam no final do trabalho
      /* Now wait for results from all worker tasks */
      for (i=0; i<N; i++){
         // tag vai conter a linha onde o Master vai inserir em G2 
         linha = (double*) malloc(M*sizeof(double));
        printf("antes do receive: %d\n",i);
        MPI_Recv(linha, M, MPI_DOUBLE, MPI_ANY_SOURCE, index, MPI_COMM_WORLD, &status);
        printf("depois do receive: %d\n",i);
      //Atualizar para G2
        for(i=0; i <M;i++){
            G2[index][i]= linha[i];
        }
      }
      //fazer Swap
      
      for(i=0;i<N;i++)
        for(j=0;j<M;j++)
            printf("G2[%d][%d]=%f \n",i,j,G2[i][j]);

   }   /* End of master code */

   /************************* workers code **********************************/
  else{
       printf("ola\n");
      /* Receive my offset, rows, neighbors and grid partition from master */
      source = MASTER;

      printf("A\n");
      MPI_Recv(&i, 1, MPI_INT, 0, index, MPI_COMM_WORLD, &status); //Receber indice para alocação de G1
      printf("Depois do 1º receive(slave) nr:%d\n",i);
      G1= (double**) malloc(sizeof(double*)*3);
      for(j=0;j<3;j++)
        G1[j]= (double*)malloc(M*sizeof(double)); 
      
      MPI_Recv(G1, 3*M, MPI_DOUBLE, source, index, MPI_COMM_WORLD, &status); // Recebe 3 linhas do Master
      printf("Depois do 2º receive(slave) nr:%d\n",i);
        linha = (double*) malloc(M*sizeof(double)); 
        printf("antes do compute\n");
        for(j=0;j<3;j++)
            printf("G1[%d][M-1]=%f\n",j,G1[j][M-1]);
        compute(G1, linha);
        printf("depois do compute\n");
        //Envia computação para o Master
        MPI_Send(linha,M, MPI_DOUBLE,MASTER,index,MPI_COMM_WORLD);
        printf("B\n");
    }

      MPI_Finalize();
}



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


int main (int argc, char *argv[]){
void inidat(), prtdat(), update();
int	taskid,                     /* this task's unique id */
	numworkers,                 /* number of worker processes */
	numtasks,                   /* number of tasks */
	averow,offset,extra,   /* for sending rows of data */
	dest, source,               /* to - from for message send-receive */
	left,right,        /* neighbor tasks */
	msgtype,                    /* for message types */
	rc,start,end,               /* misc */
	i,ix,iy,iz,it,index;              /* loop variables */
MPI_Status status;


/* First, find out my taskid and how many tasks are running */
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
   MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
   numworkers = numtasks-1;

   if (taskid == MASTER) {
      /************************* master code *******************************/
      /* Initialize grid */
      init();
      fillMatrix();
      /* Distribute work to workers.  Must first figure out how many rows to */
      /* send and what to do with extra rows.  */
      averow = 3*M;
      for (i=1; i<=N-1; i++)
      {
         /*  Now send startup information to each worker  */
         index = i;
         MPI_Send(&G1[i-1][0], averow, MPI_DOUBLE,index,MASTER,MPI_COMM_WORLD);
         printf("Sent to worker %d: rows=%d,%d,%d\n",index, index-1, index, index+1);
      }
      //Receber o que os workers mandam no final do trabalho
      /* Now wait for results from all worker tasks */
      for (i=0; i<N; i++){
         // tag vai conter a linha onde o Master vai inserir em G2 
         MPI_Recv(&linha, 1, MPI_DOUBLE, MPI_ANY_SOURCE, index, MPI_COMM_WORLD, &status);
      }
      //fazer Swap

      MPI_Finalize();
   }   /* End of master code */



   /************************* workers code **********************************/
   if (taskid != MASTER) 
   {

      /* Receive my offset, rows, neighbors and grid partition from master */
      source = MASTER;
      MPI_Recv(&G1[i-1][0], averow, MPI_DOUBLE, source, index, MPI_COMM_WORLD, &status);
        double* linha = (double*) malloc(M*sizeof(double)); 
        compute(&G1[i-1][0], linha);
        //Envia computação para o Master
        MPI_Send(linha,M, MPI_DOUBLE,MASTER,index,MPI_COMM_WORLD);
    }
      MPI_Finalize();
   }
}


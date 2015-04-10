/* File:       
 *    mpi_project.c
 *
 * Purpose:    
 *    A program for calculating PageRanks that uses MPI
 *
 * Compile:    
 *    mpicc -g -Wall -o mpi_project mpi_project.c
 * Usage:        
 *    mpiexec -n<number of processes> ./mpi_project < number of nodes>
 *    <max index of nodes>
 *
 * Input:      
 *    File with nodes and their edges called data_input
 * Output:     
 *    PageRanks of all the nodes in the input file
 *
 * Algorithm:  
 *    
 *
 * IPP:  Section 3.1 (pp. 84 and ff.)
 */
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <omp.h> 
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

const int MAX_STRING = 100;
int N,max_index;
int *l, *map;
double *r1, *r2;
int** D;

void Usage(char* prog_name);
int loadinput(int thread_count);
int saveoutput();

int main(int argc, char* argv[]) {
  int i,j,t;
  double sum;
  double d = 0.85;
  double c;
  double start, end, min_start, max_end;

  int rank, processes;
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &processes );

  if (argc != 3) Usage(argv[0]);
  N = strtol(argv[1], NULL, 10);
  max_index=strtol(argv[2],NULL,10);

  l=calloc(N,sizeof(int));
  r1=malloc(N*sizeof(double));
  r2=malloc(N*sizeof(double));
  map=calloc(max_index+1,sizeof(int));
  for (i=0; i<=max_index; i++) {
    map[i] = -1;
  }

  D=malloc(N*sizeof(int*));

  if(D == NULL) {
    fprintf(stderr, "Out of memory\n");
    exit(0);
  }
  for (i= 0; i < N; i++) {
    D[i]=calloc(N,sizeof(int));
    if(D[i] == NULL) {
      fprintf(stderr, "Out of memory\n");
      exit(0);
    }
  }

  start = MPI_Wtime();

  loadinput(processes);

  c = (1-d)/(double)N;

  for (i= 0; i < N; i++) {
    r2[i] = (double) 1/N;
  }


  int block_size;
  if ( rank == processes) {
    block_size = N - (N/processes)*(processes-1);
  } else {
    block_size = N/processes;
  }


  int recvcounts[processes];
  int displs[processes];

  for (t=0; t<processes; t++) {
    recvcounts[t]=block_size;
    displs[t]=t*block_size;
  }
  recvcounts[processes-1]=N-(N/processes)*(processes-1);
  displs[processes-1]=(processes-1)*N/processes;

  for (t=0; t<10; t++) {
    MPI_Allgatherv(r2+block_size*rank, block_size, MPI_DOUBLE, r1, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
    for (i=0; i<block_size; i++) {
      sum = 0;
      for (j=0; j<N; j++) {
	if (D[i+rank*block_size][j]) {
	  sum+= r1[j]/(double)l[j];
	}
      }
      r2[i+rank*block_size] = c + d*sum;
    }
  }
	
  end = MPI_Wtime();

  MPI_Reduce(&start,&min_start,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&end,&max_end,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
	
  if (rank == 0) {
    printf("Time taken:%f\n", max_end-min_start);
    saveoutput(r1,r2,map,N);
  }

  free(D);
  free(l);
  free(r1);
  free(r2);

  MPI_Finalize();
  return 0;
}  /* main */

/*-------------------------------------------------------------------*/
void Usage(char* prog_name) {
  fprintf(stderr, "usage: %s <number of nodes> <max node index>\n", prog_name);
  exit(0);
}  /* Usage */

/*-------------------------------------------------------------------*/
int loadinput(int thread_count)
{
  FILE* ip;
  int src, dst;

  if ((ip=fopen("data_input","r"))==NULL)
    {
      printf("error opening the input data.\n");
      return 1;
    }

  int index = 0;

  while(!feof(ip))
    {
      fscanf(ip, "%d\t%d\n", &src, &dst);
      if (map[src] == -1) {
        map[src] = index;
        src = index;
        index++;
      } else {
        src = map[src];
      }
      if (map[dst] == -1) {
        map[dst] = index;
        dst = index;
        index++;
      } else {
          dst = map[dst];
      }
      l[src]++;
      D[dst][src]=1;
    }

  int i, j;
  //#pragma omp parallel for private(i,j)
  for (i=0; i<N; i++) {
    if (l[i] == 0) {
      l[i] = N-1;
      for (j=0; j < N; j++) {
	if (j != i) {
	  D[j][i] = 1;
	}
      }
    }
  }
  fclose(ip);     
  return 0;
}

/*-------------------------------------------------------------------*/
int saveoutput(double *r1, double *r2, int *map, int N)
{
  FILE* op;

  if ((op=fopen("data_output","w"))==NULL)
    {
      printf("error opening the input data.\n");
      return 1;
    }
  int i,k,node_index;
  double sum = 0;

  for (i=0; i<N; i++) {
    e = fabs(r1[i] - r2[i])/r2[i];
    for (k=1; k<=max_index; k++) {
      if (map[k] == i) {
	node_index = k;
	break;
      }
    }
    fprintf(op, "%d\t%f\n", node_index, r1[i]);
    sum+=r2[i];
  }
  //printf("%f\n", sum);
  fclose(op);     
  return 0;
}
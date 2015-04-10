/* File:       
 *    mpi_project.c
 *
 * Purpose:    
 *    A program for calculating PageRanks that uses MPI
 *
 * Compile:    
 *    mpicc -g -Wall -std=C99 -o mpi_project mpi_project.c
 * Usage:        
 *    mpiexec -n<number of processes> ./mpi_project
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
#include "csparse.h"

const int MAX_STRING = 100;
int N,max_index;
int *l, *map;
double *r1, *r2;
int** D;
cs *ADJ;

void Usage(char* prog_name);
int loadinput(int **D, int *l, int *map, int N, int thread_count);
int saveoutput(double *r1, double *r2, int *map, int N);

int main(int argc, char* argv[]) {
  int i,j,t,k;
  double sum;
  double d = 0.85;
  double c;
  double* tmp;
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
  map=calloc(max_index,sizeof(int));
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

  loadinput(D,l,map,N,processes);

  /*if (rank==0) {
    int *Ai, *Ap, nz;
    double *Ax;
    Ai=ADJ->i;
    Ap=ADJ->p;
    Ax=ADJ->x;
    nz=ADJ->nz;
    for (i=0;i<nz;i++) {
    printf("%d %d : %g\n",Ai[i],Ap[i],Ax?Ax[i]:1);
    }
    }*/

  c = (1/(double)N)*(1-d);
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

  //printf("I'm number %d and I like %d\n",rank,recvcounts[rank]);

	
  for (t=0; t<10; t++) {
    MPI_Allgatherv(r2+block_size*rank, block_size, MPI_DOUBLE, r1, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
    for (i=0; i<block_size; i++) {
      sum = 0;
      for (j=0; j<N; j++) {
	for (k=ADJ->p[j];k<ADJ->p[j+1];k++) {
	  if (ADJ->i[k]==i+rank*block_size) {
	    //sum+=r1[j]/(double)l[j];
	    break;
	  }
	}
	if (D[i+rank*block_size][j] == 1) {
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
    printf("Time taken:%f\n", end-start);
    saveoutput(r1,r2,map,N);
  }

  free(D);
  free(l);
  free(r1);
  free(r2);
  cs_spfree(ADJ);

  MPI_Finalize();
  return 0;
}  /* main */

/*-------------------------------------------------------------------*/
void Usage(char* prog_name) {
  fprintf(stderr, "usage: %s <number of nodes> <max node index>\n", prog_name);
  exit(0);
}  /* Usage */

/*-------------------------------------------------------------------*/
int loadinput(int **D, int *l, int *map, int N, int thread_count)
{
  FILE* ip;
  int src, dst;

  if ((ip=fopen("data_input","r"))==NULL)
    {
      printf("error opening the input data.\n");
      return 1;
    }
  int k, size = 0;
  bool foundSrc, foundDst;

  ADJ=cs_spalloc(0, 0, 1, 1, 1);

  int index = 0;
  int zero;
  bool first_time = true;
  while(!feof(ip))
    {
      fscanf(ip, "%d\t%d\n", &src, &dst);
      if (first_time) {
	map[src] = index;
	zero = src;
	src = index;
	index++;
	map[dst] = index;
	dst = index;
	index++;
	first_time = false;
      } else {
	if (map[src] == 0) {
	  if (src != zero) {
	    map[src] = index;
	    src = index;
	    index++;
	  } else {
	    src = 0;
	  }
	} else {
	  src = map[src];
	}
	if (map[dst] == 0) {
	  if (dst != zero) {
	    map[dst] = index;
	    dst = index;
	    index++;
	  } else {
	    dst = 0;
	  }
	} else {
	  dst = map[dst];
	}
      }
		
	        
      l[src]++;
      D[dst][src]=1;
      cs_entry(ADJ,dst,src,1);
    }
  int i, j;
  //#pragma omp parallel for private(i,j)
  for (i=0; i<N; i++) {
    if (l[i] == 0) {
      l[i] = N-1;
      for (j=0; j < N; j++) {
	if (j != i) {
	  D[j][i] = 1;
	  cs_entry(ADJ,j,i,1);
	}
      }
    }
  }
  ADJ=cs_triplet(ADJ);
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
  double e;

  for (i=0; i<N; i++) {
    e = fabs(r1[i] - r2[i])/r2[i];
    for (k=0; k<max_index; k++) {
      if (map[k] == i) {
	node_index = k;
	break;
      }
    }
    fprintf(op, "%d\t%f\t%f\t%f\n", node_index, r1[i], r2[i], e);
    sum+=r2[i];
  }
  printf("%f\n", sum);
  fclose(op);     
  return 0;
}

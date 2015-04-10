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
#include <string.h>  /* For strlen             */
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "timer.h"

const int MAX_STRING = 100;
int N;
int *l, *map;
double *r1, *r2;
int** D;

void Usage(char* prog_name);
int loadinput(int **D, int *l, int *map, int N);
int saveoutput(double *r1, double *r2, int *map, int N);

int main(int argc, char* argv[]) {
	int i,j,t;
	double sum;
	double d = 0.85;
	double c;
	double* tmp;
	double start, end;

	if (argc != 2) Usage(argv[0]);
	N = strtol(argv[1], NULL, 10);  
	//N = 10697;

	l=calloc(N,sizeof(int));
	r1=malloc(N*sizeof(double));
	r2=malloc(N*sizeof(double));
	map=calloc(N,sizeof(int));
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

	GET_TIME(start);
	loadinput(D,l,map,N);

	for (i= 0; i < N; i++) {
		r1[i] = (double) 1/N;
	}

	c = (1/(double)N)*(1-d);

	for (t=0; t<19; t++) {
		for (i=0; i<N; i++) {
			sum = 0;
			for (j=0; j<N; j++) {
				if (D[i][j] == 1) {
					sum+= r1[j]/(double)l[j];
				}
			}
			r2[i] = c + d*sum;
		}
		tmp=r1;
		r1=r2;
		r2=tmp;
	}		
	
	GET_TIME(end);

	printf("Time taken:%f\n", end-start);

	saveoutput(r1,r2,map,N);
	
	free(D);
	free(l);
	free(r1);
	free(r2);
   return 0;
}  /* main */

/*-------------------------------------------------------------------*/
void Usage(char* prog_name) {
   fprintf(stderr, "usage: %s <number of nodes>\n", prog_name);
   exit(0);
}  /* Usage */

/*-------------------------------------------------------------------*/
int loadinput(int **D, int *l, int *map, int N)
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
	while(!feof(ip))
	{
        	fscanf(ip, "%d\t%d\n", &src, &dst);
		foundSrc = 0;
		foundDst = 0;
		for(k=0; k<size; k++) {
			if (map[k] == src) {
				src = k;
				foundSrc = 1;
			}
			if (map[k] == dst) {
				dst = k;
				foundDst = 1;
			}
		}
		if (!foundSrc) {
			map[size] = src;
			src = size;
			size++;
		}
		if (!foundDst) {
			map[size] = dst;
			dst = size;
			size++;
		}
		l[src]++;
		D[dst][src]=1;
	}
	int i, j;
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

        if ((op=fopen("serial_data_output","w"))==NULL)
        {
                printf("error opening the input data.\n");
                return 1;
        }
	int i;
	double sum = 0;
	double e;

	for (i=0; i<N; i++) {
		e = fabs(r1[i] - r2[i])/r2[i];
		fprintf(op, "%d\t%f\t%f\t%f\n", map[i], r1[i], r2[i], e);
		sum+=r1[i];
	}
	printf("%f\n", sum);
	fclose(op);     
        return 0;
}

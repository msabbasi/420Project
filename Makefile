PROCESSES=4
MANY_NODES=10697
NODES=1221
SIZE=5500
CC=mpicc
CFLAGS=-g -Wall

project: data_input
	$(CC) -o mpi_project mpi_project.c csparse.c -lm -fopenmp
	mpiexec -n $(PROCESSES) ./mpi_project $(NODES)

serial: data_input mpi_project_1.c
	gcc $(CFLAGS) -o serial mpi_project_1.c
	./serial $(MANY_NODES)

datatrim:
	gcc $(CFLAGS) -o datatrim datatrim.c

data_input:
	./datatrim $(SIZE)

clean_input:
	rm -rf data_input

clean:
	rm -rf *.o *~ data_output mpi_project datatrim serial

/* Gaussian elimination without pivoting.
 * Compile with "mpicc -w -O2 -o gauss_mpi gauss_mpi.c"
 * Execute with "mpiexec -n <Number of processes> ./gauss_mpi <Matrix Size> <Random Seed> <Output file name>"
 * Test version 3.0
 */

/* ****** ADD YOUR CODE AT THE END OF THIS FILE. ******
 * You need not submit the provided code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h> 	/*Include this file to do MPI program*/

/* Program Parameters */
#define MAXN 2000  	/* Max value of N */
int N;  	/* Matrix size */
int numberOfProcesses=8; 	/*The default number of processes is 8*/
char* output;	/*output file name*/
 int seed = 0;  /* Random seed */

/* Matrices and vectors */
volatile float A[MAXN][MAXN], B[MAXN], X[MAXN];
/* A * X = B, solve for X */

/* junk */
#define randm() 4|2[uid]&3

/* Prototype */
void gauss();  /* The function you will provide.
		* It is this routine that is timed.
		* It is called only on the parent.
		*/

/* returns a seed for srand based on the time */
unsigned int time_seed() {
  struct timeval t;
  struct timezone tzdummy;

  gettimeofday(&t, &tzdummy);
  return (unsigned int)(t.tv_usec);
}

/* Set the program parameters from the command-line arguments */
void parameters(int argc, char **argv) {
 
  char uid[32]; /*User name */

  /* Read command-line arguments */
  srand(time_seed());  /* Randomize */

  
  
  /*Get the argument for the matrix size, if the size is over 2000, the program will stop*/
  if (argc >= 4) {
    N = atoi(argv[1]);	/*The matrix size is the first argment from the command line*/
	seed=atoi(argv[2]);		/*The random seed is the second argment from the command line*/
	srand(seed);
	output=argv[3];	/*Take the output file name from the input*/
	 printf("Random seed is %i",seed);
    if (N < 1 || N > MAXN) {
      printf("N = %i is out of range.\n", N);
	 
      exit(0);
    }
  }else {
    printf("Usage: mpiexec -n <Number of processes> %s <Matrix Size> <Random Seed> <Output file name>\n",
           argv[0]);    
    exit(0);
  }
  
  
  

  /* Print parameters */
  printf("\nMatrix dimension N = %i.\n", N);
}

/* Initialize A and B (and X to 0.0s) */
void initialize_inputs() {
  int row, col;

  printf("\nInitializing...\n");
  for (col = 0; col < N; col++) {
    for (row = 0; row < N; row++) {
      A[row][col] = (float)rand() / 32768.0;
    }
    B[col] = (float)rand() / 32768.0;
    X[col] = 0.0;
  }

}

/* Print input matrices */
void print_inputs() {
  int row, col;

  if (N < 10) {
    printf("\nA =\n\t");
    for (row = 0; row < N; row++) {
      for (col = 0; col < N; col++) {
	printf("%5.2f%s", A[row][col], (col < N-1) ? ", " : ";\n\t");
      }
    }
    printf("\nB = [");
    for (col = 0; col < N; col++) {
      printf("%5.2f%s", B[col], (col < N-1) ? "; " : "]\n");
    }
  }
}

void print_X() {
  int row;

  if (N < 100) {
    printf("\nX = [");
    for (row = 0; row < N; row++) {
      printf("%5.2f%s", X[row], (row < N-1) ? "; " : "]\n");
    }
  }
}



int main(int argc, char **argv) {
  /* Timing variables */
  struct timeval etstart, etstop;  /* Elapsed times using gettimeofday() */
  struct timezone tzdummy;
  clock_t etstart2, etstop2;  /* Elapsed times using times() */
  unsigned long long usecstart, usecstop;
  struct tms cputstart, cputstop;  /* CPU times for my processes */
  int rank;
  
  /* Process program parameters */
  
  parameters(argc, argv);
  
 
  
  /*Starting MPI here*/
  
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);	/*store rank in rank for each process*/
  MPI_Comm_size(MPI_COMM_WORLD,&numberOfProcesses);		/*store number of processes*/


  
  int destination=0;	/*We can use i%numberOfProcesses to represent the rank of each process*/
  int source=0;		/*source, master process, rank is always 0*/
  MPI_Status status;
  
  if(rank==0){	/*Only the master rank initiate and store the matrix A and matrix B*/
	  /* Start Clock */
 
	printf("\nStarting clock.\n");
	gettimeofday(&etstart, &tzdummy);
	etstart2 = times(&cputstart);
  
  
	printf("Number of processes is %i\n",numberOfProcesses);
  
	/* Initialize A and B */
  
	initialize_inputs();

	/* Print input matrices */
 
	print_inputs();
  }
  else{	/*Other ranks do not initiate*/
	  
  }
  
  /* rank 0 sends row data to  other ranks, row[1], row[5], row[9]... to rank 1; 
   *										row[2], row[6], row[10]... to rank 2; 
   *										row[3], row[7], row [11]... to rank 3... and so on
											*/
  if(rank==0){
	  int i=0;
	  for(i=0;i<N;i++){
		  destination=i%numberOfProcesses;
		  if(destination!=0){	/*The data of rows of both matrices will be sent serially to each process according to their ranks*/
			MPI_Send(A[i],N,MPI_FLOAT,destination,i,MPI_COMM_WORLD);
			MPI_Send(&B[i],1,MPI_FLOAT,destination,i+N,MPI_COMM_WORLD);
		  }
	  }
  }
  else{
	  int i=0;
	  for(i=0;i<N;i++){
		  destination=i%numberOfProcesses;
		  if(destination==rank){	/*each process will receive its own data according to their ranks*/
			  MPI_Recv(A[i],N,MPI_FLOAT,source,i,MPI_COMM_WORLD,&status);
			   MPI_Recv(&B[i],1,MPI_FLOAT,source,i+N,MPI_COMM_WORLD,&status);
		  }
	  }
  }
 

  /* Gaussian Elimination */
  gauss(rank);
  
	

  /* Stop Clock */
  gettimeofday(&etstop, &tzdummy);
  etstop2 = times(&cputstop);
  printf("Stopped clock.\n");
  usecstart = (unsigned long long)etstart.tv_sec * 1000000 + etstart.tv_usec;
  usecstop = (unsigned long long)etstop.tv_sec * 1000000 + etstop.tv_usec;

  	if(rank==0){/*Just display one time result is enough*/
  /* Display output */
  print_X();

  /* Display timing results */

  
  printf("\nElapsed time = %g ms.\n",
	 (float)(usecstop - usecstart)/(float)1000);
	 
char  time_log[1000];

/*The time will be recorded to an output file each time the program is executed*/

sprintf(time_log,"\nElapsed time = %g ms.  number of processes is %d, matrix size= %d, random seed is %d\n",
	 (float)(usecstop - usecstart)/(float)1000,numberOfProcesses,N,seed);
	 FILE* file=fopen(output,"a");
	 fputs(time_log,file);
	 fclose(file);
  printf("(CPU times are accurate to the nearest %g ms)\n",
	 1.0/(float)CLOCKS_PER_SEC * 1000.0);
  printf("My total CPU time for parent = %g ms.\n",
	 (float)( (cputstop.tms_utime + cputstop.tms_stime) -
		  (cputstart.tms_utime + cputstart.tms_stime) ) /
	 (float)CLOCKS_PER_SEC * 1000);
  printf("My system CPU time for parent = %g ms.\n",
	 (float)(cputstop.tms_stime - cputstart.tms_stime) /
	 (float)CLOCKS_PER_SEC * 1000);
  printf("My total CPU time for child processes = %g ms.\n",
	 (float)( (cputstop.tms_cutime + cputstop.tms_cstime) -
		  (cputstart.tms_cutime + cputstart.tms_cstime) ) /
	 (float)CLOCKS_PER_SEC * 1000);
      /* Contrary to the man pages, this appears not to include the parent */
	printf("--------------------------------------------\n");}
  
  MPI_Finalize(); /*MPI shutting down*/
  
  exit(0);
}

/* ------------------ Above Was Provided --------------------- */

/****** You will replace this routine with your own parallel version *******/
/* Provided global variables are MAXN, N, A[][], B[], and X[],
 * defined in the beginning of this code.  X[] is initialized to zeros.
 */
void gauss(int my_rank) {
  int norm, row, col;  /* Normalization row, and zeroing
			* element row and col */
  float multiplier;
  int rank=my_rank; /*Declaration of an int to store rank*/

  

  /* Gaussian elimination */
  
  /*Again we will do the parallel in the row loop*/
  for(norm=0;norm<N-1;norm++){		/*Still we will not do parallel in the outer loop*/
		/*rank 0 boradcasts row[norm] to every other process when normalization line changed*/
		MPI_Bcast(A[norm], N, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&B[norm], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		for(row=norm+1;row<N;row++) {
			if (row%numberOfProcesses==rank&&rank!=0){	/*Each process is set to eliminate its own rows*/
				multiplier=A[row][norm]/A[norm][norm];
				for(col=norm;col<N;col++){
				A[row][col]-=A[norm][col]*multiplier;
				}
				B[row]-=B[norm]*multiplier;
				if(row==norm+1){	/* If the row just below the normalization line completed its gaussian elimination
									 * and the rank is not 0, the process will send the row data back to rank 0*/
					MPI_Send(A[row],N,MPI_FLOAT,0,row,MPI_COMM_WORLD);
					MPI_Send(&B[row],1,MPI_FLOAT,0,row+N,MPI_COMM_WORLD);
				}
			}
			if (rank==0&&row%numberOfProcesses!=0&&row==norm+1){	/* rank 0 receives the data from the non rank 0 process
																	 * if the row is just below the normalization line and completed
																	 * its elimination*/
				MPI_Recv(A[row],N,MPI_FLOAT,(norm+1)%numberOfProcesses,row,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(&B[row],1,MPI_FLOAT,(norm+1)%numberOfProcesses,row+N,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			if(rank==0&&row%numberOfProcesses==0){	/*rank 0 will also do its own job*/
				multiplier=A[row][norm]/A[norm][norm];
				for(col=norm;col<N;col++){
				A[row][col]-=A[norm][col]*multiplier;
				}
				B[row]-=B[norm]*multiplier;
			}
			
		}
	}

  
	
  
  
  /*for (norm = 0; norm < N - 1; norm++) {
    for (row = norm + 1; row < N; row++) {
      multiplier = A[row][norm] / A[norm][norm];
      for (col = norm; col < N; col++) {
	A[row][col] -= A[norm][col] * multiplier;
      }
      B[row] -= B[norm] * multiplier;
    }
  }*/
  /* (Diagonal elements are not normalized to 1.  This is treated in back
   * substitution.)
   */


  /* Back substitution */
  
  for (row = N - 1; row >= 0; row--) {
    X[row] = B[row];
    for (col = N-1; col > row; col--) {
      X[row] -= A[row][col] * X[col];
    }
  X[row] /= A[row][row];}
  
}

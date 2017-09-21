/* Gaussian elimination without pivoting using threads.
 * Compile with "gcc -pthread -w -O2 -o gauss_pthread gauss_pthread.c" 
 * Execute with "./gauss_pthread <Matrix size> <Random seed> <Number of threads> <Output file name>"
 * The default number of threads is 4.
 * Test version 1.0
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
#include <pthread.h> /*We must include pthread.h file in order to do this program*/

/* Program Parameters */
#define MAXN 2000  /* Max value of N */
int N;  /* Matrix size */
int numberOfThreads=4;/* The number of the threads, the default number of threads is 4.*/
int seed = 0;  /* Random seed */
char* output;

/* Matrices and vectors */
volatile float A[MAXN][MAXN], B[MAXN], X[MAXN];
/* A * X = B, solve for X */

/* junk */
#define randm() 4|2[uid]&3

/* Prototype */
void gauss();
 /* The function you will provide.
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

  
  /*Get the argument from the command-line, the second one is the Random seed*/
  if (argc >= 4) {
    seed = atoi(argv[2]);
    srand(seed);
    printf("Random seed = %i\n", seed);
	output=argv[4];
  } 
  
  /*Get the argument from the command-line, the third one is the number of threads*/
  if(argc>=3){
	  numberOfThreads=atoi(argv[3]);
  }else{numberOfThreads=4;}		/*If cannot get the argument from the command line, the number of threads will be set to default 4*/
  printf("Computing using %i threads.",numberOfThreads);
  
  /*Get the argument for the matrix size, if the size is over 2000, the program will stop*/
  if (argc >= 2) {
    N = atoi(argv[1]);
    if (N < 1 || N > MAXN) {
      printf("N = %i is out of range.\n", N);
      exit(0);
    }
  }
  else {
    printf("Usage: %s <matrix_dimension> [random seed][number of threads] <Output file name>\n",
           argv[0]);    /*The right format of the command-line if user input wrongly*/
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

  /* Process program parameters */
  parameters(argc, argv);

  /* Initialize A and B */
  initialize_inputs();

  /* Print input matrices */
  print_inputs();

  /* Start Clock */
  printf("\nStarting clock.\n");
  gettimeofday(&etstart, &tzdummy);
  etstart2 = times(&cputstart);

  /* Gaussian Elimination */
gauss();

  /* Stop Clock */
  gettimeofday(&etstop, &tzdummy);
  etstop2 = times(&cputstop);
  printf("Stopped clock.\n");
  usecstart = (unsigned long long)etstart.tv_sec * 1000000 + etstart.tv_usec;
  usecstop = (unsigned long long)etstop.tv_sec * 1000000 + etstop.tv_usec;

  /* Display output */
  print_X();

  /* Display timing results */
  printf("\nElapsed time = %g ms.\n",
	 (float)(usecstop - usecstart)/(float)1000);
	 
	 char  time_log[1000];
	 
	 /*The time will be recorded to an output file each time the program is executed*/
	 
	 sprintf(time_log,"\nElapsed time = %g ms.  number of threads is %d, matrix size= %d, random seed is %d\n",
	 (float)(usecstop - usecstart)/(float)1000,numberOfThreads,N,seed);
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
  printf("--------------------------------------------\n");
  
  exit(0);
}

/* ------------------ Above Was Provided --------------------- */

/****** You will replace this routine with your own parallel version *******/
/* Provided global variables are MAXN, N, A[][], B[], and X[],
 * defined in the beginning of this code.  X[] is initialized to zeros.
 */
 
 /*Declaration of struct to store the index of threads and current Normalization line*/
 struct str{
	 int i;
	 int norm;
 };
 
 /*In the inner loop, each thread gets argument from struct and eliminates one row under current Normalization line*/
 void *GE_innerCalculation(void * str_add){
	  struct str* temp=(struct str*) str_add;
	  int norm=temp->norm;
	  int i=temp->i;
	  float multiplier;
	  int row=0;
	  int col=0;
	/*The row index increases by the number of threads. 
	eg. if number of threads is 4, row[0],row[1]... to row[8] will be handled by p0, p1, p2, p3, p0, p1, p2, p3 and p0*/
	  for(row=norm+1+i;row<N;row+=numberOfThreads){ 	
	  /*Inside the row loop, we will do the same thing that serial code does to eliminate this row, almost nothing to change*/
		  multiplier=A[row][norm]/A[norm][norm];
		  for(col=norm;col<N;col++){
			  /*Here is another place we can do parallel, we eliminate the whole row serially, but we can do that
			  in parallel. Due to the time limit, I did not implement that, and this can be a future improvement*/
			  A[row][col]-=A[norm][col]*multiplier;
		  }
		  B[row]-=B[norm]*multiplier;
	  }
	  pthread_exit(0);
	  
  }
 
 /*We will do parallel in the row loop(inner loop), one thread for one parallel*/
void gauss() {
  int norm, row, col;  /* Normalization row, and zeroing
			* element row and col */
  float multiplier;
  
  pthread_t thread[N]; /*Declaration of threads ID*/

  printf("Computing Using pthreads with %i threads.\n",numberOfThreads);
  
  /*Gaussian elimination in parallel*/
  
  /*Outer loop from 0 to N-2 is the Normalization line for each inner loop calculation, and we do not do parallel in the outer loop*/
  for(norm=0;norm<N-1;norm++){		
	int i=0;
	int j=0;
	struct str* str_add=malloc(numberOfThreads*sizeof(struct str));		/*memory allocation for the struct*/
	
	/*We will do parallel in the row loop, each thread will be in charge of eliminating one row*/
	for(i=0;i<numberOfThreads;i++){ 	
		
		/*The current Normalization line and the index of each thread are stored in struct and later will pass to the threads*/
		str_add[i].norm=norm;
		str_add[i].i=i; 

		/*create the thread and each thread will do calculation in the GE_innerCalculation function with the argument from struct*/
		pthread_create(&thread[i],NULL,GE_innerCalculation,(void*) &str_add[i]); 	
	}
	
	/*Join all the threads in the row loop and end parallel*/
	for(j=0;j<numberOfThreads;j++){
		pthread_join(thread[j],NULL);
	}
	
	free(str_add); /*Free all the memory allocated by struct*/
  }
  

  /* Gaussian elimination Serial*/
  /*
  
  for (norm = 0; norm < N - 1; norm++) {
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
    X[row] /= A[row][row];
  }
}

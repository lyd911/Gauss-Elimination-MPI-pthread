The gauss.c file is a modified version of the orginal serial code, it will take another arument <file name> as the
output file.

The gauss_output, gauss_pthread_output and gauss_mpi_output are three output files I used to record the elapsed time
in matrix calcualtion from N=200 to N=2000.

For Serial code:
Compile with "gcc -O2 -o gauss gauss.c" 
Execute with "./gauss <Matrix size> <Random seed> <Output file name>"

For pthread code:
Compile with "gcc -pthread -w -O2 -o gauss_pthread gauss_pthread.c" 
Execute with "./gauss_pthread <Matrix size> <Random seed> <Number of threads> <Output file name>"

For MPI code:
Compile with "mpicc -w -O2 -o gauss_mpi gauss_mpi.c"
Execute with "mpiexec -n <Number of processes> ./gauss_mpi <Matrix Size> <Random Seed> <Output file name>"

Thank you,
Yedong Liu
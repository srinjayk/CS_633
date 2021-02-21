To run src.c

do the following

mpicc -o src src.c -lm

mpirun -np 64 -hosts csews20,csews23,csews24,csews25,csews26,csews27,csews28,csews29,csews30,csews31 ./src 1024*1024 50


output is of the form
time1 - time taken by MPI_Send
time2 - time taken by MPI_Pack
time3 - time taken by MPI_Type

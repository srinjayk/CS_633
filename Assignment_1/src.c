#include <stdio.h>
#include "mpi.h"

int main(int argc, char *argv[]) 
{
  // initialize MPI
  MPI_Init (&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank) ;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int N = atoi (argv[1]);
  int column = atoi (argv[2]);
  int count = atoi (argv[3]);
  int blocklen = atoi (argv[4]);
  int stride = atoi (argv[5]);
  int data[N][N], received[N*blocklen]; 

  for(int i=0;i<N;i++){
  	for(int j=0;j<N;j++){
  		data[i][j] = i*N + j;
  	}
  }

  // ---------------------------------------------------------------

  // Part II

  int buffer_lr[N];
  int buffer_fr[N];
  int buffer_lc[N];
  int buffer_fc[N];

  for(int i=0;i<N;i++){
  	// fill up all the parameters
  	MPI_Pack(&data[0][i]);
  }

  MPI_Send();

  for(int i=0;i<N;i++){
  	// fill up all the parameters
  	MPI_Pack(&data[i][0]);
  }

  MPI_Send();

  for(int i=0;i<N,i++){
  	// fill up all the parameters
  	MPI_Pack(&data[N-1][i])
  }

  MPI_Send();

  for(int i=0;i<N;i++){
  	// fil up all the parameters
  	MPI_Pack(&data[i][N-1]);
  }

  MPI_Send();

  MPI_Recv();

  


  // ---------------------------------------------------------------

  // ---------------------------------------------------------------
  
  // Part III
  // for the third part when we have to from the new datatypes

  // we would be defining 4 new datatypes
  // last_row => last row of any matrix
  // first_row => first row of any matrix
  // last_column => last column of any matrix
  // first_column => first column of any matrix

  MPI_Datatype last_row;
  MPI_Datatype first_row;
  MPI_Datatype last_column;
  MPI_Datatype first_column;

  int column_lr = N*(N-1);
  int column_fr = 0;
  int column_lc = N-1;
  int column_fc = 0;

  int count_lr = 1;
  int count_fr = 1;
  int count_lc = N;
  int count_fc = N;

  int blocklen_lr = N;
  int blocklen_fr = N;
  int blocklen_lc = 1;
  int blocklen_fc = 1;

  int stride_lr = 0;
  int stride_fr = 0;
  int stride_lc = N;
  int stride_fc = N;

  MPI_Type_vector (count_lr, blocklen_lr, stride_lr, MPI_INT, &last_row);
  MPI_Type_commit (&last_row);

  MPI_Type_vector (count_fr, blocklen_fr, stride_fr, MPI_INT, &first_row);
  MPI_Type_commit (&first_row);

  MPI_Type_vector (count_fc, blocklen_fc, stride_fc, MPI_INT, &first_column);
  MPI_Type_commit (&first_column);

  MPI_Type_vector (count_lc, blocklen_lc, stride_lc, MPI_INT, &last_column);
  MPI_Type_commit (&last_column);

  // ---------------------------------------------------------------------

  printf("Rank of this process is %d\n", myrank);
  printf("Totalnumber of processes is %d\n",size);

  printf ("Hello, world!\n");

  // done with MPI
  MPI_Finalize();
}


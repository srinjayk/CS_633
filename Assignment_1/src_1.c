#include <stdio.h>
#include <math.h>
#include "mpi.h"

int main(int argc, char *argv[]){

  // initialize MPI
  MPI_Init (&argc, &argv);

  // number of data points per process
  int N_sq = atoi(&argv[1]);
  // number of row and columns 
  int N = (int)(sqrt(N_sq));

  // time steps
  int time_step = atoi(&argv[2]);

  // variable to capture the position of process in the cartesian system
  int indicator = 0;

  int data[N][N];

  for(int i=0;i<N;i++){
  	for(int j=0;j<N;j++){
  		data[i][j] = i*N + j;
  	}
  }

  // get the rank and size of the current process
  int myrank, size_sq;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank) ;
  MPI_Comm_size(MPI_COMM_WORLD, &size_sq);

  int size = (int)(sqrt(size_sq));

  // row for the current process
  int row_p = myrank/size;
  // column for the current pocess
  int col_p = myrank%size;

  // upper left corner
  if((row_p == 0) && (col_p == 0)){
    indicator = 1;
  }
  // botton right corner
  else if((row_p == (size-1)) && (col_p == (size-1))){
    indicator = 2;
  }
  // upper right corner
  else if((row_p == 0) && (col_p == (size-1))){
    indicator = 3;
  }
  // bottom left corner
  else if((row_p == (size-1)) && (col_p == 0)){
    indicator = 4;
  }
  // top row
  else if(row_p == 0){
    indicator = 5;
  }
  // left column
  else if(col_p == 0){
    indicator = 6;
  }
  // bottom row
  else if(row_p == (size-1)){
    indicator = 7;
  }
  // right column
  else if(col_p == (size-1)){
    indicator = 8;
  }
  // all the interior processes
  else{
    indicator = 9;
  }


  // ---------------------------------------------------------------

  // Part I

  // Receiving the data using different MPI_Recv's 
  if(indicator == 1){
    for(int i=0;i<N;i++){
      MPI_Recv(&data[i][N-1]);
    }
    for(int i=0;i<N;i++){
      MPI_Recv(&data[N-1][i]);
    }
  }
  else if(indicator == 2){
    for(int i=0;i<N;i++){
      MPI_Recv(&data[i][0]);
    }
    for(int i=0;i<N;i++){
      MPI_Recv(&data[0][i]);
    }
  }
  else if(indicator == 3){
    for(int i=0;i<N;i++){
      MPI_Recv(&data[i][0]);
    }
    for(int i=0;i<N;i++){
      MPI_Recv(&data[N-1][i]);
    }
  }
  else if(indicator == 4){
    for(int i=0;i<N;i++){
      MPI_Recv(&data[0][i]);
    }
    for(int i=0;i<N;i++){
      MPI_Recv(&data[i][N-1]);
    }
  }
  else if(indicator == 5){
    for(int i=0;i<N;i++){
      MPI_Recv(&data[i][N-1]);
    }
    for(int i=0;i<N;i++){
      MPI_Recv(&data[N-1][i]);
    }
    for(int i=0;i<N;i++){
      MPI_Recv(&data[i][0]);
    }
  }
  else if(indicator == 6){
    for(int i=0;i<N;i++){
      MPI_Recv(&data[i][N-1]);
    }
    for(int i=0;i<N;i++){
      MPI_Recv(&data[N-1][i]);
    }
    for(int i=0;i<N;i++){
      MPI_Recv(&data[0][i]);
    }
  }
  else if(indicator == 7){
    for(int i=0;i<N;i++){
      MPI_Recv(&data[i][N-1]);
    }
    for(int i=0;i<N;i++){
      MPI_Recv(&data[i][0]);
    }
    for(int i=0;i<N;i++){
      MPI_Recv(&data[0][i]);
    }
  }
  else if(indicator == 8){
    for(int i=0;i<N;i++){
      MPI_Recv(&data[N-1][i]);
    }
    for(int i=0;i<N;i++){
      MPI_Recv(&data[i][0]);
    }
    for(int i=0;i<N;i++){
      MPI_Recv(&data[0][i]);
    }
  }
  else if(indicator == 9){
    for(int i=0;i<N;i++){
      MPI_Recv(&data[i][N-1]);
    }
    for(int i=0;i<N;i++){
      MPI_Recv(&data[N-1][i]);
    }
    for(int i=0;i<N;i++){
      MPI_Recv(&data[i][0]);
    }
    for(int i=0;i<N;i++){
      MPI_Recv(&data[0][i]);
    }
  }

  // Sending the data using different MPI_Send's
  if(indicator == 1){
    for(int i=0;i<N;i++){
      MPI_Send(&data[i][N-1]);
    }
    for(int i=0;i<N;i++){
      MPI_Send(&data[N-1][i]);
    }
  }
  else if(indicator == 2){
    for(int i=0;i<N;i++){
      MPI_Send(&data[i][0]);
    }
    for(int i=0;i<N;i++){
      MPI_Send(&data[0][i]);
    }
  }
  else if(indicator == 3){
    for(int i=0;i<N;i++){
      MPI_Send(&data[i][0]);
    }
    for(int i=0;i<N;i++){
      MPI_Send(&data[N-1][i]);
    }
  }
  else if(indicator == 4){
    for(int i=0;i<N;i++){
      MPI_Send(&data[0][i]);
    }
    for(int i=0;i<N;i++){
      MPI_Send(&data[i][N-1]);
    }
  }
  else if(indicator == 5){
    for(int i=0;i<N;i++){
      MPI_Send(&data[i][N-1]);
    }
    for(int i=0;i<N;i++){
      MPI_Send(&data[N-1][i]);
    }
    for(int i=0;i<N;i++){
      MPI_Send(&data[i][0]);
    }
  }
  else if(indicator == 6){
    for(int i=0;i<N;i++){
      MPI_Send(&data[i][N-1]);
    }
    for(int i=0;i<N;i++){
      MPI_Send(&data[N-1][i]);
    }
    for(int i=0;i<N;i++){
      MPI_Send(&data[0][i]);
    }
  }
  else if(indicator == 7){
    for(int i=0;i<N;i++){
      MPI_Send(&data[i][N-1]);
    }
    for(int i=0;i<N;i++){
      MPI_Send(&data[i][0]);
    }
    for(int i=0;i<N;i++){
      MPI_Send(&data[0][i]);
    }
  }
  else if(indicator == 8){
    for(int i=0;i<N;i++){
      MPI_Send(&data[N-1][i]);
    }
    for(int i=0;i<N;i++){
      MPI_Send(&data[i][0]);
    }
    for(int i=0;i<N;i++){
      MPI_Send(&data[0][i]);
    }
  }
  else if(indicator == 9){
    for(int i=0;i<N;i++){
      MPI_Send(&data[i][N-1]);
    }
    for(int i=0;i<N;i++){
      MPI_Send(&data[N-1][i]);
    }
    for(int i=0;i<N;i++){
      MPI_Send(&data[i][0]);
    }
    for(int i=0;i<N;i++){
      MPI_Send(&data[0][i]);
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

  // MPI_Recv();

  if(indicator == 1){
  	MPI_Recv(&buffer_lr);
  	MPI_Recv(&buffer_lc);
  }
  else if(indicator == 2){
  	MPI_Recv(&buffer_fc);
  	MPI_Recv(&buffer_fr);
  }
  else if(indicator == 3){
  	MPI_Recv(&buffer_fc);
  	MPI_Recv(&buffer_lr);
  }
  else if(indicator == 4){
  	MPI_Recv(&buffer_lc);
  	MPI_Recv(&buffer_fr);
  }
  else if(indicator == 5){
  	MPI_Recv(&buffer_lr);
  	MPI_Recv(&buffer_lc);
  	MPI_Recv(&buffer_fc);
  }
  else if(indicator == 6){
  	MPI_Recv(&buffer_lr);
  	MPI_Recv(&buffer_lc);
  	MPI_Recv(&buffer_fr);
  }
  else if(indicator == 7){
  	MPI_Recv(&buffer_lc);
  	MPI_Recv(&buffer_fc);
  	MPI_Recv(&buffer_fr);
  }
  else if(indicator == 8){
  	MPI_Recv(&buffer_lr);
  	MPI_Recv(&buffer_fc);
  	MPI_Recv(&buffer_fr);
  }
  else if(indicator == 9){
  	MPI_Recv(&buffer_lr);
  	MPI_Recv(&buffer_lc);
  	MPI_Recv(&buffer_fc);
  	MPI_Recv(&buffer_fr);
  }

  // MPI_Send();

  if(indicator == 1){
  	MPI_Send(&buffer_lr);
  	MPI_Send(&buffer_lc);
  }
  else if(indicator == 2){
  	MPI_Send(&buffer_fc);
  	MPI_Send(&buffer_fr);
  }
  else if(indicator == 3){
  	MPI_Send(&buffer_fc);
  	MPI_Send(&buffer_lr);
  }
  else if(indicator == 4){
  	MPI_Send(&buffer_lc);
  	MPI_Send(&buffer_fr);
  }
  else if(indicator == 5){
  	MPI_Send(&buffer_lr);
  	MPI_Send(&buffer_lc);
  	MPI_Send(&buffer_fc);
  }
  else if(indicator == 6){
  	MPI_Send(&buffer_lr);
  	MPI_Send(&buffer_lc);
  	MPI_Send(&buffer_fr);
  }
  else if(indicator == 7){
  	MPI_Send(&buffer_lc);
  	MPI_Send(&buffer_fc);
  	MPI_Send(&buffer_fr);
  }
  else if(indicator == 8){
  	MPI_Send(&buffer_lr);
  	MPI_Send(&buffer_fc);
  	MPI_Send(&buffer_fr);
  }
  else if(indicator == 9){
  	MPI_Send(&buffer_lr);
  	MPI_Send(&buffer_lc);
  	MPI_Send(&buffer_fc);
  	MPI_Send(&buffer_fr);
  }


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

  int recv_buffer_lr[N];
  int recv_buffer_lc[N];
  int recv_buffer_fr[N];
  int recv_buffer_fc[N];

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

  if(indicator == 1){
  	MPI_Recv(&recv_buffer_lr, 1, last_row, 1, 99, MPI_COMM_WORLD);
  	MPI_Recv(&recv_buffer_lc, 1, last_column, 1, 99, MPI_COMM_WORLD);
  }
  else if(indicator == 2){
  	MPI_Recv(&recv_buffer_fr, 1, first_row, 1, 99, MPI_COMM_WORLD);
  	MPI_Recv(&recv_buffer_fc, 1, first_column, 1, 99, MPI_COMM_WORLD);
  }
  else if(indicator == 3){
  	MPI_Recv(&recv_buffer_fc, 1, first_column, 1, 99, MPI_COMM_WORLD);
  	MPI_Recv(&recv_buffer_lr, 1, last_row, 1, 99, MPI_COMM_WORLD);
  }
  else if(indicator == 4){
  	MPI_Recv(&recv_buffer_fr, 1, first_row, 1, 99, MPI_COMM_WORLD);
  	MPI_Recv(&recv_buffer_lc, 1, last_column, 1, 99, MPI_COMM_WORLD);
  }
  else if(indicator == 5){
  	MPI_Recv(&recv_buffer_lc, 1, last_column, 1, 99, MPI_COMM_WORLD);
  	MPI_Recv(&recv_buffer_fc, 1, first_column, 1, 99, MPI_COMM_WORLD);
  	MPI_Recv(&recv_buffer_lr, 1, last_row, 1, 99, MPI_COMM_WORLD);
  }
  else if(indicator == 6){
  	MPI_Recv(&recv_buffer_fr, 1, first_row, 1, 99, MPI_COMM_WORLD);
  	MPI_Recv(&recv_buffer_lc, 1, last_column, 1, 99, MPI_COMM_WORLD);
  	MPI_Recv(&recv_buffer_lr, 1, last_row, 1, 99, MPI_COMM_WORLD);
  }
  else if(indicator == 7){
  	MPI_Recv(&recv_buffer_fr, 1, first_row, 1, 99, MPI_COMM_WORLD);
  	MPI_Recv(&recv_buffer_lc, 1, last_column, 1, 99, MPI_COMM_WORLD);
  	MPI_Recv(&recv_buffer_fc, 1, first_column, 1, 99, MPI_COMM_WORLD);
  }
  else if(indicator == 8){
  	MPI_Recv(&recv_buffer_fr, 1, first_row, 1, 99, MPI_COMM_WORLD);
  	MPI_Recv(&recv_buffer_fc, 1, first_column, 1, 99, MPI_COMM_WORLD);
  	MPI_Recv(&recv_buffer_lr, 1, last_row, 1, 99, MPI_COMM_WORLD);
  }
  else if(indicator == 9){
  	MPI_Recv(&recv_buffer_fr, 1, first_row, 1, 99, MPI_COMM_WORLD);
  	MPI_Recv(&recv_buffer_lc, 1, last_column, 1, 99, MPI_COMM_WORLD);
  	MPI_Recv(&recv_buffer_fc, 1, first_column, 1, 99, MPI_COMM_WORLD);
  	MPI_Recv(&recv_buffer_lr, 1, last_row, 1, 99, MPI_COMM_WORLD);
  }




  if(indicator == 1){
  	MPI_Send(&data[0][column_lr], 1, last_row, 1, 99, MPI_COMM_WORLD);
  	MPI_Send(&data[0][column_lc], 1, last_column, 1, 99, MPI_COMM_WORLD);
  }
  else if(indicator == 2){
  	MPI_Send(&data[0][column_fr], 1, first_row, 1, 99, MPI_COMM_WORLD);
  	MPI_Send(&data[0][column_fc], 1, first_column, 1, 99, MPI_COMM_WORLD);
  }
  else if(indicator == 3){
  	MPI_Send(&data[0][column_fc], 1, first_column, 1, 99, MPI_COMM_WORLD);
  	MPI_Send(&data[0][column_lr], 1, last_row, 1, 99, MPI_COMM_WORLD);
  }
  else if(indicator == 4){
  	MPI_Send(&data[0][column_fr], 1, first_row, 1, 99, MPI_COMM_WORLD);
  	MPI_Send(&data[0][column_lc], 1, last_column, 1, 99, MPI_COMM_WORLD);
  }
  else if(indicator == 5){
  	MPI_Send(&data[0][column_lc], 1, last_column, 1, 99, MPI_COMM_WORLD);
  	MPI_Send(&data[0][column_fc], 1, first_column, 1, 99, MPI_COMM_WORLD);
  	MPI_Send(&data[0][column_lr], 1, last_row, 1, 99, MPI_COMM_WORLD);
  }
  else if(indicator == 6){
  	MPI_Send(&data[0][column_fr], 1, first_row, 1, 99, MPI_COMM_WORLD);
  	MPI_Send(&data[0][column_lc], 1, last_column, 1, 99, MPI_COMM_WORLD);
  	MPI_Send(&data[0][column_lr], 1, last_row, 1, 99, MPI_COMM_WORLD);
  }
  else if(indicator == 7){
  	MPI_Send(&data[0][column_fr], 1, first_row, 1, 99, MPI_COMM_WORLD);
  	MPI_Send(&data[0][column_lc], 1, last_column, 1, 99, MPI_COMM_WORLD);
  	MPI_Send(&data[0][column_fc], 1, first_column, 1, 99, MPI_COMM_WORLD);
  }
  else if(indicator == 8){
  	MPI_Send(&data[0][column_fr], 1, first_row, 1, 99, MPI_COMM_WORLD);
  	MPI_Send(&data[0][column_fc], 1, first_column, 1, 99, MPI_COMM_WORLD);
  	MPI_Send(&data[0][column_lr], 1, last_row, 1, 99, MPI_COMM_WORLD);
  }
  else if(indicator == 9){
  	MPI_Send(&data[0][column_fr], 1, first_row, 1, 99, MPI_COMM_WORLD);
  	MPI_Send(&data[0][column_lc], 1, last_column, 1, 99, MPI_COMM_WORLD);
  	MPI_Send(&data[0][column_fc], 1, first_column, 1, 99, MPI_COMM_WORLD);
  	MPI_Send(&data[0][column_lr], 1, last_row, 1, 99, MPI_COMM_WORLD);
  }

  // ---------------------------------------------------------------------

  printf("Rank of this process is %d\n", myrank);
  printf("Totalnumber of processes is %d\n",size);

  printf ("Hello, world!\n");

  // done with MPI
  MPI_Finalize();
}


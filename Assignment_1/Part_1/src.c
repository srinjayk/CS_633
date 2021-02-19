#include <stdio.h>
#include <math.h>
#include "mpi.h"

int main(int argc, char *argv[]){

  // initialize MPI
  MPI_Init (&argc, &argv);

  // number of data points per process
  int N_sq = atoi(argv[1]);
  // number of row and columns 
  int N = (int)(sqrt(N_sq));

  // time steps
  int time_step = atoi(argv[2]);
  
  // for printing the given rank 
  // for running the given option
  int option = atoi(argv[3]);
  // variable to capture the position of process in the cartesian system
  int indicator = 0;

  double data[N][N];

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

  printf("Myrank:%d Row:%d Column:%d Indicator:%d\n",myrank,row_p,col_p,indicator);
  // ---------------------------------------------------------------
  
  // Part I
  // we would encoe the postion and iteration togeher in the TAG
  // in the format (length of points in row)*iteration + position

  if(1){

    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
        data[i][j] = 1+i+j;
      }
    }

    for(int i=0;i<time_step;i++){

      // for process 1
      // int MPI_Send (const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) 
      if(1){

          // lr
          if((indicator == 1) || (indicator == 3) || (indicator == 5) || (indicator == 6) || (indicator == 8) || (indicator == 9)){
            for(int j=0;j<N;j++){
              MPI_Send(&data[N-1][j], 1, MPI_DOUBLE, myrank+size, N*i+j, MPI_COMM_WORLD);
            }
          }

          // fr
          if((indicator == 2) || (indicator == 4) || (indicator == 6) || (indicator == 7) || (indicator == 8) || (indicator == 9)){
            for(int j=0;j<N;j++){
              MPI_Send(&data[0][j], 1, MPI_DOUBLE, myrank-size, N*i+j, MPI_COMM_WORLD);
            }
          }

          // lc
          if((indicator == 1) || (indicator == 4) || (indicator == 5) || (indicator == 6) || (indicator == 7) || (indicator == 9)){
            for(int j=0;j<N;j++){
              MPI_Send(&data[j][N-1], 1, MPI_DOUBLE, myrank+1, N*i+j, MPI_COMM_WORLD);
            }
          }

          // fc
          if((indicator == 2) || (indicator == 3) || (indicator == 5) || (indicator == 7) || (indicator == 8) || (indicator == 9)){
            for(int j=0;j<N;j++){
              MPI_Send(&data[j][0], 1, MPI_DOUBLE, myrank-1, N*i+j, MPI_COMM_WORLD);
            }
          }

      }


      int ind_lr = 0;
      int ind_fr = 0;
      int ind_fc = 0;
      int ind_lc = 0;

      double recv_buffer_fc[N];
      double recv_buffer_fr[N];
      double recv_buffer_lr[N];
      double recv_buffer_lc[N];

        // int MPI_Recv (void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
        MPI_Status status;
        if(1){

          // lr
          if((indicator == 1) || (indicator == 3) || (indicator == 5) || (indicator == 6) || (indicator == 8) || (indicator == 9)){
            for(int j=0;j<N;j++){
              MPI_Recv(&recv_buffer_lr[j], 1, MPI_DOUBLE, myrank+size, N*i+j, MPI_COMM_WORLD, &status);
            }
            ind_lr = 1;
          }

          // fr
          if((indicator == 2) || (indicator == 4) || (indicator == 6) || (indicator == 7) || (indicator == 8) || (indicator == 9)){
            for(int j=0;j<N;j++){
              MPI_Recv(&recv_buffer_fr[j], 1, MPI_DOUBLE, myrank-size, N*i+j, MPI_COMM_WORLD, &status);
            }
            ind_fr = 1;
          }

          // lc
          if((indicator == 1) || (indicator == 4) || (indicator == 5) || (indicator == 6) || (indicator == 7) || (indicator == 9)){
            for(int j=0;j<N;j++){
              MPI_Recv(&recv_buffer_lc[j], 1, MPI_DOUBLE, myrank+1, N*i+j, MPI_COMM_WORLD, &status);
            }
            ind_lc = 1;
          }

          // fc
          if((indicator == 2) || (indicator == 3) || (indicator == 5) || (indicator == 7) || (indicator == 8) || (indicator == 9)){
            for(int j=0;j<N;j++){
              MPI_Recv(&recv_buffer_fc[j], 1, MPI_DOUBLE, myrank-1, N*i+j, MPI_COMM_WORLD, &status);
            }
            ind_fc = 1;
          }

        }

        // compute


        double data_temp[N+2][N+2];
        data_temp[0][0] = 0;
        data_temp[N+1][N+1] = 0;
        data_temp[0][N+1] = 0;
        data_temp[N+1][0] = 0;

        for(int j=1;j<N+1;j++){
          for(int k=1;k<N+1;k++){
            data_temp[j][k] = data[j-1][k-1];
          }
        }

        
        for(int j=1;j<N+1;j++){
          //lr
          if(ind_lr != 0){
            data_temp[N+1][j] = recv_buffer_lr[j-1];
          }
          else{
            data_temp[N+1][j] = 0;
          }

          // lc
          if(ind_lc != 0){
            data_temp[j][N+1] = recv_buffer_lc[j-1];
          }
          else{
            data_temp[j][N+1] = 0;
          }

          // fc
          if(ind_fc != 0){
            data_temp[j][0] = recv_buffer_fc[j-1];
          }
          else{
            data_temp[j][0] = 0;
          }

          // fr
          if(ind_fr != 0){
            data_temp[0][j] = recv_buffer_fr[j-1];
          }
          else{
            data_temp[0][j] = 0;
          }

        }

        for(int j=1;j<N+1;j++){
          for(int k=1;k<N+1;k++){
            int count = 0;
            double sum = 0;

            if(data_temp[j+1][k] != 0){
              count++;
              sum = sum + data_temp[j+1][k];
            }

            if(data_temp[j-1][k] != 0){
              count++;
              sum = sum + data_temp[j-1][k];
            }

            if(data_temp[j][k-1] != 0){
              count++;
              sum = sum + data_temp[j][k-1];
            }

            if(data_temp[j][k+1] != 0){
              count++;
              sum = sum + data_temp[j][k+1];
            }
            
            data[j-1][k-1] = sum/count;

          }
        }

    }
  }

  // ---------------------------------------------------------------------

  if(myrank == option){
    printf("Rank of this process is %d\n", myrank);
    printf("Totalnumber of processes is %d\n",size);
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
        printf("%f ", data[i][j]);
      }
      printf("\n");
    }
  }

  // done with MPI
  MPI_Finalize();
}

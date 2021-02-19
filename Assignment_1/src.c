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
  
  // for running the given process either 1 or 2 or 3
  int option = atoi(argv[3]);

  // for printing the given rank
  int req_rank = atoi(argv[4]);

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

  if(option == 1){

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


  // ---------------------------------------------------------------

  // Part II

  if(option == 2){

    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
        data[i][j] = 1+i+j;
      }
    }

    for(int l=0;l<time_step;l++){
      // int MPI_Pack (const void *inbuf, int incount, MPI_Datatype datatype, void *outbuf, int outsize, int *position, MPI_Comm comm)
      // int MPI_Unpack (const void *inbuf, int insize, int *position, void *outbuf, int outcount, MPI_Datatype datatype, MPI_Comm comm)

      double send_buffer_lr[N];
      double send_buffer_fr[N];
      double send_buffer_lc[N];
      double send_buffer_fc[N];

      double recv_buffer_lr[N];
      double recv_buffer_fr[N];
      double recv_buffer_lc[N];
      double recv_buffer_fc[N];

      int position = 0;

      if(1){

          // lr
          if((indicator == 1) || (indicator == 3) || (indicator == 5) || (indicator == 6) || (indicator == 8) || (indicator == 9)){
            
            for(int i=0;i<N;i++){
              MPI_Pack(&data[N-1][i], 1, MPI_DOUBLE,send_buffer_lr,8*N,&position, MPI_COMM_WORLD);
            }

            MPI_Send(send_buffer_lr, position, MPI_PACKED, myrank+size, l, MPI_COMM_WORLD);
          }

          position = 0;

          // fr
          if((indicator == 2) || (indicator == 4) || (indicator == 6) || (indicator == 7) || (indicator == 8) || (indicator == 9)){
            
            // fr first row
            for(int i=0;i<N;i++){
              MPI_Pack(&data[0][i], 1, MPI_DOUBLE,send_buffer_fr,8*N,&position, MPI_COMM_WORLD);
            }

            MPI_Send(send_buffer_fr, position, MPI_PACKED, myrank-size, l, MPI_COMM_WORLD);
          }

          position = 0;

          // lc
          if((indicator == 1) || (indicator == 4) || (indicator == 5) || (indicator == 6) || (indicator == 7) || (indicator == 9)){
            
            for(int i=0;i<N;i++){
              MPI_Pack(&data[i][N-1], 1, MPI_DOUBLE,send_buffer_lc,8*N,&position, MPI_COMM_WORLD);
            }
            MPI_Send(send_buffer_lc, position, MPI_PACKED, myrank+1, l, MPI_COMM_WORLD);
          }

          position = 0;

          // fc
          if((indicator == 2) || (indicator == 3) || (indicator == 5) || (indicator == 7) || (indicator == 8) || (indicator == 9)){
            
            for(int i=0;i<N;i++){
              MPI_Pack(&data[i][0], 1, MPI_DOUBLE,send_buffer_fc,8*N,&position, MPI_COMM_WORLD);
            }

            MPI_Send(send_buffer_fc, position, MPI_PACKED, myrank-1, l, MPI_COMM_WORLD);
          }

      }


      int ind_lr = 0;
      int ind_fr = 0;
      int ind_fc = 0;
      int ind_lc = 0;

      // int MPI_Recv (void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
      MPI_Status status;

          // lr
          if((indicator == 1) || (indicator == 3) || (indicator == 5) || (indicator == 6) || (indicator == 8) || (indicator == 9)){
            MPI_Recv(recv_buffer_lr, 8*N, MPI_PACKED, myrank+size, l, MPI_COMM_WORLD, &status);
            // MPI_Get_count (&status, MPI_PACKED, &count);
            ind_lr = 1;
          }

          // fr
          if((indicator == 2) || (indicator == 4) || (indicator == 6) || (indicator == 7) || (indicator == 8) || (indicator == 9)){
            MPI_Recv(recv_buffer_fr, 8*N, MPI_PACKED, myrank-size, l, MPI_COMM_WORLD, &status);
            ind_fr = 1;
          }

          // lc
          if((indicator == 1) || (indicator == 4) || (indicator == 5) || (indicator == 6) || (indicator == 7) || (indicator == 9)){
            MPI_Recv(recv_buffer_lc, 8*N, MPI_PACKED, myrank+1, l, MPI_COMM_WORLD, &status);
            ind_lc = 1;
          }

          // fc
          if((indicator == 2) || (indicator == 3) || (indicator == 5) || (indicator == 7) || (indicator == 8) || (indicator == 9)){
            MPI_Recv(recv_buffer_fc, 8*N, MPI_PACKED, myrank-1, l, MPI_COMM_WORLD, &status);
            ind_fc = 1;
          }


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

        double buffer_lr[N];
        double buffer_fr[N];
        double buffer_lc[N];
        double buffer_fc[N];

        position = 0;

        // lr
        if(ind_lr != 0){
          for(int i=0;i<N;i++){
            MPI_Unpack(recv_buffer_lr, 8*N, &position, &buffer_lr[i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
          }
        }
        else{
          for(int i=0;i<N;i++){
            buffer_lr[i] = 0;
          }
        }

        position = 0;

        // // fr
        if(ind_fr != 0){
          for(int i=0;i<N;i++){
            MPI_Unpack(recv_buffer_fr, 8*N, &position, &buffer_fr[i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
          }
        }
        else{
          for(int i=0;i<N;i++){
            buffer_fr[i] = 0;
          }
        }

        position = 0;

        // // lc
        if(ind_lc != 0){
          for(int i=0;i<N;i++){
            MPI_Unpack(recv_buffer_lc, 8*N, &position, &buffer_lc[i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
          }
        }
        else{
          for(int i=0;i<N;i++){
            buffer_lc[i] = 0;
          }
        }

        position = 0;

        // // fc
        if(ind_fc != 0){
          for(int i=0;i<N;i++){
            MPI_Unpack(recv_buffer_fc, 8*N, &position, &buffer_fc[i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
          }
        }
        else{
          for(int i=0;i<N;i++){
            buffer_fc[i] = 0;
          }
        }


        for(int j=1;j<N+1;j++){
          //lr
          if(ind_lr != 0){
            data_temp[N+1][j] = buffer_lr[j-1];
          }
          else{
            data_temp[N+1][j] = 0;
          }

          // lc
          if(ind_lc != 0){
            data_temp[j][N+1] = buffer_lc[j-1];
          }
          else{
            data_temp[j][N+1] = 0;
          }

          // fc
          if(ind_fc != 0){
            data_temp[j][0] = buffer_fc[j-1];
          }
          else{
            data_temp[j][0] = 0;
          }

          // fr
          if(ind_fr != 0){
            data_temp[0][j] = buffer_fr[j-1];
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


  // ---------------------------------------------------------------


  // ---------------------------------------------------------------
  
  // Part III
  // for the third part when we have to from the new datatypes

  // we would be defining 4 new datatypes
  // last_row => last row of any matrix
  // first_row => first row of any matrix
  // last_column => last column of any matrix
  // first_column => first column of any matrix

  if(option == 3){

    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
        data[i][j] = 1+i+j;
      }
    }

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

    double recv_buffer_lr[N];
    double recv_buffer_lc[N];
    double recv_buffer_fr[N];
    double recv_buffer_fc[N];

  

    int blocklen_lr = N;
    int blocklen_fr = N;
    int blocklen_lc = 1;
    int blocklen_fc = 1;

    int stride_lr = 0;
    int stride_fr = 0;
    int stride_lc = N;
    int stride_fc = N;

    MPI_Type_vector (count_lr, blocklen_lr, stride_lr, MPI_DOUBLE, &last_row);
    MPI_Type_commit (&last_row);

    MPI_Type_vector (count_fr, blocklen_fr, stride_fr, MPI_DOUBLE, &first_row);
    MPI_Type_commit (&first_row);

    MPI_Type_vector (count_fc, blocklen_fc, stride_fc, MPI_DOUBLE, &first_column);
    MPI_Type_commit (&first_column);

    MPI_Type_vector (count_lc, blocklen_lc, stride_lc, MPI_DOUBLE, &last_column);
    MPI_Type_commit (&last_column);

    for(int i=0;i<time_step;i++){

      // for process 1
      // int MPI_Send (const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) 
      if(1){

          // lr
          if((indicator == 1) || (indicator == 3) || (indicator == 5) || (indicator == 6) || (indicator == 8) || (indicator == 9)){
            MPI_Send(&data[0][column_lr], 1, last_row, myrank+size, i, MPI_COMM_WORLD);
          }

          // fr
          if((indicator == 2) || (indicator == 4) || (indicator == 6) || (indicator == 7) || (indicator == 8) || (indicator == 9)){
            MPI_Send(&data[0][column_fr], 1, first_row, myrank-size, i, MPI_COMM_WORLD);
          }

          // lc
          if((indicator == 1) || (indicator == 4) || (indicator == 5) || (indicator == 6) || (indicator == 7) || (indicator == 9)){
            MPI_Send(&data[0][column_lc], 1, last_column, myrank+1, i, MPI_COMM_WORLD);
          }

          // fc
          if((indicator == 2) || (indicator == 3) || (indicator == 5) || (indicator == 7) || (indicator == 8) || (indicator == 9)){
            MPI_Send(&data[0][column_fc], 1, first_column, myrank-1, i, MPI_COMM_WORLD);
          }

      }


      int ind_lr = 0;
      int ind_fr = 0;
      int ind_fc = 0;
      int ind_lc = 0;

      // int MPI_Recv (void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
      MPI_Status status;
      if(1){

          // lr
          if((indicator == 1) || (indicator == 3) || (indicator == 5) || (indicator == 6) || (indicator == 8) || (indicator == 9)){
            MPI_Recv(recv_buffer_lr, N, MPI_DOUBLE, myrank+size, i, MPI_COMM_WORLD, &status);
            ind_lr = 1;
          }

          // fr
          if((indicator == 2) || (indicator == 4) || (indicator == 6) || (indicator == 7) || (indicator == 8) || (indicator == 9)){
            MPI_Recv(recv_buffer_fr, N, MPI_DOUBLE, myrank-size, i, MPI_COMM_WORLD, &status);
            ind_fr = 1;
          }

          // lc
          if((indicator == 1) || (indicator == 4) || (indicator == 5) || (indicator == 6) || (indicator == 7) || (indicator == 9)){
            MPI_Recv(recv_buffer_lc, N, MPI_DOUBLE, myrank+1, i, MPI_COMM_WORLD, &status);
            ind_lc = 1;
          }

          // fc
          if((indicator == 2) || (indicator == 3) || (indicator == 5) || (indicator == 7) || (indicator == 8) || (indicator == 9)){
            MPI_Recv(recv_buffer_fc, N, MPI_DOUBLE, myrank-1, i, MPI_COMM_WORLD, &status);
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


  if(myrank == req_rank){
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

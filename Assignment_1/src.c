#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "mpi.h"

double maxm(double a, double b){
  if(a>b){
    return a;
  }
  else{
    return b;
  }
}

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
  // int option = atoi(argv[3]);

  // for printing the given rank
  // int req_rank = atoi(argv[4]);

  // variable to capture the position of process in the cartesian system
  int indicator = 0;

  // double data[N][N];
  double **data = (double **)calloc(N , sizeof(double *));
  for(int i=0;i<N;i++){
    data[i] = (double *)calloc(N , sizeof(double));
  }

  // double data_temp[N+2][N+2];
  double **data_temp = (double **)calloc((N+2) , sizeof(double *));
  for(int i=0;i<(N+2);i++){
    data_temp[i] = (double *)calloc((N+2) , sizeof(double));
  }

  data_temp[0][0] = -DBL_MAX;
  data_temp[N+1][N+1] = -DBL_MAX;
  data_temp[0][N+1] = -DBL_MAX;
  data_temp[N+1][0] = -DBL_MAX;

  // get the rank and size of the current process
  int myrank, size_sq;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &size_sq);

  int size = (int)(sqrt(size_sq));

  // row for the current process
  int row_p = myrank/size;
  // column for the current pocess
  int col_p = myrank%size;

  // double time1[size_sq], time2[size_sq], time3[size_sq];
  // for(int i=0;i<size_sq;i++){
  //   time1[i] = 0;
  //   time2[i] = 0;
  //   time3[i] = 0;
  // }

  double time1, time2, time3;
  time1 = 0;
  time2 = 0;
  time3 = 0;

  double time1max, time2max, time3max;
  time1max = 0;
  time2max = 0;
  time3max = 0;

  double time1sum, time2sum, time3sum;
  time1sum = 0;
  time2sum = 0;
  time3sum = 0;

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

  // printf("Myrank:%d Row:%d Column:%d Indicator:%d\n",myrank,row_p,col_p,indicator);

  // ---------------------------------------------------------------
  
  // Part I
  // we would encoe the postion and iteration togeher in the TAG
  // in the format (length of points in row)*iteration + position

  if(1){

    time1 = 0;

    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
        data[i][j] = 1+i+j;
      }
    }

    double stime, etime;

    for(int i=0;i<time_step;i++){

      stime = MPI_Wtime();

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

        etime = MPI_Wtime();

        time1 = time1 + etime - stime;

        // compute

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
            data_temp[N+1][j] = -DBL_MAX;
          }

          // lc
          if(ind_lc != 0){
            data_temp[j][N+1] = recv_buffer_lc[j-1];
          }
          else{
            data_temp[j][N+1] = -DBL_MAX;
          }

          // fc
          if(ind_fc != 0){
            data_temp[j][0] = recv_buffer_fc[j-1];
          }
          else{
            data_temp[j][0] = -DBL_MAX;
          }

          // fr
          if(ind_fr != 0){
            data_temp[0][j] = recv_buffer_fr[j-1];
          }
          else{
            data_temp[0][j] = -DBL_MAX;
          }

        }

        for(int j=1;j<N+1;j++){
          for(int k=1;k<N+1;k++){
            int count = 0;
            double sum = 0;

            if(data_temp[j+1][k] != -DBL_MAX){
              count++;
              sum = sum + data_temp[j+1][k];
            }

            if(data_temp[j-1][k] != -DBL_MAX){
              count++;
              sum = sum + data_temp[j-1][k];
            }

            if(data_temp[j][k-1] != -DBL_MAX){
              count++;
              sum = sum + data_temp[j][k-1];
            }

            if(data_temp[j][k+1] != -DBL_MAX){
              count++;
              sum = sum + data_temp[j][k+1];
            }
            
            data[j-1][k-1] = sum/count;

          }
        }

    }

    // printf("time1 : %lf\n",time1);

  }

  // ---------------------------------------------------------------------


  // ---------------------------------------------------------------

  // Part II

  if(1){

    time2 = 0;

    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
        data[i][j] = 1+i+j;
      }
    }

    double stime, etime;

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

      stime = MPI_Wtime();

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

        etime = MPI_Wtime();

        time2 = time2 + etime - stime;

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
            buffer_lr[i] = -DBL_MAX;
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
            buffer_fr[i] = -DBL_MAX;
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
            buffer_lc[i] = -DBL_MAX;
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
            buffer_fc[i] = -DBL_MAX;
          }
        }


        for(int j=1;j<N+1;j++){
          //lr
          if(ind_lr != 0){
            data_temp[N+1][j] = buffer_lr[j-1];
          }
          else{
            data_temp[N+1][j] = -DBL_MAX;
          }

          // lc
          if(ind_lc != 0){
            data_temp[j][N+1] = buffer_lc[j-1];
          }
          else{
            data_temp[j][N+1] = -DBL_MAX;
          }

          // fc
          if(ind_fc != 0){
            data_temp[j][0] = buffer_fc[j-1];
          }
          else{
            data_temp[j][0] = -DBL_MAX;
          }

          // fr
          if(ind_fr != 0){
            data_temp[0][j] = buffer_fr[j-1];
          }
          else{
            data_temp[0][j] = -DBL_MAX;
          }

        }
        

        for(int j=1;j<N+1;j++){
          for(int k=1;k<N+1;k++){
            int count = 0;
            double sum = 0;

            if(data_temp[j+1][k] != -DBL_MAX){
              count++;
              sum = sum + data_temp[j+1][k];
            }

            if(data_temp[j-1][k] != -DBL_MAX){
              count++;
              sum = sum + data_temp[j-1][k];
            }

            if(data_temp[j][k-1] != -DBL_MAX){
              count++;
              sum = sum + data_temp[j][k-1];
            }

            if(data_temp[j][k+1] != -DBL_MAX){
              count++;
              sum = sum + data_temp[j][k+1];
            }
            
            data[j-1][k-1] = sum/count;

          }
        }


    }

    // printf("time2 : %lf\n",time2);
  }


  // ---------------------------------------------------------------


  // ---------------------------------------------------------------
  
  // Part III
  // for the third part when we have to from the new datatype

  if(1){

    time3 = 0;

    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
        data[i][j] = 1+i+j;
      }
    }

    double lrow[N];
    double lcol[N];
    double fcol[N];
    double frow[N];

    MPI_Datatype row;

    MPI_Type_vector (1, N, 0, MPI_DOUBLE, &row);
    MPI_Type_commit (&row);

    double recv_buffer_lr[N];
    double recv_buffer_lc[N];
    double recv_buffer_fr[N];
    double recv_buffer_fc[N];

    double stime, etime;

    for(int i=0;i<time_step;i++){

      for(int l=0;l<N;l++){
        lcol[l] = data[l][N-1];
        lrow[l] = data[N-1][l];
        fcol[l] = data[l][0];
        frow[l] = data[0][l];
      }

      stime = MPI_Wtime();
      // for process 1
      // int MPI_Send (const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) 
      if(1){

          // lr
          if((indicator == 1) || (indicator == 3) || (indicator == 5) || (indicator == 6) || (indicator == 8) || (indicator == 9)){
            MPI_Send(&lrow[0], 1, row, myrank+size, i, MPI_COMM_WORLD);
          }

          // fr
          if((indicator == 2) || (indicator == 4) || (indicator == 6) || (indicator == 7) || (indicator == 8) || (indicator == 9)){
            MPI_Send(&frow[0], 1, row, myrank-size, i, MPI_COMM_WORLD);
          }

          // lc
          if((indicator == 1) || (indicator == 4) || (indicator == 5) || (indicator == 6) || (indicator == 7) || (indicator == 9)){
            MPI_Send(&lcol[0], 1, row, myrank+1, i, MPI_COMM_WORLD);
          }

          // fc
          if((indicator == 2) || (indicator == 3) || (indicator == 5) || (indicator == 7) || (indicator == 8) || (indicator == 9)){
            MPI_Send(&fcol[0], 1, row, myrank-1, i, MPI_COMM_WORLD);
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

        etime = MPI_Wtime();

        time3 = time3 + etime - stime;

        // compute

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
            data_temp[N+1][j] = -DBL_MAX;
          }

          // lc
          if(ind_lc != 0){
            data_temp[j][N+1] = recv_buffer_lc[j-1];
          }
          else{
            data_temp[j][N+1] = -DBL_MAX;
          }

          // fc
          if(ind_fc != 0){
            data_temp[j][0] = recv_buffer_fc[j-1];
          }
          else{
            data_temp[j][0] = -DBL_MAX;
          }

          // fr
          if(ind_fr != 0){
            data_temp[0][j] = recv_buffer_fr[j-1];
          }
          else{
            data_temp[0][j] = -DBL_MAX;
          }

        }

        for(int j=1;j<N+1;j++){
          for(int k=1;k<N+1;k++){
            int count = 0;
            double sum = 0;

            if(data_temp[j+1][k] != -DBL_MAX){
              count++;
              sum = sum + data_temp[j+1][k];
            }

            if(data_temp[j-1][k] != -DBL_MAX){
              count++;
              sum = sum + data_temp[j-1][k];
            }

            if(data_temp[j][k-1] != -DBL_MAX){
              count++;
              sum = sum + data_temp[j][k-1];
            }

            if(data_temp[j][k+1] != -DBL_MAX){
              count++;
              sum = sum + data_temp[j][k+1];
            }
            
            data[j-1][k-1] = sum/count;

          }
        }

    }

    // printf("time3 : %lf\n",time3);

    MPI_Type_free(&row);
  }

  // ---------------------------------------------------------------------


  free(data_temp);
  free(data);

  MPI_Reduce (&time1, &time1max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (&time2, &time2max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (&time3, &time3max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  // MPI_Status status;

  // if(myrank != 0){
  //   for(int i=1;i<size_sq;i++){
  //     MPI_Send(&time1[myrank], 1, MPI_DOUBLE, 0, 10000, MPI_COMM_WORLD);
  //     MPI_Send(&time2[myrank], 1, MPI_DOUBLE, 0, 10000, MPI_COMM_WORLD);
  //     MPI_Send(&time3[myrank], 1, MPI_DOUBLE, 0, 10000, MPI_COMM_WORLD);
  //   }
  // }

  // if(myrank == 0){
  //   for(int i=1;i<size_sq;i++){
  //     MPI_Recv(&time1[i], 1, MPI_DOUBLE, i, 10000, MPI_COMM_WORLD, &status);
  //     MPI_Recv(&time2[i], 1, MPI_DOUBLE, i, 10000, MPI_COMM_WORLD, &status);
  //     MPI_Recv(&time3[i], 1, MPI_DOUBLE, i, 10000, MPI_COMM_WORLD, &status);
  //   }
  // }

  // if(myrank == 0){
  //   for(int i=0;i<size_sq;i++){
  //     time1s = time1s + time1[i];
  //     time2s = time2s + time2[i];
  //     time3s = time3s + time3[i];

  //     time1m = maxm(time1m, time1[i]);
  //     time2m = maxm(time2m, time2[i]);
  //     time3m = maxm(time3m, time3[i]);

  //   }

  //   printf("Total time 1 : %lf\n", time1s);
  //   printf("Total time 2 : %lf\n", time2s);
  //   printf("Total time 3 : %lf\n", time3s);


    // printf("Maximum time 1 : %lf\n", time1max);
    // printf("Maximum time 2 : %lf\n", time2max);
    // printf("Maximum time 3 : %lf\n", time3max);
  // }

    // if(myrank == 0){
    //   printf("Maximum time 1 : %lf\n", time1max);
    //   printf("Maximum time 2 : %lf\n", time2max);
    //   printf("Maximum time 3 : %lf\n", time3max);
    // }

    if(myrank == 0){
      printf("%lf\n", time1max);
      printf("%lf\n", time2max);
      printf("%lf\n", time3max);
    }


  // done with MPI
  MPI_Finalize();
}

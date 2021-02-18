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

  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      data[i][j] = 1+i+j;
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

  printf("Myrank:%d Row:%d Column:%d Indicator:%d\n",myrank,row_p,col_p,indicator);
  // ---------------------------------------------------------------

  // Part II

  // double send_buffer_lr[N];
  // double send_buffer_fr[N];
  // double send_buffer_lc[N];
  // double send_buffer_fc[N];

  // double recv_buffer_lr[N];
  // double recv_buffer_fr[N];
  // double recv_buffer_lc[N];
  // double recv_buffer_fc[N];

  // int position = 0;


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

    // // fr first row
    // for(int i=0;i<N;i++){
    //   MPI_Pack(&data[0][i], 1, MPI_DOUBLE,send_buffer_fr,8*N,&position, MPI_COMM_WORLD);
    //   printf("%d\n",i);
    // }

    // position = 0;

    // // lr last row
    // for(int i=0;i<N;i++){
    //   MPI_Pack(&data[N-1][i], 1, MPI_DOUBLE,send_buffer_lr,8*N,&position, MPI_COMM_WORLD);
    // }

    // position = 0;

    // // fc first column
    // for(int i=0;i<N;i++){
    //   MPI_Pack(&data[i][0], 1, MPI_DOUBLE,send_buffer_fc,8*N,&position, MPI_COMM_WORLD);
    // }

    // position = 0;

    // // lc last column
    // for(int i=0;i<N;i++){
    //   MPI_Pack(&data[i][N-1], 1, MPI_DOUBLE,send_buffer_lc,8*N,&position, MPI_COMM_WORLD);
    // }

        // if(myrank == option){

        //   for(int i=0;i<N;i++){
        //     printf("%f ", send_buffer_fc[i]);
        //   }
        //   for(int i=0;i<N;i++){
        //     printf("%f ", send_buffer_lc[i]);
        //   }
        //   for(int i=0;i<N;i++){
        //     printf("%f ", send_buffer_fr[i]);
        //   }
        //   for(int i=0;i<N;i++){
        //     printf("%f ", send_buffer_lr[i]);
        //   }

        //   printf("\n\n\n\n");


        //   // for(int i=0;i<N+2;i++){
        //   //   for(int j=0;j<N+2;j++){
        //   //     printf("%f ",data_temp[i][j]);
        //   //   }
        //   //   printf("\n");
        //   // }
        // }

    // position = 0;

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
      // if(1){

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

        // }

        // if(myrank == option){
        //   printf("%d %d %d %d \n", ind_fc, ind_lc, ind_lr, ind_fr);
        //   for(int i=0;i<N;i++){
        //     printf("%f %f %f %f\n",recv_buffer_fc[i], recv_buffer_lc[i], recv_buffer_lr[i], recv_buffer_fr[i]);
        //   }
        // }


          /// verified till here

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


        if(myrank == option){

          for(int i=0;i<N;i++){
            printf("%f ",buffer_fc[i]);
          }
          printf("\n");
          for(int i=0;i<N;i++){
            printf("%f ",buffer_lc[i]);
          }
          printf("\n");
          for(int i=0;i<N;i++){
            printf("%f ",buffer_fr[i]);
          }
          printf("\n");
          for(int i=0;i<N;i++){
            printf("%f ",buffer_lr[i]);
          }
        }


        //// verified till here srinjay


        //   // for(int i=0;i<N+2;i++){
        //   //   for(int j=0;j<N+2;j++){
        //   //     printf("%f ",data_temp[i][j]);
        //   //   }
        //   //   printf("\n");
        //   // }
        // }
        

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


  // ---------------------------------------------------------------


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


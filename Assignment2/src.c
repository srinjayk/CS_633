// Srinjay Kumar
// Ayush Soneria

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

int find_grp(char* hostname){

  const char *grp_1[16] = {"csews1","csews2","csews3","csews4","csews5","csews6","csews7","csews8","csews9","csews10","csews11","csews12","csews14","csews15","csews16","csews31"};
  const char *grp_2[16] = {"csews13","csews17","csews18","csews19","csews20","csews21","csews22","csews23","csews24","csews25","csews26","csews27","csews28","csews29","csews30","csews32"};
  const char *grp_3[13] = {"csews33","csews34","csews35","csews36","csews37","csews38","csews39","csews40","csews41","csews42","csews43","csews44","csews46"};
  const char *grp_4[14] = {"csews45","csews47","csews48","csews49","csews50","csews51","csews52","csews53","csews54","csews56","csews58","csews59","csews60","csews61"};
  const char *grp_5[17] = {"csews62","csews63","csews64","csews65","csews66","csews67","csews68","csews69","csews70","csews71","csews72","csews73","csews74","csews75","csews76","csews77","csews78"};
  const char *grp_6[14] = {"csews79","csews80","csews81","csews82","csews83","csews84","csews85","csews86","csews87","csews88","csews89","csews90","csews91","csews92"};

  for(int i=0;i<16;i++){
    if(strcmp(hostname, grp_1[i]) == 0){
      return 1;
    }
  }

  for(int i=0;i<16;i++){
    if(strcmp(hostname, grp_2[i]) == 0){
      return 2;
    }
  }

  for(int i=0;i<13;i++){
    if(strcmp(hostname, grp_3[i]) == 0){
      return 3;
    }
  }

  for(int i=0;i<14;i++){
    if(strcmp(hostname, grp_4[i]) == 0){
      return 4;
    }
  }

  for(int i=0;i<17;i++){
    if(strcmp(hostname, grp_5[i]) == 0){
      return 5;
    }
  }

  for(int i=0;i<14;i++){
    if(strcmp(hostname, grp_6[i]) == 0){
      return 6;
    }
  }

}

int find_node(char* hostname){
  const char *cpu_name[92] = {"csews1","csews2","csews3","csews4","csews5","csews6","csews7","csews8","csews9","csews10","csews11","csews12","csews13","csews14","csews15","csews16","csews17","csews18","csews19","csews20","csews21","csews22","csews23","csews24","csews25","csews26","csews27","csews28","csews29","csews30","csews31","csews32","csews33","csews34","csews35","csews36","csews37","csews38","csews39","csews40","csews41","csews42","csews43","csews44", "csews45","csews46","csews47","csews48","csews49","csews50","csews51","csews52","csews53","csews54","csews55","csews56","csews57","csews58","csews59","csews60","csews61","csews62","csews63","csews64","csews65","csews66","csews67","csews68","csews69","csews70","csews71","csews72","csews73","csews74","csews75","csews76","csews77","csews78","csews79","csews80","csews81","csews82","csews83","csews84","csews85","csews86","csews87","csews88","csews89","csews90","csews91","csews92"};
  
  for(int i=0;i<92;i++){
    if(strcmp(hostname, cpu_name[i]) == 0){
      return i+1;
    }
  }

  return 0;
}

int main(int argc, char *argv[]){

  // initialize MPI
  MPI_Init (&argc, &argv);

  // number of data points per process
  // int N = ((atoi(argv[1]))*1024)/8;

  // int N = 4;
  int N = atoi(argv[1]);
  // collective calls
  int call_type = atoi(argv[2]);

  // optimised or unoptimised
  int opt = atoi(argv[3]);

  // variable to capture the position of process in the cartesian system
  int indicator = 0;

  // get the rank and size of the current process
  int myrank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Status status;

  MPI_Group old_group;
  MPI_Comm_group(MPI_COMM_WORLD, &old_group);

  int len;
  char hostname[MPI_MAX_PROCESSOR_NAME];
  // int coreID = sched_getcpu();

  MPI_Get_processor_name (hostname, &len);

  // printf("rank %d on %s and size %d\n",myrank,hostname, sizeof(hostname));
  // printf("%d \n", find_grp(hostname));

  int grp[size];

  grp[myrank] = find_grp(hostname);

  if(myrank > 0){
    MPI_Send(&grp[myrank], 1, MPI_INT, 0, myrank, MPI_COMM_WORLD);
  }

  if(myrank == 0){
    for(int i=1;i<size;i++)
      MPI_Recv(&grp[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
  }

  MPI_Bcast(grp, size, MPI_INT, 0, MPI_COMM_WORLD);

  // if(myrank == 0){
  //   for(int i=0;i<size;i++){
  //     printf("%d\n",grp[i]);
  //   }
  // }

  int node[size];

  node[myrank] = find_node(hostname);

  if(myrank > 0){
    MPI_Send(&node[myrank], 1, MPI_INT, 0, myrank, MPI_COMM_WORLD);
  }

  if(myrank == 0){
    for(int i=1;i<size;i++)
      MPI_Recv(&node[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
  }

  MPI_Bcast(node, size, MPI_INT, 0, MPI_COMM_WORLD);

  // if(myrank == 1){
  //   for(int i=0;i<size;i++){
  //     printf("%d\n",node[i]);
  //   }
  // }

  int ppg[6];
  int ppn = 0;

  for(int i=0;i<6;i++){
    ppg[i] = 0;
  }

  for(int i=0;i<size;i++){
    ppg[grp[i]-1]++;
  }

  int a = node[0];
  for(int i=0;i<size;i++){
    if(node[i] == a){
      ppn++;
    }
  }

  int grp_ldr[6];
  for(int i=0;i<6;i++){
    grp_ldr[i] = -1;
  }

  for(int i=0;i<size;i++){
    int a = grp[i];
    if(grp_ldr[a-1] == -1){
      grp_ldr[a-1] = i;
    }
  }

  // if(myrank == 0){
  //   for(int i=0;i<6;i++){
  //     printf("%d\n",grp_ldr[i]);
  //   }
  // }

  int cnt_ldr = 0;

  for(int i=0;i<6;i++){
    if(grp_ldr[i] != -1){
      cnt_ldr++;
    }
  }

  int fin_ldr[cnt_ldr];

  int j=0;
  for(int i=0;i<6;i++){
    if(grp_ldr[i] != -1){
      fin_ldr[j] = grp_ldr[i];
      j++;
    }
    // j++;
  }

  // if(myrank == 0){
  //   for(int i=0;i<cnt_ldr;i++){
  //     printf("%d\n",fin_ldr[i]);
  //   }
  // }

  // if(myrank == 0){
  //   for(int i=0;i<6;i++){
  //     printf("%d\n", ppg[i]);
  //   }
  //   printf("%d\n",ppn);
  // }



  // call type 1
  // BCast
  if(call_type == 1){

    if(opt == 0){
      double buf[N];

      for(int i=0;i<N;i++){
        buf[i] = i+1+myrank;
      }

      double stime, etime, ttime, mtime;
      stime = MPI_Wtime();

      MPI_Bcast(buf, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      etime = MPI_Wtime();

      ttime = etime - stime;

      MPI_Reduce(&ttime, &mtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if(myrank == 0){
        printf("%lf\n", mtime);
      }
      // printf("%lf\n",ttime);

      // if(myrank == 63){
      //   for(int i=0;i<10;i++){
      //     printf("%lf\n",buf[i]);
      //   }
      // }
    }
    
    if(opt == 1){
      double buf[N];

      double stime, etime, ttime, mtime;

      for(int i=0;i<N;i++){
        buf[i] = i+1+myrank;
      }

      stime = MPI_Wtime();

      MPI_Group inter_group;

      MPI_Group_incl(old_group, cnt_ldr, fin_ldr, &inter_group);
      MPI_Comm inter_group_comm;
      MPI_Comm_create_group(MPI_COMM_WORLD, inter_group, 123, &inter_group_comm);

      for(int i=0;i<cnt_ldr;i++){
        if(myrank == fin_ldr[i])
          MPI_Bcast(buf, N, MPI_DOUBLE, 0, inter_group_comm);
      }

      // if(myrank == 1){
      //   for(int i=0;i<N;i++){
      //     printf("%lf\n",buf[i]);
      //   }
      // }

      int color = grp[myrank];

      MPI_Comm intra_group;
      MPI_Comm_split (MPI_COMM_WORLD, color, myrank, &intra_group);

      // printf("%d %d %d \n",myrank, color, grp_ldr[color - 1]);

      MPI_Bcast(buf, N,MPI_DOUBLE, 0, intra_group);
      
      etime = MPI_Wtime();

      ttime = etime - stime;
      // printf("%lf\n",ttime);

      MPI_Reduce(&ttime, &mtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if(myrank == 0){
        printf("%lf\n", mtime);
      }

      // if(myrank == 63){
      //   for(int i=0;i<10;i++){
      //     printf("%lf\n",buf[i]);
      //   }
      // }
    }

  }

  // call type 2
  // MPI Reduce 
  if(call_type == 2){

    if(opt == 0){
      double buf[N];

      for(int i=0;i<N;i++){
        buf[i] = i+1+myrank;
      }

      double stime, etime, ttime, mtime;
      stime = MPI_Wtime();

      double sol_val[N];
      MPI_Reduce(buf, sol_val, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      etime = MPI_Wtime();

      ttime = etime - stime;
      // printf("%lf\n",ttime);

      MPI_Reduce(&ttime, &mtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if(myrank == 0){
        printf("%lf\n", mtime);
      }

      // if(myrank == 0){
      //   for(int i=0;i<N;i++){
      //     printf("%lf\n", sol_val[i]);
      //   }
      // }
    }

    if(opt == 1){
      double buf[N];
      double agg_val[N];
      double fin_val[N];

      for(int i=0;i<N;i++){
        buf[i] = i+1;
        agg_val[i] = 0;
      }

      int color = grp[myrank];

      double stime1, etime1, stime2, etime2, ttime, mtime;
      // stime = MPI_Wtime();

      MPI_Comm intra_group;
      MPI_Comm_split (MPI_COMM_WORLD, color, myrank, &intra_group);

      // MPI_Bcast(buf, N,MPI_DOUBLE, 0, intra_group);
      stime1 = MPI_Wtime();
      MPI_Reduce(buf, agg_val, N, MPI_DOUBLE, MPI_SUM, 0, intra_group);
      etime1 = MPI_Wtime();

      MPI_Group inter_group;

      MPI_Group_incl(old_group, cnt_ldr, fin_ldr, &inter_group);
      MPI_Comm inter_group_comm;
      MPI_Comm_create_group(MPI_COMM_WORLD, inter_group, 123, &inter_group_comm);

      stime2 = MPI_Wtime();
      for(int i=0;i<cnt_ldr;i++){
        if(myrank == fin_ldr[i])
          MPI_Reduce(agg_val, fin_val, N, MPI_DOUBLE, MPI_SUM, 0, inter_group_comm);
      }
      etime2 = MPI_Wtime();

      ttime = etime2 + etime1 - stime2 - stime1;
      // printf("%lf\n",ttime);

      MPI_Reduce(&ttime, &mtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if(myrank == 0){
        printf("%lf\n", mtime);
      }
      // MPI_Reduce(&ttime, &mtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      // if(myrank == 0){

      // if(myrank == 0){
      //   for(int i=0;i<N;i++){
      //     printf("%lf\n", fin_val[i]);
      //   }
      // }



    }
    

  }


  // call type 3
  // MPI_Gather
  if(call_type == 3){
    N = N/size;

    if(opt == 0){
      double buf[N];

      for(int i=0;i<N;i++){
        buf[i] = 1+i+myrank;
      }

      double recvMessage[size * N];
      double stime, etime, ttime, mtime;
      stime = MPI_Wtime();

      MPI_Gather(buf, N, MPI_DOUBLE, recvMessage, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      etime = MPI_Wtime();

      ttime = etime - stime;
      // printf("%lf\n",ttime);

      MPI_Reduce(&ttime, &mtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if(myrank == 0){
        printf("%lf\n", mtime);
      }

      // if(myrank == 0){
      //   for(int i=0;i<size*N;i++){
      //     printf("%lf\n", recvMessage[i]);
      //   }
      // }
    }

    if(opt == 1){
      double buf[N+1];
      buf[0] = myrank;

      for(int i=1;i<N+1;i++){
        buf[i] = i+myrank;
      }

      int color = grp[myrank];

      double stime1, etime1, stime2, etime2, ttime, mtime;
      

      MPI_Comm intra_group;
      MPI_Comm_split (MPI_COMM_WORLD, color, myrank, &intra_group);

      double recv[(N+1)*(ppg[color-1])];

      stime1 = MPI_Wtime();
      MPI_Gather(buf, N+1, MPI_DOUBLE, recv, N+1, MPI_DOUBLE, 0, intra_group);
      etime1 = MPI_Wtime();
      // if(myrank == 0){
      //   for(int i=0;i<(N+1)*(ppg[color-1]);i++)
      //     printf("%lf\n", recv[i]);
      // }

      MPI_Group inter_group;

      MPI_Group_incl(old_group, cnt_ldr, fin_ldr, &inter_group);
      MPI_Comm inter_group_comm;
      MPI_Comm_create_group(MPI_COMM_WORLD, inter_group, 123, &inter_group_comm);

      double agg[size * (N+1)];

      // if(myrank == 0){
      //   for(int i=0;i<cnt_ldr;i++)
      //     printf("%d\n", fin_ldr[i]);
      // }

      stime2 = MPI_Wtime();
      // for(int i=0;i<cnt_ldr;i++){
      //   if(myrank == fin_ldr[i]){
      //   //   color = grp[myrank];
      //     MPI_Gather(recv, (N+1)*(ppg[color-1]), MPI_DOUBLE, agg, (N+1)*(ppg[color-1]), MPI_DOUBLE, 0, inter_group_comm);
      //   }

      //   // if(myrank == fin_ldr[i]){
      //   //   if(myrank != 0){
      //   //     MPI_Send(recv, (N+1)*(ppg[color-1]), MPI_DOUBLE, 0, i, inter_group_comm);
      //   //   }

      //   //   if(myrank == 0){
      //   //     MPI_Recv(agg)
      //   //   }
      //   // }

      // }

      int nrank;
      MPI_Comm_rank(intra_group, &nrank);

      if((nrank == 0) && (myrank > 0)){
        MPI_Send(recv, (N+1)*(ppg[color-1]), MPI_DOUBLE, 0, color-1, MPI_COMM_WORLD);
      }

      double fin[size*N];

      if(myrank == 0){
        if(ppg[1] > 0){
          double tmp[ppg[1] * (N+1)];
          MPI_Recv(tmp, ppg[1] * (N+1), MPI_DOUBLE, grp_ldr[1], 1, MPI_COMM_WORLD, &status);
          for(int i=0;i<ppg[1];i++){
            int rnk = tmp[i*(N+1)];
            for(int j=0;j<N;j++){
              fin[rnk*N + j] = tmp[i*(N+1) + j];
            }
          }
        }
        if(ppg[2] > 0){
          double tmp[ppg[2] * (N+1)];
          MPI_Recv(tmp, ppg[2] * (N+1), MPI_DOUBLE, grp_ldr[2], 2, MPI_COMM_WORLD, &status);
          for(int i=0;i<ppg[2];i++){
            int rnk = tmp[i*(N+1)];
            for(int j=0;j<N;j++){
              fin[rnk*N + j] = tmp[i*(N+1) + j];
            }
          }
        }
        if(ppg[3] > 0){
          double tmp[ppg[3] * (N+1)];
          MPI_Recv(tmp, ppg[3] * (N+1), MPI_DOUBLE, grp_ldr[3], 3, MPI_COMM_WORLD, &status);
          for(int i=0;i<ppg[3];i++){
            int rnk = tmp[i*(N+1)];
            for(int j=0;j<N;j++){
              fin[rnk*N + j] = tmp[i*(N+1) + j];
            }
          }
        }
        if(ppg[4] > 0){
          double tmp[ppg[4] * (N+1)];
          MPI_Recv(tmp, ppg[4] * (N+1), MPI_DOUBLE, grp_ldr[4], 4, MPI_COMM_WORLD, &status);
          for(int i=0;i<ppg[4];i++){
            int rnk = tmp[i*(N+1)];
            for(int j=0;j<N;j++){
              fin[rnk*N + j] = tmp[i*(N+1) + j];
            }
          }
        }
        if(ppg[5] > 0){
          double tmp[ppg[5] * (N+1)];
          MPI_Recv(tmp, ppg[5] * (N+1), MPI_DOUBLE, grp_ldr[5], 5, MPI_COMM_WORLD, &status);
          for(int i=0;i<ppg[5];i++){
            int rnk = tmp[i*(N+1)];
            for(int j=0;j<N;j++){
              fin[rnk*N + j] = tmp[i*(N+1) + j];
            }
          }
        }
        // if(ppg[0] > 0){
        //   double tmp[ppg[0] * (N+1)];
        //   MPI_Recv(tmp, ppg[0] * (N+1), MPI_DOUBLE, grp_ldr[0], 1, MPI_COMM_WORLD, &status);
        // }

        for(int i=0;i<ppg[0];i++){
            int rnk = recv[i*(N+1)];
            for(int j=0;j<N;j++){
              fin[rnk*N + j] = recv[i*(N+1) + j];
            }
          }
      }

      etime2 = MPI_Wtime();

      ttime = etime2 + etime1 - stime2 - stime1;
      // // printf("%lf\n",ttime);

      MPI_Reduce(&ttime, &mtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if(myrank == 0){
        printf("%lf\n", mtime);
      }

      // double fin[size*N];

      // for(int i=0;i<size;i++){
      //   int rnk = agg[i*(N+1)];
      //   for(int j=0;j<N;j++){
      //     fin[rnk*N + j] = agg[i*(N+1) + j];
      //   }
      // }


      // if(myrank == 0){
      //   for(int i=0;i<(N)*(size);i++)
      //     printf("%lf\n", fin[i]);
      // }


      

    }
  }
  

  // call type 4
  // MPI AlltoAllv
  if(call_type == 4){
    int displ[size];
    int cntarr[size];

    double buf[N*size];
    for(int i=0;i<N*size;i++){
      buf[i] = 1+i+myrank;
    }

    for(int i=0;i<size;i++){
      displ[i] = i*N;
      cntarr[i] = N;
    }

    double recv[size*N];

    MPI_Alltoallv(buf, cntarr, displ, MPI_DOUBLE, recv, cntarr, displ, MPI_DOUBLE, MPI_COMM_WORLD);

    // if(myrank == 0){
    //   for(int i=0;i<N*size;i++){
    //     printf("%lf\n",recv[i]);
    //   }
    // }
  }
  

  // done with MPI
  MPI_Finalize();
}


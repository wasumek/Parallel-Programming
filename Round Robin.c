#include "mpi.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc,char ** argv) {
  int dataSize = 10;//Vary this size later
  int nprocs;
  int rank;
  int source,dest,count,tag=1;
  int i;
  MPI_Status stat;
  MPI_Request request;
  double sum = 0.0;
  double * sendbuf = (double *) malloc (sizeof(double)*dataSize);
  double * sendbuf2 = (double *) malloc (sizeof(double)*dataSize); // might be needed
  /* QUESTION */
  // Is it possible to avoid the data copy from recvbuf to sendbuf?
  /* ANSWER */
  // Yes, just pass on the address of the recvbuf, hence no data copied
  double * rndData = (double *) malloc (sizeof(double)*dataSize);

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Init(&argc,&argv);
  MPI_Comm_size (comm,&nprocs);
  MPI_Comm_rank (comm,&rank);
  printf("I'm %d of %d",rank,nprocs);
  int root,last;
  root = 0;
  last = nprocs - 1;
  count = 1;

  //Every proc generates their own recvbuf vector and randomize it.
  double * recvbuf = (double *) malloc (sizeof(double)*dataSize);
  srand(time(NULL));
  for ( i=0;i<dataSize;i++) recvbuf[i] = rand();

  //initialize before the round-robin exchange
  for ( i=0;i<dataSize;i++) rndData[i] = 0;

  /* ----------------------  Round-robin Starts  ----------------------- */

  if(rank==last){ 
    dest = root; //next is first
    source = i-1;
  } else if(rank==root) {
    dest = i+1;
    source = last; //get from last
  } else {
    dest = i+1;
    source = i-1;
  }

  //Blocking sends
  MPI_Send(&sendbuf, count, MPI_INT, dest, tag, comm);
  
  //Synchronous sends
  /* MPI_Ssend(&sendbuf, count, MPI_INT, dest, tag, comm); */

  //Ready sends
  /* MPI_Rsend(&sendbuf, count, MPI_INT, dest, tag ,comm);  */

  //Non-blocking sends
  /* MPI_Isend(&sendbuf, count, MPI_INT, dest, tag, comm, &request); */

  //Normal Receive
  MPI_Recv(&recvbuf, count, MPI_INT, source, tag, comm, &stat);
  /* ----------------------  End Round-robin  ----------------------- */

  //Vector-addition of the received data to an array rndData
  for ( i=0;i<dataSize;i++) rndData[i] = rndData[i] + recvbuf[i];

  //Compute the checksum of all elements in the array rndData and adding up the entries and print if once per PE
  for ( i=0;i<dataSize;i++) sum = sum + rndData[i]; //Checksum
  printf("I'm %d with sum: %lf",rank,sum);//print once


  /* QUESTIONS */
  
  //1.To which method does MPI_Send correspond?
  /* ANSWER */
  /* For Ssend, Rsend, and Send all resulted in Deadlocking while using sendbuf. */

  //2.Does this behavior change with the message size?
  /* ANSWER */
  /* With the obsevation, it seemed that Ssend does not work with any datasize, however, Rsend and Send work when the datasize is fit within the send buffer. */

  //There will be a deadlock in some cases, for non-deadlock, plot with SCALASCA and VAMPIR
  /*by logging into cluster with ssh -Y flag and download
    module load UNITE scalasca scorep vampir
    scorep mpicc main.c -o hw3
    scan -t mpiexec -np 4 ./hw3
    ./scorep_hw3_4_trace/
    vampir scorep_hw3_4_trace/traces.otf2
  */
  MPI_Finalize();
}

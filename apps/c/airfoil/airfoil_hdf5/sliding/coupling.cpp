/*
 * Open source copyright declaration based on BSD open source template:
 * http://www.opensource.org/licenses/bsd-license.php
 *
 * This file is part of the OP2 distribution.
 *
 * Copyright (c) 2011, Mike Giles and others. Please see the AUTHORS file in
 * the main source directory for a full list of copyright holders.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The name of Mike Giles may not be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY Mike Giles ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL Mike Giles BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//
//     Nonlinear airfoil lift calculation
//
//     Written by Mike Giles, 2010-2011, based on FORTRAN code
//     by Devendra Ghate and Mike Giles, 2005
//

//
// standard headers
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

// main program

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  //Figure out who belongs in what group
  int *groups = (int *)malloc(size * sizeof(int));
  int *groups2 = (int *)malloc(size * sizeof(int));
  int *coupling_group = (int *)malloc(size * sizeof(int));
  int my_type = 0; //This is to be read from a configuration file
  MPI_Allgather(&my_type, 1, MPI_INT, groups, 1, MPI_INT, MPI_COMM_WORLD);

  //count number of groups
  int num_groups = 0;
  for (int i = 0; i < size; i++) num_groups = num_groups > groups[i] ? num_groups : groups[i];
  num_groups++;

  int iam_root = 0;

  //The global group
  MPI_Group global_grp;
  MPI_Comm_group(MPI_COMM_WORLD, &global_grp);

  //Create sub-groups and sub-communicators
  MPI_Group mpigroups[num_groups];
  MPI_Comm mpicomms[num_groups];
  int count = 0;
  for (int i = 0; i < num_groups; ++i) {
    count = 0;
    for (int j = 0; j < size; ++j) {
      if (groups[j] == i) {
        if (i==0 && count == 0 && rank == j) iam_root = 1;
        groups2[count++] = j;
      }
    }
    MPI_Group_incl(global_grp, count, groups2, &mpigroups[i]);
    MPI_Comm_create(MPI_COMM_WORLD, mpigroups[i], &mpicomms[i]);
  }

  //Figure out who is in the coupling group (id 0)
  int coupling_group_size = 0;
  int rank_in_group = -1;
  for (int j = 0; j < size; ++j) {
    if (groups[j] == 0) {
      if (j==rank) rank_in_group = coupling_group_size;
      coupling_group[coupling_group_size++] = j;
    }
  }

  //
  // op_export_init
  //

  //Setting up for export
  int node_sizes[count];
  int cell_sizes[count];

  MPI_Status statuses[count];
  MPI_Request requests[count];

  //Step 1 receive set sizes from all OP2 procs
  int recv[2];
  for (int i = 0; i < count; i++) {
    MPI_Recv(recv, 2, MPI_INT, groups2[i], 100, MPI_COMM_WORLD, &statuses[i]);
    cell_sizes[i] = recv[0];
    node_sizes[i] = recv[1];
  }

  //Step 2: receive global ndoe indices and parts of the cellsToNodes map
  int *node_indices[count];
  for (int i = 0; i < count; i++) {
    node_indices[i] = (int*)malloc(node_sizes[i]*sizeof(int));
    MPI_Recv(node_indices[i], node_sizes[i], MPI_INT, groups2[i], 101, MPI_COMM_WORLD, &statuses[i]);
  }

  int *cellsToNodes[count];
  for (int i = 0; i < count; i++) {
    cellsToNodes[i] = (int*)malloc(4*cell_sizes[i]*sizeof(int)); //Note: map dimensionality has to be known
    MPI_Recv(cellsToNodes[i], 4*cell_sizes[i], MPI_INT, groups2[i], 102, MPI_COMM_WORLD, &statuses[i]);
  }

  //Step 3: root sends list of processes that the OP2 processes will have to send to
  if (iam_root) {
    int num_target_procs = coupling_group_size;
    for (int i = 0; i < count; i++) {
      MPI_Send(&num_target_procs, 1, MPI_INT, groups2[i], 103, MPI_COMM_WORLD);
      MPI_Send(coupling_group, num_target_procs, MPI_INT, groups2[i], 104, MPI_COMM_WORLD);
    }
  }

  //
  // op_import_init
  //
  //Step 1 receive node set size from all OP2 procs
  int recv_nodesize[count];
  for (int i = 0; i < count; i++) {
    MPI_Recv(&recv_nodesize[i], 1, MPI_INT, groups2[i], 400, MPI_COMM_WORLD, &statuses[i]);
  }

  //Step 2 receive coordinate data (3D)
  char *coords[count];
  for (int i = 0; i < count; ++i) {
    coords[i] = (char*)malloc(recv_nodesize[i]*3*sizeof(double));
  }
  for (int i = 0; i < count; i++) {
    MPI_Recv(coords[i], recv_nodesize[i]*3*sizeof(double), MPI_CHAR, groups2[i], 401, MPI_COMM_WORLD, &statuses[i]);
  }

  //Divide up who sends back what - in reality this will change
  int send_begin[count];
  int send_end[count];
  for (int i = 0; i < count; i++) {
    send_begin[i] = (recv_nodesize[i]/coupling_group_size) * rank_in_group;
    send_end[i] = (recv_nodesize[i]/coupling_group_size) * (rank_in_group+1);
    if (rank_in_group == coupling_group_size-1) send_end[i] = recv_nodesize[i];
  }

  //prepare send buffers for import
  char *send_buf[count];
  for (int i = 0; i < count; ++i) {
    send_buf[i] = (char*)malloc((send_end[i]-send_begin[i])*(2*sizeof(int) + 4*sizeof(double)));
  }
  int item_size = 2*sizeof(int) + 4*sizeof(double);



  //
  // op_export and op_import - iterative loop
  //

  //  we are going to be receiving bound (1 int), coord (3 double), p (1 double)
  char *recv_buffers[count];
  for (int i = 0; i < count; ++i) {
    recv_buffers[i] = (char*)malloc(node_sizes[i]*(sizeof(int) + 4*sizeof(double)));
  }
  for (int repeat = 0; repeat < 10; repeat++) {
    //op_export - the receive side. Data is received as SoA
    for (int i = 0; i < count; ++i) {
      MPI_Irecv(recv_buffers[i], node_sizes[i]*(sizeof(int) + 4*sizeof(double)), MPI_CHAR, groups2[i], 200, MPI_COMM_WORLD, &requests[i]);
      printf("Coupling (%d) receiving %d bytes from %d\n",rank, node_sizes[i]*(sizeof(int) + 4*sizeof(double)), groups2[i]);
    }
    MPI_Waitall(count, requests, statuses);

    //op_import - send side
    //we send the data back, based on up-front division
    for (int i = 0; i < count; ++i) {
      for (int j = send_begin[i]; j < send_end[i]; j++) {
        //Pack structure: global node id and then the payload, data is sent back as AoS
        *(int*)(send_buf[i]+(j-send_begin[i])*item_size) = node_indices[i][j]; //global node id
        *(int*)(send_buf[i]+(j-send_begin[i])*item_size + sizeof(int)) = ((int*)(recv_buffers[i]))[j]; //p_bound
        *(double*)(send_buf[i]+(j-send_begin[i])*item_size + 2*sizeof(int))                  = ((double*)(recv_buffers[i] + node_sizes[i]*sizeof(int)))[3*j+0]; //coord x
        *(double*)(send_buf[i]+(j-send_begin[i])*item_size + 2*sizeof(int)+sizeof(double))   = ((double*)(recv_buffers[i] + node_sizes[i]*sizeof(int)))[3*j+1]; //coord y
        *(double*)(send_buf[i]+(j-send_begin[i])*item_size + 2*sizeof(int)+2*sizeof(double)) = ((double*)(recv_buffers[i] + node_sizes[i]*sizeof(int)))[3*j+2]; //coord z
        *(double*)(send_buf[i]+(j-send_begin[i])*item_size + 2*sizeof(int)+3*sizeof(double)) = ((double*)(recv_buffers[i] + node_sizes[i]*sizeof(int)+node_sizes[i]*3*sizeof(double)))[j]; //pressure
      }
      MPI_Isend(send_buf[i], (send_end[i]-send_begin[i])*(2*sizeof(int) + 4*sizeof(double)), MPI_CHAR, groups2[i], 500, MPI_COMM_WORLD, &requests[i]);
      printf("Coupling (%d) sending %d->%d (%d bytes) to %d\n", rank, send_begin[i], send_end[i], (send_end[i]-send_begin[i])*(2*sizeof(int) + 4*sizeof(double)), groups2[i]);
    }
    MPI_Waitall(count, requests, statuses);
    printf("Sent back\n");
  }

}

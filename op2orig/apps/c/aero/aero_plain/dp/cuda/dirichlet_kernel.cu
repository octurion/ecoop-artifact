//
// auto-generated by op2.py
//

//user function
__device__ void dirichlet_gpu( double *res) { *res = 0.0; 
}

// CUDA kernel function
__global__ void op_cuda_dirichlet(
  double *__restrict ind_arg0,
  const int *__restrict opDat0Map,
  int    block_offset,
  int   *blkmap,
  int   *offset,
  int   *nelems,
  int   *ncolors,
  int   *colors,
  int   nblocks,
  int   set_size) {

  __shared__ int    nelem, offset_b;

  extern __shared__ char shared[];

  if (blockIdx.x+blockIdx.y*gridDim.x >= nblocks) {
    return;
  }
  if (threadIdx.x==0) {

    //get sizes and shift pointers and direct-mapped data

    int blockId = blkmap[blockIdx.x + blockIdx.y*gridDim.x  + block_offset];

    nelem    = nelems[blockId];
    offset_b = offset[blockId];

  }
  __syncthreads(); // make sure all of above completed

  for ( int n=threadIdx.x; n<nelem; n+=blockDim.x ){
    int map0idx;
    map0idx = opDat0Map[n + offset_b + set_size * 0];


    //user-supplied kernel call
    dirichlet_gpu(ind_arg0+map0idx*1);
  }
}


//host stub function
void op_par_loop_dirichlet(char const *name, op_set set,
  op_arg arg0){

  int nargs = 1;
  op_arg args[1];

  args[0] = arg0;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(1);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[1].name      = name;
  OP_kernels[1].count    += 1;


  int    ninds   = 1;
  int    inds[1] = {0};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: dirichlet\n");
  }

  //get plan
  #ifdef OP_PART_SIZE_1
    int part_size = OP_PART_SIZE_1;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set->size > 0) {

    op_plan *Plan = op_plan_get(name,set,part_size,nargs,args,ninds,inds);

    //execute plan

    int block_offset = 0;
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==Plan->ncolors_core) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      #ifdef OP_BLOCK_SIZE_1
      int nthread = OP_BLOCK_SIZE_1;
      #else
      int nthread = OP_block_size;
      #endif

      dim3 nblocks = dim3(Plan->ncolblk[col] >= (1<<16) ? 65535 : Plan->ncolblk[col],
      Plan->ncolblk[col] >= (1<<16) ? (Plan->ncolblk[col]-1)/65535+1: 1, 1);
      if (Plan->ncolblk[col] > 0) {
        op_cuda_dirichlet<<<nblocks,nthread>>>(
        (double *)arg0.data_d,
        arg0.map_data_d,
        block_offset,
        Plan->blkmap,
        Plan->offset,
        Plan->nelems,
        Plan->nthrcol,
        Plan->thrcol,
        Plan->ncolblk[col],
        set->size+set->exec_size);

      }
      block_offset += Plan->ncolblk[col];
    }
    OP_kernels[1].transfer  += Plan->transfer;
    OP_kernels[1].transfer2 += Plan->transfer2;
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[1].time     += wall_t2 - wall_t1;
}

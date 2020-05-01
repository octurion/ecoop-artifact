//
// auto-generated by op2.py
//

//user function
//user function

void adt_calc_omp4_kernel(int *map0, int map0size, float *data4, int dat4size,
                          float *data5, int dat5size, float *data0,
                          int dat0size, int *col_reord, int set_size1,
                          int start, int end, int num_teams, int nthread);

// host stub function
void op_par_loop_adt_calc(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  int nargs = 6;
  op_arg args[6];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(1);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[1].name      = name;
  OP_kernels[1].count    += 1;

  int  ninds   = 1;
  int  inds[6] = {0,0,0,0,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: adt_calc\n");
  }

  // get plan
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  #ifdef OP_PART_SIZE_1
    int part_size = OP_PART_SIZE_1;
  #else
    int part_size = OP_part_size;
  #endif
  #ifdef OP_BLOCK_SIZE_1
    int nthread = OP_BLOCK_SIZE_1;
  #else
    int nthread = OP_block_size;
  #endif


  int ncolors = 0;
  int set_size1 = set->size + set->exec_size;

  if (set->size >0) {

    //Set up typed device pointers for OpenMP
    int *map0 = arg0.map_data_d;
     int map0size = arg0.map->dim * set_size1;

    float* data4 = (float*)arg4.data_d;
    int dat4size = getSetSizeFromOpArg(&arg4) * arg4.dat->dim;
    float* data5 = (float*)arg5.data_d;
    int dat5size = getSetSizeFromOpArg(&arg5) * arg5.dat->dim;
    float *data0 = (float *)arg0.data_d;
    int dat0size = getSetSizeFromOpArg(&arg0) * arg0.dat->dim;

    op_plan *Plan = op_plan_get_stage(name,set,part_size,nargs,args,ninds,inds,OP_COLOR2);
    ncolors = Plan->ncolors;
    int *col_reord = Plan->col_reord;

    // execute plan
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==1) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      int start = Plan->col_offsets[0][col];
      int end = Plan->col_offsets[0][col+1];

      adt_calc_omp4_kernel(map0, map0size, data4, dat4size, data5, dat5size,
                           data0, dat0size, col_reord, set_size1, start, end,
                           part_size != 0 ? (end - start - 1) / part_size + 1
                                          : (end - start - 1) / nthread,
                           nthread);
    }
    OP_kernels[1].transfer  += Plan->transfer;
    OP_kernels[1].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  if (OP_diags>1) deviceSync();
  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[1].time     += wall_t2 - wall_t1;
}

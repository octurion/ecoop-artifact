//
// auto-generated by op2.py
//

void res_omp4_kernel(double *data0, int dat0size, int *map1, int map1size,
                     float *arg3, float *data1, int dat1size, float *data2,
                     int dat2size, int *col_reord, int set_size1, int start,
                     int end, int num_teams, int nthread) {

  float arg3_l = *arg3;
  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size])\
    map(to:col_reord[0:set_size1],map1[0:map1size],data1[0:dat1size],data2[0:dat2size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map1idx = map1[n_op + set_size1 * 1];
    int map2idx = map1[n_op + set_size1 * 0];

    //variable mapping
    const double *A = &data0[3 * n_op];
    const float *u = &data1[2 * map1idx];
    float *du = &data2[3 * map2idx];
    const float *beta = &arg3_l;

    //inline function

    *du += (float)((*beta) * (*A) * (*u));
    //end inline func
  }

  *arg3 = arg3_l;
}

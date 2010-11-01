#include "operations.h"
#include "testing.h"


//operations on complex variables
hmc_complex complexconj(hmc_complex *in){
  hmc_complex z = *in;
  z.im = -z.im;
  return z;
}

hmc_complex complexmult(hmc_complex *a, hmc_complex *b){
  hmc_complex res;
  res.re = (*a).re*(*b).re - (*a).im*(*b).im;
  res.im = (*a).im*(*b).re + (*a).re*(*b).im;
  return res;
}

hmc_complex complexadd(hmc_complex *a, hmc_complex *b){
  hmc_complex res;
  res.re = (*a).re + (*b).re;
  res.im = (*a).im + (*b).im;
  return res;
}

hmc_error complexaccumulate(hmc_complex *inout, hmc_complex *incr){
  (*inout).re += (*incr).re;
  (*inout).im += (*incr).im;
  return HMC_SUCCESS;
}

hmc_error accumulate_su3matrix_prod(hmc_su3matrix *acc, hmc_su3matrix *mult){
  hmc_su3matrix tmp;
  multiply_su3matrices(&tmp,acc,mult);
  copy_su3matrix(acc,&tmp);
  return HMC_SUCCESS;
}


//spinor operations

hmc_error set_zero_spinor(hmc_full_spinor_field *field) {
  for (int j=0; j<NC*NSPIN; j++) {
    for (int n=0; n<VOLSPACE; n++) {
      for (int t=0; t<NTIME; t++) {
	(*field)[j][n][t].re=0;
	(*field)[j][n][t].im=0;
      }
    }
  }
  return HMC_SUCCESS;
}

hmc_float local_squarenorm(hmc_full_spinor_field *field, int spacepos, int timepos) {
  hmc_float sum=0;
  for (int j=0; j<NC*NSPIN; j++) {
    hmc_float dummy_re = (*field)[j][spacepos][timepos].re;
    hmc_float dummy_im = (*field)[j][spacepos][timepos].im;
    sum += dummy_re*dummy_re + dummy_im*dummy_im;
  }
  return sum;
}

hmc_float global_squarenorm(hmc_full_spinor_field *field) {
  hmc_float sum=0;
  for (int t=0; t<NTIME; t++) {
    for (int n=0; n<VOLSPACE; n++) {
      sum += local_squarenorm(field,n,t);
    }
  }
  return sum;
}

hmc_error fill_with_one(hmc_full_spinor_field *field, int spacepos, int timepos, int j){
  (*field)[j][spacepos][timepos].re = hmc_one_f;
  (*field)[j][spacepos][timepos].im = 0;
  return HMC_SUCCESS;
}


//gaugefield operations

hmc_error set_gaugefield_cold(hmc_gaugefield * field) {
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int mu=0; mu<NDIM; mu++) {
	hmc_su3matrix tmp;
	unit_su3matrix(&tmp);
	put_su3matrix(field, &tmp, n, t, mu);
      }
    }
  }
  return HMC_SUCCESS;
}

hmc_error adjoin_su3(hmc_gaugefield * in, hmc_gaugefield * out){
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int mu=0; mu<NDIM; mu++) {
	hmc_su3matrix tmp;
	get_su3matrix(&tmp, in, n, t, mu);
	adjoin_su3matrix(&tmp);
	put_su3matrix(out, &tmp, n, t, mu);
      }
    }
  }
  return HMC_SUCCESS;
}

hmc_complex global_trace_su3(hmc_gaugefield * field, int mu) {
  hmc_complex sum;
  sum.re = 0;
  sum.im = 0;
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      hmc_su3matrix tmp;
      get_su3matrix(&tmp, field, n, t, mu);
      sum.re += trace_su3matrix(&tmp).re;
      sum.im += trace_su3matrix(&tmp).im;
    }
  }
  return sum;
}


//operations that contain explicit SU(3) indices!!!

hmc_error copy_su3matrix(hmc_su3matrix *out, hmc_su3matrix *in){
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*out)[a][b] = (*in)[a][b];
    }
  }
  return HMC_SUCCESS;
}

hmc_error unit_su3matrix(hmc_su3matrix * u){
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      if(a!=b) {
	(*u)[a][b].re = 0;
      } else {
	(*u)[a][b].re = hmc_one_f;
      }
      (*u)[a][b].im = 0;
    }
  }
  return HMC_SUCCESS;
}

hmc_error multiply_su3matrices(hmc_su3matrix *out, hmc_su3matrix *p, hmc_su3matrix *q){
  for(int i=0; i<NC; i++) {
    for(int k=0; k<NC; k++) {
      (*out)[i][k].re=0;
      (*out)[i][k].im=0;
      for(int j=0;j<NC;j++) {
	hmc_complex tmp = complexmult(&(*p)[i][j],&(*q)[j][k]);
	complexaccumulate(&(*out)[i][k],&tmp);
      }
    }
  }
  return HMC_SUCCESS;
}


hmc_error adjoin_su3matrix(hmc_su3matrix * mat){
  hmc_su3matrix tmp;
  copy_su3matrix(&tmp, mat);
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*mat)[a][b] = complexconj(&(tmp[b][a]));
    }
  }
  return HMC_SUCCESS;
}

hmc_complex trace_su3matrix(hmc_su3matrix * mat){
  hmc_complex trace;
  trace.re=0;
  trace.im=0;
  for(int a=0; a<NC; a++) complexaccumulate(&trace,&((*mat)[a][a]));;
  return trace;
}

hmc_error get_su3matrix(hmc_su3matrix * out, hmc_gaugefield * in, int spacepos, int timepos, int mu) {
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*out)[a][b] = (*in)[a][b][mu][spacepos][timepos];
    }
  }
  return HMC_SUCCESS;
}

hmc_error put_su3matrix(hmc_gaugefield * field, hmc_su3matrix * in, int spacepos, int timepos, int mu) {
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*field)[a][b][mu][spacepos][timepos] = (*in)[a][b];
    }
  }
  return HMC_SUCCESS;
}

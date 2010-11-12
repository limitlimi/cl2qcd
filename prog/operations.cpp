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

hmc_complex complexsubtract(hmc_complex *a, hmc_complex *b){
  hmc_complex res;
  res.re = (*a).re - (*b).re;
  res.im = (*a).im - (*b).im;
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

hmc_error set_gaugefield_source(hmc_gaugefield * gaugefield, hmc_float * gaugefield_tmp, int check){
  //little check if arrays are big enough
  if (VOL4D*NDIM*NC*NC*2 != check){
    std::cout << "error in setting gaugefield to source values!! "<< std::endl << "Check global settings!!" << std::endl << std::endl;
    return HMC_STDERR;
  }

  int cter=0;
  //our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
  for (int t = 0; t<NTIME; t++){
  //if the alg is known to be correct, the next three for-loops could also be changed to one
    for (int i = 0; i<NSPACE; i++){
      for (int j = 0; j<NSPACE; j++){
        for (int k = 0; k<NSPACE; k++){
          for (int l = 0; l<NDIM; l++){
            //topology right?? What is the fastest index?? It should be colour, than dirac, than spatial, then time
            int spacepos = k + j*NSPACE + i*NSPACE*NSPACE;
            int globalpos = l + spacepos*NDIM + t*VOLSPACE*NDIM;
#ifdef _RECONSTRUCT_TWELVE_
            for (int m = 0; m<NC; m++){
              for (int n = 0; n<NC; n++){
		int ncindex = n + (NC-1)*m;
		//ildg-std: [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
		//which is stored in one single array here
		//skip NC*NC*2 cmplx numbers
		int pos = 2*n + 2*m*NC + globalpos*NC*NC*2;
		if(n<NC-1) {
		  (*gaugefield)[ncindex][(l+1)%NDIM][spacepos][t].re = gaugefield_tmp[pos];
		  (*gaugefield)[ncindex][(l+1)%NDIM][spacepos][t].im = gaugefield_tmp[pos + 1];
		}
		cter++;
	      }}
#else
            for (int m = 0; m<NC; m++){
              for (int n = 0; n<NC; n++){
                //ildg-std: [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
                //which is stored in one single array here
                //skip NC*NC*2 cmplx numbers
                int pos = 2*n + 2*m*NC + globalpos*NC*NC*2;
                (*gaugefield)[n][m][(l+1)%NDIM][spacepos][t].re = gaugefield_tmp[pos];
                (*gaugefield)[n][m][(l+1)%NDIM][spacepos][t].im = gaugefield_tmp[pos + 1];
                cter++;
	      }}
#endif
  }}}}}

  if(cter*2 != check) {
    std::cout << "error in setting gaugefield to source values! there were " << cter*2 << " vals set and not " << check << std::endl;
    return HMC_STDERR;
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

hmc_complex det_su3matrix(hmc_su3matrix * U){
#ifdef _RECONSTRUCT_TWELVE_
  hmc_complex det;
  det.re=0;
  det.im=0;
  hmc_complex subdet;
  subdet.re=0;
  subdet.im=0;
  hmc_complex tmp1;
  hmc_complex tmp2;
  hmc_complex tmp3;
  tmp1 = complexmult( &(*U)[0], &(*U)[3] );
  tmp2 = reconstruct_su3(U,2);
  tmp3 = complexmult( &tmp1, &tmp2 );
  complexaccumulate(&det,&tmp3);

  tmp1 = complexmult( &(*U)[2], &(*U)[5] );
  tmp2 = reconstruct_su3(U,0);
  tmp3 = complexmult( &tmp1, &tmp2 );
  complexaccumulate(&det,&tmp3);

  tmp1 = complexmult( &(*U)[4], &(*U)[1] );
  tmp2 = reconstruct_su3(U,1);
  tmp3 = complexmult( &tmp1, &tmp2 );
  complexaccumulate(&det,&tmp3);

  tmp1 = complexmult( &(*U)[3], &(*U)[4] );
  tmp2 = reconstruct_su3(U,0);
  tmp3 = complexmult( &tmp1, &tmp2 );
  complexaccumulate(&subdet,&tmp3);

  tmp1 = complexmult( &(*U)[5], &(*U)[0] );
  tmp2 = reconstruct_su3(U,1);
  tmp3 = complexmult( &tmp1, &tmp2 );
  complexaccumulate(&subdet,&tmp3);

  tmp1 = complexmult( &(*U)[1], &(*U)[2] );
  tmp2 = reconstruct_su3(U,2);
  tmp3 = complexmult( &tmp1, &tmp2 );
  complexaccumulate(&subdet,&tmp3);

  det.re -= subdet.re;
  det.im -= subdet.im;

#else
  hmc_complex det, det1, det2, det3, det4, det5, det6, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
  det.re=0;
  det.im=0;
  tmp1 = complexmult( &(*U)[1][1], &(*U)[2][2] );
  det1 = complexmult( &(*U)[0][0] , &tmp1);
  tmp2 = complexmult( &(*U)[1][2], &(*U)[2][0] );
  det2 = complexmult( &(*U)[0][1] , &tmp2);
  tmp3 = complexmult( &(*U)[1][0], &(*U)[2][1] );
  det3 = complexmult( &(*U)[0][2] , &tmp3);
  tmp4 = complexmult( &(*U)[1][1], &(*U)[2][0] );
  det4 = complexmult( &(*U)[0][2] , &tmp4);
  tmp5 = complexmult( &(*U)[1][0], &(*U)[2][2] );
  det5 = complexmult( &(*U)[0][1] , &tmp5);
  tmp6 = complexmult( &(*U)[1][2], &(*U)[2][1] );
  det6 = complexmult( &(*U)[0][0] , &tmp6);

det.re = det1.re + det2.re + det3.re - det4.re - det5.re - det6.re;
det.im = det1.im + det2.im + det3.im - det4.im - det5.im - det6.im;


//LZ: what is that?
 /*
U[0]*U[4]*U[8] + U[1]*U[5]*U[6] + U[2]*U[3]*U[7] - U[2]*U[4]*U[6] - U[1]*U[3]*U[8] - U[0]*U[5]*U[7]

(*U)[0][0]*(*U)[1][1]*(*U)[2][2]
(*U)[0][1]*(*U)[1][2]*(*U)[2][0] + 
(*U)[0][2]*(*U)[1][0]*(*U)[2][1] - 
(*U)[0][2]*(*U)[1][1]*(*U)[2][0] - 
(*U)[0][1]*(*U)[1][0]*(*U)[2][2] - 
(*U)[0][0]*(*U)[1][2]*(*U)[2][1];*/

#endif
  return det;
}


hmc_error copy_su3matrix(hmc_su3matrix *out, hmc_su3matrix *in){
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) {
    (*out)[n] = (*in)[n];
  }
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*out)[a][b] = (*in)[a][b];
    }
  }
#endif
  return HMC_SUCCESS;
}

hmc_error unit_su3matrix(hmc_su3matrix * u){
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) {
    if( n%NC == 0) {
      (*u)[n].re = hmc_one_f;
    } else {
      (*u)[n].re = 0;
    }
    (*u)[n].im = 0;
  }
#else
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
#endif
  return HMC_SUCCESS;
}

#ifdef _RECONSTRUCT_TWELVE_
hmc_complex reconstruct_su3(hmc_su3matrix *in, int ncomp){
  int jplusone = (ncomp+1)%NC;
  int jplustwo = (ncomp+2)%NC;
  hmc_complex first = complexmult(&((*in)[(NC-1)*jplusone]),&((*in)[1+(NC-1)*jplustwo]));
  hmc_complex second = complexmult(&((*in)[(NC-1)*jplustwo]),&((*in)[1+(NC-1)*jplusone]));
  hmc_complex result = complexsubtract(&first,&second);
  return complexconj(&result);
}
#endif


hmc_error multiply_su3matrices(hmc_su3matrix *out, hmc_su3matrix *p, hmc_su3matrix *q){
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) {
      (*out)[n].re=0;
      (*out)[n].im=0;
      for(int j=0;j<NC;j++) {
	int k = (int)(n/(NC-1));
	int i = n - (NC-1)*k;
	int np = i + (NC-1)*j;
	hmc_complex qcomponent;
	if(j==2) {
	  int nq = j + (NC-1)*k;
	  qcomponent = (*q)[nq];
	} else {
	  qcomponent = reconstruct_su3(q,k);
	}
	hmc_complex tmp = complexmult(&(*p)[np],&qcomponent);
	complexaccumulate(&(*out)[n],&tmp);
      }
    }
#else
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
#endif
  return HMC_SUCCESS;
}


hmc_error adjoin_su3matrix(hmc_su3matrix * mat){
#ifdef _RECONSTRUCT_TWELVE_
  hmc_su3matrix tmp;
  copy_su3matrix(&tmp, mat);
  for(int n=0; n<NC*(NC-1); n++) {
    int j = (int)(n/(NC-1));
    int i = n - (NC-1)*j;
    hmc_complex element;
    if ( j==2 ) {
      element = reconstruct_su3(&tmp,i);
    } else {
      int nnew = j + (NC-1)*i;
      element = tmp[nnew];
    }
    (*mat)[n] = complexconj(&element);
  }
#else
  hmc_su3matrix tmp;
  copy_su3matrix(&tmp, mat);
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*mat)[a][b] = complexconj(&(tmp[b][a]));
    }
  }
#endif
  return HMC_SUCCESS;
}

hmc_complex trace_su3matrix(hmc_su3matrix * mat){
  hmc_complex trace;
#ifdef _RECONSTRUCT_TWELVE_
  trace = reconstruct_su3(mat,NC-1);
  for(int n=0; n<(NC-1); n++) complexaccumulate(&trace,&((*mat)[NC*n]));
#else
  trace.re=0;
  trace.im=0;
  for(int a=0; a<NC; a++) complexaccumulate(&trace,&((*mat)[a][a]));;
#endif
  return trace;
}

hmc_error get_su3matrix(hmc_su3matrix * out, hmc_gaugefield * in, int spacepos, int timepos, int mu) {
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) (*out)[n] = (*in)[n][mu][spacepos][timepos];
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*out)[a][b] = (*in)[a][b][mu][spacepos][timepos];
    }
  }
#endif
  return HMC_SUCCESS;
}

hmc_error put_su3matrix(hmc_gaugefield * field, hmc_su3matrix * in, int spacepos, int timepos, int mu) {
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) (*field)[n][mu][spacepos][timepos] = (*in)[n];
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*field)[a][b][mu][spacepos][timepos] = (*in)[a][b];
    }
  }
#endif
  return HMC_SUCCESS;
}

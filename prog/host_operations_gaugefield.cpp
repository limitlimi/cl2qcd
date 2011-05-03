#include "host_operations_gaugefield.h"

hmc_error copy_to_ocl_format(hmc_ocl_gaugefield* host_gaugefield,hmc_gaugefield* gaugefield){
  for(int spacepos=0; spacepos<NSPACE*NSPACE*NSPACE; spacepos++) {
    for(int t=0; t<NTIME; t++) {
      for(int mu=0; mu<NDIM; mu++) {
	hmc_su3matrix tmp;
	get_su3matrix(&tmp, gaugefield, spacepos, t, mu);
	for(int b=0; b<NC; b++) {
#ifdef _RECONSTRUCT_TWELVE_
	  for(int a=0; a<NC-1; a++) {
	    int n = a + (NC-1)*b;
	    host_gaugefield[ocl_gaugefield_element(0,a,b,mu,spacepos,t)] = tmp[n].re;
	    host_gaugefield[ocl_gaugefield_element(1,a,b,mu,spacepos,t)] = tmp[n].im;
	  }
#else
	  for(int a=0; a<NC; a++) {
	    host_gaugefield[ocl_gaugefield_element(0,a,b,mu,spacepos,t)] = tmp[a][b].re;
	    host_gaugefield[ocl_gaugefield_element(1,a,b,mu,spacepos,t)] = tmp[a][b].im;
	  }
#endif
	}
      }
    }
  }
  return HMC_SUCCESS;
}

hmc_error copy_from_ocl_format(hmc_gaugefield* gaugefield,hmc_ocl_gaugefield* host_gaugefield){
  for(int spacepos=0; spacepos<NSPACE*NSPACE*NSPACE; spacepos++) {
    for(int t=0; t<NTIME; t++) {
      for(int mu=0; mu<NDIM; mu++) {
				hmc_su3matrix tmp;
				for(int b=0; b<NC; b++) {
#ifdef _RECONSTRUCT_TWELVE_
	  			for(int a=0; a<NC-1; a++) {
	    			int n = a + (NC-1)*b;
	    			tmp[n].re = host_gaugefield[ocl_gaugefield_element(0,a,b,mu,spacepos,t)];
	    			tmp[n].im = host_gaugefield[ocl_gaugefield_element(1,a,b,mu,spacepos,t)];
	  			}
#else
	  			for(int a=0; a<NC; a++) {
						tmp[a][b].re = host_gaugefield[ocl_gaugefield_element(0,a,b,mu,spacepos,t)];
						tmp[a][b].im = host_gaugefield[ocl_gaugefield_element(1,a,b,mu,spacepos,t)];
	  			}
#endif
	  	put_su3matrix(gaugefield, &tmp, spacepos, t, mu);
				}}}}
  return HMC_SUCCESS;
}
 
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

hmc_error set_gaugefield_hot(hmc_gaugefield * field) {
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int mu=0; mu<NDIM; mu++) {
	hmc_su3matrix tmp;
	random_su3matrix(&tmp);
	put_su3matrix(field, &tmp, n, t, mu);
      }
    }
  }
  return HMC_SUCCESS;
}

hmc_error copy_gaugefield_from_ildg_format(hmc_gaugefield * gaugefield, hmc_float * gaugefield_tmp, int check){
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
            int spacepos = k + j*NSPACE + i*NSPACE*NSPACE;
            int globalpos = l + spacepos*NDIM + t*VOLSPACE*NDIM;
#ifdef _RECONSTRUCT_TWELVE_
            for (int m = 0; m<NC; m++){
              for (int n = 0; n<NC; n++){
		int ncindex = m + (NC-1)*n;
		//ildg-std: [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
		//which is stored in one single array here
		//skip NC*NC*2 cmplx numbers
		int pos = 2*n + 2*m*NC + globalpos*NC*NC*2;
		if(m<NC-1) {
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
                (*gaugefield)[m][n][(l+1)%NDIM][spacepos][t].re = gaugefield_tmp[pos];
                (*gaugefield)[m][n][(l+1)%NDIM][spacepos][t].im = gaugefield_tmp[pos + 1];
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

hmc_error copy_gaugefield_to_ildg_format(ildg_gaugefield * dest, hmc_gaugefield * source){
  
  int cter=0;
  //our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
  for (int t = 0; t<NTIME; t++){
  //if the alg is known to be correct, the next three for-loops could also be changed to one
    for (int i = 0; i<NSPACE; i++){
      for (int j = 0; j<NSPACE; j++){
        for (int k = 0; k<NSPACE; k++){
          for (int l = 0; l<NDIM; l++){
            int spacepos = k + j*NSPACE + i*NSPACE*NSPACE;
            int globalpos = l + spacepos*NDIM + t*VOLSPACE*NDIM;
#ifdef _RECONSTRUCT_TWELVE_
            for (int m = 0; m<NC; m++){
              for (int n = 0; n<NC; n++){
		int ncindex = m + (NC-1)*n;
		//ildg-std: [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
		//which is stored in one single array here
		//skip NC*NC*2 cmplx numbers
		int pos = 2*n + 2*m*NC + globalpos*NC*NC*2;
		if(m<NC-1) {
		  (dest[0])[pos]     = (*source)[ncindex][(l+1)%NDIM][spacepos][t].re;
		  (dest[0])[pos + 1] = (*source)[ncindex][(l+1)%NDIM][spacepos][t].im;
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
                (dest[0])[pos]     = ((*source)[m][n][(l+1)%NDIM][spacepos][t]).re;
                (dest[0])[pos + 1] = ((*source)[m][n][(l+1)%NDIM][spacepos][t]).im;
                cter++;
	      }}
#endif
  }}}}}
  
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

#endif
  return det;
}


hmc_error project_su3(hmc_su3matrix *U){

  //Extract initial vectors
  hmc_complex a[NC];
  hmc_complex b[NC];
  hmc_complex c[NC];
#ifdef _RECONSTRUCT_TWELVE_
  a[0] = (*U)[0];
  a[1] = (*U)[1];
  a[2] = reconstruct_su3(U,0);
  b[0] = (*U)[2];
  b[1] = (*U)[3];
  b[2] = reconstruct_su3(U,1);
  c[0] = (*U)[4];
  c[1] = (*U)[5];
  c[2] = reconstruct_su3(U,2);
#else
    for (int i = 0; i<NC; i++){
     a[i] = (*U)[i][0];
     b[i] = (*U)[i][1];
     c[i] = (*U)[i][2];
    }
#endif

  //New SU3-Matrix
  //first vector
  //norm
  hmc_float norm = 0.;
  for (int i=0; i<NC; i++){
    hmc_complex tmp = complexconj(&(a[i]));
    tmp = complexmult (& a[i], & tmp);
    norm += tmp.re;
  }
  norm = 1./sqrt(norm);
  //rescale
  for (int i=0; i<NC; i++){
    //perhaps define a new complex-function for multiplying with a real number
    a[i].re *= norm;
    a[i].im *= norm;
  }
  
  //second vector
  //orthogonal vector
  hmc_complex factor;
  factor.re = 0.0;
  factor.im = 0.0;
  for (int i=0; i<NC; i++){
    hmc_complex tmp;
    tmp = complexconj (&(b[i]));
    tmp = complexmult (&(a[i]), &tmp);
    factor = complexadd (&factor, &tmp);
  }
  for (int i=0; i<NC; i++){
    hmc_complex tmp;
    tmp = complexmult(&factor, &(a[i]));
    b[i] = complexsubtract(&(b[i]), &tmp); 
  }
  
//norm
  norm = 0;
  for (int i=0; i<NC; i++)
  {
    hmc_complex tmp;
    tmp = complexconj(&(b[i]));
    tmp = complexmult (&(b[i]), &tmp);
    norm +=  tmp.re;
  }
  norm = 1./sqrt(norm);
  //rescale
  for  (int i=0; i<NC; i++){
    b[i].re *= norm;
    b[i].im *= norm;
  }
  
  //third vector 
  //orthogonal vector
  hmc_complex tmp;
  hmc_complex tmp2;
  tmp = complexmult(&(a[1]), &(b[2]));
  tmp = complexconj(&tmp);
  tmp2 = complexmult(&(a[2]), &(b[1]));
  tmp2 = complexconj(&tmp2);
  c[0] = complexsubtract(&tmp, &tmp2);
  tmp = complexmult(&(a[2]), &(b[0]));
  tmp = complexconj(&tmp);
  tmp2 = complexmult(&(a[0]), &(b[2]));
  tmp2 = complexconj(&tmp2);
  c[1] = complexsubtract(&tmp, &tmp2);
#ifndef _RECONSTRUCT_TWELVE_
  tmp = complexmult(&(a[0]), &(b[1]));
  tmp = complexconj(&tmp);
  tmp2 = complexmult(&(a[1]), &(b[0]));
  tmp2 = complexconj(&tmp2);
  c[2] = complexsubtract(&tmp, &tmp2);
#endif

#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC-1; n++) { 
    (*U)[n] = a[n];
    (*U)[n+NC-1] = b[n];
    (*U)[n+NC+NC-2] = c[n];
  }
  #else
  for(int i=0; i<NC; i++) {
      (*U)[i][0] = a[i];
      (*U)[i][1] = b[i];
      (*U)[i][2] = c[i];
  }
#endif
  return HMC_SUCCESS;
}

hmc_error project_su3_old(hmc_su3matrix *U){
//old code
  hmc_complex det = det_su3matrix(U);
  hmc_float detsqunorm = det.re*det.re + det.im*det.im;

  hmc_float phi;
  if(det.re*det.re<projectioneps) { 
  phi = PI/2.;
  } else {
  phi = atan(det.im/det.re);
  if(det.re<0) phi += PI;
  }

  hmc_complex norm;
  norm.re = pow(detsqunorm,hmc_one_f/6.)*cos(phi/3.);
  norm.im = pow(detsqunorm,hmc_one_f/6.)*sin(phi/3.);
  
  hmc_float normsqunorm = norm.re*norm.re+norm.im*norm.im;

  #ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) { 
    hmc_complex tmp = (*U)[n];
    (*U)[n].re = (tmp.re*norm.re+tmp.im*norm.im)/normsqunorm;
    (*U)[n].im = (tmp.im*norm.re-tmp.re*norm.im)/normsqunorm;
  }
  #else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      hmc_complex tmp = (*U)[a][b]; 
      (*U)[a][b].re = (tmp.re*norm.re+tmp.im*norm.im)/normsqunorm;
      (*U)[a][b].im = (tmp.im*norm.re-tmp.re*norm.im)/normsqunorm;
    }
  }
  #endif
  return HMC_SUCCESS;
}

/** @todo memcpy ... */
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

/** @todo memcpy ... */
hmc_error copy_staplematrix(hmc_staplematrix *out, hmc_staplematrix *in){
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*NC; n++) {
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

/** @todo memset ... */
hmc_error zero_su3matrix(hmc_su3matrix * u){
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) {
    (*u)[n].re = 0;
    (*u)[n].im = 0;
  }
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*u)[a][b].re = 0;
      (*u)[a][b].im = 0;
    }
  }
#endif
  return HMC_SUCCESS;
}

/** @todo memset ... */
hmc_error zero_staplematrix(hmc_staplematrix * u){
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*NC; n++) {
    (*u)[n].re = 0;
    (*u)[n].im = 0;
  }
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*u)[a][b].re = 0;
      (*u)[a][b].im = 0;
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

hmc_error random_su3matrix(hmc_su3matrix * u){
  printf("random su3matrix needs to be implemented...\n");
  exit(HMC_UNDEFINEDERROR);
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
	  qcomponent = reconstruct_su3(q,k);
	} else {
	  int nq = j + (NC-1)*k;
	  qcomponent = (*q)[nq];
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

// wird wohl doch nicht gebraucht...
// hmc_error add_su3matrices(hmc_su3matrix *p, hmc_su3matrix *q){
// #ifdef _RECONSTRUCT_TWELVE
//   for(int n=0; n<NC*(NC-1); n++) {
//       (*p)[n].re = (*p)[n].re + (*q)[n].re;
//       (*p)[n].im = (*p)[n].im + (*q)[n].im;
// #else
//   for(int i=0; i<NC; i++) {
//     for(int k=0; k<NC; k++) {
//       (*p)[i][k].re = (*p)[i][k].re + (*q)[i][k].re;
//       (*p)[i][k].im = (*p)[i][k].im + (*q)[i][k].im;
//     }
//   }
// #endif
// }


hmc_error multiply_staplematrix(hmc_staplematrix *out, hmc_su3matrix *p, hmc_staplematrix *q){
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
// 	  qcomponent = reconstruct_su3(q,k);
          qcomponent = (*q)[NC*(NC-1)+k];
	} else {
	  int nq = j + (NC-1)*k;
	  qcomponent = (*q)[nq];
	}
	hmc_complex tmp = complexmult(&(*p)[np],&qcomponent);
	complexaccumulate(&(*out)[n],&tmp);
      }
    }
    //the left components:
    hmc_complex X = reconstruct_su3(p,0);
    hmc_complex Y = reconstruct_su3(p,1);
    hmc_complex Z = reconstruct_su3(p,2);
    hmc_complex tmp;
    ((*out)[6]).re=0;
    ((*out)[6]).im=0;
    ((*out)[7]).re=0;
    ((*out)[7]).im=0;
    ((*out)[8]).re=0;
    ((*out)[8]).im=0;
    
    tmp = complexmult(&X,&(*q)[0]);
    complexaccumulate(&(*out)[6],&tmp);
    tmp = complexmult(&Y,&(*q)[1]);
    complexaccumulate(&(*out)[6],&tmp);
    tmp = complexmult(&Z,&(*q)[6]);
    complexaccumulate(&(*out)[6],&tmp);

    tmp = complexmult(&X,&(*q)[2]);
    complexaccumulate(&(*out)[7],&tmp);
    tmp = complexmult(&Y,&(*q)[3]);
    complexaccumulate(&(*out)[7],&tmp);
    tmp = complexmult(&Z,&(*q)[7]);
    complexaccumulate(&(*out)[7],&tmp);

    tmp = complexmult(&X,&(*q)[4]);
    complexaccumulate(&(*out)[8],&tmp);
    tmp = complexmult(&Y,&(*q)[5]);
    complexaccumulate(&(*out)[8],&tmp);
    tmp = complexmult(&Z,&(*q)[8]);
    complexaccumulate(&(*out)[8],&tmp);
    
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

hmc_error accumulate_su3matrices_add(hmc_staplematrix *p, hmc_su3matrix *q){
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) {
    complexaccumulate(&(*p)[n], &(*q)[n]);
  }
  for(int n=NC*(NC-1);  n<NC*NC; n++) {
    hmc_complex tmp = reconstruct_su3(q, n-NC*(NC-1)); 
    complexaccumulate(&(*p)[n], &(tmp));
  }  
#else

  for(int i=0; i<NC; i++) {
    for(int k=0; k<NC; k++) {
      complexaccumulate(&(*p)[i][k],&(*q)[i][k]);
    }
  }
#endif
  return HMC_SUCCESS;
}

hmc_error accumulate_su3matrix_prod(hmc_su3matrix *acc, hmc_su3matrix *mult){
  hmc_su3matrix tmp;
  multiply_su3matrices(&tmp,acc,mult);
  copy_su3matrix(acc,&tmp);
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

void reduction (hmc_complex dest[su2_entries], hmc_staplematrix src, const int rand){
#ifdef _RECONSTRUCT_TWELVE_
  if(rand == 1)
  { 
    dest[0] = src[0]; 
    dest[1] = src[2]; 
    dest[2] = src[1]; 
    dest[3] = src[3]; 
  } 
  else if (rand==2) 
  { 
    dest[0] = src[3]; 
    dest[1] = src[5]; 
    dest[2] = src[7];
    dest[3] = src[8];
  } 
  else if (rand==3)
  { 
    dest[0] = src[0]; 
    dest[1] = src[4]; 
    dest[2] = src[6];
    dest[3] = src[8];
  } 
  else std::cout<<"error at reduction, rand not 1,2,3"<<std::endl; 
#else
  if(rand == 1)
  {
    dest[0] = src[0][0];
    dest[1] = src[0][1];
    dest[2] = src[1][0];
    dest[3] = src[1][1];
  }
  else if (rand==2)
  {
    dest[0] = src[1][1];
    dest[1] = src[1][2];
    dest[2] = src[2][1];
    dest[3] = src[2][2];
  }
  else if (rand==3)
  {
    dest[0] = src[0][0];
    dest[1] = src[0][2];
    dest[2] = src[2][0];
    dest[3] = src[2][2];
  }
  else
    std::cout<<"error at reduction, rand not 1,2,3"<<std::endl;
#endif
}

// return an SU2 matrix (std basis) extended to SU3 (std basis)
void extend (hmc_su3matrix * dest, const int random, hmc_complex src[su2_entries]){
#ifdef _RECONSTRUCT_TWELVE_
  if (random == 1){
    (*dest)[0] = src[0];
    (*dest)[2] = src[1];
    (*dest)[4] = hmc_complex_zero;
    (*dest)[1] = src[2];
    (*dest)[3] = src[3];
    (*dest)[5] = hmc_complex_zero;
  }
  else if (random == 2){
    (*dest)[0] = hmc_complex_one;
    (*dest)[2] = hmc_complex_zero;
    (*dest)[4] = hmc_complex_zero;
    (*dest)[1] = hmc_complex_zero;
    (*dest)[3] = src[0];
    (*dest)[5] = src[1];
  }
  else if (random == 3){
    (*dest)[0] = src[0];
    (*dest)[2] = hmc_complex_zero;
    (*dest)[4] = src[1];
    (*dest)[1] = hmc_complex_zero;
    (*dest)[3] = hmc_complex_one;
    (*dest)[5] = hmc_complex_zero;
  }
  else
    std::cout<<"error at extend, random not 1,2,3"<<std::endl;

#else
  if (random == 1){
    (*dest)[0][0] = src[0];
    (*dest)[0][1] = src[1];
    (*dest)[0][2] = hmc_complex_zero;
    (*dest)[1][0] = src[2];
    (*dest)[1][1] = src[3];
    (*dest)[1][2] = hmc_complex_zero;
    (*dest)[2][0] = hmc_complex_zero;
    (*dest)[2][1] = hmc_complex_zero;
    (*dest)[2][2] = hmc_complex_one;
  }
  else if (random == 2){
    (*dest)[0][0] = hmc_complex_one;
    (*dest)[0][1] = hmc_complex_zero;
    (*dest)[0][2] = hmc_complex_zero;
    (*dest)[1][0] = hmc_complex_zero;
    (*dest)[1][1] = src[0];
    (*dest)[1][2] = src[1];
    (*dest)[2][0] = hmc_complex_zero;
    (*dest)[2][1] = src[2];
    (*dest)[2][2] = src[3];
  }
  else if (random == 3){
    (*dest)[0][0] = src[0];
    (*dest)[0][1] = hmc_complex_zero;
    (*dest)[0][2] = src[1];
    (*dest)[1][0] = hmc_complex_zero;
    (*dest)[1][1] = hmc_complex_one;
    (*dest)[1][2] = hmc_complex_zero;
    (*dest)[2][0] = src[2];
    (*dest)[2][1] = hmc_complex_zero;
    (*dest)[2][2] = src[3];
  }
  else
    std::cout<<"error at extend, random not 1,2,3"<<std::endl;
#endif
}

void gaugefield_apply_bc(hmc_su3matrix * in, hmc_float theta){
  hmc_float tmp1,tmp2;
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) {
    tmp1 = ((*in)[n]).re;
    tmp2 = ((*in)[n]).im;
    ((*in)[n]).re = cos(theta)*tmp1 - sin(theta)*tmp2;
    ((*in)[n]).im = sin(theta)*tmp1 + cos(theta)*tmp2;
  }
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      tmp1 = ((*in)[a][b]).re;
      tmp2 = ((*in)[a][b]).im;
      ((*in)[a][b]).re = cos(theta)*tmp1 - sin(theta)*tmp2;
      ((*in)[a][b]).im = sin(theta)*tmp1 + cos(theta)*tmp2;
    }
  }
#endif
  return;
}

// replace link in with e^(mu.re)*(cos(mu.im) + i*sin(mu.im))
void gaugefield_apply_chem_pot(hmc_su3matrix * u, hmc_su3matrix * udagger, hmc_float chem_pot_re, hmc_float chem_pot_im){
  hmc_float tmp1,tmp2;
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) {
    tmp1 = ((*u)[n]).re;
    tmp2 = ((*u)[n]).im;
    ((*u)[n]).re = exp(chem_pot_re)*( cos(chem_pot_im)*tmp1 - sin(chem_pot_im)*tmp2 );
    ((*u)[n]).im = exp(chem_pot_re)*( sin(chem_pot_im)*tmp1 + cos(chem_pot_im)*tmp2 );
    tmp1 = ((*udagger)[n]).re;
    tmp2 = ((*udagger)[n]).im;
    ((*udagger)[n]).re = exp(-chem_pot_re)*( cos(chem_pot_im)*tmp1 + sin(chem_pot_im)*tmp2 );
    ((*udagger)[n]).im = exp(-chem_pot_re)*( -sin(chem_pot_im)*tmp1 + cos(chem_pot_im)*tmp2 );
  }
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      tmp1 = ((*u)[a][b]).re;
      tmp2 = ((*u)[a][b]).im;
      ((*u)[a][b]).re = exp(chem_pot_re)*( cos(chem_pot_im)*tmp1 - sin(chem_pot_im)*tmp2 );
      ((*u)[a][b]).im = exp(chem_pot_re)*( sin(chem_pot_im)*tmp1 + cos(chem_pot_im)*tmp2 );
      tmp1 = ((*udagger)[a][b]).re;
      tmp2 = ((*udagger)[a][b]).im;
      ((*udagger)[a][b]).re = exp(-chem_pot_re)*( cos(chem_pot_im)*tmp1 + sin(chem_pot_im)*tmp2 );
      ((*udagger)[a][b]).im = exp(-chem_pot_re)*( -sin(chem_pot_im)*tmp1 + cos(chem_pot_im)*tmp2 );
    }
  }
#endif
  return;
}

void local_polyakov(hmc_gaugefield * field, hmc_su3matrix * prod, int n){
	unit_su3matrix(prod);
	for(int t=0; t<NTIME; t++) {
		hmc_su3matrix tmp;
		get_su3matrix(&tmp,field,n,t,0);
		accumulate_su3matrix_prod(prod,&tmp);
	}
	return;
}

void local_plaquette(hmc_gaugefield * field, hmc_su3matrix * prod, int n, int t, int mu, int nu ){
	hmc_su3matrix tmp;
	//u_mu(x)
	get_su3matrix(prod,field,n,t,mu);
	//u_nu(x+mu)
	if(mu==0) {
	  int newt = (t+1)%NTIME; //(haha)
	  get_su3matrix(&tmp,field,n,newt,nu);
	} else {
	  get_su3matrix(&tmp,field,get_neighbor(n,mu),t,nu);
	}
	accumulate_su3matrix_prod(prod,&tmp);
	//adjoint(u_mu(x+nu))
	if(nu==0) {
	  int newt = (t+1)%NTIME;
	  get_su3matrix(&tmp,field,n,newt,mu);
	} else {
	  get_su3matrix(&tmp,field,get_neighbor(n,nu),t,mu);
	}
	adjoin_su3matrix(&tmp);
	accumulate_su3matrix_prod(prod,&tmp);
	//adjoint(u_nu(x))
	get_su3matrix(&tmp,field,n,t,nu);
	adjoin_su3matrix(&tmp);
	accumulate_su3matrix_prod(prod,&tmp);
		
	return;
}

/** @todo memcpy ... */
hmc_error copy_gaugefield(hmc_gaugefield * source, hmc_gaugefield * dest){
	// copies source to destination within cpu memory, layer for gaugefield array
	return complexcopy((hmc_complex *)source, (hmc_complex *)dest, GAUGEFIELDSIZE); // SL: not tested
}


//gauge-momenta operations
//TODO CP: these should go into a seperate file like host_operations_gaugemomenta.cpp

/** @todo memcpy ... */
hmc_error copy_gaugemomenta(hmc_gauge_momentum * source, hmc_gauge_momentum * dest){
	// copies source to destination within cpu memory, layer for momentum array
	return complexcopy((hmc_complex *)source, (hmc_complex *)dest, GAUGEMOMENTASIZE); // SL: not tested
}

//gaugemomentum is just a hmc_float vector of length GAUGEMOMENTASIZE
hmc_error gaugemomenta_squarenorm(hmc_gauge_momentum * in, hmc_float * result){
	//make sure result is zero
	(*result) = 0.;
	hmc_float sum = 0.;
	for(int i = 0; i<GAUGEMOMENTASIZE; i++){
		sum += ((in)[i])*((in)[i]);
	}
	(*result) = sum;
	return HMC_SUCCESS;
}

/** @todo memset... */
hmc_error set_zero_gaugemomenta(hmc_gauge_momentum * in){
	for(int i = 0; i<GAUGEMOMENTASIZE; i++){
		(in[i]) = 0.;
	}
	return HMC_SUCCESS;
}

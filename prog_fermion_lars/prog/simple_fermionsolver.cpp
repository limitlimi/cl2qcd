#include "simple_fermionsolver.h"

hmc_error simple_correlator(hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, int cgmax){

  hmc_spinor_field in[SPINORFIELDSIZE];
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int a=0; a<NSPIN; a++) {
	for(int j=0; j<NC; j++) {
	  fill_with_one(in, n, t, a, j);
	}
      }
    }
  }
  hmc_float norm=global_squarenorm(in);
  norm = sqrt(norm);
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    in[n].re /= norm;
    in[n].im /=norm;
  }


  //pseudo scalar, flavour multiplet
  hmc_complex correlator_ps[NSPACE];
  for(int z=0; z<NSPACE; z++) {
    correlator_ps[z].re = 0;
    correlator_ps[z].im = 0;
  }

  hmc_spinor_field b[SPINORFIELDSIZE];
  hmc_spinor_field phi[SPINORFIELDSIZE];

  for(int k=0; k<NC*NSPIN; k++) {
    create_point_source(b,k,0,0,kappa,mu,gaugefield);
    simple_solver(in, phi, b, gaugefield, kappa, mu, cgmax);

    for(int timepos = 0; timepos<NTIME; timepos++) {
      for(int spacepos = 0; spacepos<VOLSPACE; spacepos++) {
	for(int alpha = 0; alpha<NSPIN; alpha++) {
	  for(int c = 0; c<NC; c++) {
	    //	    int j = spinor_element(alpha,c);
	    int n = spinor_field_element(alpha, c, spacepos, timepos);
	    int z = get_spacecoord(spacepos, 3);
	    hmc_complex tmp = phi[n];
	    hmc_complex ctmp = complexconj(&tmp);
	    hmc_complex incr = complexmult(&ctmp,&tmp);
	    correlator_ps[z].re += incr.re;
	    correlator_ps[z].im += incr.im;
	  }
	}
      }
    }
  }

  printf("pseudo scalar correlator:\n");
  for(int z=0; z<NSPACE; z++) {
    printf("%d\t(%e,%e)\n",z,correlator_ps[z].re,correlator_ps[z].im);
  }
  return HMC_SUCCESS;
}



hmc_error simple_solver(hmc_spinor_field* in, hmc_spinor_field* out, hmc_spinor_field* b, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, int cgmax){

  //convert to kappa format
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    in[n].re *= sqrt(2.*kappa);
    in[n].im *= sqrt(2.*kappa);
  }

  //  create_point_source(in,4,1,0,kappa,mu,gaugefield);
  
  simple_solve_bicgstab(in, b, gaugefield, kappa, mu, cgmax);

  //convert from kappa format
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    out[n].re = in[n].re/sqrt(2.*kappa);
    out[n].im = in[n].im/sqrt(2.*kappa);
  }

  return HMC_SUCCESS;
}


hmc_error simple_solve_bicgstab(hmc_spinor_field* inout, hmc_eoprec_spinor_field* source, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, int cgmax){

  //BiCGStab according to hep-lat/9404013

  hmc_spinor_field* rn = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* rhat = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* v = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* p = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* s = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* t = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_complex rho;
  hmc_complex rho_next;
  hmc_complex alpha;
  hmc_complex omega;
  hmc_complex beta;

  hmc_complex tmp1;
  hmc_complex tmp2;
  hmc_complex tmp3;


  for(int iter=0; iter<cgmax; iter++){

    if(iter%iter_refresh==0) {
      //fresh start
      M(inout,rn,gaugefield,kappa,mu);
      for (int n=0; n<SPINORFIELDSIZE; n++) {
	rn[n].re = source[n].re - rn[n].re;
	rn[n].im = source[n].im - rn[n].im;
	rhat[n].re = rn[n].re;
	rhat[n].im = rn[n].im;
      }

      //      hmc_complex check_scalarprod = scalar_product(rn,rhat);
      //      printf("check rn: (%e,%e)\n",check_scalarprod.re,check_scalarprod.im);
      
      printf("true residue squared: %e\n",global_squarenorm(rn));

      rho.re = hmc_one_f;
      alpha.re = hmc_one_f;
      omega.re = hmc_one_f;
      rho.im = 0;
      alpha.im = 0;
      omega.im = 0;

      set_zero_spinor(v);
      set_zero_spinor(p);
    }

    //    printf("rhat-norm: %e\n",global_squarenorm(rhat));
    //    printf("rnnorm: %e\n",global_squarenorm(rn));

    rho_next = scalar_product(rhat,rn);
    //    printf("rho: (%f,%f)\n",rho.re,rho.im);
    //    printf("rho_next: (%f,%f)\n",rho_next.re,rho_next.im);
    tmp1 = complexdivide(&rho_next,&rho);
    rho = rho_next;
    //    printf("tmp1: (%e,%e)\n",tmp1.re,tmp1.im);
    tmp2 = complexdivide(&alpha,&omega);
    beta = complexmult(&tmp1,&tmp2);

    //    printf("tmp1 -- tmp2: (%e,%e) -- (%e,%e)\n",tmp1.re,tmp1.im,tmp2.re,tmp2.im);
    //    printf("beta-rhonext/rho: (%e,%e)\n",beta.re-tmp1.re,beta.im-tmp1.im);

    //    printf("beta = (%e,%e)\n",beta.re,beta.im);

    for (int n=0; n<SPINORFIELDSIZE; n++) {
      tmp1 = complexmult(&omega,&v[n]);
      tmp2.re = p[n].re - tmp1.re;
      tmp2.im = p[n].im - tmp1.im;
      tmp3 = complexmult(&beta,&tmp2);
      p[n].re = rn[n].re + tmp3.re;
      p[n].im = rn[n].im + tmp3.im;
    }

    //    hmc_float sum=0;
    //    for(int n=0; n<SPINORFIELDSIZE; n++) {
    //      hmc_float tmpre = p[n].re-rn[n].re;
    //      hmc_float tmpim = p[n].im-rn[n].im;
    //      sum += tmpre*tmpre+tmpim*tmpim;
    //    }
    //    printf("check sum: %e\n",sum);


    M(p,v,gaugefield,kappa,mu);

    tmp1 = scalar_product(rhat,v);
    alpha = complexdivide(&rho,&tmp1);

    //   printf("alpha = (%e,%e)\n",alpha.re,alpha.im);

    for (int n=0; n<SPINORFIELDSIZE; n++) {
      tmp1 = complexmult(&alpha,&v[n]);
      s[n].re = rn[n].re - tmp1.re;
      s[n].im = rn[n].im - tmp1.im;
    }

    M(s,t,gaugefield,kappa,mu);

    tmp1 = scalar_product(t,s);
    tmp2 = scalar_product(t,t);

    //    hmc_float mycheck = global_squarenorm(t);
    //    printf("check norm: (%e,%e)\n",tmp2.re/mycheck,tmp2.im/mycheck);

    omega = complexdivide(&tmp1,&tmp2);

    for (int n=0; n<SPINORFIELDSIZE; n++) {
      tmp1 = complexmult(&omega,&s[n]);
      tmp2 = complexmult(&alpha,&p[n]);
      tmp3 = complexmult(&omega,&t[n]);
      inout[n].re += tmp1.re + tmp2.re;
      inout[n].im += tmp1.im + tmp2.im;
      rn[n].re = s[n].re - tmp3.re;
      rn[n].im = s[n].im - tmp3.im;
    }
    //    hmc_complex check2 = scalar_product(rn,rhat);
    //    printf("check rn: (%e,%e)\n",check2.re,check2.im);

    hmc_float resid = global_squarenorm(rn);
    //    printf("%d\t%e\n",iter,resid);

    if(resid<epssquare) {
      hmc_spinor_field* aux = new hmc_spinor_field[SPINORFIELDSIZE];
      M(inout,aux,gaugefield,kappa,mu);
      for (int n=0; n<SPINORFIELDSIZE; n++) {
	aux[n].re = source[n].re - aux[n].re;
	aux[n].im = source[n].im - aux[n].im;
      }
      hmc_float trueresid = global_squarenorm(aux);
      //      printf("true residue squared: %e\n",trueresid);
      printf("%d\t%e\t%e\n",iter,resid,trueresid);

      if(trueresid<epssquare) return HMC_SUCCESS;

      delete [] aux;
    }

      //    if((iter+1)%10==0) return 0;

  }

  delete [] rn;
  delete [] rhat;
  delete [] v;
  delete [] p;
  delete [] s;
  delete [] t;
  return HMC_SUCCESS;
}



hmc_error create_point_source(hmc_spinor_field* b, int i, int spacepos, int timepos, hmc_float kappa, hmc_float mu, hmc_gaugefield* gaugefield){
  set_zero_spinor(b);

  int color = spinor_color(i);
  int spin = spinor_spin(i,color);

  b[spinor_field_element(spin,color,spacepos,timepos)].re = sqrt(2.*kappa);

  return HMC_SUCCESS;
}


/*
hmc_error M(hmc_spinor_field* in, hmc_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu){
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    out[n].re = n*.001*in[n].re + in[n].im*2412;
    out[n].im = in[n].im;
  }

  out[0].im = out[3].re - 7*in[3].re;


  return HMC_SUCCESS;
}
*/


hmc_error M(hmc_spinor_field* in, hmc_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu){
  
  for(int spacepos=0; spacepos<VOLSPACE; spacepos++) {
    for(int timepos=0; timepos<NTIME; timepos++) {
      hmc_spinor spinout[SPINORSIZE];
      hmc_spinor spintmp[SPINORSIZE];
      get_spinor_from_field(in,spinout,spacepos,timepos);
      hmc_float twistfactor = 2*kappa*mu;
      multiply_spinor_factor_gamma5(spinout,spintmp,twistfactor);
      spinors_accumulate(spinout,spintmp);
      put_spinor_to_field(spinout,out,spacepos,timepos);
    }
  }

    
  hmc_spinor_field tmp[SPINORFIELDSIZE];
  dslash(in,tmp,gaugefield);


  for(int n=0; n<SPINORFIELDSIZE; n++) {
    out[n].re -= kappa*tmp[n].re;
    out[n].im -= kappa*tmp[n].im;
  }

  return HMC_SUCCESS;
}



hmc_error dslash(hmc_spinor_field* in, hmc_spinor_field* out, hmc_gaugefield* gaugefield){

  hmc_spinor spinout[SPINORSIZE];

  for(int spacepos=0; spacepos<VOLSPACE; spacepos++) {
    for(int timepos=0; timepos<NTIME; timepos++) {

      hmc_spinor spinnext[SPINORSIZE];
      hmc_spinor spinprev[SPINORSIZE];
      hmc_spinor tmp[SPINORSIZE];
      hmc_su3matrix u;
      hmc_su3matrix udagger;
      int next;
      int prev;

      set_local_zero_spinor(spinout);    

      // spinout = U_0*(r-gamma_0)*spinnext + U^dagger_0(x-hat0) * (r+gamma_0)*spinprev
      next = (timepos+1)%NTIME;
      prev = (timepos-1+NTIME)%NTIME;
      //implement APBC!!!
      get_spinor_from_field(in, spinnext, spacepos, next);
      get_spinor_from_field(in, spinprev, spacepos, prev);

      get_su3matrix(&u,gaugefield,spacepos,timepos,0);
      get_su3matrix(&udagger,gaugefield,spacepos,prev,0);
      adjoin_su3matrix(&udagger);

      multiply_spinor_gamma0(spinnext,tmp);
      real_multiply_spinor(tmp,-hmc_one_f);
      spinors_accumulate(spinnext,tmp);
      su3matrix_times_spinor(&u,spinnext,tmp);
      spinors_accumulate(spinout,tmp);

      multiply_spinor_gamma0(spinprev,tmp);
      spinors_accumulate(spinprev,tmp);
      su3matrix_times_spinor(&udagger,spinprev,tmp);
      spinors_accumulate(spinout,tmp);

    // spinout += U_1*(r-gamma_1)*spinnext + U^dagger_1(x-hat1) * (r+gamma_1)*spinprev
      next = get_neighbor(spacepos,1);
      prev = get_lower_neighbor(spacepos,1);
      get_spinor_from_field(in, spinnext, next, timepos);
      get_spinor_from_field(in, spinprev, prev, timepos);
      get_su3matrix(&u,gaugefield,spacepos,timepos,1);
      get_su3matrix(&udagger,gaugefield,prev,timepos,1);
      adjoin_su3matrix(&udagger);

      multiply_spinor_gamma1(spinnext,tmp);
      real_multiply_spinor(tmp,-hmc_one_f);
      spinors_accumulate(spinnext,tmp);
      su3matrix_times_spinor(&u,spinnext,tmp);
      spinors_accumulate(spinout,tmp);

      multiply_spinor_gamma1(spinprev,tmp);
      spinors_accumulate(spinprev,tmp);
      su3matrix_times_spinor(&u,spinprev,tmp);
      spinors_accumulate(spinout,tmp);


    // spinout += U_2*(r-gamma_2)*spinnext + U^dagger_2(x-hat2) * (r+gamma_2)*spinprev
      next = get_neighbor(spacepos,2);
      prev = get_lower_neighbor(spacepos,2);
      get_spinor_from_field(in, spinnext, next, timepos);
      get_spinor_from_field(in, spinprev, prev, timepos);
      get_su3matrix(&u,gaugefield,spacepos,timepos,2);
      get_su3matrix(&udagger,gaugefield,prev,timepos,2);
      adjoin_su3matrix(&udagger);

      multiply_spinor_gamma2(spinnext,tmp);
      real_multiply_spinor(tmp,-hmc_one_f);
      spinors_accumulate(spinnext,tmp);
      su3matrix_times_spinor(&u,spinnext,tmp);
      spinors_accumulate(spinout,tmp);

      multiply_spinor_gamma2(spinprev,tmp);
      spinors_accumulate(spinprev,tmp);
      su3matrix_times_spinor(&u,spinprev,tmp);
      spinors_accumulate(spinout,tmp);    


    // spinout += U_3*(r-gamma_3)*spinnext + U^dagger_3(x-hat3) * (r+gamma_3)*spinprev
      next = get_neighbor(spacepos,3);
      prev = get_lower_neighbor(spacepos,3);
      get_spinor_from_field(in, spinnext, next, timepos);
      get_spinor_from_field(in, spinprev, prev, timepos);
      get_su3matrix(&u,gaugefield,spacepos,timepos,3);
      get_su3matrix(&udagger,gaugefield,prev,timepos,3);
      adjoin_su3matrix(&udagger);

      multiply_spinor_gamma3(spinnext,tmp);
      real_multiply_spinor(tmp,-hmc_one_f);
      spinors_accumulate(spinnext,tmp);
      su3matrix_times_spinor(&u,spinnext,tmp);
      spinors_accumulate(spinout,tmp);

      multiply_spinor_gamma3(spinprev,tmp);
      spinors_accumulate(spinprev,tmp);
      su3matrix_times_spinor(&u,spinprev,tmp);
      spinors_accumulate(spinout,tmp);


      put_spinor_to_field(spinout,out,spacepos,timepos);

    }
  }


  return HMC_SUCCESS;
}



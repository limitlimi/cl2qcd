/** @file
 * Device code implementing 3x3 matrices
 */

//opencl_operations_matrix.cl

int inline ocl_3x3matrix_element(int a, int b)
{
	return a + 3*b;
}

void multiply_3x3matrix (__private hmc_ocl_3x3matrix *out, __private hmc_ocl_3x3matrix *p,__private hmc_ocl_3x3matrix *q)
{
	out[0].re = p[0].re*q[0].re + p[3].re*q[1].re + p[6].re*q[2].re -
	            (p[0].im*q[0].im + p[3].im*q[1].im + p[6].im*q[2].im );
	out[3].re = p[0].re*q[3].re + p[3].re*q[4].re + p[6].re*q[5].re-
	            (p[0].im*q[3].im + p[3].im*q[4].im + p[6].im*q[5].im);
	out[6].re = p[0].re*q[6].re + p[3].re*q[7].re + p[6].re*q[8].re -
	            (p[0].im*q[6].im + p[3].im*q[7].im + p[6].im*q[8].im);
	out[1].re = p[1].re*q[0].re + p[4].re*q[1].re + p[7].re*q[2].re -
	            (p[1].im*q[0].im + p[4].im*q[1].im + p[7].im*q[2].im);
	out[4].re = p[1].re*q[3].re + p[4].re*q[4].re + p[7].re*q[5].re -
	            (p[1].im*q[3].im + p[4].im*q[4].im + p[7].im*q[5].im);
	out[7].re = p[1].re*q[6].re + p[4].re*q[7].re + p[7].re*q[8].re -
	            (p[1].im*q[6].im + p[4].im*q[7].im + p[7].im*q[8].im);
	out[2].re = p[2].re*q[0].re + p[5].re*q[1].re + p[8].re*q[2].re -
	            (p[2].im*q[0].im + p[5].im*q[1].im + p[8].im*q[2].im);
	out[5].re = p[2].re*q[3].re + p[5].re*q[4].re + p[8].re*q[5].re -
	            (p[2].im*q[3].im + p[5].im*q[4].im + p[8].im*q[5].im);
	out[8].re = p[2].re*q[6].re + p[5].re*q[7].re + p[8].re*q[8].re -
	            (p[2].im*q[6].im + p[5].im*q[7].im + p[8].im*q[8].im);

	out[0].im = p[0].re*q[0].im + p[3].re*q[1].im + p[6].re*q[2].im +
	            p[0].im*q[0].re + p[3].im*q[1].re + p[6].im*q[2].re;
	out[3].im = p[0].re*q[3].im + p[3].re*q[4].im + p[6].re*q[5].im +
	            p[0].im*q[3].re + p[3].im*q[4].re + p[6].im*q[5].re;
	out[6].im = p[0].re*q[6].im + p[3].re*q[7].im + p[6].re*q[8].im +
	            p[0].im*q[6].re + p[3].im*q[7].re + p[6].im*q[8].re;
	out[1].im = p[1].re*q[0].im + p[4].re*q[1].im + p[7].re*q[2].im +
	            p[1].im*q[0].re + p[4].im*q[1].re + p[7].im*q[2].re;
	out[4].im = p[1].re*q[3].im + p[4].re*q[4].im + p[7].re*q[5].im +
	            p[1].im*q[3].re + p[4].im*q[4].re + p[7].im*q[5].re ;
	out[7].im = p[1].re*q[6].im + p[4].re*q[7].im + p[7].re*q[8].im +
	            p[1].im*q[6].re + p[4].im*q[7].re + p[7].im*q[8].re;
	out[2].im = p[2].re*q[0].im + p[5].re*q[1].im + p[8].re*q[2].im +
	            p[2].im*q[0].re + p[5].im*q[1].re + p[8].im*q[2].re;
	out[5].im = p[2].re*q[3].im + p[5].re*q[4].im + p[8].re*q[5].im +
	            p[2].im*q[3].re + p[5].im*q[4].re + p[8].im*q[5].re;
	out[8].im = p[2].re*q[6].im + p[5].re*q[7].im + p[8].re*q[8].im +
	            p[2].im*q[6].re + p[5].im*q[7].re + p[8].im*q[8].re;
}


void add_3x3matrix (__private hmc_ocl_3x3matrix *out, __private hmc_ocl_3x3matrix *p, __private hmc_ocl_3x3matrix *q)
{
	for (int i=0; i<9; i++){
	  out[i].re = p[i].re + q[i].re; 
	  out[i].im = p[i].im + q[i].im;
	}
}


void subtract_3x3matrix (__private hmc_ocl_3x3matrix *out, __private hmc_ocl_3x3matrix *p,__private hmc_ocl_3x3matrix *q)
{
	for (int i=0; i<9; i++){
	  out[i].re = p[i].re - q[i].re; 
	  out[i].im = p[i].im - q[i].im;
	}
}
  


void copy_3x3_matrix(__private hmc_ocl_3x3matrix *out, __private hmc_ocl_3x3matrix *in)
{
	for (int i=0; i<9; i++){
	  out[i] = in[i]; 
	}
}



void accumulate_su3matrix_3x3_add(__private hmc_ocl_3x3matrix *out, __private hmc_ocl_su3matrix *q){
#ifdef _RECONSTRUCT_TWELVE_
	for(int n=0; n<NC*(NC-1); n++) {
		out[n] = complexadd(out[n] , q[n]);
	}
	for(int n=NC*(NC-1);  n<NC*NC; n++) {
		hmc_complex tmp = reconstruct_su3(q, n-NC*(NC-1));
		out[n] = complexadd(out[n] , tmp);
	}
#else
	for(int k=0; k<NC*NC; k++) {
		out[k] = complexadd(out[k] , q[k]);
	}
#endif
  
} 



void trace_3x3matrix (__private hmc_complex * out, __private hmc_ocl_3x3matrix *q){
  
  out[0].re = q[0].re;
  out[0].im = q[0].im;
  for (int i=1; i< 3; i++){
    out[i].re += q[i].re;
    out[i].im += q[i].im;
  }
}



void adjoint_3x3matrix (__private hmc_ocl_3x3matrix * out, __private hmc_ocl_3x3matrix *q){
	for(int a=0; a<3; a++) {
		for(int b=0; b<3; b++) {
			out[ocl_3x3matrix_element(a,b)] = complexconj(q[ocl_3x3matrix_element(b,a)]);
		}
	}
}


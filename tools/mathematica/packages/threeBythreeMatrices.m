(* ::Package:: *)

BeginPackage["ThreeBythreeMatrices`"]

cold3x3mat::usage = 
	"coldSU3mat gives the 3-by-3 Identity matrix"

nonTrivial3x3mat::usage =
	"nonTrivialSU3mat gives a non-trivially filled 3-by-3 matrix"

ascending3x3mat::usage =
	"ascending3x3mat gives a 3-by-3 matrix filled with ascending complex numbers"

stapleSum::usage = 
	"stapleSum[s] takes an SU(3) matrix as argument and compute the sum of staples."

rectStapleSum::usage = 
	"rectStapleSum[s] takes an SU(3) matrix as argument and compute the sum of rectangles staples."

tracelessAntihermitianPart::usage =
	"tracelessAntihermitianPart takes a 3-by-3 matrix as argument and computes its traceless
	 anti-hermitian part which is another 3-by-3 matrix."

mat3x3FromKroneckerProductOf3ComponentsVectors::usage=
	"mat3x3FromKroneckerProductOf3ComponentsVectors calculates the Dirac-Trace of the matrix resulting
	from multiplying U*V^dagger =  u*v^dagger + w*x^dagger, where u, v, w, x are SU(3)-vectors
	(using spinprojection) see tr_v_times_u_dagger in operations_su3vec.cl"

mat3x3FromKroneckerProduct::usage=
	"mat3x3FromKroneckerProduct calculates the Dirac-Trace of the matrix resulting
	from multiplying X*Y^dagger =  x1*y1^dagger + x2*y2^dagger + x3*y3^dagger + x4*y4^dagger , where x1,x2,x3,x4,y1,y2,y3,y4 are SU(3)-vectors
	(using spinprojection) see tr_dirac_x_times_y_dagger in operations_su3vec.cl"

Begin["Private`"]

cold3x3mat:=
	Module[ {coldSU3mat=IdentityMatrix[3]},
	coldSU3mat
	]

nonTrivial3x3mat:=
	Module[ {nonTrivialSU3mat={{0.130189 + I*0.260378, 0.260378 + I*0.390567, 
  0.520756 + I*0.650945}, {.572742 + I*.403041, .371222 + 
   I*.321726, -.449002 + I*(-.258088)}, {0.0000000000000000111022 + 
   I*.651751, .0271563 + I*(-.733219), -.0271563 + I*.190094}}},
	nonTrivialSU3mat
	]

ascending3x3mat:=
	Module[ {ascending3x3mat={{1 + I*2, 3 + I*4, 5 + I*6}, {7 + I*8, 9 + I*10, 11 + I*12}, {13 + I*14, 15 + I*16, 17 + I*18}}},
	ascending3x3mat
	]

stapleSum[u_]:=
	Module[{stapleSum=3*(u.ConjugateTranspose[u].ConjugateTranspose[u] + 
   ConjugateTranspose[u].ConjugateTranspose[u].u)},
	stapleSum
	]

rectStapleSum[u_]:=
	Module[{rectStapleSum=3*(u.u.ConjugateTranspose[u].ConjugateTranspose[u].ConjugateTranspose[u] +
							 u.ConjugateTranspose[u].ConjugateTranspose[u].ConjugateTranspose[u].u +
							 ConjugateTranspose[u].ConjugateTranspose[u].ConjugateTranspose[u].u.u +
							 ConjugateTranspose[u].ConjugateTranspose[u].ConjugateTranspose[u].u.u +
							 u.ConjugateTranspose[u].ConjugateTranspose[u].ConjugateTranspose[u].u +
							 u.u.ConjugateTranspose[u].ConjugateTranspose[u].ConjugateTranspose[u])},
	rectStapleSum
	]

tracelessAntihermitianPart[u_]:=
	Module[{tracelessAntihermitianPart= (1/2)*(u - ConjugateTranspose[u]) -
										(1/6)*Tr[u - ConjugateTranspose[u]]*IdentityMatrix[3]},
	tracelessAntihermitianPart
	]

mat3x3FromKroneckerProductOf3ComponentsVectors[s1_,s2_,s3_,s4_]:=
	Module[{mat3x3FromKroneckerProductOf3ComponentsVectors=KroneckerProduct[s1, ConjugateTranspose[s2]]
															+ KroneckerProduct[s3, ConjugateTranspose[s4]]},
	mat3x3FromKroneckerProductOf3ComponentsVectors
	]

mat3x3FromKroneckerProduct[x_,y_]:=
	Module[{mat3x3FromKroneckerProduct=KroneckerProduct[{x[[1]], x[[2]], x[[3]]}, ConjugateTranspose[{y[[1]], y[[2]], y[[3]]}]]
															+ KroneckerProduct[{x[[4]], x[[5]], x[[6]]}, ConjugateTranspose[{y[[4]], y[[5]], y[[6]]}]]
															+ KroneckerProduct[{x[[7]], x[[8]], x[[9]]}, ConjugateTranspose[{y[[7]], y[[8]], y[[9]]}]]
															+ KroneckerProduct[{x[[10]], x[[11]], x[[12]]}, ConjugateTranspose[{y[[10]], y[[11]], y[[12]]}]]},
	mat3x3FromKroneckerProduct
	]

End[]

EndPackage[]







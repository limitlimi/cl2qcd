(* ::Package:: *)

SetOptions[SelectedNotebook[], PrintPrecision -> 20]


Get["threeBythreeMatrices.m",Path->{NotebookDirectory[]}]
Get["twelveComponentsVectors.m",Path->{NotebookDirectory[]}]
Get["threeComponentsVectors.m",Path->{NotebookDirectory[]}]
Get["eightComponentsVectors.m",Path->{NotebookDirectory[]}]
Get["sigmaMuNuMatricesTensorColorIdentity.m",Path->{NotebookDirectory[]}]
Get["real.m",Path->{NotebookDirectory[]}]
Get["vectors.m",Path->{NotebookDirectory[]}]


Square41[s1_,s2_]:=I*(mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma41[cold3x3mat].s1,s2]+mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma41[cold3x3mat].s2,s1]) 
Square14[s1_,s2_]:=I*(mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma14[cold3x3mat].s1,s2]+mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma14[cold3x3mat].s2,s1])
Square42[s1_,s2_]:=I*(mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma42[cold3x3mat].s1,s2]+mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma42[cold3x3mat].s2,s1])
Square24[s1_,s2_]:=I*(mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma24[cold3x3mat].s1,s2]+mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma24[cold3x3mat].s2,s1])
Square43[s1_,s2_]:=I*(mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma43[cold3x3mat].s1,s2]+mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma43[cold3x3mat].s2,s1])
Square34[s1_,s2_]:=I*(mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma34[cold3x3mat].s1,s2]+mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma34[cold3x3mat].s2,s1])
Square32[s1_,s2_]:=I*(mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma32[cold3x3mat].s1,s2]+mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma32[cold3x3mat].s2,s1])
Square23[s1_,s2_]:=I*(mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma23[cold3x3mat].s1,s2]+mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma23[cold3x3mat].s2,s1])
Square31[s1_,s2_]:=I*(mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma31[cold3x3mat].s1,s2]+mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma31[cold3x3mat].s2,s1])
Square13[s1_,s2_]:=I*(mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma13[cold3x3mat].s1,s2]+mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma13[cold3x3mat].s2,s1])
Square21[s1_,s2_]:=I*(mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma21[cold3x3mat].s1,s2]+mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma21[cold3x3mat].s2,s1])
Square12[s1_,s2_]:=I*(mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma12[cold3x3mat].s1,s2]+mat3x3FromKroneckerProduct[Gamma5[cold3x3mat].Sigma12[cold3x3mat].s2,s1])
Square41[spinorCold, spinorCold]
Square31[spinorCold, spinorCold]
Square41[spinorCold, spinorCold].IdentityMatrix[3]


(*let us assume a filling for the gaugefield with ascending3x3Matrix filling in TDIR=4 and unit_matrixsu3 (to be changed in ->nonTrivialSu3Matrix) filling in all other DIRs*)
diagram1aup41[s1_,s2_]:=nonTrivial3x3mat.ConjugateTranspose[ascending3x3mat].Square41[s1, s2].ConjugateTranspose[nonTrivial3x3mat]
diagram1aup14[s1_,s2_]:=ascending3x3mat.ConjugateTranspose[nonTrivial3x3mat].Square14[s1, s2].ConjugateTranspose[ascending3x3mat]
diagram1aup42[s1_,s2_]:=nonTrivial3x3mat.ConjugateTranspose[ascending3x3mat].Square42[s1, s2].ConjugateTranspose[nonTrivial3x3mat]
diagram1aup24[s1_,s2_]:=ascending3x3mat.ConjugateTranspose[nonTrivial3x3mat].Square24[s1, s2].ConjugateTranspose[ascending3x3mat]
diagram1aup43[s1_,s2_]:=nonTrivial3x3mat.ConjugateTranspose[ascending3x3mat].Square43[s1, s2].ConjugateTranspose[nonTrivial3x3mat]
diagram1aup34[s1_,s2_]:=ascending3x3mat.ConjugateTranspose[nonTrivial3x3mat].Square34[s1, s2].ConjugateTranspose[ascending3x3mat]
diagram1aup32[s1_,s2_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].Square32[s1, s2].ConjugateTranspose[nonTrivial3x3mat]
diagram1aup23[s1_,s2_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].Square23[s1, s2].ConjugateTranspose[nonTrivial3x3mat]
diagram1aup31[s1_,s2_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].Square31[s1, s2].ConjugateTranspose[nonTrivial3x3mat]
diagram1aup13[s1_,s2_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].Square13[s1, s2].ConjugateTranspose[nonTrivial3x3mat]
diagram1aup21[s1_,s2_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].Square21[s1, s2].ConjugateTranspose[nonTrivial3x3mat]
diagram1aup12[s1_,s2_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].Square12[s1, s2].ConjugateTranspose[nonTrivial3x3mat]
diagram1aup41[spinorCold,spinorAscendingComplex]


(*let us assume a filling for the gaugefield with ascending3x3Matrix filling in TDIR=4 and unit_matrixsu3 (to be changed in ->nonTrivialSu3Matrix) filling in all other DIRs*)
diagram1adown41[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].Square41[s1, s2].nonTrivial3x3mat
diagram1adown14[s1_,s2_]:=-ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].Square14[s1, s2].ascending3x3mat
diagram1adown42[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].Square42[s1, s2].nonTrivial3x3mat
diagram1adown24[s1_,s2_]:=-ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].Square24[s1, s2].ascending3x3mat
diagram1adown43[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].Square43[s1, s2].nonTrivial3x3mat
diagram1adown34[s1_,s2_]:=-ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].Square34[s1, s2].ascending3x3mat
diagram1adown32[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Square32[s1, s2].nonTrivial3x3mat
diagram1adown23[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Square23[s1, s2].nonTrivial3x3mat
diagram1adown31[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Square31[s1, s2].nonTrivial3x3mat
diagram1adown13[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Square13[s1, s2].nonTrivial3x3mat
diagram1adown21[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Square21[s1, s2].nonTrivial3x3mat
diagram1adown12[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Square12[s1, s2].nonTrivial3x3mat
diagram1adown41[spinorCold,spinorAscendingComplex]


(*let us assume a filling for the gaugefield with ascending3x3Matrix filling in TDIR=4 and unit_matrixsu3 (to be changed in ->nonTrivialSu3Matrix) filling in all other DIRs*)
diagram1bup41[s1_,s2_]:=nonTrivial3x3mat.ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].Square41[s1, s2]
diagram1bup14[s1_,s2_]:=ascending3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].Square14[s1, s2]
diagram1bup42[s1_,s2_]:=nonTrivial3x3mat.ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].Square42[s1, s2]
diagram1bup24[s1_,s2_]:=ascending3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].Square24[s1, s2]
diagram1bup43[s1_,s2_]:=nonTrivial3x3mat.ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].Square43[s1, s2]
diagram1bup34[s1_,s2_]:=ascending3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].Square34[s1, s2]
diagram1bup32[s1_,s2_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Square32[s1, s2]
diagram1bup23[s1_,s2_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Square23[s1, s2]
diagram1bup31[s1_,s2_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Square31[s1, s2]
diagram1bup13[s1_,s2_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Square13[s1, s2]
diagram1bup21[s1_,s2_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Square21[s1, s2]
diagram1bup12[s1_,s2_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Square12[s1, s2]
diagram1bup41[spinorCold,spinorAscendingComplex]


(*let us assume a filling for the gaugefield with ascending3x3Matrix filling in TDIR=4 and unit_matrixsu3 (to be changed in ->nonTrivialSu3Matrix) filling in all other DIRs*)
diagram1bdown41[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].nonTrivial3x3mat.Square41[s1, s2]
diagram1bdown14[s1_,s2_]:=-ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].ascending3x3mat.Square14[s1, s2]
diagram1bdown42[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].nonTrivial3x3mat.Square42[s1, s2]
diagram1bdown24[s1_,s2_]:=-ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].ascending3x3mat.Square24[s1, s2]
diagram1bdown43[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].nonTrivial3x3mat.Square43[s1, s2]
diagram1bdown34[s1_,s2_]:=-ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].ascending3x3mat.Square34[s1, s2]
diagram1bdown32[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat.Square32[s1, s2]
diagram1bdown23[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat.Square23[s1, s2]
diagram1bdown31[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat.Square31[s1, s2]
diagram1bdown13[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat.Square13[s1, s2]
diagram1bdown21[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat.Square21[s1, s2]
diagram1bdown12[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat.Square12[s1, s2]
diagram1bdown41[spinorCold,spinorAscendingComplex]


(*let us assume a filling for the gaugefield with ascending3x3Matrix filling in TDIR=4 and unit_matrixsu3 (to be changed in ->nonTrivialSu3Matrix) filling in all other DIRs*)
diagram1cup41[s1_,s2_]:=Square41[s1, s2].nonTrivial3x3mat.ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram1cup14[s1_,s2_]:=Square14[s1, s2].ascending3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat]
diagram1cup42[s1_,s2_]:=Square42[s1, s2].nonTrivial3x3mat.ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram1cup24[s1_,s2_]:=Square24[s1, s2].ascending3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat]
diagram1cup43[s1_,s2_]:=Square43[s1, s2].nonTrivial3x3mat.ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram1cup34[s1_,s2_]:=Square34[s1, s2].ascending3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat]
diagram1cup32[s1_,s2_]:=Square32[s1, s2].nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram1cup23[s1_,s2_]:=Square23[s1, s2].nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram1cup31[s1_,s2_]:=Square31[s1, s2].nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram1cup13[s1_,s2_]:=Square13[s1, s2].nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram1cup21[s1_,s2_]:=Square21[s1, s2].nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram1cup12[s1_,s2_]:=Square12[s1, s2].nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram1cup41[spinorCold,spinorAscendingComplex]


(*let us assume a filling for the gaugefield with ascending3x3Matrix filling in TDIR=4 and unit_matrixsu3 (to be changed in ->nonTrivialSu3Matrix) filling in all other DIRs*)
diagram1cdown41[s1_,s2_]:=-Square41[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].nonTrivial3x3mat
diagram1cdown14[s1_,s2_]:=-Square14[s1, s2].ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].ascending3x3mat
diagram1cdown42[s1_,s2_]:=-Square42[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].nonTrivial3x3mat
diagram1cdown24[s1_,s2_]:=-Square24[s1, s2].ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].ascending3x3mat
diagram1cdown43[s1_,s2_]:=-Square43[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].nonTrivial3x3mat
diagram1cdown34[s1_,s2_]:=-Square34[s1, s2].ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].ascending3x3mat
diagram1cdown32[s1_,s2_]:=-Square32[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram1cdown23[s1_,s2_]:=-Square23[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram1cdown31[s1_,s2_]:=-Square31[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram1cdown13[s1_,s2_]:=-Square13[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram1cdown21[s1_,s2_]:=-Square21[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram1cdown12[s1_,s2_]:=-Square12[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram1cdown41[spinorCold,spinorAscendingComplex]


(*let us assume a filling for the gaugefield with ascending3x3Matrix filling in TDIR=4 and unit_matrixsu3 (to be changed in ->nonTrivialSu3Matrix) filling in all other DIRs*)
diagram1dup41[s1_,s2_]:=nonTrivial3x3mat.Square41[s1, s2].ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram1dup14[s1_,s2_]:=ascending3x3mat.Square14[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat]
diagram1dup42[s1_,s2_]:=nonTrivial3x3mat.Square42[s1, s2].ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram1dup24[s1_,s2_]:=ascending3x3mat.Square24[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat]
diagram1dup43[s1_,s2_]:=nonTrivial3x3mat.Square43[s1, s2].ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram1dup34[s1_,s2_]:=ascending3x3mat.Square34[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat]
diagram1dup32[s1_,s2_]:=nonTrivial3x3mat.Square32[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram1dup23[s1_,s2_]:=nonTrivial3x3mat.Square23[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram1dup31[s1_,s2_]:=nonTrivial3x3mat.Square31[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram1dup13[s1_,s2_]:=nonTrivial3x3mat.Square13[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram1dup21[s1_,s2_]:=nonTrivial3x3mat.Square21[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram1dup12[s1_,s2_]:=nonTrivial3x3mat.Square12[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram1dup41[spinorCold,spinorAscendingComplex]


(*let us assume a filling for the gaugefield with ascending3x3Matrix filling in TDIR=4 and unit_matrixsu3 (to be changed in ->nonTrivialSu3Matrix) filling in all other DIRs*)
diagram1ddown41[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].Square41[s1, s2].ConjugateTranspose[ascending3x3mat].nonTrivial3x3mat
diagram1ddown14[s1_,s2_]:=-ConjugateTranspose[ascending3x3mat].Square14[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ascending3x3mat
diagram1ddown42[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].Square42[s1, s2].ConjugateTranspose[ascending3x3mat].nonTrivial3x3mat
diagram1ddown24[s1_,s2_]:=-ConjugateTranspose[ascending3x3mat].Square24[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ascending3x3mat
diagram1ddown43[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].Square43[s1, s2].ConjugateTranspose[ascending3x3mat].nonTrivial3x3mat
diagram1ddown34[s1_,s2_]:=-ConjugateTranspose[ascending3x3mat].Square34[s1, s2].ConjugateTranspose[nonTrivial3x3mat].ascending3x3mat
diagram1ddown32[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].Square32[s1, s2].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram1ddown23[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].Square23[s1, s2].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram1ddown31[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].Square31[s1, s2].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram1ddown13[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].Square13[s1, s2].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram1ddown21[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].Square21[s1, s2].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram1ddown12[s1_,s2_]:=-ConjugateTranspose[nonTrivial3x3mat].Square12[s1, s2].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram1ddown41[spinorCold,spinorAscendingComplex]


summedUpDiagrams41[s1_,s2_]:=diagram1aup41[s1,s2]+diagram1adown41[s1,s2]+diagram1bup41[s1,s2]+diagram1bdown41[s1,s2]+diagram1cup41[s1,s2]+diagram1cdown41[s1,s2]+diagram1dup41[s1,s2]+diagram1ddown41[s1,s2]
summedUpDiagrams14[s1_,s2_]:=diagram1aup14[s1,s2]+diagram1adown14[s1,s2]+diagram1bup14[s1,s2]+diagram1bdown14[s1,s2]+diagram1cup14[s1,s2]+diagram1cdown14[s1,s2]+diagram1dup14[s1,s2]+diagram1ddown14[s1,s2]
summedUpDiagrams42[s1_,s2_]:=diagram1aup42[s1,s2]+diagram1adown42[s1,s2]+diagram1bup42[s1,s2]+diagram1bdown42[s1,s2]+diagram1cup42[s1,s2]+diagram1cdown42[s1,s2]+diagram1dup42[s1,s2]+diagram1ddown42[s1,s2]
summedUpDiagrams24[s1_,s2_]:=diagram1aup24[s1,s2]+diagram1adown24[s1,s2]+diagram1bup24[s1,s2]+diagram1bdown24[s1,s2]+diagram1cup24[s1,s2]+diagram1cdown24[s1,s2]+diagram1dup24[s1,s2]+diagram1ddown24[s1,s2]
summedUpDiagrams43[s1_,s2_]:=diagram1aup43[s1,s2]+diagram1adown43[s1,s2]+diagram1bup43[s1,s2]+diagram1bdown43[s1,s2]+diagram1cup43[s1,s2]+diagram1cdown43[s1,s2]+diagram1dup43[s1,s2]+diagram1ddown43[s1,s2]
summedUpDiagrams34[s1_,s2_]:=diagram1aup34[s1,s2]+diagram1adown34[s1,s2]+diagram1bup34[s1,s2]+diagram1bdown34[s1,s2]+diagram1cup34[s1,s2]+diagram1cdown34[s1,s2]+diagram1dup34[s1,s2]+diagram1ddown34[s1,s2]
summedUpDiagrams32[s1_,s2_]:=diagram1aup32[s1,s2]+diagram1adown32[s1,s2]+diagram1bup32[s1,s2]+diagram1bdown32[s1,s2]+diagram1cup32[s1,s2]+diagram1cdown32[s1,s2]+diagram1dup32[s1,s2]+diagram1ddown32[s1,s2]
summedUpDiagrams23[s1_,s2_]:=diagram1aup23[s1,s2]+diagram1adown23[s1,s2]+diagram1bup23[s1,s2]+diagram1bdown23[s1,s2]+diagram1cup23[s1,s2]+diagram1cdown23[s1,s2]+diagram1dup23[s1,s2]+diagram1ddown23[s1,s2]
summedUpDiagrams31[s1_,s2_]:=diagram1aup31[s1,s2]+diagram1adown31[s1,s2]+diagram1bup31[s1,s2]+diagram1bdown31[s1,s2]+diagram1cup31[s1,s2]+diagram1cdown31[s1,s2]+diagram1dup31[s1,s2]+diagram1ddown31[s1,s2]
summedUpDiagrams13[s1_,s2_]:=diagram1aup13[s1,s2]+diagram1adown13[s1,s2]+diagram1bup13[s1,s2]+diagram1bdown13[s1,s2]+diagram1cup13[s1,s2]+diagram1cdown13[s1,s2]+diagram1dup13[s1,s2]+diagram1ddown13[s1,s2]
summedUpDiagrams21[s1_,s2_]:=diagram1aup21[s1,s2]+diagram1adown21[s1,s2]+diagram1bup21[s1,s2]+diagram1bdown21[s1,s2]+diagram1cup21[s1,s2]+diagram1cdown21[s1,s2]+diagram1dup21[s1,s2]+diagram1ddown21[s1,s2]
summedUpDiagrams12[s1_,s2_]:=diagram1aup12[s1,s2]+diagram1adown12[s1,s2]+diagram1bup12[s1,s2]+diagram1bdown12[s1,s2]+diagram1cup12[s1,s2]+diagram1cdown12[s1,s2]+diagram1dup12[s1,s2]+diagram1ddown12[s1,s2]
summedUpDiagrams41[spinorCold,spinorAscendingComplex]
summedUpDiagrams42[spinorCold,spinorAscendingComplex]
summedUpDiagrams43[spinorCold,spinorAscendingComplex]


c0hat[k_]:=1./(1+64*k*k)
factor[k_,csw_]:=0.125 * k * c0hat[k] * csw
factor[nonTrivialRealPar,nonTrivialRealPar]
summedUpDiagrams4WithFactor[k_,csw_,s1_,s2_]:=factor[k,csw]*(summedUpDiagrams41[s1,s2]+summedUpDiagrams42[s1,s2]+summedUpDiagrams43[s1,s2])
summedUpDiagrams3WithFactor[k_,csw_,s1_,s2_]:=factor[k,csw]*(summedUpDiagrams31[s1,s2]+summedUpDiagrams32[s1,s2]+summedUpDiagrams34[s1,s2])
summedUpDiagrams2WithFactor[k_,csw_,s1_,s2_]:=factor[k,csw]*(summedUpDiagrams21[s1,s2]+summedUpDiagrams23[s1,s2]+summedUpDiagrams24[s1,s2])
summedUpDiagrams1WithFactor[k_,csw_,s1_,s2_]:=factor[k,csw]*(summedUpDiagrams12[s1,s2]+summedUpDiagrams13[s1,s2]+summedUpDiagrams14[s1,s2])
summedUpDiagrams4WithFactor[nonTrivialRealPar,nonTrivialRealPar,spinorCold,spinorAscendingComplex]
summedUpDiagrams3WithFactor[nonTrivialRealPar,nonTrivialRealPar,spinorCold,spinorAscendingComplex]
summedUpDiagrams2WithFactor[nonTrivialRealPar,nonTrivialRealPar,spinorCold,spinorAscendingComplex]
summedUpDiagrams1WithFactor[nonTrivialRealPar,nonTrivialRealPar,spinorCold,spinorAscendingComplex]


fermionForceClover1Eo0[k_,csw_,s1_,s2_]:=algebraElement[summedUpDiagrams4WithFactor[k,csw,s1,s2]]
fermionForceClover1Eo1[k_,csw_,s1_,s2_]:=algebraElement[summedUpDiagrams1WithFactor[k,csw,s1,s2]]
fermionForceClover1Eo2[k_,csw_,s1_,s2_]:=algebraElement[summedUpDiagrams2WithFactor[k,csw,s1,s2]]
fermionForceClover1Eo3[k_,csw_,s1_,s2_]:=algebraElement[summedUpDiagrams3WithFactor[k,csw,s1,s2]]
fermionForceClover1Eo0[nonTrivialRealPar,nonTrivialRealPar,spinorCold,spinorAscendingComplex]
fermionForceClover1Eo1[nonTrivialRealPar,nonTrivialRealPar,spinorCold,spinorAscendingComplex]
fermionForceClover1Eo2[nonTrivialRealPar,nonTrivialRealPar,spinorCold,spinorAscendingComplex]
fermionForceClover1Eo3[nonTrivialRealPar,nonTrivialRealPar,spinorCold,spinorAscendingComplex]


(*Here we take into account that the calculation of the force with the corresponding gaugemomentum update concerns only ODD/EVEN sites. For all other sites (and directions) the gaugemomentum is unchanged.*)

FermionForceClover[gm_,k_,csw_,s1_,s2_]:= squareNorm[gm]*4 +squareNorm[gm + fermionForceClover1Eo0[k,csw,s1,s2]] +
											squareNorm[gm + fermionForceClover1Eo1[k,csw,s1,s2]] +
											squareNorm[gm + fermionForceClover1Eo2[k,csw,s1,s2]] +
											squareNorm[gm + fermionForceClover1Eo3[k,csw,s1,s2]]
FermionForceClover[gaugeMomOne,nonTrivialRealPar,nonTrivialRealPar,spinorAscendingComplex,spinorAscendingComplex]*4^4/2
FermionForceClover[gaugeMomAscending,nonTrivialRealPar,nonTrivialRealPar,spinorAscendingComplex,spinorAscendingComplex]


(*triangle insertion*)


N[(Pi/3*90/100)^2]
N[9/10*Pi/4]

N[(9/40)*Pi]
N[(9/10)*Pi-(2/3)*Pi]




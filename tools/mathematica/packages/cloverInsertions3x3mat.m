(* ::Package:: *)

SetOptions[SelectedNotebook[], PrintPrecision -> 20]


Get["threeBythreeMatrices.m",Path->{NotebookDirectory[]}]
Get["twelveComponentsVectors.m",Path->{NotebookDirectory[]}]
Get["threeComponentsVectors.m",Path->{NotebookDirectory[]}]
Get["eightComponentsVectors.m",Path->{NotebookDirectory[]}]
Get["sigmaMuNuMatricesTensorColorIdentity.m",Path->{NotebookDirectory[]}]
Get["sixBysixMatrices.m", Path->{NotebookDirectory[]}]
Get["twelveBytwelveMatrices.m", Path->{NotebookDirectory[]}]
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


c0hat[k_]:=1./(1+64*k*k)
factor[k_,csw_]:=0.125 * k * c0hat[k] * csw
summedUpDiagrams4WithFactor[k_,csw_,s1_,s2_]:=factor[k,csw]*(summedUpDiagrams41[s1,s2]+summedUpDiagrams42[s1,s2]+summedUpDiagrams43[s1,s2])
summedUpDiagrams3WithFactor[k_,csw_,s1_,s2_]:=factor[k,csw]*(summedUpDiagrams31[s1,s2]+summedUpDiagrams32[s1,s2]+summedUpDiagrams34[s1,s2])
summedUpDiagrams2WithFactor[k_,csw_,s1_,s2_]:=factor[k,csw]*(summedUpDiagrams21[s1,s2]+summedUpDiagrams23[s1,s2]+summedUpDiagrams24[s1,s2])
summedUpDiagrams1WithFactor[k_,csw_,s1_,s2_]:=factor[k,csw]*(summedUpDiagrams12[s1,s2]+summedUpDiagrams13[s1,s2]+summedUpDiagrams14[s1,s2])


fermionForceClover1Eo0[k_,csw_,s1_,s2_]:=algebraElement[summedUpDiagrams4WithFactor[k,csw,s1,s2]]
fermionForceClover1Eo1[k_,csw_,s1_,s2_]:=algebraElement[summedUpDiagrams1WithFactor[k,csw,s1,s2]]
fermionForceClover1Eo2[k_,csw_,s1_,s2_]:=algebraElement[summedUpDiagrams2WithFactor[k,csw,s1,s2]]
fermionForceClover1Eo3[k_,csw_,s1_,s2_]:=algebraElement[summedUpDiagrams3WithFactor[k,csw,s1,s2]]


(*Here we take into account that the calculation of the force with the corresponding gaugemomentum update concerns only ODD/EVEN sites. For all other sites (and directions) the gaugemomentum is unchanged.*)
FermionForceClover[gm_,k_,csw_,s1_,s2_]:= squareNorm[gm]*4 +squareNorm[gm + fermionForceClover1Eo0[k,csw,s1,s2]] +
											squareNorm[gm + fermionForceClover1Eo1[k,csw,s1,s2]] +
											squareNorm[gm + fermionForceClover1Eo2[k,csw,s1,s2]] +
											squareNorm[gm + fermionForceClover1Eo3[k,csw,s1,s2]]
FermionForceClover[gaugeMomOne,nonTrivialRealPar,nonTrivialRealPar,spinorAscendingComplex,spinorAscendingComplex]*4^4/2
FermionForceClover[gaugeMomAscending,nonTrivialRealPar,nonTrivialRealPar,spinorAscendingComplex,spinorAscendingComplex]


(*triangle insertion*)
Triangle41[m_]:=I*diracTraceOf12x12Matrix[Sigma41[cold3x3mat].m]
Triangle14[m_]:=I*diracTraceOf12x12Matrix[Sigma14[cold3x3mat].m]
Triangle31[m_]:=I*diracTraceOf12x12Matrix[Sigma31[cold3x3mat].m]
Triangle13[m_]:=I*diracTraceOf12x12Matrix[Sigma13[cold3x3mat].m]
Triangle21[m_]:=I*diracTraceOf12x12Matrix[Sigma21[cold3x3mat].m]
Triangle12[m_]:=I*diracTraceOf12x12Matrix[Sigma12[cold3x3mat].m]
Triangle23[m_]:=I*diracTraceOf12x12Matrix[Sigma23[cold3x3mat].m]
Triangle32[m_]:=I*diracTraceOf12x12Matrix[Sigma32[cold3x3mat].m]
Triangle24[m_]:=I*diracTraceOf12x12Matrix[Sigma24[cold3x3mat].m]
Triangle42[m_]:=I*diracTraceOf12x12Matrix[Sigma42[cold3x3mat].m]
Triangle34[m_]:=I*diracTraceOf12x12Matrix[Sigma34[cold3x3mat].m]
Triangle43[m_]:=I*diracTraceOf12x12Matrix[Sigma43[cold3x3mat].m]


(*let us assume a filling for the gaugefield with ascending3x3Matrix filling in TDIR=4 and unit_matrixsu3 (to be changed in ->nonTrivialSu3Matrix) filling in all other DIRs*)
diagram2aup41[m_]:=nonTrivial3x3mat.ConjugateTranspose[ascending3x3mat].Triangle41[m].ConjugateTranspose[nonTrivial3x3mat]
diagram2aup14[m_]:=ascending3x3mat.ConjugateTranspose[nonTrivial3x3mat].Triangle14[m].ConjugateTranspose[ascending3x3mat]
diagram2aup42[m_]:=nonTrivial3x3mat.ConjugateTranspose[ascending3x3mat].Triangle42[m].ConjugateTranspose[nonTrivial3x3mat]
diagram2aup24[m_]:=ascending3x3mat.ConjugateTranspose[nonTrivial3x3mat].Triangle24[m].ConjugateTranspose[ascending3x3mat]
diagram2aup43[m_]:=nonTrivial3x3mat.ConjugateTranspose[ascending3x3mat].Triangle43[m].ConjugateTranspose[nonTrivial3x3mat]
diagram2aup34[m_]:=ascending3x3mat.ConjugateTranspose[nonTrivial3x3mat].Triangle34[m].ConjugateTranspose[ascending3x3mat]
diagram2aup32[m_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].Triangle32[m].ConjugateTranspose[nonTrivial3x3mat]
diagram2aup23[m_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].Triangle23[m].ConjugateTranspose[nonTrivial3x3mat]
diagram2aup31[m_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].Triangle31[m].ConjugateTranspose[nonTrivial3x3mat]
diagram2aup13[m_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].Triangle13[m].ConjugateTranspose[nonTrivial3x3mat]
diagram2aup21[m_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].Triangle21[m].ConjugateTranspose[nonTrivial3x3mat]
diagram2aup12[m_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].Triangle12[m].ConjugateTranspose[nonTrivial3x3mat]


(*let us assume a filling for the gaugefield with ascending3x3Matrix filling in TDIR=4 and unit_matrixsu3 (to be changed in ->nonTrivialSu3Matrix) filling in all other DIRs*)
diagram2adown41[m_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].Triangle41[m].nonTrivial3x3mat
diagram2adown14[m_]:=-ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].Triangle14[m].ascending3x3mat
diagram2adown42[m_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].Triangle42[m].nonTrivial3x3mat
diagram2adown24[m_]:=-ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].Triangle24[m].ascending3x3mat
diagram2adown43[m_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].Triangle43[m].nonTrivial3x3mat
diagram2adown34[m_]:=-ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].Triangle34[m].ascending3x3mat
diagram2adown32[m_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Triangle32[m].nonTrivial3x3mat
diagram2adown23[m_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Triangle23[m].nonTrivial3x3mat
diagram2adown31[m_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Triangle31[m].nonTrivial3x3mat
diagram2adown13[m_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Triangle13[m].nonTrivial3x3mat
diagram2adown21[m_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Triangle21[m].nonTrivial3x3mat
diagram2adown12[m_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Triangle12[m].nonTrivial3x3mat


(*let us assume a filling for the gaugefield with ascending3x3Matrix filling in TDIR=4 and unit_matrixsu3 (to be changed in ->nonTrivialSu3Matrix) filling in all other DIRs*)
diagram2bup41[m_]:=nonTrivial3x3mat.ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].Triangle41[m]
diagram2bup14[m_]:=ascending3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].Triangle14[m]
diagram2bup42[m_]:=nonTrivial3x3mat.ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].Triangle42[m]
diagram2bup24[m_]:=ascending3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].Triangle24[m]
diagram2bup43[m_]:=nonTrivial3x3mat.ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].Triangle43[m]
diagram2bup34[m_]:=ascending3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].Triangle34[m]
diagram2bup32[m_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Triangle32[m]
diagram2bup23[m_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Triangle23[m]
diagram2bup31[m_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Triangle31[m]
diagram2bup13[m_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Triangle13[m]
diagram2bup21[m_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Triangle21[m]
diagram2bup12[m_]:=nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].Triangle12[m]


(*let us assume a filling for the gaugefield with ascending3x3Matrix filling in TDIR=4 and unit_matrixsu3 (to be changed in ->nonTrivialSu3Matrix) filling in all other DIRs*)
diagram2bdown41[m_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].nonTrivial3x3mat.Triangle41[m]
diagram2bdown14[m_]:=-ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].ascending3x3mat.Triangle14[m]
diagram2bdown42[m_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].nonTrivial3x3mat.Triangle42[m]
diagram2bdown24[m_]:=-ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].ascending3x3mat.Triangle24[m]
diagram2bdown43[m_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].nonTrivial3x3mat.Triangle43[m]
diagram2bdown34[m_]:=-ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].ascending3x3mat.Triangle34[m]
diagram2bdown32[m_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat.Triangle32[m]
diagram2bdown23[m_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat.Triangle23[m]
diagram2bdown31[m_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat.Triangle31[m]
diagram2bdown13[m_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat.Triangle13[m]
diagram2bdown21[m_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat.Triangle21[m]
diagram2bdown12[m_]:=-ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat.Triangle12[m]


(*let us assume a filling for the gaugefield with ascending3x3Matrix filling in TDIR=4 and unit_matrixsu3 (to be changed in ->nonTrivialSu3Matrix) filling in all other DIRs*)
diagram2cup41[m_]:=Triangle41[m].nonTrivial3x3mat.ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram2cup14[m_]:=Triangle14[m].ascending3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat]
diagram2cup42[m_]:=Triangle42[m].nonTrivial3x3mat.ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram2cup24[m_]:=Triangle24[m].ascending3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat]
diagram2cup43[m_]:=Triangle43[m].nonTrivial3x3mat.ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram2cup34[m_]:=Triangle34[m].ascending3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat]
diagram2cup32[m_]:=Triangle32[m].nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram2cup23[m_]:=Triangle23[m].nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram2cup31[m_]:=Triangle31[m].nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram2cup13[m_]:=Triangle13[m].nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram2cup21[m_]:=Triangle21[m].nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram2cup12[m_]:=Triangle12[m].nonTrivial3x3mat.ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]


(*let us assume a filling for the gaugefield with ascending3x3Matrix filling in TDIR=4 and unit_matrixsu3 (to be changed in ->nonTrivialSu3Matrix) filling in all other DIRs*)
diagram2cdown41[m_]:=-Triangle41[m].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].nonTrivial3x3mat
diagram2cdown14[m_]:=-Triangle14[m].ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].ascending3x3mat
diagram2cdown42[m_]:=-Triangle42[m].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].nonTrivial3x3mat
diagram2cdown24[m_]:=-Triangle24[m].ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].ascending3x3mat
diagram2cdown43[m_]:=-Triangle43[m].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat].nonTrivial3x3mat
diagram2cdown34[m_]:=-Triangle34[m].ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat].ascending3x3mat
diagram2cdown32[m_]:=-Triangle32[m].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram2cdown23[m_]:=-Triangle23[m].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram2cdown31[m_]:=-Triangle31[m].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram2cdown13[m_]:=-Triangle13[m].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram2cdown21[m_]:=-Triangle21[m].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram2cdown12[m_]:=-Triangle12[m].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat


(*let us assume a filling for the gaugefield with ascending3x3Matrix filling in TDIR=4 and unit_matrixsu3 (to be changed in ->nonTrivialSu3Matrix) filling in all other DIRs*)
diagram2dup41[m_]:=nonTrivial3x3mat.Triangle41[m].ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram2dup14[m_]:=ascending3x3mat.Triangle14[m].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat]
diagram2dup42[m_]:=nonTrivial3x3mat.Triangle42[m].ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram2dup24[m_]:=ascending3x3mat.Triangle24[m].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat]
diagram2dup43[m_]:=nonTrivial3x3mat.Triangle43[m].ConjugateTranspose[ascending3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram2dup34[m_]:=ascending3x3mat.Triangle34[m].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[ascending3x3mat]
diagram2dup32[m_]:=nonTrivial3x3mat.Triangle32[m].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram2dup23[m_]:=nonTrivial3x3mat.Triangle23[m].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram2dup31[m_]:=nonTrivial3x3mat.Triangle31[m].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram2dup13[m_]:=nonTrivial3x3mat.Triangle13[m].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram2dup21[m_]:=nonTrivial3x3mat.Triangle21[m].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]
diagram2dup12[m_]:=nonTrivial3x3mat.Triangle12[m].ConjugateTranspose[nonTrivial3x3mat].ConjugateTranspose[nonTrivial3x3mat]


(*let us assume a filling for the gaugefield with ascending3x3Matrix filling in TDIR=4 and unit_matrixsu3 (to be changed in ->nonTrivialSu3Matrix) filling in all other DIRs*)
diagram2ddown41[m_]:=-ConjugateTranspose[nonTrivial3x3mat].Triangle41[m].ConjugateTranspose[ascending3x3mat].nonTrivial3x3mat
diagram2ddown14[m_]:=-ConjugateTranspose[ascending3x3mat].Triangle14[m].ConjugateTranspose[nonTrivial3x3mat].ascending3x3mat
diagram2ddown42[m_]:=-ConjugateTranspose[nonTrivial3x3mat].Triangle42[m].ConjugateTranspose[ascending3x3mat].nonTrivial3x3mat
diagram2ddown24[m_]:=-ConjugateTranspose[ascending3x3mat].Triangle24[m].ConjugateTranspose[nonTrivial3x3mat].ascending3x3mat
diagram2ddown43[m_]:=-ConjugateTranspose[nonTrivial3x3mat].Triangle43[m].ConjugateTranspose[ascending3x3mat].nonTrivial3x3mat
diagram2ddown34[m_]:=-ConjugateTranspose[ascending3x3mat].Triangle34[m].ConjugateTranspose[nonTrivial3x3mat].ascending3x3mat
diagram2ddown32[m_]:=-ConjugateTranspose[nonTrivial3x3mat].Triangle32[m].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram2ddown23[m_]:=-ConjugateTranspose[nonTrivial3x3mat].Triangle23[m].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram2ddown31[m_]:=-ConjugateTranspose[nonTrivial3x3mat].Triangle31[m].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram2ddown13[m_]:=-ConjugateTranspose[nonTrivial3x3mat].Triangle13[m].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram2ddown21[m_]:=-ConjugateTranspose[nonTrivial3x3mat].Triangle21[m].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat
diagram2ddown12[m_]:=-ConjugateTranspose[nonTrivial3x3mat].Triangle12[m].ConjugateTranspose[nonTrivial3x3mat].nonTrivial3x3mat


summedUpDiagrams41E[m_]:=diagram2bup41[m]+diagram2bdown41[m]+diagram2dup41[m]+diagram2ddown41[m]
summedUpDiagrams14E[m_]:=diagram2bup14[m]+diagram2bdown14[m]+diagram2dup14[m]+diagram2ddown14[m]
summedUpDiagrams42E[m_]:=diagram2bup42[m]+diagram2bdown42[m]+diagram2dup42[m]+diagram2ddown42[m]
summedUpDiagrams24E[m_]:=diagram2bup24[m]+diagram2bdown24[m]+diagram2dup24[m]+diagram2ddown24[m]
summedUpDiagrams43E[m_]:=diagram2bup43[m]+diagram2bdown43[m]+diagram2dup43[m]+diagram2ddown43[m]
summedUpDiagrams34E[m_]:=diagram2bup34[m]+diagram2bdown34[m]+diagram2dup34[m]+diagram2ddown34[m]
summedUpDiagrams32E[m_]:=diagram2bup32[m]+diagram2bdown32[m]+diagram2dup32[m]+diagram2ddown32[m]
summedUpDiagrams23E[m_]:=diagram2bup23[m]+diagram2bdown23[m]+diagram2dup23[m]+diagram2ddown23[m]
summedUpDiagrams31E[m_]:=diagram2bup31[m]+diagram2bdown31[m]+diagram2dup31[m]+diagram2ddown31[m]
summedUpDiagrams13E[m_]:=diagram2bup13[m]+diagram2bdown13[m]+diagram2dup13[m]+diagram2ddown13[m]
summedUpDiagrams21E[m_]:=diagram2bup21[m]+diagram2bdown21[m]+diagram2dup21[m]+diagram2ddown21[m]
summedUpDiagrams12E[m_]:=diagram2bup12[m]+diagram2bdown12[m]+diagram2dup12[m]+diagram2ddown12[m]
summedUpDiagrams41O[m_]:=diagram2aup41[m]+diagram2adown41[m]+diagram2cup41[m]+diagram2cdown41[m]
summedUpDiagrams14O[m_]:=diagram2aup14[m]+diagram2adown14[m]+diagram2cup14[m]+diagram2cdown14[m]
summedUpDiagrams42O[m_]:=diagram2aup42[m]+diagram2adown42[m]+diagram2cup42[m]+diagram2cdown42[m]
summedUpDiagrams24O[m_]:=diagram2aup24[m]+diagram2adown24[m]+diagram2cup24[m]+diagram2cdown24[m]
summedUpDiagrams43O[m_]:=diagram2aup43[m]+diagram2adown43[m]+diagram2cup43[m]+diagram2cdown43[m]
summedUpDiagrams34O[m_]:=diagram2aup34[m]+diagram2adown34[m]+diagram2cup34[m]+diagram2cdown34[m]
summedUpDiagrams32O[m_]:=diagram2aup32[m]+diagram2adown32[m]+diagram2cup32[m]+diagram2cdown32[m]
summedUpDiagrams23O[m_]:=diagram2aup23[m]+diagram2adown23[m]+diagram2cup23[m]+diagram2cdown23[m]
summedUpDiagrams31O[m_]:=diagram2aup31[m]+diagram2adown31[m]+diagram2cup31[m]+diagram2cdown31[m]
summedUpDiagrams13O[m_]:=diagram2aup13[m]+diagram2adown13[m]+diagram2cup13[m]+diagram2cdown13[m]
summedUpDiagrams21O[m_]:=diagram2aup21[m]+diagram2adown21[m]+diagram2cup21[m]+diagram2cdown21[m]
summedUpDiagrams12O[m_]:=diagram2aup12[m]+diagram2adown12[m]+diagram2cup12[m]+diagram2cdown12[m]


factor[k_,csw_]:=0.5 * k * csw
summedUpDiagrams4WithFactorE[k_,csw_,m_]:=factor[k,csw]*(summedUpDiagrams41E[m]+summedUpDiagrams42E[m]+summedUpDiagrams43E[m])
summedUpDiagrams3WithFactorE[k_,csw_,m_]:=factor[k,csw]*(summedUpDiagrams31E[m]+summedUpDiagrams32E[m]+summedUpDiagrams34E[m])
summedUpDiagrams2WithFactorE[k_,csw_,m_]:=factor[k,csw]*(summedUpDiagrams21E[m]+summedUpDiagrams23E[m]+summedUpDiagrams24E[m])
summedUpDiagrams1WithFactorE[k_,csw_,m_]:=factor[k,csw]*(summedUpDiagrams12E[m]+summedUpDiagrams13E[m]+summedUpDiagrams14E[m])
summedUpDiagrams4WithFactorO[k_,csw_,m_]:=factor[k,csw]*(summedUpDiagrams41O[m]+summedUpDiagrams42O[m]+summedUpDiagrams43O[m])
summedUpDiagrams3WithFactorO[k_,csw_,m_]:=factor[k,csw]*(summedUpDiagrams31O[m]+summedUpDiagrams32O[m]+summedUpDiagrams34O[m])
summedUpDiagrams2WithFactorO[k_,csw_,m_]:=factor[k,csw]*(summedUpDiagrams21O[m]+summedUpDiagrams23O[m]+summedUpDiagrams24O[m])
summedUpDiagrams1WithFactorO[k_,csw_,m_]:=factor[k,csw]*(summedUpDiagrams12O[m]+summedUpDiagrams13O[m]+summedUpDiagrams14O[m])


fermionForceClover2Eo0E[k_,csw_,m_]:=algebraElement[summedUpDiagrams4WithFactorE[k,csw,m]]
fermionForceClover2Eo1E[k_,csw_,m_]:=algebraElement[summedUpDiagrams1WithFactorE[k,csw,m]]
fermionForceClover2Eo2E[k_,csw_,m_]:=algebraElement[summedUpDiagrams2WithFactorE[k,csw,m]]
fermionForceClover2Eo3E[k_,csw_,m_]:=algebraElement[summedUpDiagrams3WithFactorE[k,csw,m]]
fermionForceClover2Eo0O[k_,csw_,m_]:=algebraElement[summedUpDiagrams4WithFactorO[k,csw,m]]
fermionForceClover2Eo1O[k_,csw_,m_]:=algebraElement[summedUpDiagrams1WithFactorO[k,csw,m]]
fermionForceClover2Eo2O[k_,csw_,m_]:=algebraElement[summedUpDiagrams2WithFactorO[k,csw,m]]
fermionForceClover2Eo3O[k_,csw_,m_]:=algebraElement[summedUpDiagrams3WithFactorO[k,csw,m]]


(*Here we take into account that the calculation of the force with the corresponding gaugemomentum update concerns only ODD/EVEN sites. For all other sites (and directions) the gaugemomentum is unchanged.*)
FermionForceClover2E[gm_,k_,csw_,m_]:= squareNorm[gm]*4 +squareNorm[gm + fermionForceClover2Eo0E[k,csw,m]] +
											squareNorm[gm + fermionForceClover2Eo1E[k,csw,m]] +
											squareNorm[gm + fermionForceClover2Eo2E[k,csw,m]] +
											squareNorm[gm + fermionForceClover2Eo3E[k,csw,m]]
FermionForceClover2O[gm_,k_,csw_,m_]:= squareNorm[gm]*4 +squareNorm[gm + fermionForceClover2Eo0O[k,csw,m]] +
											squareNorm[gm + fermionForceClover2Eo1O[k,csw,m]] +
											squareNorm[gm + fermionForceClover2Eo2O[k,csw,m]] +
											squareNorm[gm + fermionForceClover2Eo3O[k,csw,m]]
FermionForceClover2E[gaugeMomOne,nonTrivialRealPar,nonTrivialRealPar,cold12x12mat]*4^4/2
FermionForceClover2O[gaugeMomAscending,nonTrivialRealPar,nonTrivialRealPar,cold12x12mat]*4^4/2
FermionForceClover2O[gaugeMomOne,nonTrivialRealPar,nonTrivialRealPar,nonTrivial12x12mat]*4^4/2
FermionForceClover2E[gaugeMomAscending,nonTrivialRealPar,nonTrivialRealPar,nonTrivial12x12mat]*4^4/2


FermionForceClover2E[gaugeMomOne,nonTrivialRealPar,nonTrivialRealPar,cold12x12mat]
FermionForceClover2O[gaugeMomAscending,nonTrivialRealPar,nonTrivialRealPar,cold12x12mat]
FermionForceClover2O[gaugeMomOne,nonTrivialRealPar,nonTrivialRealPar,nonTrivial12x12mat]
FermionForceClover2E[gaugeMomAscending,nonTrivialRealPar,nonTrivialRealPar,nonTrivial12x12mat]




(* ::Package:: *)

SetOptions[SelectedNotebook[], PrintPrecision -> 16]

(*Get the needed packages for these tests*)
Get[FileNameJoin[{ParentDirectory[ParentDirectory[ParentDirectory[NotebookDirectory[]]]], "packages/threeBythreeMatrices.m"}], Path -> {NotebookDirectory[]}]
Get[FileNameJoin[{ParentDirectory[ParentDirectory[ParentDirectory[NotebookDirectory[]]]], "packages/pauliMatricesTensorColorIdentity.m"}], Path -> {NotebookDirectory[]}]


(*Field-Strenght-Tensor
FieldStrengthTensor[u_] := 4*u.u.u.u - 4*ConjugateTranspose[u].ConjugateTranspose[u].ConjugateTranspose[u].ConjugateTranspose[u]*)
FieldStrengthTensor[u_] := u.u.ConjugateTranspose[u].ConjugateTranspose[u]+u.ConjugateTranspose[u].ConjugateTranspose[u].u+ConjugateTranspose[u].ConjugateTranspose[u].u.u+ConjugateTranspose[u].u.u.ConjugateTranspose[u]-u.u.ConjugateTranspose[u].ConjugateTranspose[u]-ConjugateTranspose[u].u.u.ConjugateTranspose[u]-ConjugateTranspose[u].ConjugateTranspose[u].u.u-u.ConjugateTranspose[u].ConjugateTranspose[u].u
(*Field-Strength-Tensor tests*)
FieldStrengthTensor[cold3x3mat];
FieldStrengthTensor[nonTrivial3x3mat];


(*Build upperleft 6x6 block of 1 + T*)
csw = 0.5;
factor = 1/16 * csw * I;
F = FieldStrengthTensor[nonTrivial3x3mat];
EPlusB = 8*F + ( 4*F - 4*F );
EMinusB = 8*F - ( 4*F - 4*F );
UpperLeftBlock = IdentityMatrix[6] + factor * (Pauli1[EMinusB] + Pauli2[EMinusB] + Pauli3[EMinusB]); 
LowerRightBlock = IdentityMatrix[6] - factor * (Pauli1[EPlusB] + Pauli2[EPlusB] + Pauli3[EPlusB]);


(*Invert 6x6 blocks*)
UpperLeftBlockInverse = Inverse[UpperLeftBlock, Method -> "CofactorExpansion"];
LowerRightBlockInverse = Inverse[LowerRightBlock];


(*Sum up all entries of 6x6 blocks*)
Ns = 4; Nt = 4;
UpperLeftBlockInverseSum = Total [UpperLeftBlockInverse, 2];
UpperLeftBlockInverseKernelResult = Ns^3 * Nt * UpperLeftBlockInverseSum;

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
FieldStrengthTensor0i[u_,v_]:=u.v.ConjugateTranspose[u].ConjugateTranspose[v]+v.ConjugateTranspose[u].ConjugateTranspose[v].u+ConjugateTranspose[u].ConjugateTranspose[v].u.v+ConjugateTranspose[v].u.v.ConjugateTranspose[u]-
v.u.ConjugateTranspose[v].ConjugateTranspose[u]-ConjugateTranspose[u].v.u.ConjugateTranspose[v]-ConjugateTranspose[v].ConjugateTranspose[u].v.u-u.ConjugateTranspose[v].ConjugateTranspose[u].v


(*Build upperleft 6x6 block of 1 + T for uniformly and non-uniformly filled gaugefield*)
csw = 0.5;
factor = 1/16 * csw * I;
Ns = 4; Nt = 4;
(*Build upperleft 6x6 block of 1 + T for uniformly filled gaugefield*)
F = FieldStrengthTensor[nonTrivial3x3mat];
EPlusB = 8*F + ( 4*F - 4*F );
EMinusB = 8*F - ( 4*F - 4*F );
UpperLeftBlock = IdentityMatrix[6] + factor * (Pauli1[EMinusB] + Pauli2[EMinusB] + Pauli3[EMinusB]); 
LowerRightBlock = IdentityMatrix[6] - factor * (Pauli1[EPlusB] + Pauli2[EPlusB] + Pauli3[EPlusB]);

(*Invert 6x6 blocks*)
UpperLeftBlockInverse = Inverse[UpperLeftBlock, Method -> "CofactorExpansion"];
LowerRightBlockInverse = Inverse[LowerRightBlock];

(*Sum up all entries of 6x6 blocks*)
UpperLeftBlockInverseSum = Total [UpperLeftBlockInverse, 2];
UpperLeftBlockInverseKernelResult = Ns^3 * Nt * UpperLeftBlockInverseSum

(*Build upperleft 6x6 block of 1 + T for non-uniformly filled gaugefield (es. ascending in one dir, nonTrivial in the other) the NU suffix stands for Not Uniform*)
F0i = FieldStrengthTensor0i[ascending3x3mat,nonTrivial3x3mat];
EPlusBNU = 8*F0i + ( 4*F - 4*F );
EMinusBNU = 8*F0i - ( 4*F - 4*F );
UpperLeftBlockNU = IdentityMatrix[6] + factor * (Pauli1[EMinusBNU] + Pauli2[EMinusBNU] + Pauli3[EMinusBNU]); 
LowerRightBlockNU = IdentityMatrix[6] - factor * (Pauli1[EPlusBNU] + Pauli2[EPlusBNU] + Pauli3[EPlusBNU]);

(*Invert 6x6 blocks*)
UpperLeftBlockInverseNU = Inverse[UpperLeftBlockNU, Method -> "CofactorExpansion"];
LowerRightBlockInverseNU = Inverse[LowerRightBlockNU];

(*Sum up all entries of 6x6 blocks*)
UpperLeftBlockInverseSumNU = Total [UpperLeftBlockInverseNU, 2];
UpperLeftBlockInverseKernelResultNU = Ns^3 * Nt * UpperLeftBlockInverseSumNU




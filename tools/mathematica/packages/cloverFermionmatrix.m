(* ::Package:: *)

SetOptions[SelectedNotebook[], PrintPrecision -> 16]
Get["cloverPackage.m",Path->{NotebookDirectory[]}]


(*Clover Fermionmatrix for a specific site, for constant gauge- and \
spinorfield*)


BeginPackage["clover`"]

Clover::usage =
	"compute Clover Matrix 1+T taking kappa, csw,  a spinor filling and a 3x3matrix filling."
CloverInverse::usage = 
	"compute inverse Clover Matrix (1+T)^(-1) taking kappa, csw,  a spinor filling and a 3x3matrix filling."


Begin["Private`"]
Needs["cloverPackage`"]
Needs["PauliMatricesTensorColorIdentity`"]

Clover[k_, csw_, s_, u_, v_] :=
	Module[{clover = ArrayFlatten[{{UpperLeftBlock[csw,u,v],0},{0,LowerRightBlock[csw,u,v]}}].s},
	clover
	]

CloverInverse[k_, csw_, s_, u_, v_] :=
	Module[{cloverInv = ArrayFlatten[{{Inverse[UpperLeftBlock[csw,u,v], Method -> "CofactorExpansion"],0},{0,Inverse[LowerRightBlock[csw,u,v], Method -> "CofactorExpansion"]}}].s},
	cloverInv
	]

End[]

EndPackage[]




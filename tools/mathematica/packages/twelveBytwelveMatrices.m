(* ::Package:: *)

BeginPackage["TwelveByTwelveMatrices`"]
cold12x12mat::usage =
	"cold12x12mat gives the 12-by-12 block diagonal matrix with two 6-by-6 identity matrices as diagonal blocks"
nonTrivial12x12mat::usage = 
	"nonTrivial12x12mat gives a non trivial 12-by-12 block diagonal matrix with two 6-by-6 non trivial matrices as diagonal blocks"

Begin["Private`"]

cold12x12mat:=
	Module[{cold6x6mat=ArrayFlatten[{{IdentityMatrix[6],ConstantArray[0, {6, 6}]},{ConstantArray[0, {6, 6}],IdentityMatrix[6]}}]},
	cold6x6mat
	]

nonTrivial12x12mat:=
	Module[{nonTrivial6x6mat=ArrayFlatten[{{Table[6*(x-1)+y,{x,6},{y,6}],ConstantArray[0, {6, 6}]},{ConstantArray[0, {6, 6}],Table[6*(x-1)+y,{x,6},{y,6}]}}]},
	nonTrivial6x6mat
	]

End[]

EndPackage[]







(* ::Package:: *)

SetOptions[SelectedNotebook[], PrintPrecision -> 16]


Get["threeBythreeMatrices.m",Path->{NotebookDirectory[]}]
Get["twelveComponentsVectors.m",Path->{NotebookDirectory[]}]


gamma1:={{0, 0, 0, -I}, {0, 0, -I, 0}, {0, I, 0, 0}, {I, 0, 0, 0}}
gamma2:={{0, 0, 0, -1.}, {0, 0, 1., 0}, {0, 1., 0, 0}, {-1., 0, 0, 0}}
gamma3:={{0, 0, -I, 0}, {0, 0, 0, I}, {I, 0, 0, 0}, {0, -I, 0, 0}}
gamma4:={{0, 0, -1., 0}, {0, 0, 0, -1.}, {-1., 0, 0, 0}, {0, -1., 0, 0}}


BeginPackage["SigmaMuNuMatricesTensorColorIdentity`"]

Sigma14::usage = 
	"Sigma14[x] gives the sigma_{14} = (I/2.)*[gamma_1,gamma_4] matrix, tensor multiplied
	with the 3-by-3 x matrix in color space."

Sigma41::usage = 
	"Sigma41[x] gives the sigma_{41} = (I/2.)*[gamma_4,gamma_1] matrix, tensor multiplied
	with the 3-by-3 x matrix in color space."

Sigma24::usage = 
	"Sigma24[x] gives the sigma_{24} = (I/2.)*[gamma_2,gamma_4] matrix, tensor multiplied
	with the 3-by-3 x matrix in color space."

Sigma42::usage = 
	"Sigma42[x] gives the sigma_{42} = (I/2.)*[gamma_4,gamma_2] matrix, tensor multiplied
	with the 3-by-3 x matrix in color space."

Sigma34::usage = 
	"Sigma34[x] gives the sigma_{34} = (I/2.)*[gamma_3,gamma_4] matrix, tensor multiplied
	with the 3-by-3 x matrix in color space."

Sigma43::usage = 
	"Sigma43[x] gives the sigma_{43} = (I/2.)*[gamma_4,gamma_3] matrix, tensor multiplied
	with the 3-by-3 x matrix in color space."

Sigma12::usage = 
	"Sigma14[x] gives the sigma_{12} = (I/2.)*[gamma_1,gamma_2] matrix, tensor multiplied
	with the 3-by-3 x matrix in color space."

Sigma21::usage = 
	"Sigma21[x] gives the sigma_{21} = (I/2.)*[gamma_2,gamma_1] matrix, tensor multiplied
	with the 3-by-3 x matrix in color space."

Sigma13::usage = 
	"Sigma13[x] gives the sigma_{13} = (I/2.)*[gamma_1,gamma_3] matrix, tensor multiplied
	with the 3-by-3 x matrix in color space."

Sigma31::usage = 
	"Sigma31[x] gives the sigma_{31} = (I/2.)*[gamma_3,gamma_1] matrix, tensor multiplied
	with the 3-by-3 x matrix in color space."

Sigma23::usage = 
	"Sigma23[x] gives the sigma_{23} = (I/2.)*[gamma_2,gamma_3] matrix, tensor multiplied
	with the 3-by-3 x matrix in color space."

Sigma32::usage = 
	"Sigma32[x] gives the sigma_{32} = (I/2.)*[gamma_3,gamma_2] matrix, tensor multiplied
	with the 3-by-3 x matrix in color space."

Begin["Private`"]

Sigma14[x_] := 
	Module[ {sigma14=ArrayFlatten[{{0, -1.*x, 0, 0}, {-1.*x, 0, 0, 0}, {0, 0, 0, 1.*x}, {0, 0, 1.*x,0}}]},
	sigma14
	]

Sigma41[x_] := 
	Module[ {sigma41=ArrayFlatten[{{0, 1.*x, 0, 0}, {1.*x, 0, 0, 0}, {0, 0, 0, -1.*x}, {0, 0, -1.*x,0}}]},
	sigma41
	]

Sigma24[x_] := 
	Module[ {sigma24=ArrayFlatten[{{0, I*x, 0, 0}, {-I*x, 0, 0, 0}, {0, 0, 0, -I*x}, {0, 0, I*x, 0}}]},
	sigma24
	]

Sigma42[x_] := 
	Module[ {sigma42=ArrayFlatten[{{0, -I*x, 0, 0}, {I*x, 0, 0, 0}, {0, 0, 0, I*x}, {0, 0, -I*x, 0}}]},
	sigma42
	]

Sigma34[x_] := 
	Module[ {sigma34=ArrayFlatten[{{-1.*x, 0, 0, 0}, {0, 1.*x, 0, 0}, {0, 0, 1.*x, 0}, {0, 0, 0, -1.*x}}]},
	sigma34
	]

Sigma43[x_] := 
	Module[ {sigma43=ArrayFlatten[{{1.*x, 0, 0, 0}, {0, -1.*x, 0, 0}, {0, 0, -1.*x, 0}, {0, 0, 0, 1.*x}}]},
	sigma43
	]

Sigma12[x_] := 
	Module[ {sigma12=ArrayFlatten[{{-1.*x, 0, 0, 0}, {0, 1.*x, 0, 0}, {0, 0, -1.*x, 0}, {0, 0, 0, 1.*x}}]},
	sigma12
	]

Sigma21[x_] := 
	Module[ {sigma21=ArrayFlatten[{{1.*x, 0, 0, 0}, {0, -1.*x, 0, 0}, {0, 0, 1.*x, 0}, {0, 0, 0, -1.*x}}]},
	sigma21
	]

Sigma13[x_] := 
	Module[ {sigma13=ArrayFlatten[{{0, -I*x, 0, 0}, {I*x, 0, 0, 0}, {0, 0, 0, -I*x}, {0, 0, I*x,0}}]},
	sigma13
	]

Sigma31[x_] := 
	Module[ {sigma31=ArrayFlatten[{{0, I*x, 0, 0}, {-I*x, 0, 0, 0}, {0, 0, 0, I*x}, {0, 0, -I*x,0}}]},
	sigma31
	]

Sigma23[x_] := 
	Module[ {sigma23=ArrayFlatten[{{0, -1.*x, 0, 0}, {-1.*x, 0, 0, 0}, {0, 0, 0, -1.*x}, {0, 0, -1.*x, 0}}]},
	sigma23
	]

Sigma32[x_] := 
	Module[ {sigma32=ArrayFlatten[{{0, 1.*x, 0, 0}, {1.*x, 0, 0, 0}, {0, 0, 0, 1.*x}, {0, 0, 1.*x, 0}}]},
	sigma32
	]

End[]

EndPackage[]


Sigma14[cold3x3mat].spinorCold;


(*Sigma14[x_]:=ArrayFlatten[{{0, -1.*x, 0, 0}, {-1.*x, 0, 0, 0}, {0, 0, 0, 1.*x}, {0, 0, 1.*x,0}}]
Sigma41[x_]:=ArrayFlatten[{{0, 1.*x, 0, 0}, {1.*x, 0, 0, 0}, {0, 0, 0, -1.*x}, {0, 0, -1.*x,0}}]
Sigma13[x_]:=ArrayFlatten[{{0, -I*x, 0, 0}, {I*x, 0, 0, 0}, {0, 0, 0, -I*x}, {0, 0, I*x,0}}]
Sigma31[x_]:=ArrayFlatten[{{0, I*x, 0, 0}, {-I*x, 0, 0, 0}, {0, 0, 0, I*x}, {0, 0, -I*x,0}}]
Sigma12[x_]:=ArrayFlatten[{{-1.*x, 0, 0, 0}, {0, 1.*x, 0, 0}, {0, 0, -1.*x, 0}, {0, 0, 0, 1.*x}}]
Sigma21[x_]:=ArrayFlatten[{{1.*x, 0, 0, 0}, {0, -1.*x, 0, 0}, {0, 0, 1.*x, 0}, {0, 0, 0, -1.*x}}]
Sigma23[x_]:=ArrayFlatten[{{0, -1.*x, 0, 0}, {-1.*x, 0, 0, 0}, {0, 0, 0, -1.*x}, {0, 0, -1.*x, 0}}]
Sigma32[x_]:=ArrayFlatten[{{0, 1.*x, 0, 0}, {1.*x, 0, 0, 0}, {0, 0, 0, 1.*x}, {0, 0, 1.*x, 0}}]
Sigma24[x_]:=ArrayFlatten[{{0, I*x, 0, 0}, {-I*x, 0, 0, 0}, {0, 0, 0, -I*x}, {0, 0, I*x, 0}}]
Sigma42[x_]:=ArrayFlatten[{{0, -I*x, 0, 0}, {I*x, 0, 0, 0}, {0, 0, 0, I*x}, {0, 0, -I*x, 0}}]
Sigma34[x_]:=ArrayFlatten[{{-1.*x, 0, 0, 0}, {0, 1.*x, 0, 0}, {0, 0, 1.*x, 0}, {0, 0, 0, -1.*x}}]
Sigma43[x_]:=ArrayFlatten[{{1.*x, 0, 0, 0}, {0, -1.*x, 0, 0}, {0, 0, -1.*x, 0}, {0, 0, 0, 1.*x}}]*)

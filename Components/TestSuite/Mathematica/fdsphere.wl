(* ::Package:: *)

BeginPackage["FDSphere`",{"NumericalDifferentialEquationAnalysis`"}];


$FDSphereOrder::usage="Order of finite differences schemes";
fdsphgrid::usage="fdsphgrid[n] spherical finite difference grid";
fdsphintegral::usage="fdintegral[rg, p, k] compute integral weights for given grid rg and order p";
fdsphdiff::usage="fdsphdiff[rg, p, k] compute differentiation matrix of order p";
fdsphlapl::usage="fdsphlapl[rg, k] Compute the finite differences spherical laplacian";


Begin["`Private`"];


$FDSphereOrder=4;
fdsphgrid[n_]=Subdivide[0,1,n-1];
fdsphdiff[rg_,p_,k_:$FDSphereOrder]:=Module[{ds},
ds=NDSolve`FiniteDifferenceDerivative[p,rg,"DifferenceOrder"->k]["DifferentiationMatrix"];
ds]
fdsphintegral[rg_,p_,k_:$FDSphereOrder]:=Module[{drm,drp,ds,wgts},
drm=Flatten[{rg[[1]],rg[[1;;-2]]}]-rg;
drp=Flatten[{rg[[2;;]],rg[[-1]]}]-rg;
ds=Table[NDSolve`FiniteDifferenceDerivative[i,rg,"DifferenceOrder"->k]["DifferentiationMatrix"],{i,0,p}];
wgts=Sum[1/(2 i!) (drp^i-drm^i) . ds[[i]],{i,1,Length[ds]}];
wgts
]
fdsphlapl[rg_,l_,k_:$FDSphereOrder]:=Module[{ds,slapl,invrg},
ds=Table[NDSolve`FiniteDifferenceDerivative[i,rg,"DifferenceOrder"->k]["DifferentiationMatrix"],{i,0,2}];
invrg=Flatten[{0,1/rg[[2;;]]}];
slapl=ds[[3]]+2 invrg ds[[2]]-l(l+1)invrg^2 ds[[1]];
slapl]


End[];


EndPackage[];

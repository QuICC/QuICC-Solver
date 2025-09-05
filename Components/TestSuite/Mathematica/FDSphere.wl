(* ::Package:: *)

BeginPackage["FDSphere`",{"NumericalDifferentialEquationAnalysis`"}];


$mpprec::usage="Precision for Worland computations, default is 100";
$FDSphereOrder::usage="Order of finite differences schemes";
$FDSphereTaylorOrder::usage="Order of Taylor expansion for integration";
fdsphgrid::usage="fdsphgrid[n] default spherical finite difference grid";
fdsphchegrid::usage="fdsphchegrid[n] Chebyshev spherical finite difference grid";
fdsphunigrid::usage="fdsphunigrid[n] Uniform spherical finite difference grid";

fdsphintegral::usage="fdintegral[rg, l, t, k] compute integral weights for given grid rg and order p";

fdsphzt::usage="fdsphzt[n] zero top";
fdsphzb::usage="fdsphzb[n] zero bottom";
fdsphztb::usage="fdsphztb[n] zero top and bottom";
fdsphid::usage="fdsphid[rg, l] identity";
fdsphproj::usage="fdsphidz[rg, l] identity with zeros are r=0, r=1";
fdsphr1::usage="fdsphr1[rg, l] identity";
fdsphdivr1::usage="fdsphidivr1[rg, l] identity";
fdsphdivr1d1r1::usage="fdsphdivr1d1[rg, l, k] identity";
fdsphdiff::usage="fdsphdiff[rg, l, p, k] compute differentiation matrix of order p";
fdsphd1::usage="fdsphd1[rg, l, k] first derivative";
fdsphd2::usage="fdsphd2[rg, l, k] second derivative";
fdsphd3::usage="fdsphd3[rg, l, k] third derivative";
fdsphd4::usage="fdsphd[rg, l, k] fourth derivative";
fdsphlapl::usage="fdsphlapl[rg, l, k] Compute the finite differences spherical laplacian";

fdsphzerol0::usage="fdsphzerol0[f] zero if l = 0";


Begin["`Private`"];


If[!ValueQ[$mpprec], $mpprec=100];


If[!ValueQ[$FDSphereOrder], $FDSphereOrder=2];
If[!ValueQ[$FDSphereTaylorOrder], $FDSphereTaylorOrder=4];


fdsphunigrid[n_]:=N[Subdivide[0,1,n-1],$mpprec];
fdsphchegrid[n_]:=N[Flatten[{0,Table[Cos[(2k-1)/(4n) \[Pi]],{k,n,1,-1}],1}],$mpprec];
fdsphgrid[n_]:=fdsphunigrid[n];


fdsphzt[n_]:=DiagonalMatrix[Flatten[{0,Table[1,n-1]}]]
fdsphzb[n_]:=DiagonalMatrix[Flatten[{Table[1,n-1],0}]]
fdsphztb[n_]:=DiagonalMatrix[Flatten[{0,Table[1,n-2],0}]]


fdsphdiff[rg_,l_,p_,k_:$FDSphereOrder]:=Module[{ds},
ds=NDSolve`FiniteDifferenceDerivative[p,rg,"DifferenceOrder"->k]["DifferentiationMatrix"];
ds]
fdsphr1[rg_,l_]:=fdsphztb[Length@rg] . DiagonalMatrix[rg];
fdsphdivr1[rg_,l_]:=fdsphztb[Length@rg] . DiagonalMatrix[Flatten[{0,Table[1/rg[[i]],{i,2,Length@rg}]}]];
fdsphdivr1d1r1[rg_,l_,k_:$FDSphereOrder]:=fdsphdivr1[rg,l] . fdsphdiff[rg,l,1,k] . DiagonalMatrix[rg];
fdsphid[rg_,l_]:=fdsphdiff[rg,l,0,$FDSphereOrder];
fdsphproj[rg_,l_]:=If[l>0,fdsphzt[Length@rg] . fdsphid[rg,l],fdsphid[rg,l]];
fdsphd1[rg_,l_,k_:$FDSphereOrder]:=fdsphztb[Length@rg] . fdsphdiff[rg,l,1,k];
fdsphd2[rg_,l_,k_:$FDSphereOrder]:=fdsphztb[Length@rg] . fdsphdiff[rg,l,2,k];
fdsphd3[rg_,l_,k_:$FDSphereOrder]:=fdsphztb[Length@rg] . fdsphdiff[rg,l,3,k];
fdsphd4[rg_,l_,k_:$FDSphereOrder]:=fdsphztb[Length@rg] . fdsphdiff[rg,l,4,k];
fdsphlapl[rg_,l_,k_:$FDSphereOrder]:=Module[{ds,slapl,invrg},
ds=Table[NDSolve`FiniteDifferenceDerivative[i,rg,"DifferenceOrder"->k]["DifferentiationMatrix"],{i,0,2}];
invrg=Flatten[{0,1/rg[[2;;]]}];
slapl=ds[[3]]+2 invrg ds[[2]]-l(l+1)invrg^2 ds[[1]];
fdsphztb[Length@rg] . slapl]


fdsphintegral[rg_,l_,t_:$FDSphereTaylorOrder,k_:$FDSphereOrder]:=Module[{drm,drp,ds,wgts},
drm=Flatten[{rg[[1]],rg[[1;;-2]]}]-rg;
drp=Flatten[{rg[[2;;]],rg[[-1]]}]-rg;
ds=Table[NDSolve`FiniteDifferenceDerivative[i,rg,"DifferenceOrder"->k]["DifferentiationMatrix"],{i,0,t}];
wgts=Sum[1/(2 i!) (drp^i-drm^i) . ds[[i]],{i,1,Length[ds]}];
wgts
]


fdsphzerol0[f_]:=Module[{g},g=Function[{rg,l},If[l>0,f[rg,l],SparseArray[{},{Length[rg],Length[rg]}]]];g]


End[];


EndPackage[];

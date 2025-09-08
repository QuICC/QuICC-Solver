(* ::Package:: *)

BeginPackage["FDSphere`",{"NumericalDifferentialEquationAnalysis`"}];


$mpprec::usage="Precision for Worland computations, default is 100";
$FDSphereOrder::usage="Order of finite differences schemes";
$FDSphereTaylorOrder::usage="Order of Taylor expansion for integration";
fdsphgrid::usage="fdsphgrid[n] default spherical finite difference grid";
fdsphchegrid::usage="fdsphchegrid[n] Chebyshev spherical finite difference grid";
fdsphunigrid::usage="fdsphunigrid[n] Uniform spherical finite difference grid";

fdsphintegral::usage="fdintegral[rg, l, t, k] compute integral weights for given grid rg and order p";

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

fdsphzt::usage="fdsphzt[f] zero top";
fdsphzb::usage="fdsphzb[f] zero bottom";
fdsphztb::usage="fdsphztb[f] zero top and bottom";
fdsphzerol0::usage="fdsphzerol0[f] zero if l = 0";


Begin["`Private`"];


If[!ValueQ[$mpprec], $mpprec=100];


If[!ValueQ[$FDSphereOrder], $FDSphereOrder=2];
If[!ValueQ[$FDSphereTaylorOrder], $FDSphereTaylorOrder=4];


fdsphunigrid[n_]:=N[Subdivide[0,1,n-1],$mpprec];
fdsphchegrid[n_]:=N[Flatten[{0,Table[Cos[(2k-1)/(4n) \[Pi]],{k,n,1,-1}],1}],$mpprec];
fdsphgrid[n_]:=fdsphunigrid[n];


fdsphdiff[rg_,l_,p_,k_:$FDSphereOrder]:=Module[{ds},
ds=NDSolve`FiniteDifferenceDerivative[p,rg,"DifferenceOrder"->k]["DifferentiationMatrix"];
ds]
fdsphr1[rg_,l_]:=DiagonalMatrix[rg];
fdsphdivr1[rg_,l_]:=Module[{dvir1},divr1=DiagonalMatrix[Flatten[{0,Table[1/rg[[i]],{i,2,Length@rg}]}]]+DiagonalMatrix[Table[If[l==1&&n==1,1,0],{n,1,Length@rg}]] . fdsphd1[rg,l];divr1]
fdsphdivr1d1r1[rg_,l_,k_:$FDSphereOrder]:=Module[{divr1d1r1},divr1d1r1=fdsphdivr1[rg,l] + fdsphdiff[rg,l,1,k];If[l==0||l>1,divr1d1r1=fdsphzt[fdsphid][rg,l] . divr1d1r1];divr1d1r1]
fdsphid[rg_,l_]:=fdsphdiff[rg,l,0,$FDSphereOrder];
fdsphproj[rg_,l_]:=If[l>0,fdsphid[rg,l],fdsphid[rg,l]];
fdsphd1[rg_,l_,k_:$FDSphereOrder]:=Module[{d1},d1=fdsphdiff[rg,l,1,k];If[l==0||l>1,d1=fdsphzt[fdsphid][rg,l] . d1];d1]
fdsphd2[rg_,l_,k_:$FDSphereOrder]:=Module[{d2},d2=fdsphdiff[rg,l,2,k];If[l==1||l>2,d2=fdsphzt[fdsphid][rg,l] . d2];d2]
fdsphd3[rg_,l_,k_:$FDSphereOrder]:=Module[{d3},d3=fdsphdiff[rg,l,3,k];If[l==0||l==2||l>3,d3=fdsphzt[fdsphid][rg,l] . d3];d3]
fdsphd4[rg_,l_,k_:$FDSphereOrder]:=Module[{d4},d4=fdsphdiff[rg,l,4,k];If[l==1||l==3||l>4,d4=fdsphzt[fdsphid][rg,l] . d4];d4]
fdsphlapl[rg_,l_,k_:$FDSphereOrder]:=Module[{id,d1,d2,slapl,invr,invr2},
id=fdsphid[rg,l];
d1=fdsphd1[rg,l,k];
d2=fdsphd2[rg,l,k];
invr=fdsphdivr1[rg,l];
invr2=invr . invr;
If[l==0,
slapl=d2+2 invr . d1;,
slapl=d2+ 2 fdsphdivr1[rg,l-1] . d1-l(l+1)invr2;
slapl=fdsphzt[fdsphid][rg,l] . slapl;
];
slapl]


fdsphintegral[rg_,l_,t_:$FDSphereTaylorOrder,k_:$FDSphereOrder]:=Module[{drm,drp,ds,wgts},
drm=Flatten[{rg[[1]],rg[[1;;-2]]}]-rg;
drp=Flatten[{rg[[2;;]],rg[[-1]]}]-rg;
ds=Table[NDSolve`FiniteDifferenceDerivative[i,rg,"DifferenceOrder"->k]["DifferentiationMatrix"],{i,0,t}];
wgts=Sum[1/(2 i!) (drp^i-drm^i) . ds[[i]],{i,1,Length[ds]}];
wgts
]


fdsphzt[f_]:=Module[{g},g=Function[{rg,l},DiagonalMatrix[Flatten[{0,Table[1,Length@rg-1]}] . f[rg,l]]];g]
fdsphzb[f_]:=Module[{g},g=Function[{rg,l},DiagonalMatrix[Flatten[{Table[1,Length@rg-1],0}] . f[rg,l]]];g]
fdsphztb[f_]:=Module[{g},g=Function[{rg,l},DiagonalMatrix[Flatten[{0,Table[1,Length@rg-2],0}] . f[rg,l]]];g]
fdsphzerol0[f_]:=Module[{g},g=Function[{rg,l},If[l>0,f[rg,l],SparseArray[{},{Length[rg],Length[rg]}]]];g]


End[];


EndPackage[];

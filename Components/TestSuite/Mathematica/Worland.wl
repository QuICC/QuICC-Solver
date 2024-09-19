(* ::Package:: *)

BeginPackage["Worland`",{"NumericalDifferentialEquationAnalysis`"}];


$WorlandType::usage="Type of Jones-Worland basis to use: Chebyshev, Legendre, CylEnergy or SphEnergy";
$mpprec::usage="Precision for Worland computations, default is 100";
$wche\[Alpha]::usage="\[Alpha] parameter of Chebyshev type Jacobi";
$wched\[Beta]::usage="\[Beta] = l + $wd\[Beta] parameter of Chebyshev type Jacobi";
wchegrid::usage="wchegrid[n] computes Chebyshev type quadrature grid";
wcheweights::usage="wcheweights[n] computes Chebyshev type quadrature weights";
$wleg\[Alpha]::usage="\[Alpha] parameter of Legendre type Jacobi";
$wlegd\[Beta]::usage="\[Beta] = l + $wd\[Beta] parameter of Legendre type Jacobi";
wleggrid::usage="wleggrid[n] computes Legendre type quadrature grid";
wlegweights::usage="wlegweights[n] computes Legendre type quadrature weights";
$wcyl\[Alpha]::usage="\[Alpha] parameter of cylindrical energy type Jacobi";
$wcyld\[Beta]::usage="\[Beta] = l + $wd\[Beta] parameter of cyindrical energy type Jacobi";
wcylgrid::usage="wcylgrid[n] computes cylindrical energy type quadrature grid";
wcylweights::usage="wcylweights[n] computes cylindrical energy type quadrature weights";
$wsph\[Alpha]::usage="\[Alpha] parameter of spherical energy type Jacobi";
$wsphd\[Beta]::usage="\[Beta] = l + $wd\[Beta] parameter of spherical energy type Jacobi";
wsphgrid::usage="wsphgrid[n] computes spherical energy type quadrature grid";
wsphweights::usage="wsphweights[n] computes spherical energy type quadrature weights";

(*Generic interface*)
$w\[Alpha]::usage="\[Alpha] parameter of Jacobi polynomial";
$wd\[Beta]::usage="\[Beta] = l + $wd\[Beta] parameter of Jacobi polynomial";
wgrid::usage="wgrid[n] computes the Worland quadrature grid";
wweights::usage="wweights[n] computes the Worland quadrature weights";
wnorm::usage="wnorm[k,\[Alpha],\[Beta]] norm of Worland polynomial k";

(* Operators to work on grid*)
Wnl::usage="Wnl[n,l,r]";
rWnl::usage="rWnl[n,l,r]";
dWnl::usage="dWnl[n,l,r]";
d2Wnl::usage="d2Wnl[n,l,r]";
slaplWnl::usage="slaplWnl[n,l,r]";
divrWnl::usage="divrWnl[n,l,r]";
drWnl::usage="drWnl[n,l,r]";
divrdrWnl::usage="divrdrWnl[n,l,r]";
rddivrWnl::usage="rddivrWnl[n,l,r]";
insulatingWnl::usage="insulatingWnl[n,l,r]";

ddivrdrWnl::usage="ddivrdrWnl[n,l,r]";
claplhWnl::usage="claplhWnl[n,l,r]";
dclaplhWnl::usage="dclaplhWnl[n,l,r]";
divrclaplhWnl::usage="divrclaplhWnl[n,l,r]";
(* More complex operators *)
divrdrWnlExplicit::usage="divrdrWnlExplicit[n,l,r]";
divrdrWnlImplicit::usage="divrdrWnlImplicit[n,l,r]";


Begin["`Private`"];


If[!ValueQ[$WorlandType], $WorlandType="Chebyshev"];
Print["Worland:: Using " <> ToString[$WorlandType] <> " Worland type"];
If[!ValueQ[$mpprec], $mpprec=100];
Print["Worland:: Number of digits is set to " <> ToString[$mpprec]];


(*Gaussian quadratures*)
(*Chebyshev type*)
$wche\[Alpha]=-1/2;$wched\[Beta]=-1/2;
wchegrid[n_]:=N[Table[Cos[(2k-1)/(4n) \[Pi]],{k,n,1,-1}],$mpprec];
wcheweights[n_]:=Table[\[Pi]/(2n),{k,1,n}];
(*Legendre quadrature*)
lgrid[n_]:=GaussianQuadratureWeights[n,-1,1,$mpprec][[;;,1]];
lweights[n_]:=GaussianQuadratureWeights[n,-1,1,$mpprec][[;;,2]];
(*Legendre type*)
$wleg\[Alpha]=0;$wlegd\[Beta]=-1/2;
wleggrid[n_]:=lgrid[2n][[n+1;;]];
wlegweights[n_]:=lweights[2n][[n+1;;]];
(*Cylindrical energy type*)
$wcyl\[Alpha]=0;$wcyld\[Beta]=0;
wcylgrid[n_]:=(lgrid[n]+1)/2;
wcylweights[n_]:=lweights[n] wcylgrid[n]/2;
(*Spherical energy type*)
$wsph\[Alpha]=0;$wsphd\[Beta]=1/2;
wsphgrid[n_]:=lgrid[2n][[n+1;;]];
wsphweights[n_]:=lweights[2n][[n+1;;]] wsphgrid[n]^2;
(*Generic interface*)
Switch[$WorlandType,
"Chebyshev",
$w\[Alpha]=$wche\[Alpha];$wd\[Beta]=$wched\[Beta];
wgrid[n_]:=wchegrid[n];
wweights[n_]:=wcheweights[n];,
"Legendre",
$w\[Alpha]=$wleg\[Alpha];$wd\[Beta]=$wlegd\[Beta];
wgrid[n_]:=wleggrid[n];
wweights[n_]:=wlegweights[n];,
"CylEnergy",
$w\[Alpha]=$wcyl\[Alpha];$wd\[Beta]=$wcyld\[Beta];
wgrid[n_]:=wcylgrid[n];
wweights[n_]:=wcylweights[n];,
"SphEnergy",
$w\[Alpha]=$wsph\[Alpha];$wd\[Beta]=$wsphd\[Beta];
wgrid[n_]:=wsphgrid[n];
wweights[n_]:=wsphweights[n];
];


wnorm[k_,\[Alpha]_,\[Beta]_]:=If[k==0,Sqrt[1/2 (Gamma[\[Alpha]+1]Gamma[\[Beta]+1])/Gamma[\[Alpha]+\[Beta]+2]],Sqrt[1/(2(2k+\[Alpha]+\[Beta]+1)) (Gamma[k+\[Alpha]+1]Gamma[k+\[Beta]+1])/(Gamma[k+\[Alpha]+\[Beta]+1]Gamma[k+1])]];


(*dWorland*)
dWorland0[l_,t_]=Simplify[D[t^l JacobiP[0,$w\[Alpha],l+$wd\[Beta],2t^2-1],{t,1}]];
dWorland[n_,l_,t_]=Simplify[D[t^l JacobiP[n,$w\[Alpha],l+$wd\[Beta],2t^2-1],{t,1}]];
(*d2Worland*)
d2Worland0[l_,t_]=Simplify[D[t^l JacobiP[0,$w\[Alpha],l+$wd\[Beta],2t^2-1],{t,2}]];
d2Worland1[l_,t_]=Simplify[D[t^l JacobiP[1,$w\[Alpha],l+$wd\[Beta],2t^2-1],{t,2}]];
d2Worland[n_,l_,t_]=Simplify[D[t^l JacobiP[n,$w\[Alpha],l+$wd\[Beta],2t^2-1],{t,2}]];
(*slaplWorland*)
slaplWorland0[l_,t_]:=If[Length[t]==0,0,Table[0,{i,1,Length[t]}]];
slaplWorland1[l_,t_]=Simplify[1/t^2 D[t^2 D[t^l JacobiP[1,$w\[Alpha],l+$wd\[Beta],2t^2-1],t],t]-(l(l+1))/t^2 t^l JacobiP[1,$w\[Alpha],l+$wd\[Beta],2t^2-1]];
slaplWorland[n_,l_,t_]=Simplify[1/t^2 D[t^2 D[t^l JacobiP[n,$w\[Alpha],l+$wd\[Beta],2t^2-1],t],t]-(l(l+1))/t^2 t^l JacobiP[n,$w\[Alpha],l+$wd\[Beta],2t^2-1]];
(*claplhWorland*)
claplhWorland0[l_,t_]:=If[Length[t]==0,0,Table[0,{i,1,Length[t]}]];
claplhWorland1[l_,t_]=Simplify[1/t D[t D[t^l JacobiP[1,$w\[Alpha],l+$wd\[Beta],2t^2-1],t],t]-l^2/t^2 t^l JacobiP[1,$w\[Alpha],l+$wd\[Beta],2t^2-1]];
claplhWorland[n_,l_,t_]=Simplify[1/t D[t D[t^l JacobiP[n,$w\[Alpha],l+$wd\[Beta],2t^2-1],t],t]-l^2/t^2 t^l JacobiP[n,$w\[Alpha],l+$wd\[Beta],2t^2-1]];
(*dclaplhWorland*)
dclaplhWorland0[l_,t_]:=If[Length[t]==0,0,Table[0,{i,1,Length[t]}]];
dclaplhWorland1[l_,t_]=Simplify[D[1/t D[t D[t^l JacobiP[1,$w\[Alpha],l+$wd\[Beta],2t^2-1],t],t]-l^2/t^2 t^l JacobiP[1,$w\[Alpha],l+$wd\[Beta],2t^2-1],t]];
dclaplhWorland2[l_,t_]=Simplify[D[1/t D[t D[t^l JacobiP[2,$w\[Alpha],l+$wd\[Beta],2t^2-1],t],t]-l^2/t^2 t^l JacobiP[2,$w\[Alpha],l+$wd\[Beta],2t^2-1],t]];
dclaplhWorland[n_,l_,t_]=Simplify[D[1/t D[t D[t^l JacobiP[n,$w\[Alpha],l+$wd\[Beta],2t^2-1],t],t]-l^2/t^2 t^l JacobiP[n,$w\[Alpha],l+$wd\[Beta],2t^2-1],t]];
(*divrclaplhWorland*)
divrclaplhWorland0[l_,t_]:=If[Length[t]==0,0,Table[0,{i,1,Length[t]}]];
divrclaplhWorland1[l_,t_]=Simplify[1/t (1/t D[t D[t^l JacobiP[1,$w\[Alpha],l+$wd\[Beta],2t^2-1],t],t]-l^2/t^2 t^l JacobiP[1,$w\[Alpha],l+$wd\[Beta],2t^2-1])];
divrclaplhWorland[n_,l_,t_]=Simplify[1/t (1/t D[t D[t^l JacobiP[n,$w\[Alpha],l+$wd\[Beta],2t^2-1],t],t]-l^2/t^2 t^l JacobiP[n,$w\[Alpha],l+$wd\[Beta],2t^2-1])];
(*drWorland*)
drWorland0[l_,t_]=Simplify[D[t^(l+1) JacobiP[0,$w\[Alpha],l+$wd\[Beta],2t^2-1],t]];
drWorland[n_,l_,t_]=Simplify[D[t^(l+1) JacobiP[n,$w\[Alpha],l+$wd\[Beta],2t^2-1],t]];
(*divrdrWorland*)
divrdrWorland0[l_,t_]=Simplify[1/t D[t^(l+1) JacobiP[0,$w\[Alpha],l+$wd\[Beta],2t^2-1],t]];
divrdrWorland[n_,l_,t_]=Simplify[1/t D[t^(l+1) JacobiP[n,$w\[Alpha],l+$wd\[Beta],2t^2-1],t]];
(*rddivrWorland*)
rddivrWorland0[l_,t_]=Simplify[t D[t^(l-1) JacobiP[0,$w\[Alpha],l+$wd\[Beta],2t^2-1],{t,1}]];
rddivrWorland[n_,l_,t_]=Simplify[t D[t^(l-1) JacobiP[n,$w\[Alpha],l+$wd\[Beta],2t^2-1],{t,1}]];
(*insulating BC*)
insulatingWorland0[l_,t_]=Simplify[D[t^l JacobiP[0,$w\[Alpha],l+$wd\[Beta],2t^2-1],{t,1}]+(l+1)t^(l-1) JacobiP[0,$w\[Alpha],l+$wd\[Beta],2t^2-1]];
insulatingWorland[n_,l_,t_]=Simplify[D[t^l JacobiP[n,$w\[Alpha],l+$wd\[Beta],2t^2-1],{t,1}]+(l+1)t^(l-1) JacobiP[n,$w\[Alpha],l+$wd\[Beta],2t^2-1]];
(*ddivrdrWorland*)
ddivrdrWorland0[l_,t_]=Simplify[D[1/t D[t^(l+1) JacobiP[0,$w\[Alpha],l+$wd\[Beta],2t^2-1],t],t]];
ddivrdrWorland1[l_,t_]=Simplify[D[1/t D[t^(l+1) JacobiP[1,$w\[Alpha],l+$wd\[Beta],2t^2-1],t],t]];
ddivrdrWorland[n_,l_,t_]=Simplify[D[1/t D[t^(l+1) JacobiP[n,$w\[Alpha],l+$wd\[Beta],2t^2-1],t],t]];


(* Operators to work on grid*)
Wnl[n_,l_,r_]:=r^l JacobiP[n,$w\[Alpha],l+$wd\[Beta],2r^2-1]/wnorm[n,$w\[Alpha],l+$wd\[Beta]]
rWnl[n_,l_,r_]:=r^(l+1) JacobiP[n,$w\[Alpha],l+$wd\[Beta],2r^2-1]/wnorm[n,$w\[Alpha],l+$wd\[Beta]]
dWnl[n_,l_,r_]:=If[n>0,dWorland[n,l,r],dWorland0[l,r]]/wnorm[n,$w\[Alpha],l+$wd\[Beta]];
d2Wnl[n_,l_,r_]:=If[n>0,If[n>1,d2Worland[n,l,r],d2Worland1[l,r]],d2Worland0[l,r]]/wnorm[n,$w\[Alpha],l+$wd\[Beta]];
slaplWnl[n_,l_,r_]:=If[n>0,If[n>1,slaplWorland[n,l,r],slaplWorland1[l,r]],slaplWorland0[l,r]]/wnorm[n,$w\[Alpha],l+$wd\[Beta]];
divrWnl[n_,l_,r_]:=r^(l-1) JacobiP[n,$w\[Alpha],l+$wd\[Beta],2r^2-1]/wnorm[n,$w\[Alpha],l+$wd\[Beta]];
drWnl[n_,l_,r_]:=If[n>0,drWorland[n,l,r],drWorland0[l,r]]/wnorm[n,$w\[Alpha],l+$wd\[Beta]];
divrdrWnl[n_,l_,r_]:=If[n>0,divrdrWorland[n,l,r],divrdrWorland0[l,r]]/wnorm[n,$w\[Alpha],l+$wd\[Beta]];
rddivrWnl[n_,l_,r_]:=If[n>0,rddivrWorland[n,l,r],rddivrWorland0[l,r]]/wnorm[n,$w\[Alpha],l+$wd\[Beta]];
insulatingWnl[n_,l_,r_]:=If[n>0,insulatingWorland[n,l,r],insulatingWorland0[l,r]]/wnorm[n,$w\[Alpha],l+$wd\[Beta]];

ddivrdrWnl[n_,l_,r_]:=If[n>0,If[n>1,ddivrdrWorland[n,l,r],ddivrdrWorland1[l,r]],ddivrdrWorland0[l,r]]/wnorm[n,$w\[Alpha],l+$wd\[Beta]];
claplhWnl[n_,l_,r_]:=If[n>0,If[n>1,claplhWorland[n,l,r],claplhWorland1[l,r]],claplhWorland0[l,r]]/wnorm[n,$w\[Alpha],l+$wd\[Beta]];
dclaplhWnl[n_,l_,r_]:=If[n>0,If[n>1,If[n>2,dclaplhWorland[n,l,r],dclaplhWorland2[l,r]],dclaplhWorland1[l,r]],dclaplhWorland0[l,r]]/wnorm[n,$w\[Alpha],l+$wd\[Beta]];
divrclaplhWnl[n_,l_,r_]:=If[n>0,If[n>1,divrclaplhWorland[n,l,r],divrclaplhWorland1[l,r]],divrclaplhWorland0[l,r]]/wnorm[n,$w\[Alpha],l+$wd\[Beta]];


(* More complex operators *)
divrdrWnlExplicit[maxN_,l_,r_,w_]:=Module[{opA,opB,opC,op},
opA=ParallelTable[Wnl[i,l,r],{i,0,maxN+1}];
opB=ParallelTable[dWnl[i,l,r],{i,0,maxN+1}];
opC=ParallelTable[Wnl[i,l,r],{i,0,maxN}];
op=(opC . DiagonalMatrix[w]) . DiagonalMatrix[1/r] . Transpose[opB] . (opA . DiagonalMatrix[w]) . DiagonalMatrix[r];
Transpose[op]
]
divrdrWnlImplicit[maxN_,l_,r_,w_]:=Module[{opA,opB,opC,op},
If[l==0,
op =divrdrWnlExplicit[maxN,l,r,w];
,
opA=ParallelTable[Wnl[i,l-1,r],{i,0,maxN+1}];
opB=ParallelTable[divrdrWnl[i,l-1,r],{i,0,maxN+1}];
opC=ParallelTable[Wnl[i,l,r],{i,0,maxN}];
op=Transpose[((opC . DiagonalMatrix[w]) . Transpose[opB] . (opA . DiagonalMatrix[w]))];
];
op
]


End[];


EndPackage[];

(* ::Package:: *)

BeginPackage["Bessel`",{"NumericalDifferentialEquationAnalysis`"}];


$mpprec::usage="Precision for Worland computations, default is 100";


bgrid::usage="wgrid[n] computes the Worland quadrature grid";
bweights::usage="wweights[n] computes the Worland quadrature weights";
norm::usage="norm[k,\[Alpha],\[Beta]] norm of Worland polynomial k";

(* Zeros *)
torZero::usage="torZero[n,l]";
polZero::usage="polZero[n,l]";
torZeros::usage="torZeros[n,l]";
polZeros::usage="polZeros[n,l]";

(* Operators to work on grid*)
TorSphJnl::usage="TorSphJnl[n,l,r]";
PolSphJnl::usage="TorSphJnl[n,l,r]";
rTorSphJnl::usage="rTorSphJnl[n,l,r]";
rPolSphJnl::usage="rTorSphJnl[n,l,r]";
dTorSphJnl::usage="dTorSphJnl[n,l,r]";
dPolSphJnl::usage="dTorSphJnl[n,l,r]";
slaplTorSphJnl::usage="slaplTorSphJnl[n,l,r]";
slaplPolSphJnl::usage="slaplTorSphJnl[n,l,r]";
divrTorSphJnl::usage="divrTorSphJnl[n,l,r]";
divrPolSphJnl::usage="divrTorSphJnl[n,l,r]";
drTorSphJnl::usage="drTorSphJnl[n,l,r]";
drPolSphJnl::usage="drTorSphJnl[n,l,r]";
divrdrTorSphJnl::usage="divrdrTorSphJnl[n,l,r]";
divrdrPolSphJnl::usage="divrdrTorSphJnl[n,l,r]";
rddivrTorSphJnl::usage="rddivrTorSphJnl[n,l,r]";
rddivrPolSphJnl::usage="rddivrTorSphJnl[n,l,r]";


Begin["`Private`"];


(*Gaussian quadratures*)
(*Legendre quadrature*)
lgrid[n_]:=GaussianQuadratureWeights[n,-1,1,$mpprec][[;;,1]];
lweights[n_]:=GaussianQuadratureWeights[n,-1,1,$mpprec][[;;,2]];
(*Spherical energy type*)
wsphgrid[n_]:=lgrid[2n][[n+1;;]];
wsphweights[n_]:=lweights[2n][[n+1;;]] wsphgrid[n]^2;
(*Generic interface*)
bgrid[n_]:=wsphgrid[n];
bweights[n_]:=wsphweights[n];


(*Roots*)
torZero[n_,l_]:=BesselJZero[l+1/2,n+1];
polZero[n_,l_]:=BesselJZero[l-1/2,n+1];
torZeros[n_,l_]:=Table[torZero[i,l],{i,0,n}];
polZeros[n_,l_]:=Table[polZero[i,l] ,{i,0,n}];


(*Jnl*)
Jnl[k_,l_,t_]:=SphericalBesselJ[l, k t]
(*rJnl*)
rJnl[k_,l_,t_]=Simplify[t SphericalBesselJ[l, k t]];
(*dJnl*)
dJnl[k_,l_,t_]=FullSimplify[D[Jnl[k,l,t],{t,1}]];
(*slaplJnl*)
slaplJnl[k_,l_,t_]=Simplify[1/t^2 D[t^2 D[Jnl[k,l,t],t],t]-(l(l+1))/t^2 Jnl[k,l,t]];
(*divrJnl*)
divrJnl[k_,l_,t_]=Simplify[1/t Jnl[k,l,t]];
(*drJnl*)
drJnl[k_,l_,t_]=Simplify[D[t Jnl[k,l,t],t]];
(*divrdrJnl*)
divrdrJnl[k_,l_,t_]=Simplify[1/t D[t Jnl[k,l,t],t]];


(* Operators to work on grid*)
TorSphJnl[n_,l_,r_]:=Jnl[torZero[n,l],l,r]
PolSphJnl[n_,l_,r_]:=Jnl[polZero[n,l],l,r]
rTorSphJnl[n_,l_,r_]:=rJnl[torZero[n,l],l,r]
rPolSphJnl[n_,l_,r_]:=rJnl[polZero[n,l],l,r]
dTorSphJnl[n_,l_,r_]:=dJnl[torZero[n,l],l,r]
dPolSphJnl[n_,l_,r_]:=dJnl[polZero[n,l],l,r]
slaplTorSphJnl[n_,l_,r_]:=slaplJnl[torZero[n,l],l,r]
slaplPolSphJnl[n_,l_,r_]:=slaplJnl[polZero[n,l],l,r]
divrTorSphJnl[n_,l_,r_]:=divrJnl[torZero[n,l],l,r]
divrPolSphJnl[n_,l_,r_]:=divrJnl[polZero[n,l],l,r]
drTorSphJnl[n_,l_,r_]:=drJnl[torZero[n,l],l,r]
drPolSphJnl[n_,l_,r_]:=drJnl[polZero[n,l],l,r]
divrdrTorSphJnl[n_,l_,r_]:=divrdrJnl[torZero[n,l],l,r]
divrdrPolSphJnl[n_,l_,r_]:=divrdrJnl[polZero[n,l],l,r]


End[];


EndPackage[];

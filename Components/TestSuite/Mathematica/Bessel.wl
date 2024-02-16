(* ::Package:: *)

BeginPackage["Bessel`",{"NumericalDifferentialEquationAnalysis`"}];


$mpprec::usage="Precision for Worland computations, default is 100";


bgrid::usage="wgrid[n] computes the Worland quadrature grid";
bweights::usage="wweights[n] computes the Worland quadrature weights";
norm::usage="norm[k,l] norm of Bessel k";

(* Zeros *)
torZero::usage="torZero[n,l]";
polZero::usage="polZero[n,l]";
torZeros::usage="torZeros[n,l]";
polZeros::usage="polZeros[n,l]";

(* Operators to work on grid*)
TorSphJnl::usage="TorSphJnl[n,l,r]";
PolSphJnl::usage="PolSphJnl[n,l,r]";
rTorSphJnl::usage="rTorSphJnl[n,l,r]";
rPolSphJnl::usage="rPolSphJnl[n,l,r]";
dTorSphJnl::usage="dTorSphJnl[n,l,r]";
dPolSphJnl::usage="dPolSphJnl[n,l,r]";
slaplTorSphJnl::usage="slaplTorSphJnl[n,l,r]";
slaplPolSphJnl::usage="slaplPolSphJnl[n,l,r]";
divrTorSphJnl::usage="divrTorSphJnl[n,l,r]";
divrPolSphJnl::usage="divrPolSphJnl[n,l,r]";
drTorSphJnl::usage="drTorSphJnl[n,l,r]";
drPolSphJnl::usage="drPolSphJnl[n,l,r]";
(*divrdrJnl operators*)
divrdrTorSphJnl::usage="divrdrTorSphJnl[n,l,r]";
divrdrTorSphJnlExplicit::usage="divrdrTorSphJnlExplicit[n,l,r]";
divrdrTorSphJnlImplicit::usage="divrdrTorSphJnlImplicit[n,l,r]";
divrdrPolSphJnl::usage="divrdrPolSphJnl[n,l,r]";
divrdrPolSphJnlExplicit::usage="divrdrPolSphJnlExplicit[n,l,r]";
divrdrPolSphJnlImplicit::usage="divrdrPolSphJnlImplicit[n,l,r]";


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


norm[k_,l_]:=Abs[SphericalBesselJ[l,k]/Sqrt[2]];


(*Jnl*)
Jnl[k_,l_,t_]:=SphericalBesselJ[l, k t]/norm[k,l]
(*rJnl*)
rJnl[k_,l_,t_]=Simplify[t Jnl[k,l,t]];
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


(*divrdrJnl operators*)
divrdrJnlExplicit[ks_,l_,r_,w_]:=Module[{opA,opB,opC,op},
opA=Table[Jnl[k,l,r],{k,ks}];
opB=Table[dJnl[k,l,r],{k,ks}];
opC=Table[Jnl[k,l,r],{k,ks[[;;-2]]}];
op=(opC . DiagonalMatrix[w]) . DiagonalMatrix[1/r] . Transpose[opB] . (opA . DiagonalMatrix[w]) . DiagonalMatrix[r];
Transpose[op]
]
divrdrJnlImplicit[ks_,ksm1_,l_,r_,w_]:=Module[{opA,opB,opC,op},
If[l==0,
op =divrdrJnlExplicit[ks,l,r,w];
,
opA=Table[Jnl[k,l-1,r],{k,ksm1}];
opB=Table[divrdrJnl[k,l-1,r],{k,ksm1}];
opC=Table[Jnl[k,l,r],{k,ks[[;;-2]]}];
op=Transpose[((opC . DiagonalMatrix[w]) . Transpose[opB] . (opA . DiagonalMatrix[w]))];
];
op
]


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
divrdrTorSphJnlExplicit[maxN_,l_,r_,w_]:=Module[{ks},
	ks = Table[torZero[n,l],{n,0,maxN+1}];
	divrdrJnlExplicit[ks,l,r,w]
	]
divrdrTorSphJnlImplicit[maxN_,l_,r_,w_]:=Module[{ks,ksm1},
	ks = Table[torZero[n,l],{n,0,maxN+1}];
	ksm1 = Table[torZero[n,l-1],{n,0,maxN+1}];
	divrdrJnlImplicit[ks,ksm1,l,r,w]
	]
divrdrPolSphJnlExplicit[maxN_,l_,r_,w_]:=Module[{ks},
	ks = Table[polZero[n,l],{n,0,maxN+1}];
	divrdrJnlExplicit[ks,l,r,w]
	]
divrdrPolSphJnlImplicit[maxN_,l_,r_,w_]:=Module[{ks,ksm1},
	ks = Table[polZero[n,l],{n,0,maxN+1}];
	ksm1 = Table[polZero[n,l-1],{n,0,maxN+1}];
	divrdrJnlImplicit[ks,ksm1,l,r,w]
	]


End[];


EndPackage[];

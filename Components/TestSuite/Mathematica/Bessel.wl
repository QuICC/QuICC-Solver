(* ::Package:: *)

BeginPackage["Bessel`",{"NumericalDifferentialEquationAnalysis`"}];


$mpprec::usage="Precision for Worland computations, default is 100";


bgrid::usage="wgrid[n] computes the Worland quadrature grid";
bweights::usage="wweights[n] computes the Worland quadrature weights";
norm::usage="norm[k,l,dv] norm of Bessel k";

$valueD\[Nu]::usage="nu = l + d\[Nu] for value boundary condition";
$insulatingD\[Nu]::usage="nu = l + d\[Nu] for insulating boundary condition";

(* Zeros *)
valueZero::usage="valueZero[n,l]";
insulatingZero::usage="insulatingZero[n,l]";
valueZeros::usage="valueZeros[n,l]";
insulatingZeros::usage="insulatingZeros[n,l]";

(* Operators to work on grid*)
ValueSphJnl::usage="ValueSphJnl[n,l,r]";
InsulatingSphJnl::usage="InsulatingSphJnl[n,l,r]";
ValueRSphJnl::usage="ValueRSphJnl[n,l,r]";
InsulatingRSphJnl::usage="InsulatingRSphJnl[n,l,r]";
ValueDSphJnl::usage="ValueDSphJnl[n,l,r]";
InsulatingDSphJnl::usage="InsulatingDSphJnl[n,l,r]";
ValueSlaplSphJnl::usage="ValueSlaplSphJnl[n,l,r]";
InsulatingSlaplSphJnl::usage="InsulatingSlaplSphJnl[n,l,r]";
ValueDivrSphJnl::usage="ValueDivrSphJnl[n,l,r]";
InsulatingDivrSphJnl::usage="InsulatingDivrSphJnl[n,l,r]";
ValueDrSphJnl::usage="ValueDrSphJnl[n,l,r]";
InsulatingDrSphJnl::usage="InsulatingDrSphJnl[n,l,r]";
(*divrdrJnl operators*)
ValueDivrdrSphJnl::usage="ValueDivrdrSphJnl[n,l,r]";
ValueDivrdrSphJnlExplicit::usage="ValueDivrdrSphJnlExplicit[n,l,r]";
ValueDivrdrSphJnlImplicit::usage="ValueDivrdrSphJnlImplicit[n,l,r]";
InsulatingDivrdrSphJnl::usage="InsulatingDivrdrSphJnl[n,l,r]";
InsulatingDivrdrSphJnlExplicit::usage="InsulatingDivrdrSphJnlExplicit[n,l,r]";
InsulatingDivrdrSphJnlImplicit::usage="InsulatingDivrdrSphJnlImplicit[n,l,r]";


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


$valueD\[Nu] =1/2;
$insulatingD\[Nu] =-(1/2);
(*Roots*)
getZero[n_,l_,d\[Nu]_]:=BesselJZero[l+d\[Nu],n+1]
valueZero[n_,l_]:=getZero[n,l,$valueD\[Nu]];
insulatingZero[n_,l_]:=getZero[n,l,$insulatingD\[Nu]];
valueZeros[n_,l_]:=Table[valueZero[i,l],{i,0,n}];
insulatingZeros[n_,l_]:=Table[insulatingZero[i,l] ,{i,0,n}];


norm[k_,l_,d\[Nu]_]:=Module[{},
	If[d\[Nu]==$valueD\[Nu],
		Abs[SphericalBesselJ[l+1,k]/Sqrt[2]],
		If[d\[Nu]==$insulatingD\[Nu],
			Abs[SphericalBesselJ[l,k]/Sqrt[2]],
			Abs[SphericalBesselJ[l,k]/Sqrt[2]]
			]
		]
]


(*Jnl*)
Jnl[k_,l_,t_,d\[Nu]_]:=SphericalBesselJ[l, k t]/norm[k,l,d\[Nu]]
(*rJnl*)
rJnl[k_,l_,t_,d\[Nu]_]=Simplify[t Jnl[k,l,t,d\[Nu]]];
(*dJnl*)
dJnl[k_,l_,t_,d\[Nu]_]=FullSimplify[D[Jnl[k,l,t,d\[Nu]],{t,1}]];
(*slaplJnl*)
slaplJnl[k_,l_,t_,d\[Nu]_]=Simplify[1/t^2 D[t^2 D[Jnl[k,l,t,d\[Nu]],t],t]-(l(l+1))/t^2 Jnl[k,l,t,d\[Nu]]];
(*divrJnl*)
divrJnl[k_,l_,t_,d\[Nu]_]=Simplify[1/t Jnl[k,l,t,d\[Nu]]];
(*drJnl*)
drJnl[k_,l_,t_,d\[Nu]_]=Simplify[D[t Jnl[k,l,t,d\[Nu]],t]];
(*divrdrJnl*)
divrdrJnl[k_,l_,t_,d\[Nu]_]=Simplify[1/t D[t Jnl[k,l,t,d\[Nu]],t]];


(*divrdrJnl operators*)
divrdrJnlExplicit[ks_,l_,r_,d\[Nu]_,w_]:=Module[{opA,opB,opC,op},
opA=Table[Jnl[k,l,r,d\[Nu]],{k,ks}];
opB=Table[dJnl[k,l,r,d\[Nu]],{k,ks}];
opC=Table[Jnl[k,l,r,d\[Nu]],{k,ks[[;;-2]]}];
op=(opC . DiagonalMatrix[w]) . DiagonalMatrix[1/r] . Transpose[opB] . (opA . DiagonalMatrix[w]) . DiagonalMatrix[r];
Transpose[op]
]
divrdrJnlImplicit[ks_,ksm1_,l_,r_,d\[Nu]_,w_]:=Module[{opA,opB,opC,op},
If[l==0,
op =divrdrJnlExplicit[ks,l,r,d\[Nu],w];
,
opA=Table[Jnl[k,l-1,r,d\[Nu]],{k,ksm1}];
opB=Table[divrdrJnl[k,l-1,r,d\[Nu]],{k,ksm1}];
opC=Table[Jnl[k,l,r,d\[Nu]],{k,ks[[;;-2]]}];
op=Transpose[((opC . DiagonalMatrix[w]) . Transpose[opB] . (opA . DiagonalMatrix[w]))];
];
op
]


(* Operators to work on grid*)
ValueSphJnl[n_,l_,r_]:=Jnl[valueZero[n,l],l,r,$valueD\[Nu]]
InsulatingSphJnl[n_,l_,r_]:=Jnl[insulatingZero[n,l],l,r,$insulatingD\[Nu]]
ValueRSphJnl[n_,l_,r_]:=rJnl[valueZero[n,l],l,r,$valueD\[Nu]]
InsulatingRSphJnl[n_,l_,r_]:=rJnl[insulatingZero[n,l],l,r,$insulatingD\[Nu]]
ValueDSphJnl[n_,l_,r_]:=dJnl[valueZero[n,l],l,r,$valueD\[Nu]]
InsulatingDSphJnl[n_,l_,r_]:=dJnl[insulatingZero[n,l],l,r,$insulatingD\[Nu]]
ValueSlaplSphJnl[n_,l_,r_]:=slaplJnl[valueZero[n,l],l,r,$valueD\[Nu]]
InsulatingSlaplSphJnl[n_,l_,r_]:=slaplJnl[insulatingZero[n,l],l,r,$insulatingD\[Nu]]
ValueDivrSphJnl[n_,l_,r_]:=divrJnl[valueZero[n,l],l,r,$valueD\[Nu]]
InsulatingDivrSphJnl[n_,l_,r_]:=divrJnl[insulatingZero[n,l],l,r,$insulatingD\[Nu]]
ValueDrSphJnl[n_,l_,r_]:=drJnl[valueZero[n,l],l,r,$valueD\[Nu]]
InsulatingDrSphJnl[n_,l_,r_]:=drJnl[insulatingZero[n,l],l,r,$insulatingD\[Nu]]
ValueDivrdrSphJnl[n_,l_,r_]:=divrdrJnl[valueZero[n,l],l,r,$valueD\[Nu]]
InsulatingDivrdrSphJnl[n_,l_,r_]:=divrdrJnl[insulatingZero[n,l],l,r,$insulatingD\[Nu]]
ValueDivrdrSphJnlExplicit[maxN_,l_,r_,w_]:=Module[{ks},
	ks = Table[valueZero[n,l],{n,0,maxN+1}];
	divrdrJnlExplicit[ks,l,r,$valueD\[Nu],w]
	]
InsulatingDivrdrSphJnlExplicit[maxN_,l_,r_,w_]:=Module[{ks},
	ks = Table[insulatingZero[n,l],{n,0,maxN+1}];
	divrdrJnlExplicit[ks,l,r,$insulatingD\[Nu],w]
	]
ValueDivrdrSphJnlImplicit[maxN_,l_,r_,w_]:=Module[{ks,ksm1},
	ks = Table[valueZero[n,l],{n,0,maxN+1}];
	ksm1 = Table[valueZero[n,l-1],{n,0,maxN+1}];
	divrdrJnlImplicit[ks,ksm1,l,r,$valueD\[Nu],w]
	]
InsulatingDivrdrSphJnlImplicit[maxN_,l_,r_,w_]:=Module[{ks,ksm1},
	ks = Table[insulatingZero[n,l],{n,0,maxN+1}];
	ksm1 = Table[insulatingZero[n,l-1],{n,0,maxN+1}];
	divrdrJnlImplicit[ks,ksm1,l,r,$insulatingD\[Nu],w]
	]


End[];


EndPackage[];

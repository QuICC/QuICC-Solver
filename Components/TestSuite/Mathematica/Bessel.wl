(* ::Package:: *)

BeginPackage["Bessel`",{"NumericalDifferentialEquationAnalysis`"}];


$mpprec::usage="Precision for Worland computations, default is 100";


bgrid::usage="wgrid[n] computes the Worland quadrature grid";
bweights::usage="wweights[n] computes the Worland quadrature weights";
bnorm::usage="bnorm[k,l,dv] norm of Bessel k";

$valueD\[Nu]::usage="nu = l + d\[Nu] for value boundary condition";
$insulatingD\[Nu]::usage="nu = l + d\[Nu] for insulating boundary condition";
$nsD\[Nu]::usage="nu = l + d\[Nu] for no-slip boundary condition";

(* Zeros *)
getZero::usage="getZero[n,l,d\[Nu]]"
valueZero::usage="valueZero[n,l]";
insulatingZero::usage="insulatingZero[n,l]";
nsZero::usage="nsZero[n,l]";
valueZeros::usage="valueZeros[n,l]";
insulatingZeros::usage="insulatingZeros[n,l]";
nsZeros::usage="nsZeros[n,l]";
sfZero::usage="sfZero[n,l]";
sfZeros::usage="sfZeros[n,l]";

(* Operators to work on grid*)
Jnl::usage="Jnl[k,l,t,d\[Nu]]";
lowerJnl::usage="lowerJnl[k,l,t,d\[Nu]]";
raiseJnl::usage="raiseJnl[k,l,t,d\[Nu]]";
ValueSphJnl::usage="ValueSphJnl[n,l,r]";
InsulatingSphJnl::usage="InsulatingSphJnl[n,l,r]";
StressFreeSphJnl::usage="StressFreeSphJnl[n,l,r]";
NoSlipSphJnl::usage="NoSlipSphJnl[n,l,r]";
ValueRSphJnl::usage="ValueRSphJnl[n,l,r]";
InsulatingRSphJnl::usage="InsulatingRSphJnl[n,l,r]";
NoSlipRSphJnl::usage="NoSlipRSphJnl[n,l,r]";
ValueDSphJnl::usage="ValueDSphJnl[n,l,r]";
InsulatingDSphJnl::usage="InsulatingDSphJnl[n,l,r]";
NoSlipDSphJnl::usage="NoSlipDSphJnl[n,l,r]";
ValueD2SphJnl::usage="ValueD2SphJnl[n,l,r]";
InsulatingD2SphJnl::usage="InsulatingD2SphJnl[n,l,r]";
NoSlipD2SphJnl::usage="NoSlipD2SphJnl[n,l,r]";
ValueSlaplSphJnl::usage="ValueSlaplSphJnl[n,l,r]";
InsulatingSlaplSphJnl::usage="InsulatingSlaplSphJnl[n,l,r]";
NoSlipSlaplSphJnl::usage="NoSlipSlaplSphJnl[n,l,r]";
ValueDivrSphJnl::usage="ValueDivrSphJnl[n,l,r]";
InsulatingDivrSphJnl::usage="InsulatingDivrSphJnl[n,l,r]";
NoSlipDivrSphJnl::usage="NoSlipDivrSphJnl[n,l,r]";
ValueDrSphJnl::usage="ValueDrSphJnl[n,l,r]";
InsulatingDrSphJnl::usage="InsulatingDrSphJnl[n,l,r]";
NoSlipDrSphJnl::usage="NoSlipDrSphJnl[n,l,r]";
ValueRaiseSphJnl::usage="ValueRaiseSphJnl[n,l,r]";
InsulatingRaiseSphJnl::usage="InsulatingRaiseSphJnl[n,l,r]";
ValueLowerSphJnl::usage="ValueLowerSphJnl[n,l,r]";
InsulatingLowerSphJnl::usage="InsulatingLowerSphJnl[n,l,r]";
(*divrdrJnl operators*)
ValueDivrdrSphJnl::usage="ValueDivrdrSphJnl[n,l,r]";
ValueDivrdrSphJnlExplicit::usage="ValueDivrdrSphJnlExplicit[n,l,r]";
ValueDivrdrSphJnlImplicit::usage="ValueDivrdrSphJnlImplicit[n,l,r]";
InsulatingDivrdrSphJnl::usage="InsulatingDivrdrSphJnl[n,l,r]";
InsulatingDivrdrSphJnlExplicit::usage="InsulatingDivrdrSphJnlExplicit[n,l,r]";
InsulatingDivrdrSphJnlImplicit::usage="InsulatingDivrdrSphJnlImplicit[n,l,r]";
NoSlipDivrdrSphJnl::usage="NoSlipDivrdrSphJnl[n,l,r]";
NoSlipDivrdrSphJnlExplicit::usage="NoSlipDivrdrSphJnlExplicit[n,l,r]";
NoSlipDivrdrSphJnlImplicit::usage="NoSlipDivrdrSphJnlImplicit[n,l,r]";
(*Boundary operators*)
ValueRDDivrSphJnl::usage="ValueRDDivrSphJnl[n,l,r]";
InsulatingRDDivrSphJnl::usage="ValueRDDivrSphJnl[n,l,r]";
ValueInsulatingSphereSphJnl::usage="ValueInsulatinSphereSphJnl[n,l,r]";
InsulatingInsulatingSphereSphJnl::usage="ValueInsulatingSphereSphJnl[n,l,r]";
rddivrJnl::usage="rddivrJnl[k,l,t,d\[Nu]]";


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
$nsD\[Nu]=3/2;
(*Roots*)
getZero[n_,l_,d\[Nu]_]:=BesselJZero[l+d\[Nu],n+1]
valueZero[n_,l_]:=getZero[n,l,$valueD\[Nu]];
insulatingZero[n_,l_]:=getZero[n,l,$insulatingD\[Nu]];
nsZero[n_,l_]:=If[n==0,0,getZero[n-1,l,$nsD\[Nu]]];
valueZeros[n_,l_]:=Table[valueZero[i,l],{i,0,n}];
insulatingZeros[n_,l_]:=Table[insulatingZero[i,l] ,{i,0,n}];
nsZeros[n_,l_]:=Table[nsZero[i,l] ,{i,0,n}];
(* Stress-free Roots *)
scanRoot[l_,z0_,iN_,dk_,zf_]:=Module[{k0,k1,y0,y1,i,z=-1},
For[i=0,i<iN,i++,
k0=z0+dk/iN i;
y0=zf[l,k0];
k1=z0+dk/iN (i+1);
y1=zf[l,k1];
If[y0 y1 <0,z=k1;Break[];]
];
z
]
sfZeros[n_,l_]:=Module[{dk,iN,z1,z2,zs,zf,fp=32,wp=100},
zf[p_,k_]=(D[SphericalBesselJ[p,k r],r]-1/r SphericalBesselJ[p,k])/.{r->1};
(* Scan for starting roots *)
dk=Max[2l,10];
iN=2dk;
z1=scanRoot[l,2,iN,dk,zf];
If[z1>0,
z2=scanRoot[l,z1,iN,dk,zf];
If[z2>0,
zs={k/.FindRoot[zf[l,k],{k,z1},WorkingPrecision->wp,PrecisionGoal->fp],k/.FindRoot[zf[l,k],{k,z2},WorkingPrecision->wp,PrecisionGoal->fp]};
For[i=0,i<n-2,i++,
dz=zs[[-1]]-zs[[-2]];
zn=k/.FindRoot[zf[l,k],{k,zs[[-1]]+dz},WorkingPrecision->wp,PrecisionGoal->fp];
AppendTo[zs,zn];
];
,Abort[];];
,Abort[];];
zs];
sfZero[n_,l_]:=sfZeros[n+1,l][[n+1]];


bnorm[k_,l_,d\[Nu]_]:=Module[{},
	If[d\[Nu]==$valueD\[Nu],
		Abs[SphericalBesselJ[l+1,k]/Sqrt[2]],
		If[d\[Nu]==$insulatingD\[Nu],
			Abs[SphericalBesselJ[l,k]/Sqrt[2]],
			If[d\[Nu]==$nsD\[Nu],
				Abs[SphericalBesselJ[l+2,k]/Sqrt[2]],
				1]
			]
		]
]


(*Jnl*)
Jnl[k_,l_,t_,d\[Nu]_]:=If[k==0,Sqrt[3+2l]t^l,SphericalBesselJ[l, k t]/bnorm[k,l,d\[Nu]]]
(*rJnl*)
rJnl[k_,l_,t_,d\[Nu]_]=Simplify[t Jnl[k,l,t,d\[Nu]]];
(*dJnl*)
dJnl[k_,l_,t_,d\[Nu]_]=FullSimplify[D[Jnl[k,l,t,d\[Nu]],{t,1}]];
(*d2Jnl*)
d2Jnl[k_,l_,t_,d\[Nu]_]=FullSimplify[D[Jnl[k,l,t,d\[Nu]],{t,2}]];
(*slaplJnl*)
slaplJnl[k_,l_,t_,d\[Nu]_]=Simplify[1/t^2 D[t^2 D[Jnl[k,l,t,d\[Nu]],t],t]-(l(l+1))/t^2 Jnl[k,l,t,d\[Nu]]];
(*divrJnl*)
divrJnl[k_,l_,t_,d\[Nu]_]=Simplify[1/t Jnl[k,l,t,d\[Nu]]];
(*drJnl*)
drJnl[k_,l_,t_,d\[Nu]_]=Simplify[D[t Jnl[k,l,t,d\[Nu]],t]];
(*divrdrJnl*)
divrdrJnl[k_,l_,t_,d\[Nu]_]=Simplify[1/t D[t Jnl[k,l,t,d\[Nu]],t]];
(*raiseJnl*)
raiseJnl[k_,l_,t_,d\[Nu]_]:=Simplify[k SphericalBesselJ[l+1, k t]/bnorm[k,l,d\[Nu]]];
(*lowerJnl*)
lowerJnl[k_,l_,t_,d\[Nu]_]:=Simplify[k SphericalBesselJ[l-1, k t]/bnorm[k,l,d\[Nu]]];
(* Stres-free toroidal *)
rddivrJnl[k_,l_,t_,d\[Nu]_]:=dJnl[k,l,t,d\[Nu]]- 1/t Jnl[k,l,t,d\[Nu]]


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
StressFreeSphJnl[n_,l_,r_]:=Jnl[sfZero[n,l],l,r,-42]
NoSlipSphJnl[n_,l_,r_]:=Jnl[nsZero[n,l],l,r,$nsD\[Nu]]

ValueRSphJnl[n_,l_,r_]:=rJnl[valueZero[n,l],l,r,$valueD\[Nu]]
InsulatingRSphJnl[n_,l_,r_]:=rJnl[insulatingZero[n,l],l,r,$insulatingD\[Nu]]
NoSlipRSphJnl[n_,l_,r_]:=rJnl[nsZero[n,l],l,r,$nsD\[Nu]]

ValueDSphJnl[n_,l_,r_]:=dJnl[valueZero[n,l],l,r,$valueD\[Nu]]
InsulatingDSphJnl[n_,l_,r_]:=dJnl[insulatingZero[n,l],l,r,$insulatingD\[Nu]]
NoSlipDSphJnl[n_,l_,r_]:=dJnl[nsZero[n,l],l,r,$nsD\[Nu]]

ValueD2SphJnl[n_,l_,r_]:=d2Jnl[valueZero[n,l],l,r,$valueD\[Nu]]
InsulatingD2SphJnl[n_,l_,r_]:=d2Jnl[insulatingZero[n,l],l,r,$insulatingD\[Nu]]
NoSlipD2SphJnl[n_,l_,r_]:=d2Jnl[nsZero[n,l],l,r,$nsD\[Nu]]

ValueSlaplSphJnl[n_,l_,r_]:=slaplJnl[valueZero[n,l],l,r,$valueD\[Nu]]
InsulatingSlaplSphJnl[n_,l_,r_]:=slaplJnl[insulatingZero[n,l],l,r,$insulatingD\[Nu]]
NoSlipSlaplSphJnl[n_,l_,r_]:=slaplJnl[nsZero[n,l],l,r,$nsD\[Nu]]

ValueDivrSphJnl[n_,l_,r_]:=divrJnl[valueZero[n,l],l,r,$valueD\[Nu]]
InsulatingDivrSphJnl[n_,l_,r_]:=divrJnl[insulatingZero[n,l],l,r,$insulatingD\[Nu]]
NoSlipDivrSphJnl[n_,l_,r_]:=divrJnl[nsZero[n,l],l,r,$nsD\[Nu]]

ValueDrSphJnl[n_,l_,r_]:=drJnl[valueZero[n,l],l,r,$valueD\[Nu]]
InsulatingDrSphJnl[n_,l_,r_]:=drJnl[insulatingZero[n,l],l,r,$insulatingD\[Nu]]
NoSlipDrSphJnl[n_,l_,r_]:=drJnl[nsZero[n,l],l,r,$nsD\[Nu]]

ValueDivrdrSphJnl[n_,l_,r_]:=divrdrJnl[valueZero[n,l],l,r,$valueD\[Nu]]
InsulatingDivrdrSphJnl[n_,l_,r_]:=divrdrJnl[insulatingZero[n,l],l,r,$insulatingD\[Nu]]
NoSlipDivrdrSphJnl[n_,l_,r_]:=divrdrJnl[nsZero[n,l],l,r,$nsD\[Nu]]

ValueRaiseSphJnl[n_,l_,r_]:=raiseJnl[valueZero[n,l],l,r,$valueD\[Nu]]
InsulatingRaiseSphJnl[n_,l_,r_]:=raiseJnl[insulatingZero[n,l],l,r,$insulatingD\[Nu]]

ValueLowerSphJnl[n_,l_,r_]:=lowerJnl[valueZero[n,l],l,r,$valueD\[Nu]]
InsulatingLowerSphJnl[n_,l_,r_]:=lowerJnl[insulatingZero[n,l],l,r,$insulatingD\[Nu]]

ValueDivrdrSphJnlExplicit[maxN_,l_,r_,w_]:=Module[{ks},
	ks = Table[valueZero[n,l],{n,0,maxN+1}];
	divrdrJnlExplicit[ks,l,r,$valueD\[Nu],w]
	]
InsulatingDivrdrSphJnlExplicit[maxN_,l_,r_,w_]:=Module[{ks},
	ks = Table[insulatingZero[n,l],{n,0,maxN+1}];
	divrdrJnlExplicit[ks,l,r,$insulatingD\[Nu],w]
	]
NoSlipDivrdrSphJnlExplicit[maxN_,l_,r_,w_]:=Module[{ks},
	ks = Table[nsZero[n,l],{n,0,maxN+1}];
	divrdrJnlExplicit[ks,l,r,$nsD\[Nu],w]
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
NoSlipDivrdrSphJnlImplicit[maxN_,l_,r_,w_]:=Module[{ks,ksm1},
	ks = Table[nsZero[n,l],{n,0,maxN+1}];
	ksm1 = Table[nsZero[n,l-1],{n,0,maxN+1}];
	divrdrJnlImplicit[ks,ksm1,l,r,$nsD\[Nu],w]
	]


ValueInsulatingSphereSphJnl[n_,l_,r_]:=ValueDSphJnl[n,l,r]+(l+1)/r ValueSphJnl[n,l,r]
InsulatingInsulatingSphereSphJnl[n_,l_,r_]:=InsulatingDSphJnl[n,l,r]+(l+1)/r InsulatingSphJnl[n,l,r]
ValueRDDivrSphJnl[n_,l_,r_]:=ValueDSphJnl[n,l,r]-1/r ValueSphJnl[n,l,r]
InsulatingRDDivrSphJnl[n_,l_,r_]:=InsulatingDSphJnl[n,l,r]-1/r InsulatingSphJnl[n,l,r]


End[];


EndPackage[];

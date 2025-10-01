(* ::Package:: *)

BeginPackage["ALegendre`",{"NumericalDifferentialEquationAnalysis`"}];


$mpprec::usage="Precision for ALegendre computations, default is 100";

(*Generic interface*)
lgrid::usage="lgrid[n] computes the ALegendre quadrature grid";
lweights::usage="lweights[n] computes the ALegendre quadrature weights";
lnorm::usage="lnorm[k,\[Alpha],\[Beta]] norm of ALegendre polynomial k";

(* Operators to work on grid*)
Plm::usage="Plm[l,m,t]";
dPlm::usage="dPlm[l,m,t]";
divsinPlm::usage="divsinPlm[l,m,r]";


Begin["`Private`"];


lgrid[n_]:=If[ToExpression["ValueQ[lgrid"<>ToString[n]<>"]"],ToExpression["lgrid"<>ToString[n]],ToExpression["Set[lgrid"<>ToString[n]<>",GaussianQuadratureWeights["<>ToString[n]<>",-1,1,mpprec][[;;,1]]]"]];
lweights[n_]:=If[ToExpression["ValueQ[lweights"<>ToString[n]<>"]"],ToExpression["lweights"<>ToString[n]],ToExpression["Set[lweights"<>ToString[n]<>",GaussianQuadratureWeights["<>ToString[n]<>",-1,1,mpprec][[;;,2]]]"]];


lnorm[l_,m_]:=Sqrt[2\[Pi] (2(l+m)!)/((2l+1)(l-m)!)]


dLegendreP[l_,m_,t_]=D[LegendreP[l,m,t],t];


Plm[l_,m_,t_]:=LegendreP[l,m,t]/lnorm[l,m]


dPlm[l_,m_,t_]:=-Sqrt[1-t^2]dLegendreP[l,m,t]/lnorm[l,m]
divsinPlm[l_,m_,t_]:=If[m==0,0t,1/Sqrt[1-t^2] LegendreP[l,m,t]/lnorm[l,m]]


End[];


EndPackage[];

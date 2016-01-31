(* Mathematica package *)

BeginPackage["QuantumMechanics`Derivative`"]
(* Exported symbols added here with SymbolName::usage *) 


d::usage :=
"Used like the built-in function D, 
but implements a better solution for embbeded Conjugate calls,
provided the function depends on real parameters";


Begin["`Private`"]
(* Implementation of the package *)


(* Substitute of the built in function D *)
d[ex_,vars__]:=
	Module[{res,lis={},prepareVars},
		SetAttributes[expSubstitute,HoldRest];
		SetAttributes[replaceConjugate,HoldRest];
		prepareVars[vars2_]:=Module[{},
				If[Depth[vars2]<=2,
					vars2,
					Map[#[[1]]&,vars2]
				]
			];
		res=expSubstitute[ex,{Sequence[prepareVars[{vars}]]},replaceConjugate,lis];(* transformation *)
		res=D[res,Sequence[vars]]; (* perform the derivation *)
		res=expSubstitute[res,{vars},replacePlaceHolder,lis]; (* back transformation *)
		res
	]


(* Common Substitution Function *)
expSubstitute[exp_,{vars__},getReplacement_:Function[{ex,lis,vars},ex],lis_,level_:0]:=
	Module[{return,i},
		If[
			AtomQ[exp],
			Return[exp];
		];
		return=Level[exp,1,Heads->True];
		For[
			i=2,
			i<=Length[return],
			i=i+1,
			If[!AtomQ[return[[i]]],
				return[[i]]=expSubstitute[return[[i]],{vars},getReplacement,lis,level+1];
			];
		];
		return=getReplacement[return,lis,{vars}];
		return=Apply[First[return],Drop[return,1]];
		return
	]


(* Replacement of Conjugate with Placeholder *)
replaceConjugate[x_,lis_,{vars__}]:=
	Module[{},
		If[SameQ[x[[1]],Conjugate],
			AppendTo[lis,x[[2]]];
			{\[Psi][Length[lis]],vars},
			x
		]
	]


(* Replacement of Placeholder back to Original Expression *)
replacePlaceHolder[x_,lis_,{vars__}]:=
	Module[{head,n1,m1,pat1=Derivative[n__][\[Psi][m_]],pat2=\[Psi][m_]},
		head=x[[1]];
		Which[Count[{head},pat1]>0,
			n1=head/.pat1->n;
			m1=head/.pat1->m;
			{Conjugate,D[lis[[m1]],vars]},
			Count[{head},pat2]>0,
			m1=head/.pat2->m;
			{Conjugate,lis[[m1]]},
			True,
			x
		]
	]


End[]

EndPackage[]

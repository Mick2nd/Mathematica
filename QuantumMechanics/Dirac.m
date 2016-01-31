(* Mathematica package *)

BeginPackage[
	"QuantumMechanics`Dirac`",
	{
		"QuantumMechanics`Physical`"
	}
	]
(* Exported symbols added here with SymbolName::usage *) 

spin::usage :=
"Dirac Matrix for x,y,z direction, use Subscript notation";

spinRot::usage :=
"Dirac Matrix for rotated direction";

\[Theta]::usage :=
"Angle between z axis and direction vector";

\[Phi]::usage :=
"Angle between x axis and projection of the direction vector";



Begin["`Private`"]
(* Implementation of the package *)

Subscript[spin,3]=\[HBar]/2 {{1,0},{0,-1}};

Subscript[spin,1]=\[HBar]/2 {{0,1},{1,0}};

Subscript[spin,2]=\[HBar]/2 {{0,-I},{I,0}};

spinRot=Sum[Subscript[spin,i] Subscript[n,i],{i,3}]//FullSimplify


Subscript[n,1]=Sin[\[Theta]] Cos[\[Phi]];

Subscript[n,2]=Sin[\[Theta]] Sin[\[Phi]];

Subscript[n,3]=Cos[\[Theta]];


End[]

EndPackage[]

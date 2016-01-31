(* Mathematica package *)

BeginPackage[
	"QuantumMechanics`Schroedinger`",
	{
		"QuantumMechanics`Physical`"
	}
	]
(* Exported symbols added here with SymbolName::usage *) 

Schroedinger::usage :=
"Gives an expression which, set to Zero, equals the Schroedinger equation of wave mechanics";

Begin["`Private`"]
(* Implementation of the package *)

Schroedinger[\[Psi]_,x_List,v_,chart_:"Cartesian"]:=I \[HBar] D[\[Psi][x,t],t]-(-\[HBar]^2/2/m Laplacian[\[Psi][x,t],x]+v[x] \[Psi][x,t])

End[]

EndPackage[]

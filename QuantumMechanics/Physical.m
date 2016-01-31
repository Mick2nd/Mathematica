(* Mathematica package *)


BeginPackage["QuantumMechanics`Physical`"]
(* Exported symbols added here with SymbolName::usage *)

\[HBar]::usage :=
"Planck Constant";

m::usage :=
"Mass";

t::usage :=
"Time variable";

Begin["`Private`"]
(* Implementation of the package *)

dummy := \[HBar] m t

End[]

EndPackage[]

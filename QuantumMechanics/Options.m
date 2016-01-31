(* Mathematica package *)

BeginPackage["QuantumMechanics`Options`"]
(* Exported symbols added here with SymbolName::usage *) 


conjugateTranspose::usage :=
"Like the built-in ConjugateTranspose, but with additional Assumptions option";

hermitianMatrixQ::usage :=
"Like the built-in HermitianMatrixQ, but with additional Assumptions option";

Begin["`Private`"]
(* Implementation of the package *)

Options[conjugateTranspose]={Assumptions->True}

conjugateTranspose[op_List,opt:OptionsPattern[]]:=
	Map[Simplify[#,Assumptions->OptionValue[Assumptions]]&,ConjugateTranspose[op]]

Options[hermitianMatrixQ]={Assumptions->True};

hermitianMatrixQ[op_List,opt:OptionsPattern[]]:=
	Module[{ct=conjugateTranspose[op,Assumptions->OptionValue[Assumptions]]},
		MatrixQ[op]&& ct===op
	]

End[]

EndPackage[]

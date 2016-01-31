(* Mathematica Package *)

(* Created by the Wolfram Workbench 26.07.2015 *)

BeginPackage[
	"QuantumMechanics`Quantum`", 
	{
		"QuantumMechanics`Derivative`",
		"QuantumMechanics`Options`",
		"QuantumMechanics`Dirac`",
		"QuantumMechanics`Schroedinger`",
		"QuantumMechanics`Physical`"
	}]
(* Exported symbols added here with SymbolName::usage *) 


eigenValues::usage :=
"Calculates the Eigen Values of an Operator given in matrix form";

eigenKets::usage :=
"Calculates the Eigen Kets of an Operator given in matrix form";

eigenBras::usage :=
"Calculates the Eigen Bras of an Operator given in matrix form";

braOf::usage :=
"Calculates the Eigen Bra of a Ket given as vector";

expectationValue::usage :=
"Calculates the Expectation Value of an Operator related to a given State Ket";

dispersion::usage :=
"Calculates the Dispersion of an Operator related to a given State Ket";

commutator::usage :=
"Calculates the Commutator of 2 Operators";

anticommutator::usage :=
"Calculates the Anti Commutator of 2 Operators";

uncertaintyRelationLeft::usage :=
"Calculates the Left Hand Side of the Uncertainty Relation";

uncertaintyRelationRight::usage :=
"Calculates the Right Hand Side of the Uncertainty Relation";

proofUncertaintyRelation::usage :=
"Checks the Uncertainty Relation for 2 given Operators";

uncertaintyProduct::usage :=
"Calculates the Uncertainty Product of x and p operator related to a given Wave Function";


Begin["`Private`"]
(* Implementation of the package *)

(* Calculates the Eigen Values of an Operator given in matrix form *)
eigenValues[op_List]:=Eigenvalues[op]

(* Calculates the Eigen Kets of an Operator given in matrix form *)
eigenKets[op_List]:=FullSimplify/@Normalize/@Eigenvectors[op]

(* Calculates the Eigen Bras of an Operator given in matrix form *)
eigenBras[op_List]:=Conjugate[eigenKets[op]];

(* Calculates the Eigen Bra of a Ket given as vector  *)
braOf[ket_List]:=Conjugate[ket];

(* Calculates the Expectation Value of an Operator related to a given State Ket *)
expectationValue[op_List,ket_List]:=
	Module[{bra=Conjugate[ket]},
		bra.(op.ket)
	]

(* Calculates the Dispersion of an Operator related to a given State Ket *)
dispersion[op_List,ket_List]:=
	Module[{exp=expectationValue[op,ket],id=IdentityMatrix[Length[op]]},
		expectationValue[op.op-exp exp id,ket]
	]

(* Calculates the Expectation Value of an Operator related to a given Wave Function *)
expectationValue[op_,\[Phi]_,area_List]:=
	Integrate[Conjugate[\[Phi]] op[\[Phi]],area];

(* Calculates the Dispersion of an Operator related to a given Wave Function *)
dispersion[op_,\[Phi]_,area_List]:=
	Module[{ev},
		ev=expectationValue[op,\[Phi],area];
		expectationValue[Nest[op,#,2]&,\[Phi],area]-ev^2
	]


(* Calculates the Commutator of 2 Operators *)
commutator[op1_List,op2_List]:=op1.op2-op2.op1

(* Calculates the Anti Commutator of 2 Operators *)
anticommutator[op1_List,op2_List]:=op1.op2+op2.op1

(* Calculates the Left Hand Side of the Uncertainty Relation *)
uncertaintyRelationLeft[op1_List,op2_List,ket_List]:=dispersion[op1,ket] dispersion[op2,ket]

(* Calculates the Right Hand Side of the Uncertainty Relation *)
uncertaintyRelationRight[op1_List,op2_List,ket_List]:=1/4 Abs[expectationValue[commutator[op1,op2],ket]]^2

(* Checks the Uncertainty Relation for 2 given Operators *)
proofUncertaintyRelation[op1_List,op2_List,ket_List]:=uncertaintyRelationLeft[op1,op2,ket]>=uncertaintyRelationRight[op1,op2,ket]

(* Calculates the Uncertainty Product of x and p operator related to a given Wave Function *)
uncertaintyProduct[\[Phi]_,area_List]:=
	Module[{coord},
		coord=area[[1]];
		dispersion[coord #&,\[Phi],area] dispersion[-I \[HBar] D[#,coord]&,\[Phi],area]
	]


End[]

EndPackage[]

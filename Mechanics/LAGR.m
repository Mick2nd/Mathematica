(* ::Package:: *)

(*
    Package LAGR for calculation of LAGRANGE equations  
    This is a port of the Diploma work from 1978  
    Usage:  
    - Define all functions  
    - Call chooseLocations for input and output  
    - Choose the locations by clicking the buttons  
    - Define input data  
    - Invoke lagr to read, calculate, write  
    - Read the output  
    
    Differences to the Diploma work  
    - use PC instead of Mainframe IBM  
    - use Mathematica CAS instead of PL1/FORMAC on Job Basis  
    - work interactively on PC  
    - Calculate one system per invocation instead of N  
    - Work with any count of files  
*)

BeginPackage["Mechanics`LAGR`"]


chooseLocations::usage:=
"Used to locate the input / output files.
First call this function, then use the buttons to invoke the File Open / File Save dialogs.";

lagr::usage:=
"Reads input, calculates results and writes output.
Create input file, then call lagr[], then read output.
For the format of input files look at a sample file.";

q::usage:=
	"Generalized coordinates, can be used in input.";

En::usage:=
	"Rotation Matrix, can be used in input, f.i. to specify forces";

xn::usage:=
	"Distance Vector, can be used in input, f.i. to specify forces";

m::usage:=
	"Mass. Can be substituted in output";

\[Theta]::usage:=
	"Inertia Tensor. Can be substituted in output";

c::usage:=
	"Spring constants. Can be substituted in output";

lambda::usage:=
	"Spring Lengths. Can be substituted in output";

xs::usage:=
	"Spring coordinates. Can be substituted in output";

k::usage:=
	"Damper constants. Can be substituted in output";
	
xd::usage:=
	"Damper coordinates. Can be substituted in output";

qd::usage:=
	"Generalized velocity. Generalized forces of Dampers depend on it";

	
Begin["`Private`"]


(*
	The function used to select input and output file
	- invoke the function
	- press the buttons to show the File Open / File Save dialog
	- select the input / output files
*)
  chooseLocations[]:=Module[{},
  inputFn = "Programmieren/Mathematica/IO/LagrInput-Diplom.txt";
  outputFn = "Programmieren/Mathematica/IO/LagrOutput-Diplom.txt";
  Print[{FileNameSetter[Dynamic[inputFn], "Open"], Dynamic[inputFn]}];
  Print[{FileNameSetter[Dynamic[outputFn], "Save"], Dynamic[outputFn]}];
  ]


(*
	The Main calculation function
	- reads the input as selected by chooseLocations
	- performs calculations
	- writes the output data either to the screen and the output file
*)
  lagr[] := Module[{},
  readLAGR[];
  calcLAGR[];
  writeLAGRAll[];
  ]


(*
	Reads the input file and stores the data
*)
  readLAGR[] := Module[{stm, l,k, Rlk,zlk, ck, cl, gF},
  If[IntegerQ[ degreesOfFreedom], init[degreesOfFreedom],l=0];
  stm = OpenRead[inputFn];
  numBodies = read[stm];
  degreesOfFreedom = read[stm];
  relativeMovement = read[stm];
  printAssignment["Number of Bodies", {}, numBodies];
  printAssignment["Degrees of Freedom", {}, degreesOfFreedom];
  printAssignment["Relative Movement", {}, relativeMovement];
  En[0] = read[stm];
  xn[0] = read[stm];
  For[l = 1, l <= numBodies, l++,
  If[relativeMovement,
  (k = read[stm];
  Rlk = read[stm];
  zlk = read[stm];
  ck =read[stm];
  cl =read[stm];
  En[l]= Rlk . En[k];
  xn[l] = xn[k] + (ck + zlk) . En[k] - cl.En[l];),
  (En[l]= read[stm];
  xn[l] = read[stm];)
  ];
  printAssignment["Rotation Matrix ", {l}, En[l]];
  printAssignment["Distance Vector ", {l}, xn[l]];
  \[Theta][l] = getInertiaTensor[l];
  gF = read[stm];
  If[gF,
  K[l]= read[stm];
  M[l]=  readMomentTensor[stm];,
  K[l]= {0., 0., 0.};
  M[l]=  Table[0.,{3},{3}];
  ];
  printAssignment["Force Vector ", {l}, K[l]];
  printAssignment["Moment Tensor ", {l}, M[l]];
  ]
  readSprings[stm];
  readDampers[stm];
  readCases[stm];
  Close[stm];
  ]


(*
	Initializes the data structures
	Mainly unsets defined symbols as set by a previous session
*)
  init[f_Integer]:=Module[{p, i, j, k},
  p = Flatten[Table[{i, j},{i,f},{j,f}], 1];
  Unset[g[#[[1]],#[[2]]]]& /@ p;
  p =Flatten[Table[{i, j, k},{i,f},{j,0,f},{k,0,j}], 2];
  Unset[\[CapitalGamma][#[[1]],#[[2]],#[[3]]]]& /@ p;
  p = Table[{i},{i,f}];
  Unset[Q[#[[1]]]]& /@ p;
  Unset[Qd[#[[1]]]]& /@ p;
  ]


(*
	Reads a 3x3 matrix defined by 9 single elements
	No longer used
*)
  readMatrix[stm_] := Module[{E},
  E = Table[read[stm],{3},{3}];
  Return [E];
  ]
  

(*
	Reads a 3 vector defined by 3 single elements
	No longer used as this can be done by a single read from a list
*)
  readVector[stm_] := Module[{x},
  x = Table[read[stm],{3}];
  Return[x];
  ]
  

(*
	Reads the Momentum Tensor
	Reads the Momentum Vector and converts it to the Momentum Tensor
*)
  readMomentTensor[stm_] :=Module[{x},
  x = read[stm];
  Return[{
  {0., -x[[3]], x[[2]]},
  {x[[3]], 0., -x[[1]]},
  {-x[[2]], x[[1]], 0.}
  }];
  ]
  

(*
	Returns the Inertia Tensor
	The Tensor (3x3 Matrix) consists of 9 symmetric elements of the form \[Theta][l, i, j]
	These symbols represent the Inertia properties of a body related to rotation
*)
  getInertiaTensor[l_Integer] := Module[{i,j},    (* symmetric tensor !! *)
  Return[Table[If[i<j,\[Theta][l, i, j],\[Theta][l, j, i]], {i,3},{j,3}]];
  ]


(*
	Reads the Springs definitions
	A Spring may be placed at 2 fixed points of 2 bodies
	The Spring force is proportional to the distance of the 2 points
	The generic spring point coordinates and spring constants may be replaced as defined by a substitution list
*)
  readSprings[stm_] := Module[{no},
  no = read[stm];
  printAssignment["Number of Springs",{}, no];
  Springs = Table[read[stm],{no},{4}];
  printAssignment["Springs' data",{}, Springs];
  ]
  

(*
	Reads the Dampers definitions
	A Damper may be placed at 2 fixed points of 2 bodies
	The Damper force is proportional to the relative velocity of the 2 points
	The generic damper point coordinates and damper constants may be replaced as defined by a substitution list
*)
  readDampers[stm_] := Module[{no},
  no=read[stm];
  printAssignment["Number of Dampers",{}, no];
  Dampers = Table[read[stm],{no},{3}];
  printAssignment["Dampers' data",{}, Dampers];
  ]

(*
	Reads the parameter specification cases of the given system
	This is a list of parameter substitutions
	An empty list is prepended to handle the common case
*)
  readCases[stm_] := Module[{no},
  no = read[stm];
  printAssignment["Number of Special Cases", {}, no];
  SystemCases = Table[read[stm],{no}];
  SystemCases = Join[{{}}, SystemCases];
  printAssignment["System Cases", {}, SystemCases];
  ]
  

(*
	Reads an item from the input stream
	overreads empty items
*)
  read[stm_] := Module[{exp},
  While[(exp=Read[stm])== Null];
  Return[exp];
  ]



(*
	Writes the output of the calculation for all parameter specification cases
*)
  writeLAGRAll[] := Module[{stm, case, no},
  Print["\nResults of Calculation"];
  stm = OpenWrite[outputFn];
  writeLAGR[stm,#] & /@ SystemCases;
  Close[stm];
  ]

(*
	Writes the output of the calculation as are:
	- the Metrik
	- the CHRISTOFFELsymbols
	- the generalized forces
	The symmetry is taken into account and only non-disappearing results are written
*)
  writeLAGR[stm_, case_List] := Module[{p},
  If[Length[case] == 0,
  Print["\nCommon Case"];
  WriteString[stm, "\nCommon Case\n"];,
  writeAssignment[stm, "\nSpecial Case", {}, Position[SystemCases,case][[1,1]]-1];
  ];
  p = Flatten[Table[{i, j},{i,degreesOfFreedom},{j,i}], 1];
  writeAssignment[stm,"g",{#[[1]],#[[2]]}, g[#[[1]], #[[2]]] /.case] & /@ p;
  p =Flatten[Table[{i, j, k},{i,degreesOfFreedom},{j,0,degreesOfFreedom},{k, 0,j}], 2];
  writeAssignment[stm,"Gam",{#[[1]],#[[2]],#[[3]]}, \[CapitalGamma][#[[1]], #[[2]],#[[3]]] /.case] & /@ p;
  p = Table[{i},{i,degreesOfFreedom}];
  writeAssignment[stm,"Q",{#[[1]]}, Q[#[[1]]] /.case] & /@ p;
  printAssignment["Potential of Springs",{}, U[] /.case];
  printAssignment["Generalized Force of Dampers ",{#[[1]]}, Qd[#[[1]]] /.case] & /@ p;
  ]


(*
	Writes a result in the form <name> = <value>
	Writes to screen and the output stream
*)
  writeAssignment[stm_,name_String, list_List, exp_]:=Module[{nameStr},
  nameStr = buildName[name,list];
  Print[nameStr, " = ", exp];
  If[(!IntegerQ[exp] ||!exp== 0) && (!NumericQ[exp] || ! exp == 0.),
  WriteString[stm,nameStr, " = "];
  Write[stm, exp];,
  Print["is 0"];
  ]
  ]


(*
	Build a name aa an output identifier using optionally an index list
*)
  buildName[name_String, indices_List]:= Module[{str},
  str = name;
  str = StringJoin @@ Prepend[Table[ToString[indices[[i]]],{i,Length[indices]}], str];
  Return[str];
  ]


(*
	Writes a result in the form <name> = <value>
	Writes to screen only
*)
  printAssignment[name_String, list_List, exp_]:= Module[{},
  Print[buildName[name,list], " = ", exp]]


(*
    Calculation of results:  
    - CHRISTOFFELsymbols  
    - Metrik  
    - generalized Forces  
*)  
  calcLAGR[] := Module[{},
  Length[En]]


(*
    Calculation of results:  
    - CHRISTOFFELsymbols  
*)
  \[CapitalGamma][\[Rho]_Integer,\[Mu]_Integer,\[Nu]_Integer]:= 
  \[CapitalGamma][\[Rho],\[Mu],\[Nu]]=1/2(\!\(
\*SubscriptBox[\(\[PartialD]\), \(q[\[Nu]]\)]\(g[\[Rho], \[Mu]]\)\) + \!\(
\*SubscriptBox[\(\[PartialD]\), \(q[\[Mu]]\)]\(g[\[Nu], \ \[Rho]]\)\)-\!\(
\*SubscriptBox[\(\[PartialD]\), \(q[\[Rho]]\)]\(g[\[Mu], \ \[Nu]]\)\)) // FullSimplify // Expand // Chop;


(*
    Calculation of results:  
    - Metrik  
*)
  g[\[Mu]_Integer,\[Nu]_Integer] :=
  g[\[Mu],\[Nu]] = Module[{i,j,k,l1},
  Return[Sum[m[l1]*D[xn[l1],q[\[Mu]]]. D[xn[l1],q[\[Nu]]]+ Sum[\[Theta][l1][[ j, k]]*D[En[l1][[j, i]],q[\[Mu]]]*D[En[l1][[k, i]],q[\[Nu]]],{i,3}, {j,3}, {k,3}],{l1,numBodies}]  // FullSimplify // Expand // Chop];
  ]


(*
    Calculation of results:  
    - generalized Forces  
*)  
  Q[\[Rho]_Integer] :=
  Q[\[Rho]] = Module[{i, j,k, l1},
  Return[
  (Sum[D[xn[l1],q[\[Rho]]].K[l1]+
  1/2Sum[D[En[l1][[j,i]],q[\[Rho]]]*En[l1][[j,k]]*M[l1][[i,k]],{j,3},{i,3},{k,3}],{l1, numBodies}] // FullSimplify // Expand // Chop)
  - D[U[],q[\[Rho]]]+Qd[\[Rho]]];
  ]


(*
	Calculation of the Potential of Springs
*)
  U[] :=Module[{},
  Return[(Plus @@ (calcU[#]& /@ Springs)) // Chop];
  ]


(*
	Calculation of the Potential of a single Spring
*)
  calcU[spring_List]:= Module[{k1, k2, rule, options, sq, dx, u,i},
  {k1, k2, rule, options} = spring;
  dx = calcF[spring];
  sq = dx.dx;
  u = c[k1,k2]*(1/2 sq - lambda[k1,k2]*Sqrt[sq]);
  u = u /.rule // Chop;
  u = FullSimplify[u,options] // Expand // Chop;
  Return[u ];
  ]

(*
	Calculates the coordinates differences of a single Spring
*)
  calcF[spring_List] := Module[{k1,k2,rule, options,dx,i},
  {k1, k2, rule, options} = spring;
  xs[k1,k2] = Table[xs[k1, k2, i],{i,3}];
  xs[k2,k1] = Table[xs[k2, k1, i],{i,3}];
  dx = Evaluate[(xn[k1] + xs[k1,k2].En[k1]) - (xn[k2] + xs[k2, k1].En[k2])];
  Return[dx];
  ]
  

(*
	Calculation of the generalized forces of the Dampers
*)
  Qd[\[Rho]_Integer]:=
  Qd[\[Rho]]=Module[{k1,k2,rule,dx, gf,\[Nu]},
  Return[(Plus @@ 
  (Function[damper,
  {k1, k2, rule} = damper;
  dx = calcD[damper];
  gf = Plus @@ Table[-k[k1,k2]*D[dx,q[\[Rho]]].D[dx,q[\[Nu]]]*qd[\[Nu]],{\[Nu],degreesOfFreedom}];
  gf = (gf /. rule // Chop);
  Return[gf];
  ] /@ Dampers)) // FullSimplify // Expand // Chop
  ];
  ]
  

(*
	Calculates the coordinates differences of a single Damper
*)
  calcD[damper_List] := Module[{k1,k2,rule,dx,i},
  {k1, k2, rule} = damper;
  xd[k1,k2] = Table[xd[k1, k2, i],{i,3}];
  xd[k2,k1] = Table[xd[k2, k1, i],{i,3}];
  dx = (xn[k1] + xd[k1,k2].En[k1]) - (xn[k2] + xd[k2, k1].En[k2]);
  Return[dx];
  ]


End[]

EndPackage[]

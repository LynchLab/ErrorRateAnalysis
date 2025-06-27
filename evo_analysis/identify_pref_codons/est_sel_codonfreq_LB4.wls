#!/usr/bin/env wolframscript
(* ::Package:: *)

## Convert observed frequencies of 4-nucleotides in different codon groups to population-scaled selection coefficients (S) under Li-Bulmer 4-allele model 

# arguments
# input table of counts of codons across codon groups
fin = $ScriptCommandLine[[2]]
# mutation-rate matrix
fmut = $ScriptCommandLine[[3]]
# output selection coefficients corresponding to 4 bases in each group
fsel = $ScriptCommandLine[[4]] 
# initial vector of selection coefficients
Si = {0.1, 0.2, 0.3};

# load required functions
Clear[P1,P2,P3,P4,Px1,Px2,Px3,Px4]
<<"identify_pref_codons/functions_codonfreq_to_sel_LB4.wls";
## Functions to calculate difference in observed and expected frequencies under given selection coefficients
f1[S2_,S3_,S4_]:=Evaluate[Px1[S2,S3,S4]-P1]
f2[S2_,S3_,S4_]:=Evaluate[Px2[S2,S3,S4]-P2]
f3[S2_,S3_,S4_]:=Evaluate[Px3[S2,S3,S4]-P3]
f4[S2_,S3_,S4_]:=Evaluate[Px4[S2,S3,S4]-P4]
## Jacobian Matrix for Newton-Raphson Approximation
JM[S2_,S3_,S4_]:=Evaluate[D[{f2[S2,S3,S4],f3[S2,S3,S4],f4[S2,S3,S4]},{{S2,S3,S4}}]]

# load data and mutation matrix
df = Import[fin, "Dataset", "HeaderLines" -> {1, 1}];
U = Import[fmut, "Dataset", "HeaderLines" -> {1, 1}];
u12 = U[1, 2];
u13 = U[1, 3];
u14 = U[1, 4];
u21 = U[2, 1];
u23 = U[2, 3];
u24 = U[2, 4];
# number of codon groups
nf = Length[df];
# provide any value for Ne
Ne = 1000;
# initialize output array
selout = ConstantArray[0, {nf, 4}];

# For each codon group
# convert counts to frequencies
# apply newton-raphson iterations
For[ i = 1, i <= nf, i++, 
	Clear[P1, P2, P3, P4];

	Ps = Values[Normal[df[[i, All]]]]/Total[df[[i, All]]];
	P1 = Ps[[1]]; P2 = Ps[[2]]; P3 = Ps[[3]]; P4 = Ps[[4]];

	X = Si;

	For[j = 0, j < 20, j++,

		J = JM[X[[1]], X[[2]], X[[3]]]; 

		F = {f2[X[[1]], X[[2]], X[[3]]], f3[X[[1]], X[[2]], X[[3]]], f4[X[[1]], X[[2]], X[[3]]]}; 

		If[FailureQ[Check[Y = LinearSolve[J, F], $Failed]], Break[]]; 
  		X = X - Y
	];

 	Sf = X;
 	selout[[i, 2 ;; 4]] = Sf;
]

# save output
Export[fsel, selout];

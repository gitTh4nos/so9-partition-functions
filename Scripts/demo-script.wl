(* ::Package:: *)

#!/usr/bin/env wolframscript
(* ::Package arguments::*)

Print["[Checkpoint] Script started."];

args = Rest[$ScriptCommandLine];
If[Length[args] < 2,
  Print["Usage: wolframscript -script demo-script.wl <order> <output_file>"];
  Exit[]
];

ThetaOrder = ToExpression[args[[1]]];
OutputPath = args[[2]];

Print["--- B4 Multiplicity Demo ---"];
Print["Target q-order: ", ThetaOrder];
Print["Output Path:    ", OutputPath];
Print["----------------------------"];

Print["[Checkpoint] Arguments parsed. Loading Section 1..."];

(*Section 1: B4 data*)

rho = {7, 5, 3, 1};
etaStar = {4, 0, 0, 0};

WeightAdd[aList_, bList_] := MapThread[Plus, {aList, bList}];
WeightSubtract[aList_, bList_] := MapThread[Subtract, {aList, bList}];

Mu[{l_, m_, r_, s_}] := Module[
{lambda},
  lambda = {2*l + 2*m + 2*r + s, 2*m + 2*r + s, 2*r + s, s};
  WeightAdd[lambda, rho]
];

RepFromMu[muList_] := Module[
{lam, dl, dm, dr, s},

  lam = WeightSubtract[muList, rho];
  dl = lam[[1]] - lam[[2]];
  dm = lam[[2]] - lam[[3]];
  dr = lam[[3]] - lam[[4]];
  s = lam[[4]];
  
  If[
  Min[dl, dm, dr, s] < 0 || OddQ[dl] || OddQ[dm] || OddQ[dr],
   $Failed,
   {dl/2, dm/2, dr/2, s}
  ]
];

(*Section 2: Weyl group*)

weylRhoShifts = Flatten[
  Table[
   Module[
   {base, pSign, weight, totalSign},
    pSign = Signature[perm];
    base = rho[[perm]];
    weight = signs*base;
    totalSign = pSign*Times @@ signs;
    {weight, totalSign}
   ], 
   {perm, Permutations[Range[4]]}, 
   {signs, Tuples[{-1, 1}, 4]}
  ], 
 1];

zero = {0, 0, 0, 0};
z1 = {2, 0, 0, 0};
z2 = {2, 2, 0, 0};
z3 = {1, 1, 1, -1};
z4 = {1, 1, 1, 1};

factors = {
  {zero, z3},
  {z1, z3},
  {z2, z3},
  {z2, WeightAdd[z1, z3]},
  {zero, z4}, 
  {z1, z4}, 
  {z2, z4},
  {z2, WeightAdd[z1, z4]}
};

denominator = {12, 8, 4, 0};
allChoices = Tuples[{1, 2}, 8];

allWeights = Map[Total[MapThread[Part, {factors, #}]] - denominator &, allChoices];

zq = Counts[allWeights];
zqTail = KeyDrop[zq, {etaStar}];

(*Section 3: Ansatz generation*)

AnsatzReps[n_] := Module[
{reps = {}},
  Do[
   If[
   a + 3 b + 2 d <= n - 1 
   && a + 3 b + 6 c + 3 d <= n + 1
   && a + 3 b + 6 c + 4 d <= n + 2 
   && a + 3 b + 6 c + 5 d <= n + 6,
    AppendTo[reps, {a, b, c, d}]
   ], 
   {a, 0, n - 1}, {b, 0, Floor[(n - 1)/3]}, {c, 0, Floor[(n + 1)/6]}, {d, 0, Floor[(n - 1)/2]}
  ];
  Sort[reps]
];

SignedDominantRep[weight_] := Module[
{absWeight, dominant, sign, rep},
  absWeight = Abs[weight];
  
  If[
  MemberQ[absWeight, 0] || Length[DeleteDuplicates[absWeight]] < 4,
   {$Failed, 0},
   
   dominant = Sort[absWeight, Greater];
   sign = Signature[Ordering[absWeight, All, Greater]] * Product[Sign[w], {w, weight}];
   rep = RepFromMu[dominant];
   {rep, sign}
  ]
];

Print["[Checkpoint] Generating Theta Polynomials..."];

(*Section 4: Theta Functions Evaluation*)

maxQ = 2*ThetaOrder + 3;

Hpoly = 1 + Sum[(-1)^n*Q^(n*(n + 1))*Sum[x^(2*k), {k, -n, n}], {n, 1, Ceiling[Sqrt[maxQ]]}];
Cpoly = 1 + Sum[Q^(n^2)*(x^(2*n) + x^(-2*n)), {n, 1, Ceiling[Sqrt[maxQ]]}];
Dpoly = 1 + Sum[(-1)^n*Q^(n^2)*(x^(2*n) + x^(-2*n)), {n, 1, Ceiling[Sqrt[maxQ]]}];
Bpoly = Sum[Q^(n*(n + 1))*(x^(2*n + 1) + x^(-(2*n + 1))), {n, 0, Ceiling[Sqrt[maxQ]]}];

Fpoly = Normal[Series[1/Hpoly, {Q, 0, maxQ}]];
G3poly = Expand[Normal[Series[Fpoly*Cpoly, {Q, 0, maxQ}]]];
G4poly = Expand[Normal[Series[Fpoly*Dpoly, {Q, 0, maxQ}]]];
GBpoly = Expand[Normal[Series[Fpoly*Bpoly, {Q, 0, maxQ}]]];

Clear[c3, c4, cB];
c3[e_] := c3[e] = Coefficient[G3poly, x, e];
c4[e_] := c4[e] = Coefficient[G4poly, x, e];
cB[e_] := cB[e] = Coefficient[GBpoly, x, e];

ThetaWeight[weight_, K_] := Module[{e1, e2, e3, e4, ta, tb, val3, val4, valB},
  {e1, e2, e3, e4} = weight;
  
  ta = 2*K + 1;
  tb = 2*K;
  
  val3 = Coefficient[Expand[c3[e1]*c3[e2]*c3[e3]*c3[e4]], Q, ta];
  val4 = Coefficient[Expand[c4[e1]*c4[e2]*c4[e3]*c4[e4]], Q, ta];
  valB = Coefficient[Expand[cB[e1]*cB[e2]*cB[e3]*cB[e4]], Q, tb];
  
  (val3 - val4 + valB)/2
];

Print["[Checkpoint] Entering Main Loop..."];

(*Section 5: Main Loop*)

Clear[a];
a[_] = 0;

ansatz = AnsatzReps[ThetaOrder];
descendingAnsatz = SortBy[ansatz, -Mu[#] &];

Print["Computing multiplicities for ", Length[descendingAnsatz], " representations at order ", ThetaOrder, "..."];

Do[
  Module[
  {mu, gamma, rhs, correction, shiftedWeight, alpha, neighborRep, sg, domResult},
   mu = Mu[rep];
   gamma = WeightAdd[mu, etaStar];
   rhs = 0;
   
   Do[
    shiftedWeight = WeightSubtract[gamma, shift[[1]]];
    rhs += shift[[2]]*ThetaWeight[shiftedWeight, ThetaOrder],
    {shift, weylRhoShifts}
   ];
   
   correction = 0;
   KeyValueMap[
    Function[{eta, coeff},
     alpha = WeightSubtract[gamma, eta];
     domResult = SignedDominantRep[alpha];
     neighborRep = domResult[[1]];
     sg = domResult[[2]];
     
     If[neighborRep =!= $Failed,
      correction += coeff*sg*a[neighborRep]
     ];
    ], 
    zqTail
   ];
   
   a[rep] = rhs - correction;
  ], 
  {rep, descendingAnsatz}
];

nonzeroMultiplicities = Select[AssociationMap[a, descendingAnsatz], # != 0 &];
Print["Done. Nonzero representations: ", Length[nonzeroMultiplicities]];

(*Section 6: Saving data*)

finalOutput = <|ThetaOrder -> nonzeroMultiplicities|>;
Export[OutputPath, finalOutput];

Print["Saved output successfully to ", OutputPath];

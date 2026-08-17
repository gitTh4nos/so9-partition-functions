(* ::Package:: *)

(* SO(9) / B4 partition-function decomposition over GF(2^61-1). *)

(*==========================Instructions for Batch mode========================*)

(* On the bash command line type:
wolfram -script math_Code_simplified_batch.wl -Maxlevel N -kernels K -output outputfile -Validpts V
where the defaults are N=40, K=$ProcessorNumber, outputfile= results.m, V=3*)

(*Make sure to have a path set to Mathematica by inserting the following in  .zprofile :
PATH=$PATH:/Applications/Wolfram.app/Contents/MacOS *)

(*A logfile is written to math_Code_simplified_batch.log *)

(*To load the results into Mathematica put in the Input line
data=Uncompress[Get["results.m"]]; *)

(*Example commands:*)
(*
data["Multiplicities", 40]  (*Give all multiplicities at level 40*)
data["Solutions", 40, "InheritedAssociation"] (*Inherited multiplicities*)
data["Solutions", 40, "Stats"]
*)
(* All choices for run["Solutions",level,$string]*)
(*
$string="Level", "Labels","Coefficients","Association","NonzeroAssociation",
"InheritedAssociation","InheritedLabels","SolvedLabels","Stats","Validation"
*)



(* ============ BATCH MODE SETUP ============ *)

(*If you want to put the output file in a specific directory*)
SetDirectory["../DataFiles"];


$HistoryLength = 0;
$Epilog = Null;
Off[General::spell];
Off[General::shdw];
Off[General::stop];

(* Logging *)
logFile = "math_Code.log";
logStream = OpenWrite[logFile];

batchLog[msg___] := Module[{output},
  output = ToString[StringJoin @@ Map[ToString, {msg}]];
  Print[output];
  WriteString[logStream, output <> "\n"];
  Flush[logStream];
];

batchLog["========================================================"];
batchLog["math_Code_simplified_batch.wl - Generate SO(9) Multiplicities"];
batchLog["Start: ", DateString[]];
batchLog["Process ID: ", $ProcessID];
batchLog["========================================================"];

(* ============ COMMAND-LINE PARAMETERS ============ *)

batchLog["Parsing arguments..."];

cmdArgs = Rest[$CommandLine];
params = <|
  "Maxlevel" -> 40,
  "output" -> "results.m",
  "kernels" -> $ProcessorCount,
  "Validpts" -> 3
|>;
Do[
  If[i < Length[cmdArgs],
    Switch[cmdArgs[[i]],
      "-Maxlevel", params["Maxlevel"] = ToExpression[cmdArgs[[i + 1]]],
      "-output", params["output"] = cmdArgs[[i + 1]],
      "-kernels", params["kernels"] = ToExpression[cmdArgs[[i + 1]]],
      "-Validpts", params["Validpts"]=ToExpression[cmdArgs[[i + 1]]],
       _, Null
    ]
  ],
  {i, 1, Length[cmdArgs], 2}
];

(* Parse representation list *)
Maxlevel=params["Maxlevel"];kernels=params["kernels"];
output=params["output"];Valid=params["Validpts"];

batchLog["Parameters:"];
batchLog["  Maxlevel: ", params["Maxlevel"]];
batchLog["  Output: ", params["output"]];
batchLog["  Kernels: ", params["kernels"]];
batchLog["  Number of validation points: ", params["Validpts"]];



ClearAll[
  levelLabels, muB4, characterData, powerDifferenceTable,
  alternantFromTable, characterValuesMod, primePoint, primePoints,
  seriesMulMod, seriesInvMod, zqxValMod, thetaABCoeffsMod,
  partitionRHSValueMod, partitionProjectedRHSValuesMod,
  pMax, offsetS, knownCoefficient, inheritedValue, splitInheritedLabels,
  coefficientAt, parallelMap, initializeParallelSolver,
  unknownLayoutsMod, buildCharacterCacheMod,
  prepareParityBlock, prepareParityBlockCached,
  solvePreparedBlock, solvePreparedBlocks,
  validateCoefficientsMod, validateLevelSolutionMod,
  solveLevelMod, solveLevelsMod
];

$B4Modulus = 2^61 - 1;
$B4Rho = {7, 5, 3, 1};

(* ---------- B4 labels and modular characters ---------- *)

levelLabels[n_Integer?NonNegative] := levelLabels[n] = Sort @ Reap[
    Do[
      If[
        ell + 3 m + 2 s <= n &&
        ell + 3 m + 3 s <= n + 1 &&
        ell + 3 m + 6 r + 3 s <= n + 2 &&
        ell + 3 m + 6 r + 4 s <= n + 3 &&
        ell + 3 m + 6 r + 5 s <= n + 7,
        Sow[{ell, m, r, s}]
      ],
      {ell, 0, n},
      {m, 0, Floor[n/3]},
      {r, 0, Floor[(n + 2)/6]},
      {s, 0, Floor[n/2]}
    ]
  ][[2, 1]];

muB4[{ell_, m_, r_, s_}] := {
  2 ell + 2 m + 2 r + s + 7,
  2 m + 2 r + s + 5,
  2 r + s + 3,
  s + 1
};

characterData[labels_List] := <|
  "Labels" -> labels,
  "Mu" -> Developer`ToPackedArray[muB4 /@ labels]
|>;

powerDifferenceTable[y_List, maxExponent_Integer, modulus_Integer] := Module[
  {yy = Mod[y, modulus], inverse, positive, negative},
  If[Length[yy] =!= 4 || MemberQ[yy, 0], Return[$Failed]];
  inverse = PowerMod[#, -1, modulus] & /@ yy;
  positive = Table[
    NestList[Mod[# yy[[i]], modulus] &, 1, maxExponent],
    {i, 4}
  ];
  negative = Table[
    NestList[Mod[# inverse[[i]], modulus] &, 1, maxExponent],
    {i, 4}
  ];
  Mod[positive - negative, modulus]
];

alternantFromTable[mu_List, table_, modulus_Integer] :=
  Det[table[[All, mu + 1]], Modulus -> modulus];

characterValuesMod::badpoint =
  "The Weyl denominator vanishes at finite-field point `1`.";

characterValuesMod[
  data_Association,
  y_List,
  modulus_Integer : $B4Modulus
] := Module[{mus = data["Mu"], table, denominator, numerators},
  table = powerDifferenceTable[y, Max[Flatten[mus]], modulus];
  If[table === $Failed, Return[$Failed]];
  denominator = alternantFromTable[$B4Rho, table, modulus];
  If[denominator == 0,
    Message[characterValuesMod::badpoint, y];
    Return[$Failed]
  ];
  numerators = alternantFromTable[#, table, modulus] & /@ mus;
  Mod[PowerMod[denominator, -1, modulus] numerators, modulus]
];

primePoint[i_Integer?Positive] := Prime /@ Range[4 i - 3, 4 i];
primePoints[n_Integer?NonNegative, start_Integer : 1] :=
  primePoint /@ Range[start, start + n - 1];

(* ---------- Normalized partition-function value ---------- *)

seriesMulMod[
  a_List,
  b_List,
  length_Integer,
  modulus_Integer : $B4Modulus
] := Module[{out = ConstantArray[0, length]},
  Do[
    If[a[[i + 1]] =!= 0,
      Do[
        out[[i + j + 1]] = Mod[
          out[[i + j + 1]] + a[[i + 1]] b[[j + 1]],
          modulus
        ],
        {j, 0, length - 1 - i}
      ]
    ],
    {i, 0, length - 1}
  ];
  out
];

seriesInvMod[
  a_List,
  length_Integer,
  modulus_Integer : $B4Modulus
] := Module[{out = ConstantArray[0, length], sum},
  out[[1]] = PowerMod[a[[1]], -1, modulus];
  Do[
    sum = Sum[a[[j + 1]] out[[degree - j + 1]], {j, 1, degree}];
    out[[degree + 1]] = Mod[-out[[1]] sum, modulus],
    {degree, 1, length - 1}
  ];
  out
];

zqxValMod[xs_List, modulus_Integer : $B4Modulus] := Module[
  {x1, x2, x3, x4, z1, z2, z3, z4, numerator, denominator},
  {x1, x2, x3, x4} = Mod[xs, modulus];
  z1 = Mod[x1^2, modulus];
  z2 = Mod[z1 x2^2, modulus];
  z3 = Mod[x1 x2 x3 PowerMod[x4, -1, modulus], modulus];
  z4 = Mod[x1 x2 x3 x4, modulus];
  numerator = Times @@ Mod[{
    1 + z3, z1 + z3, z2 + z3, z2 + z1 z3,
    1 + z4, z1 + z4, z2 + z4, z2 + z1 z4
  }, modulus];
  denominator = Mod[z1^2 z2^2 z3^2 z4^2, modulus];
  Mod[numerator PowerMod[denominator, -1, modulus], modulus]
];

thetaABCoeffsMod[
  xs_List,
  k_Integer?NonNegative,
  modulus_Integer : $B4Modulus
] := Module[
  {kexp = k + 1, length, f, p3, p4, pb, x, xi, diff,
   ai, aiInv, bv, cv, dv, xpow, xipow, x2, xi2, degree, n,
   apart, bpart},

  length = 2 kexp + 3;
  f = p3 = p4 = pb = ReplacePart[ConstantArray[0, length], 1 -> 1];

  Do[
    x = Mod[x0, modulus];
    xi = PowerMod[x, -1, modulus];
    diff = Mod[x - xi, modulus];
    x2 = Mod[x^2, modulus];
    xi2 = Mod[xi^2, modulus];

    ai = ConstantArray[0, length];
    xpow = x;
    xipow = xi;
    For[n = 0, n (n + 1) < length, n++,
      degree = n (n + 1);
      ai[[degree + 1]] = Mod[(-1)^n (xpow - xipow), modulus];
      xpow = Mod[xpow x2, modulus];
      xipow = Mod[xipow xi2, modulus];
    ];
    aiInv = Mod[diff seriesInvMod[ai, length, modulus], modulus];
    f = seriesMulMod[f, aiInv, length, modulus];

    bv = ConstantArray[0, length];
    xpow = x;
    xipow = xi;
    For[n = 0, n (n + 1) < length, n++,
      degree = n (n + 1);
      bv[[degree + 1]] = Mod[xpow + xipow, modulus];
      xpow = Mod[xpow x2, modulus];
      xipow = Mod[xipow xi2, modulus];
    ];
    pb = seriesMulMod[pb, bv, length, modulus];

    cv = dv = ReplacePart[ConstantArray[0, length], 1 -> 1];
    xpow = x2;
    xipow = xi2;
    For[n = 1, n^2 < length, n++,
      degree = n^2;
      cv[[degree + 1]] = Mod[xpow + xipow, modulus];
      dv[[degree + 1]] = Mod[(-1)^n (xpow + xipow), modulus];
      xpow = Mod[xpow x2, modulus];
      xipow = Mod[xipow xi2, modulus];
    ];
    p3 = seriesMulMod[p3, cv, length, modulus];
    p4 = seriesMulMod[p4, dv, length, modulus],
    {x0, xs}
  ];

  apart = seriesMulMod[f, Mod[p3 - p4, modulus], length, modulus];
  bpart = seriesMulMod[f, pb, length, modulus];
  {apart[[2 kexp + 2]], bpart[[2 kexp + 1]]}
];

partitionRHSValueMod[
  k_Integer?NonNegative,
  y_List,
  modulus_Integer : $B4Modulus
] := Module[{a, b, inv2 = PowerMod[2, -1, modulus]},
  {a, b} = thetaABCoeffsMod[y, k, modulus];
  Mod[
    (a + b) inv2 PowerMod[zqxValMod[y, modulus], -1, modulus],
    modulus
  ]
];

(*
  Central-parity projection.

  For a B4 highest weight {ell,m,r,s},

    chi_{ell,m,r,s}(-y1,y2,y3,y4) = (-1)^s chi_{ell,m,r,s}(y).

  Consequently

    Z_even(y) = (Z(y) + Z(-y1,y2,y3,y4))/2,
    Z_odd (y) = (Z(y) - Z(-y1,y2,y3,y4))/2

  contain only even-s and odd-s characters, respectively.  In the theta
  representation A is unchanged and B changes sign, so both projections are
  obtained from one evaluation of {A,B}.  This replaces one dense system by
  two smaller independent systems.
*)
partitionProjectedRHSValuesMod[
  k_Integer?NonNegative,
  y_List,
  modulus_Integer : $B4Modulus
] := Module[{a, b, inv2, yFlip, rhs, rhsFlip},
  {a, b} = thetaABCoeffsMod[y, k, modulus];
  inv2 = PowerMod[2, -1, modulus];
  yFlip = ReplacePart[Mod[y, modulus], 1 -> Mod[-y[[1]], modulus]];
  rhs = Mod[
    (a + b) inv2 PowerMod[zqxValMod[y, modulus], -1, modulus],
    modulus
  ];
  rhsFlip = Mod[
    (a - b) inv2 PowerMod[zqxValMod[yFlip, modulus], -1, modulus],
    modulus
  ];
  Mod[inv2 {rhs + rhsFlip, rhs - rhsFlip}, modulus]
];

(* ---------- Inherited multiplicities ---------- *)

pMax[n_Integer, q_Integer, r_Integer, s_Integer] := Which[
  s >= 4, n - 3 q - 6 r - 5 s + 7,
  s == 3, n - 3 q - 6 r - 9,
  s == 2, n - 3 q - 6 r - 5,
  s == 1 && r >= 1, n - 3 q - 6 r - 1,
  s == 1, n - 3 q - 2,
  r >= 1, n - 3 q - 6 r + 2,
  True, n - 3 q
];

offsetS[s_Integer] := Which[s > 3, 0, s >= 1, 1, True, 2];

knownCoefficient[
  results_Association,
  level_Integer,
  rep_List
] := Module[{previous = Lookup[results, level, Missing["UnknownLevel"]]},
  If[MissingQ[previous],
    Missing["UnknownLevel"],
    If[KeyExistsQ[previous, rep], previous[rep], 0]
  ]
];

inheritedValue[
  {p_Integer, q_Integer, r_Integer, s_Integer},
  n_Integer?Positive,
  results_Association
] := Catch @ Module[{pm = pMax[n, q, r, s], value, shift},

  Do[
    If[p > Floor[pm/2] + Floor[shift/2],
      value = knownCoefficient[results, n - shift, {p - shift, q, r, s}];
      If[!MissingQ[value], Throw[value]]
    ],
    {shift, Min[p, n - 1], 1, -1}
  ];

  Do[
    If[p > pm - q + shift - 1,
      value = knownCoefficient[
        results, n - 3 shift, {p, q - shift, r, s}
      ];
      If[!MissingQ[value], Throw[value]]
    ],
    {shift, Min[q, Floor[(n - 1)/3]], 1, -1}
  ];

  Do[
    If[p > pm - r + shift - 1 + offsetS[s],
      value = knownCoefficient[
        results, n - 6 shift, {p, q, r - shift, s}
      ];
      If[!MissingQ[value], Throw[value]]
    ],
    {shift, Min[r, Floor[(n - 1)/6]], 1, -1}
  ];

  Missing["NotInherited"]
];

splitInheritedLabels[
  n_Integer?Positive,
  labels_List,
  results_Association
] := Module[{classified, inheritedPairs, unknown, nonzeroPairs},
  classified = {#, inheritedValue[#, n, results]} & /@ labels;
  inheritedPairs = Select[classified, !MissingQ[Last[#]] &];
  unknown = First /@ Select[classified, MissingQ[Last[#]] &];
  nonzeroPairs = Select[inheritedPairs, Last[#] =!= 0 &];
  <|
    "Association" -> Association[Rule @@@ nonzeroPairs],
    "Labels" -> (First /@ inheritedPairs),
    "UnknownLabels" -> unknown
  |>
];

coefficientAt[association_Association, rep_List] :=
  If[KeyExistsQ[association, rep], association[rep], 0];

(* ---------- Parallel recurrence-reduced solver ---------- *)

parallelMap[function_, items_List, useParallel_] := If[
  TrueQ[useParallel] && $KernelCount > 1 && Length[items] > 1,
  ParallelMap[function, items, Method -> "CoarsestGrained"],
  Map[function, items]
];

initializeParallelSolver[] := Module[{},
  If[$KernelCount == 0, Return[0]];
  If[
    !TrueQ[$B4ParallelInitialized] ||
    $B4ParallelKernelCount =!= $KernelCount,
    DistributeDefinitions[
      $B4Modulus, $B4Rho,
      powerDifferenceTable, alternantFromTable, characterValuesMod,
      seriesMulMod, seriesInvMod, zqxValMod, thetaABCoeffsMod,
      partitionRHSValueMod, partitionProjectedRHSValuesMod,
      solvePreparedBlock
    ];
    $B4ParallelInitialized = True;
    $B4ParallelKernelCount = $KernelCount;
  ];
  $KernelCount
];

(*
  Determine the unknown labels using only level availability.  The numerical
  values of lower-level multiplicities do not affect whether a recurrence
  applies, so empty Associations are sufficient for this planning pass.
*)
unknownLayoutsMod[levels_List, initialResults_Association : <||>] := Module[
  {available, layouts = <||>, split},
  available = AssociationThread[
    Keys[initialResults] -> ConstantArray[<||>, Length[initialResults]]
  ];
  Do[
    split = splitInheritedLabels[level, levelLabels[level], available];
    AssociateTo[layouts, level -> split["UnknownLabels"]];
    AssociateTo[available, level -> <||>],
    {level, Sort @ DeleteDuplicates[levels]}
  ];
  layouts
];

(*
  Build one character table for a complete multi-level run.  Rows correspond
  to prime torus points and columns to B4 labels.  Every parity block at every
  level is then obtained by Part extraction; no character determinant is
  recomputed inside the level loop.
*)
buildCharacterCacheMod[
  levels_List,
  initialResults_Association : <||>,
  modulus_Integer : $B4Modulus,
  startPoint_Integer : 1,
  useParallel_ : True
] := Module[
  {t0 = AbsoluteTime[], orderedLevels, layouts, labels, maxPoints, points,
   data, table},

  orderedLevels = Sort @ DeleteDuplicates[levels];
  layouts = unknownLayoutsMod[orderedLevels, initialResults];
  labels = Sort @ DeleteDuplicates @ Flatten[levelLabels /@ orderedLevels, 1];
  maxPoints = Max[
    0,
    Sequence @@ Flatten[
      Table[
        Count[layouts[level], rep_ /; Mod[Last[rep], 2] == parity],
        {level, orderedLevels}, {parity, 0, 1}
      ]
    ]
  ];
  points = primePoints[maxPoints, startPoint];
  data = characterData[labels];
  table = Developer`ToPackedArray @ parallelMap[
    characterValuesMod[data, #, modulus] &,
    points,
    useParallel
  ];

  <|
    "Labels" -> labels,
    "Index" -> AssociationThread[labels -> Range[Length[labels]]],
    "Points" -> points,
    "Table" -> table,
    "Layouts" -> layouts,
    "Modulus" -> modulus,
    "StartPoint" -> startPoint,
    "BuildSeconds" -> N[AbsoluteTime[] - t0]
  |>
];

prepareParityBlock[
  parity_Integer,
  unknownLabels_List,
  inherited_Association,
  points_List,
  projectedRHS_List,
  modulus_Integer,
  useParallel_
] := Module[
  {blockLabels, inheritedLabels, blockData, inheritedData, evaluations,
   matrix, rhs, inheritedCoefficients},

  blockLabels = Select[unknownLabels, Mod[Last[#], 2] == parity &];
  If[blockLabels === {}, Return[Nothing]];

  inheritedLabels = Select[Keys[inherited], Mod[Last[#], 2] == parity &];
  blockData = characterData[blockLabels];
  inheritedData = If[inheritedLabels === {}, None, characterData[inheritedLabels]];

  evaluations = parallelMap[
    Function[point,
      {
        characterValuesMod[blockData, point, modulus],
        If[inheritedData === None,
          {},
          characterValuesMod[inheritedData, point, modulus]
        ]
      }
    ],
    points,
    useParallel
  ];

  matrix = evaluations[[All, 1]];
  rhs = projectedRHS;
  If[inheritedLabels =!= {},
    inheritedCoefficients = coefficientAt[inherited, #] & /@ inheritedLabels;
    rhs = Mod[rhs - evaluations[[All, 2]] . inheritedCoefficients, modulus]
  ];

  <|
    "Parity" -> parity,
    "Labels" -> blockLabels,
    "Matrix" -> matrix,
    "RHS" -> rhs
  |>
];

prepareParityBlockCached[
  parity_Integer,
  unknownLabels_List,
  inherited_Association,
  projectedRHS_List,
  cache_Association,
  modulus_Integer
] := Module[
  {blockLabels, inheritedLabels, rows, blockColumns, inheritedColumns,
   matrix, rhs, inheritedCoefficients},

  blockLabels = Select[unknownLabels, Mod[Last[#], 2] == parity &];
  If[blockLabels === {}, Return[Nothing]];
  inheritedLabels = Select[Keys[inherited], Mod[Last[#], 2] == parity &];

  rows = Range[Length[blockLabels]];
  blockColumns = cache["Index"][#] & /@ blockLabels;
  matrix = cache["Table"][[rows, blockColumns]];
  rhs = projectedRHS;

  If[inheritedLabels =!= {},
    inheritedColumns = cache["Index"][#] & /@ inheritedLabels;
    inheritedCoefficients = coefficientAt[inherited, #] & /@ inheritedLabels;
    rhs = Mod[
      rhs - cache["Table"][[rows, inheritedColumns]] . inheritedCoefficients,
      modulus
    ]
  ];

  <|
    "Parity" -> parity,
    "Labels" -> blockLabels,
    "Matrix" -> matrix,
    "RHS" -> rhs
  |>
];

solvePreparedBlock[block_Association, modulus_Integer] := Module[
  {coefficients, nonzero},
  coefficients = Mod[
    LinearSolve[block["Matrix"], block["RHS"], Modulus -> modulus],
    modulus
  ];
  nonzero = AssociationThread[
    Pick[block["Labels"], Unitize[coefficients], 1] ->
    Pick[coefficients, Unitize[coefficients], 1]
  ];
  <|
    "Parity" -> block["Parity"],
    "Association" -> nonzero
  |>
];

solvePreparedBlocks[blocks_List, modulus_Integer, useParallel_] := If[
  TrueQ[useParallel] && $KernelCount > 1 && Length[blocks] > 1,
  WaitAll[
    With[{block = #, prime = modulus},
      ParallelSubmit[solvePreparedBlock[block, prime]]
    ] & /@ blocks
  ],
  solvePreparedBlock[#, modulus] & /@ blocks
];

validateCoefficientsMod[
  k_Integer?Positive,
  coefficients_Association,
  points_List,
  modulus_Integer,
  useParallel_
] := Module[{labels = Keys[coefficients], data, values, lhs, rhs},
  If[points === {},
    Return[<|
      "Points" -> {},
      "Checked" -> False,
      "Valid" -> Missing["NotChecked"]
    |>]
  ];
  data = If[labels === {}, None, characterData[labels]];
  values = parallelMap[
    Function[point,
      {
        If[data === None, {}, characterValuesMod[data, point, modulus]],
        partitionRHSValueMod[k, point, modulus]
      }
    ],
    points,
    useParallel
  ];
  lhs = If[
    labels === {},
    ConstantArray[0, Length[points]],
    Mod[values[[All, 1]] . (coefficientAt[coefficients, #] & /@ labels), modulus]
  ];
  rhs = values[[All, 2]];
  <|
    "Points" -> points,
    "LHS" -> lhs,
    "RHS" -> rhs,
    "Checked" -> True,
    "Valid" -> (lhs === rhs)
  |>
];

validateLevelSolutionMod[
  solution_Association,
  nPoints_Integer : 3,
  modulus_Integer : $B4Modulus,
  useParallel_ : True
] := Module[{maxBlock, points},
  If[TrueQ[useParallel], initializeParallelSolver[]];
  maxBlock = Max[
    solution["Stats", "ParityBlocks", "Even", "Solved"],
    solution["Stats", "ParityBlocks", "Odd", "Solved"]
  ];
  points = primePoints[nPoints, maxBlock + 1];
  validateCoefficientsMod[
    solution["Level"],
    solution["NonzeroAssociation"],
    points,
    modulus,
    useParallel
  ]
];

Options[solveLevelMod] = {
  "Modulus" -> $B4Modulus,
  "StartPoint" -> 1,
  "ValidationPoints" -> 3,
  "Parallel" -> True,
  "CharacterCache" -> None
};

solveLevelMod[
  k_Integer?Positive,
  previousResults_Association : <||>,
  OptionsPattern[]
] := Module[
  {modulus = OptionValue["Modulus"], start = OptionValue["StartPoint"],
   nValidation = OptionValue["ValidationPoints"],
   useParallel = TrueQ[OptionValue["Parallel"]], allLabels, split,
   characterCache = OptionValue["CharacterCache"], t0 = AbsoluteTime[],
   inherited, unknown, unknownByParity, maxBlock, points, projectedRHS,
   blocks, solvedBlocks, solved, merged, orderedNonzero, fullCoefficients,
   validationPoints, validation, parityStats},

  If[useParallel, initializeParallelSolver[]];

  allLabels = levelLabels[k];
  split = splitInheritedLabels[k, allLabels, previousResults];
  inherited = split["Association"];
  unknown = split["UnknownLabels"];
  unknownByParity = Table[Select[unknown, Mod[Last[#], 2] == parity &], {parity, 0, 1}];
  maxBlock = Max[Length /@ unknownByParity];

  If[maxBlock > 0,
    points = If[
      AssociationQ[characterCache],
      Take[characterCache["Points"], maxBlock],
      primePoints[maxBlock, start]
    ];
    projectedRHS = parallelMap[
      partitionProjectedRHSValuesMod[k, #, modulus] &,
      points,
      useParallel
    ];
    blocks = DeleteCases[
      Table[
        With[{count = Length[unknownByParity[[parity + 1]]]},
          If[count == 0,
            Nothing,
            If[AssociationQ[characterCache],
              prepareParityBlockCached[
                parity,
                unknown,
                inherited,
                Take[projectedRHS[[All, parity + 1]], count],
                characterCache,
                modulus
              ],
              prepareParityBlock[
                parity,
                unknown,
                inherited,
                Take[points, count],
                Take[projectedRHS[[All, parity + 1]], count],
                modulus,
                useParallel
              ]
            ]
          ]
        ],
        {parity, 0, 1}
      ],
      Nothing
    ];
    solvedBlocks = solvePreparedBlocks[blocks, modulus, useParallel];
    solved = Join @@ (#["Association"] & /@ solvedBlocks),
    solved = <||>
  ];

  merged = Join[inherited, solved];
  orderedNonzero = Association @ Cases[
    allLabels,
    rep_ /; coefficientAt[merged, rep] =!= 0 :>
      (rep -> coefficientAt[merged, rep])
  ];
  fullCoefficients = coefficientAt[orderedNonzero, #] & /@ allLabels;

  validationPoints = primePoints[nValidation, start + maxBlock];
  validation = validateCoefficientsMod[
    k, orderedNonzero, validationPoints, modulus, useParallel
  ];
  parityStats = <|
    "Even" -> <|
      "InheritedNonzero" -> Count[Keys[inherited], rep_ /; EvenQ[Last[rep]]],
      "Solved" -> Length[unknownByParity[[1]]]
    |>,
    "Odd" -> <|
      "InheritedNonzero" -> Count[Keys[inherited], rep_ /; OddQ[Last[rep]]],
      "Solved" -> Length[unknownByParity[[2]]]
    |>
  |>;

  <|
    "Level" -> k,
    "Labels" -> allLabels,
    "Coefficients" -> fullCoefficients,
    "Association" -> AssociationThread[allLabels -> fullCoefficients],
    "NonzeroAssociation" -> orderedNonzero,
    "InheritedAssociation" -> inherited,
    "InheritedLabels" -> split["Labels"],
    "SolvedLabels" -> unknown,
    "Stats" -> <|
      "Total" -> Length[allLabels],
      "Inherited" -> Length[split["Labels"]],
      "InheritedNonzero" -> Length[inherited],
      "Solved" -> Length[unknown],
      "ParityBlocks" -> parityStats,
      "ElapsedSeconds" -> N[AbsoluteTime[] - t0]
    |>,
    "Validation" -> validation
  |>
];

Options[solveLevelsMod] = Join[
  Options[solveLevelMod],
  {
    "InitialMultiplicities" -> <||>,
    "PrecomputeCharacters" -> True,
    "Progress" -> True
  }
];

solveLevelsMod[maxLevel_Integer?Positive, opts : OptionsPattern[]] :=
  solveLevelsMod[Range[maxLevel], opts];

solveLevelsMod[levels_List, OptionsPattern[]] := Module[
  {orderedLevels = Sort @ DeleteDuplicates[levels],
   results = OptionValue["InitialMultiplicities"], solutions = <||>, solution,
   levelOptions, useParallel = TrueQ[OptionValue["Parallel"]],
   cache = None, cacheStats, t0 = AbsoluteTime[]},

  If[useParallel, initializeParallelSolver[]];
  If[TrueQ[OptionValue["PrecomputeCharacters"]],
    cache = buildCharacterCacheMod[
      orderedLevels,
      results,
      OptionValue["Modulus"],
      OptionValue["StartPoint"],
      useParallel
    ]
  ];

  levelOptions = FilterRules[
    {
      "Modulus" -> OptionValue["Modulus"],
      "StartPoint" -> OptionValue["StartPoint"],
      "ValidationPoints" -> OptionValue["ValidationPoints"],
      "Parallel" -> OptionValue["Parallel"],
      "CharacterCache" -> cache
    },
    Options[solveLevelMod]
  ];

  Do[
    solution = solveLevelMod[level, results, Sequence @@ levelOptions];
    AssociateTo[solutions, level -> solution];
    AssociateTo[results, level -> solution["NonzeroAssociation"]];
    If[TrueQ[OptionValue["Progress"]],
      batchLog[
        "level ", level,
        ": total=", solution["Stats", "Total"],
        ", inherited=", solution["Stats", "Inherited"],
        ", solved=", solution["Stats", "Solved"],
        ", valid=", solution["Validation", "Valid"]
      ]
    ],
    {level, orderedLevels}
  ];

  cacheStats = If[
    AssociationQ[cache],
    <|
      "Rows" -> Length[cache["Points"]],
      "Columns" -> Length[cache["Labels"]],
      "BuildSeconds" -> cache["BuildSeconds"]
    |>,
    None
  ];
  <|
    "Solutions" -> solutions,
    "Multiplicities" -> results,
    "CharacterCacheStats" -> cacheStats,
    "ElapsedSeconds" -> N[AbsoluteTime[] - t0]
  |>
];


(* Main run:*)
batchLog["Launching kernels"]
   LaunchKernels[kernels];
   run = solveLevelsMod[
     Maxlevel,
     "Parallel" -> True,
     "PrecomputeCharacters" -> True,
     "ValidationPoints" -> 0
   ];
   (*run["Multiplicities", 40]
   run["Solutions", 40, "InheritedAssociation"]
   run["Solutions", 40, "Stats"]*)

 (* Validate only the final level after the timed production run. *)
   batchLog["Validate solution:",
   validateLevelSolutionMod[run["Solutions", Maxlevel], Valid][["Valid"]]];

 (*  A direct level solve is also valid, but inheritance is available only when
   previous levels are supplied:
*)
  (* sol40 = solveLevelMod[40, run["Multiplicities"]];*)



Export[output,Compress[run]];


batchLog["Total time:",run[["ElapsedSeconds"]]];


(*run[["Multiplicities",2]][{2,0,0,0}]*)

(* ::Package:: *)

(*Batch Routine to compute  tables of coefficients in appendix A *)

(*==========================Instructions for Batch mode========================*)

(* On the bash command line type:
wolfram -script coeff_gen_1.wl -reps "{{p1,q1,r1,s1},{p2,q2,r2,s2},...}" -Nmax Nmax -kernels K
 -output outputfile -coeffFile ***.m
where the defaults are "{{p1,q1,r1,s1},{p2,q2,r2,s2},...}"= "{{0,0,0,0}}", Nmax=100, 
K=$ProcessorNumber, outputfile= results.m, "coeffFile"="Blank" so it generates the coefficients*)

(*Make sure to have a path set to Mathematica by inserting the following in  .zprofile :
PATH=$PATH:/Applications/Wolfram.app/Contents/MacOS *)

(*A logfile is written to coeff_gen_1.log *)


SetDirectory["../DataFiles/Super"];

(* ============ BATCH MODE SETUP ============ *)

$HistoryLength = 0;
$Epilog = Null;
Off[General::spell];
Off[General::shdw];
Off[General::stop];

(* Logging *)
logFile = "coeff_gen_1.log";
logStream = OpenWrite[logFile];

batchLog[msg___] := Module[{output},
  output = ToString[StringJoin @@ Map[ToString, {msg}]];
  Print[output];
  WriteString[logStream, output <> "\n"];
  Flush[logStream];
];

batchLog["========================================================"];
batchLog["coeff_gen_1.wl - Generate leading Rademacher"];
batchLog["Start: ", DateString[]];
batchLog["Process ID: ", $ProcessID];
batchLog["========================================================"];

(* ============ COMMAND-LINE PARAMETERS ============ *)

batchLog["Parsing arguments..."];

cmdArgs = Rest[$CommandLine];
params = <|
  "reps" -> "{{0,0,0,0}}",
  "Nmax" -> 100,
  "output" -> "results.m",
  "kernels" -> $ProcessorCount,
  "coeffFile" -> "Blank"
|>;
Do[
  If[i < Length[cmdArgs],
    Switch[cmdArgs[[i]],
      "-reps", params["reps"] = cmdArgs[[i + 1]],
      "-Nmax", params["Nmax"] = ToExpression[cmdArgs[[i + 1]]],
      "-output", params["output"] = cmdArgs[[i + 1]],
      "-kernels", params["kernels"] = ToExpression[cmdArgs[[i + 1]]],
      "-coeffFile", params["coeffFile"] = cmdArgs[[i + 1]],
      _, Null
    ]
  ],
  {i, 1, Length[cmdArgs], 2}
];

(* Parse representation list *)
replist = ToExpression[params["reps"]];Nmax=params["Nmax"];coeffFile=params["coeffFile"];kernels=params["kernels"];output=params["output"];

batchLog["Parameters:"];
batchLog["  Representations: ", Length[replist]];
batchLog["  Nmax: ", params["Nmax"]];
batchLog["  Output: ", params["output"]];
batchLog["  Kernels: ", params["kernels"]];
batchLog["  coeffFile: ", params["coeffFile"]];

(*--------------------------------Definitions---------------------------------------*)

(*Routine to find the expansion of the measure divided by ZQ to order t^2NMax*)ClearAll[expandF];
expandF[NMax_Integer]:=Module[{w,Sval,Tval,gcoef,glist,hlist,f0,n,k,j},w={\[Omega][1],\[Omega][2],\[Omega][3],\[Omega][4]};
(*Symmetric power-sum quantities S_n,T_n*)Sval[n_]:=3 Sum[w[[i]]^(2 n),{i,4}]+2 Sum[(w[[i]]-w[[j]])^(2 n)+(w[[i]]+w[[j]])^(2 n),{i,4},{j,i+1,4}];
Tval[n_]:=Sum[(w[[1]]+e2 w[[2]]+e3 w[[3]]+e4 w[[4]])^(2 n),{e2,{-1,1}},{e3,{-1,1}},{e4,{-1,1}}];
(*g_n:coefficients of log(f/f(0))*)gcoef[n_]:=((-1)^n BernoulliB[2 n]/(2 n (2 n)!))*(2^(2 n) Sval[n]-(2^(2 n)-1) Tval[n]);
glist=Table[Expand[gcoef[n]],{n,1,NMax}];
(*Exponentiate via recurrence:h_k=(1/k) Sum_{j=1..k} j g_j h_{k-j}*)hlist=ConstantArray[0,NMax+1];
hlist[[1]]=1;
Do[hlist[[k+1]]=Expand[(1/k) Sum[j glist[[j]] hlist[[k-j+1]],{j,1,k}]],{k,1,NMax}];
(*f(0) prefactor*)f0=Expand[Product[w[[i]]^3,{i,4}]*Product[(w[[i]]^2-w[[j]]^2)^2,{i,4},{j,i+1,4}]/256];
(*Returns coefficient of t^(2k),k=0,...,NMax*)Table[Expand[f0 hlist[[k+1]]],{k,0,NMax}]];
om={\[Omega][1],\[Omega][2],\[Omega][3],\[Omega][4]};

(*----1. Convolution of two even series in t----*)
multiplyEvenSeries[a_List,b_List]:=Module[{Nm,k,j},Nm=Min[Length[a],Length[b]]-1;
Table[Expand@Sum[a[[j+1]] b[[k-j+1]],{j,0,k}],{k,0,Nm}]];

(*Single integral*)Isinh[n_Integer]:=If[OddQ[n],4 n! (1-2^(-(n+1))) Zeta[n+1],0];

(*Integrate one polynomial coefficient over all 4 omega_i*)
integrateCoef[poly_,vars_]:=Module[{rules},rules=CoefficientRules[Expand[poly],vars];
Sum[r[[2]] (Times@@(Isinh/@r[[1]])),{r,rules}]];
Routines to compute characters and their series expansion
ClearAll[weylOrbitB4,weightMultsB4,weightListB4];

weylOrbitB4[mu_List]:=DeleteDuplicates@Flatten[Table[eps perm,{perm,Permutations[mu]},{eps,Tuples[{-1,1},4]}],1];

weightMultsB4[dynkin_List]:=Module[{simpleRoots,posRoots,rho,Lambda,LpSq,LSq,dominantQ,toDom,c1max,c2max,c3max,c4max,candidates,sorted,mults,alpha,k,muCand,muDom,denom,num,mu},simpleRoots={{1,-1,0,0},{0,1,-1,0},{0,0,1,-1},{0,0,0,1}};
posRoots=Join[Flatten[Table[UnitVector[4,i]-UnitVector[4,j],{i,3},{j,i+1,4}],1],Flatten[Table[UnitVector[4,i]+UnitVector[4,j],{i,3},{j,i+1,4}],1],Table[UnitVector[4,i],{i,4}]];
rho={7/2,5/2,3/2,1/2};
Lambda=dynkin[[1]]{1,0,0,0}+dynkin[[2]]{1,1,0,0}+dynkin[[3]]{1,1,1,0}+dynkin[[4]]{1/2,1/2,1/2,1/2};
LpSq=(Lambda+rho) . (Lambda+rho);
LSq=Lambda . Lambda;
dominantQ[m_]:=m[[1]]>=m[[2]]>=m[[3]]>=m[[4]]>=0;
toDom[m_]:=Sort[Abs[m],Greater];
(*Bounds:\[CapitalLambda]\[Minus]\[Mu]=\[CapitalSigma] d_i \[Alpha]_i with 0\[LessEqual]d_i\[LessEqual](\[Alpha]-coord of \[CapitalLambda])_i since for dominant \[Mu] every \[Alpha]-coord of \[Mu] is\[GreaterEqual]0. For \[CapitalLambda]=(\[CapitalLambda]1,\[CapitalLambda]2,\[CapitalLambda]3,\[CapitalLambda]4) the \[Alpha]-coords are the partial sums.*)c1max=Floor[Lambda[[1]]];
c2max=Floor[Lambda[[1]]+Lambda[[2]]];
c3max=Floor[Lambda[[1]]+Lambda[[2]]+Lambda[[3]]];
c4max=Floor[Total[Lambda]];
candidates=DeleteDuplicates@Flatten[Table[With[{m=Lambda-c1 simpleRoots[[1]]-c2 simpleRoots[[2]]-c3 simpleRoots[[3]]-c4 simpleRoots[[4]]},If[dominantQ[m],{m},{}]],{c1,0,c1max},{c2,0,c2max},{c3,0,c3max},{c4,0,c4max}],4];
sorted=SortBy[candidates,-(#+rho) . (#+rho)&];
mults=<||>;mults[Lambda]=1;
Do[If[mu=!=Lambda,denom=LpSq-(mu+rho) . (mu+rho);
num=0;
Do[k=1;
While[muCand=mu+k alpha;muCand . muCand<=LSq,muDom=toDom[muCand];
num+=2 Lookup[mults,Key[muDom],0] muCand . alpha;
k++],{alpha,posRoots}];
mults[mu]=num/denom],{mu,sorted}];
Select[mults,#>0&]];

weightListB4[dynkin_List]:=Flatten[KeyValueMap[Function[{dom,m},{m,#}&/@weylOrbitB4[dom]],weightMultsB4[dynkin]],1];
characterSeries[dynkin_List,KMax_Integer]:=Module[{om=Array[\[Omega],4],weights,coefs,mult,mu,mudotom,k,kmax},weights=weightListB4[dynkin];kmax=Floor[(KMax+1)/2];
coefs=ConstantArray[0,kmax+1];
Do[{mult,mu}=wt;
mudotom=2 mu . om;
(*k=0:(i*mudotom)^0=1 including for the zero weight*)coefs[[1]]+=mult ;
Do[coefs[[k+1]]+=mult  (I mudotom)^(2k)/(2k)!,{k,1,kmax}],{wt,weights}];
Expand/@(ComplexExpand[Re[#]]&/@coefs)];

(*Routine to put together Weyl expansion with the other expansion*)
Weylcoeffs[rep_List,coefs_List,NMax_Integer]:=Module[{weylreps,prod},
weylreps=characterSeries[rep,NMax];
prod=multiplyEvenSeries[coefs,weylreps];(8 2^32)/(384(2\[Pi])^4) integrateCoef[#,om]&/@prod]

(*-----------------------------------Main---------------------------------------------*)

If[coeffFile=="Blank",batchLog["Generating coefficient file coefs"<>ToString[Nmax/2]<>".m"];
coefs=expandF[Nmax];Export["coefs"<>ToString[Nmax/2]<>".m",Compress[coefs]],
coefs=Uncompress[Get[coeffFile]];batchLog["coefficient file loaded"]];

(*Start output file *)Export[output,Compress[{}]];

LaunchKernels[Min[Length[replist],kernels]];
batchLog["Distributing definitions"];
DistributeDefinitions[coefs,replist,Weylcoeffs,Nmax,stream];

batchLog["Starting computation"];
Parallelize[
Do[
temp[i]={{replist[[i]],Weylcoeffs[replist[[i]],coefs,Nmax]}};
rfile=Uncompress[Get[output]];rfile=Join[rfile,temp[i]];
Export[output,Compress[rfile]];
out=ToString[replist[[i]]]<>" completed";
Save[logFile,out];,
{i,Length[replist]}]
];

batchLog["Done!"];
Exit[];

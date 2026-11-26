## Mathematica code

```mathematica
(* ================================================================ *)
(* Cantilever beam topology optimization (SIMP + OC update)         *)
(* Aspect ratio 5:2                                                 *)
(* Marcus Aßmus                                                     *)
(* 2019                                                             *)
(* ================================================================ *)

ClearAll["Global`*"];
tStart = AbsoluteTime[];

(* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *)
(* Part 1: model, assembly, optimization loop                       *)
(* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *)

(* ---------- 1. Parameters ---------- *)
nelx = 100;                    (* elements along x *)
nely = Round[2 nelx/5];        (* elements along y (aspect ratio 5:2) *)
volfrac = 0.4;                 (* target volume fraction *)
penal = 3.;                    (* SIMP penalization *)
rmin = 1.6;                    (* sensitivity filter radius, in elements *)
maxIter = 60;                  (* maximum iterations *)
E0 = 1.;                       (* Young's modulus, solid *)
Emin = 10.^-9;                 (* Young's modulus, void *)
nu = 0.3;                      (* Poisson's ratio *)
move = 0.2;                    (* OC move limit *)
nel = nelx nely;
ndof = 2 (nelx + 1) (nely + 1);

(* ---------- 2. Element stiffness matrix (Q4, plane stress, unit square) ---------- *)
(* Node order (counterclockwise, y up): lower-left, lower-right, upper-right, upper-left *)
kk = {1/2 - nu/6, 1/8 + nu/8, -1/4 - nu/12, -1/8 + 3 nu/8, -1/4 + nu/12, -1/8 - nu/8, nu/6, 1/8 - 3 nu/8};
kkIdx = {{1, 2, 3, 4, 5, 6, 7, 8}, {2, 1, 8, 7, 6, 5, 4, 3}, {3, 8, 1, 6, 7, 4, 5, 2}, {4, 7, 6, 1, 8, 3, 2, 5},
   {5, 6, 7, 8, 1, 2, 3, 4}, {6, 5, 4, 3, 2, 1, 8, 7}, {7, 4, 5, 2, 3, 8, 1, 6}, {8, 3, 2, 5, 4, 7, 6, 1}};
KE = Map[kk[[#]] &, kkIdx]/(1 - nu^2);
KEflat = Flatten[KE];

(* independent check of KE: 2x2 Gauss integration of B^T D B over the unit square *)
Dmat = 1/(1 - nu^2) {{1, nu, 0}, {nu, 1, 0}, {0, 0, (1 - nu)/2}};
gp = 1./Sqrt[3.];
KEgauss = Sum[
   Module[{dN, Bm},
    dN = 1/2 {{-(1 - et), 1 - et, 1 + et, -(1 + et)}, {-(1 - xi), -(1 + xi), 1 + xi, 1 - xi}};
    Bm = {Flatten[Table[{dN[[1, n]], 0}, {n, 4}]],
          Flatten[Table[{0, dN[[2, n]]}, {n, 4}]],
          Flatten[Table[{dN[[2, n]], dN[[1, n]]}, {n, 4}]]};
    Transpose[Bm].Dmat.Bm/4],
   {xi, {-gp, gp}}, {et, {-gp, gp}}];
Print["KE symmetric: ", Max[Abs[KE - Transpose[KE]]] < 10.^-12];
Print["KE vs. Gauss integration, max deviation (should be ~1e-16): ", Max[Abs[KE - KEgauss]]];
Print["KE rigid-rotation residual (should be ~1e-16): ", Max[Abs[KE.{0, 0, 0, 1, -1, 1, -1, 0}]]];

(* ---------- 3. Mesh connectivity and boundary conditions ---------- *)
(* node number = (nely+1)(column-1) + row, row 1 = top edge; element index = (elx-1) nely + ely *)
elIdx = Range[nel];
elx = Quotient[elIdx - 1, nely] + 1;
ely = Mod[elIdx - 1, nely] + 1;
n1 = (nely + 1) (elx - 1) + ely;      (* upper-left node of each element *)
n2 = (nely + 1) elx + ely;            (* upper-right node of each element *)
(* element DOFs in the order lower-left, lower-right, upper-right, upper-left (x, y per node) *)
edofList = Transpose[{2 n1 + 1, 2 n1 + 2, 2 n2 + 1, 2 n2 + 2, 2 n2 - 1, 2 n2, 2 n1 - 1, 2 n1}];
edofFlat = Flatten[edofList];

(* clamped left edge, downward unit load at the right edge, mid-height *)
fixedDofs = Range[2 (nely + 1)];
freeDofs = Complement[Range[ndof], fixedDofs];
Fvec = ConstantArray[0., ndof];
midRows = If[EvenQ[nely], {nely/2 + 1}, {(nely + 1)/2, (nely + 3)/2}];   (* odd nely: split load over 2 nodes *)
loadNodes = (nely + 1) nelx + midRows;
Do[Fvec[[2 nd]] = -0.1/Length[loadNodes], {nd, loadNodes}];

(* rigid-body rotation field of the whole mesh, used only for a check *)
rotVec = Flatten[Table[{-(nely + 1 - r), c}, {c, 0, nelx}, {r, 1, nely + 1}]];

(* ---------- 4. Assembly operators: K = T . blockdiag(E_e KE) . T^T ---------- *)
Tfull = SparseArray[Transpose[{edofFlat, Range[8 nel]}] -> ConstantArray[1., 8 nel], {ndof, 8 nel}];
TfullT = Transpose[Tfull];
Tfree = Tfull[[freeDofs]];
TfreeT = Transpose[Tfree];
bdLoc = Flatten[Table[{a, b}, {a, 8}, {b, 8}], 1];
bdRow = Flatten[Outer[Plus, 8 (Range[nel] - 1), bdLoc[[All, 1]]]];
bdCol = Flatten[Outer[Plus, 8 (Range[nel] - 1), bdLoc[[All, 2]]]];
bdPos = Transpose[{bdRow, bdCol}];
bdOf[Ee_] := SparseArray[bdPos -> Flatten[Outer[Times, Ee, KEflat]], {8 nel, 8 nel}];

feSolve[xd_] := Module[{Ee, Kff, Uv, Ue, ce},
   Ee = Emin + xd^penal (E0 - Emin);
   Kff = Tfree.bdOf[Ee].TfreeT;
   Uv = ConstantArray[0., ndof];
   Uv[[freeDofs]] = LinearSolve[Kff, Fvec[[freeDofs]]];
   Ue = Partition[Uv[[edofFlat]], 8];
   ce = Total[(Ue.KE) Ue, {2}];                (* u_e^T KE u_e per element *)
   {Ee, Uv, ce, Total[Ee ce]}];

(* ---------- 5. Sensitivity filter as a sparse matrix ---------- *)
rOff = Ceiling[rmin] - 1;
offsets = Select[Flatten[Table[{di, dj}, {di, -rOff, rOff}, {dj, -rOff, rOff}], 1],
   Sqrt[#[[1]]^2 + #[[2]]^2] < rmin &];
Hmat = Plus @@ Map[
   Function[o, Module[{i2, j2, ok, ids},
     i2 = elx + o[[1]]; j2 = ely + o[[2]];
     ok = UnitStep[i2 - 1] UnitStep[nelx - i2] UnitStep[j2 - 1] UnitStep[nely - j2];
     ids = Flatten[Position[ok, 1]];
     SparseArray[Transpose[{ids, (i2[[ids]] - 1) nely + j2[[ids]]}] ->
       ConstantArray[rmin - Sqrt[1. o[[1]]^2 + o[[2]]^2], Length[ids]], {nel, nel}]]],
   offsets];
Hsum = Normal[Hmat.ConstantArray[1., nel]];

(* ---------- 6. OC update (bisection on the Lagrange multiplier), exact min/max ---------- *)
ocUpdate = Compile[{{x, _Real, 1}, {dcn, _Real, 1}, {vf, _Real}, {mv, _Real}},
   Module[{l1 = 0., l2 = 1.*^9, lmid = 0., xn = x, n = Length[x]},
    While[(l2 - l1)/(l1 + l2) > 1.*^-3,
     lmid = 0.5 (l2 + l1);
     xn = Table[Max[0.001, Max[x[[i]] - mv, Min[1., Min[x[[i]] + mv, x[[i]] Sqrt[Abs[dcn[[i]]]/lmid]]]]], {i, n}];
     If[Total[xn]/n > vf, l1 = lmid, l2 = lmid]];
    xn]];

(* ---------- 7. Optimization loop ---------- *)
x = ConstantArray[volfrac, nel];
complianceHistory = {};
converged = False;
Do[
  {Ee, Uvec, ce, comp} = feSolve[x];
  If[iter == 1,
   Print["Global rigid-rotation residual (should be ~1e-12 or smaller): ",
    Max[Abs[Normal[Tfull.Normal[bdOf[Ee].Normal[TfullT.rotVec]]]]]];
   Print["Energy check |sum(E_e c_e) - F.U| / compliance (should be ~1e-8 or smaller): ",
    Abs[comp - Fvec.Uvec]/comp]];
  AppendTo[complianceHistory, comp];
  dc = -penal x^(penal - 1) (E0 - Emin) ce;
  dcn = Normal[Hmat.(x dc)]/(Clip[x, {10.^-3, 1.}] Hsum);
  xnew = ocUpdate[x, dcn, volfrac, move];
  change = Max[Abs[xnew - x]];
  x = xnew;
  If[iter <= 3 || Mod[iter, 10] == 0,
   Print["iter ", iter, "  compliance = ", comp, "  volume = ", Mean[x], "  change = ", change]];
  If[change < 0.01 && iter > 15, converged = True; Print["Converged at iteration ", iter]; Break[]],
  {iter, 1, maxIter}];
If[! converged, Print["Stopped at maxIter = ", maxIter, " (change = ", change, "); increase maxIter if needed."]];

(* analysis of the FINAL design: the history values belong to the design before each update *)
{EeFinal, UvecFinal, ceFinal, compFinal} = feSolve[x];
Print["Final design: compliance = ", compFinal, "  volume fraction = ", Mean[x]];
Print["Optimization run time (s): ", AbsoluteTime[] - tStart];

(* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *)
(* Part 2: local bending stiffness EI, plots                        *)
(* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *)

xMat = Partition[x, nely];                 (* nelx x nely; rows = x, columns = y (1 = top) *)
yc = Range[nely];                          (* element-row coordinate, element height 1 *)
eEffMat = Emin + xMat^penal (E0 - Emin);   (* SIMP stiffness, as used in the FE model *)

(* neutral axis of every cross-section = modulus-weighted centroid (not simply mid-height) *)
yNeutral = Normal[eEffMat.yc]/Total[eEffMat, {2}];
(* own inertia 1/12 of a unit square + Steiner term about the neutral axis *)
iLocalMat = 1./12. + (ConstantArray[yc, nelx] - yNeutral)^2;
eiLocal = eEffMat iLocalMat;               (* EI contribution of every element *)
eiTotal = Total[eiLocal, {2}];             (* total cross-sectional EI(x) *)
eiTotalSpan = Transpose[{Range[nelx], eiTotal}];

(* solid reference: sum of (1/12 + y^2) over nely rows = nely^3/12 *)
eiSolidRef = E0 nely^3/12.;
Print["max EI(x) / EI of solid beam: ", Max[eiTotal]/eiSolidRef];
Print["max neutral-axis shift vs. mid-height (0 for a symmetric design): ",
  Max[Abs[yNeutral - (nely + 1)/2.]]];

(* void (rho < 0.2) is masked in the colour plot only *)
mask = UnitStep[xMat - 0.2];
eiField = mask eiLocal;
maxEI = Max[eiField];
eiPlotData = eiField/maxEI - (1. - mask);  (* void cells flagged with -1 *)
eiColorFunc[val_] := If[val < 0, White, ColorData["Rainbow"][0.15 + 0.85 val]];

topologyPlot = Legended[
   ArrayPlot[Transpose[xMat],
    ColorFunction -> (GrayLevel[1 - #] &), ColorFunctionScaling -> False,
    PlotRange -> {0, 1}, Frame -> True, FrameLabel -> {"y", "x"},
    ImageSize -> 580],
   BarLegend[{GrayLevel[1 - #] &, {0, 1}}, ColorFunctionScaling -> False,
    LegendLabel -> Style["\[Rho]"]]];

stiffnessSpanPlot = ListLinePlot[eiTotalSpan,
   PlotRange -> {{1, nelx}, {0, eiSolidRef 1.05}}, Filling -> Axis,
   FillingStyle -> Directive[Opacity[0.2], RGBColor[0.85, 0.25, 0.15]],
   PlotStyle -> Directive[Thick, RGBColor[0.85, 0.25, 0.15]],
   GridLines -> {Automatic, {eiSolidRef}}, GridLinesStyle -> Directive[Dashed, Gray],
   Epilog -> {Gray, Text[Style["Solid beam reference EI", Italic, 10], {nelx 0.78, eiSolidRef 0.94}]},
   Frame -> True, FrameLabel -> {"x", "EI(x)"},
    ImageSize -> 520];
   
convergencePlot = ListLinePlot[complianceHistory, Frame -> True, GridLines -> Automatic,
   FrameLabel -> {"iteration", "compliance (strain energy)"},
   PlotStyle -> Directive[Thick, RGBColor[0.1, 0.35, 0.7]], ImageSize -> 520];

Print["Total run time (s): ", AbsoluteTime[] - tStart];
Column[{topologyPlot, stiffnessSpanPlot, convergencePlot}, Spacings -> 2, Alignment -> Center]

(* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *)
(* Part 3: 3D-STL-Export                                            *)
(* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *)

xMat = Partition[x, nely];   (* density matrix nelx x nely (also defined in part 2) *)
rhoThreshold = 0.5;     (* Schwellwert fuer Festmaterial *)
smoothRadius = 3;       (* Glaettungsradius in Elementeinheiten; 0 = roh *)
beamThickness = 10.;    (* Bauteilstaerke / Extrusionstiefe *)
stlName = "cantilever_volume_3d.stl";
stlDir = If[StringQ[NotebookDirectory[]], NotebookDirectory[], $HomeDirectory];

(* orientierung test: in ImageMesh, is the FIRST row of pixels at y = max (default) or at y = 0?*)
imgTopIsYmax = TrueQ[Quiet[RegionCentroid[ImageMesh[Binarize[Image[{{1., 0.}, {0., 0.}}], 0.5]]][[2]] > 1.]];
Print["ImageMesh: erste Bildzeile bei y = max: ", imgTopIsYmax];

(* 1. create a smooth 2D contour and triangulate it cleanly (row 1 of xMat^T = top edge, y upward) *)
densData = If[imgTopIsYmax, Transpose[xMat], Reverse[Transpose[xMat]]];
densImg = Image[densData, "Real"];
binImg = If[smoothRadius > 0, Binarize[GaussianFilter[densImg, smoothRadius], rhoThreshold], Binarize[densImg, rhoThreshold]];
mesh2D = ImageMesh[binImg, Method -> "Exact"];
tri2D = TriangulateMesh[DiscretizeRegion[mesh2D]];

(* 2. extrude to 3d volumen *)
vol3D = BoundaryDiscretizeRegion[RegionProduct[tri2D, Line[{{0.}, {beamThickness}}]]];

(* 3. proof and export *)
Print["BoundaryMeshRegion: ", Head[vol3D] === BoundaryMeshRegion, " | Volumen STL: ", Volume[vol3D],
  " | Volumen aus Dichten: ", Total[x] beamThickness];
stlPath = FileNameJoin[{stlDir, stlName}];
Export[stlPath, vol3D];
Print["3d-volume-STL saved: ", stlPath];

(* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *)
(* Part 4: FEM load test of the final design                        *)
(* deformed shape colour-coded by bending normal stress             *)
(* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *)

testForce = -0.1;      (* total vertical force; displacements and stresses are linear in it *)
dispScale = 1;     (* magnification of the displacements in the plot *)
rhoCutoff = 0.2;      (* only elements with density >= rhoCutoff are drawn *)

(* 1. Displacements: FE solution of the final design (Part 1), scaled to testForce *)
loadTotal = Total[Fvec[[Range[2, ndof, 2]]]];
Utest = UvecFinal testForce/loadTotal;

(* 2. Strains/stresses at the element centre (unit square, plane stress) *)
(* DOF order per element: LL, LR, UR, UL; B0 = strain-displacement matrix at xi = eta = 0 *)
ueMat = Partition[Utest[[edofFlat]], 8];
B0 = 0.5 {{-1., 0., 1., 0., 1., 0., -1., 0.},
    {0., -1., 0., -1., 0., 1., 0., 1.},
    {-1., -1., -1., 1., 1., 1., 1., -1.}};
epsXX = ueMat.B0[[1]];
epsYY = ueMat.B0[[2]];
eeTest = Emin + x^penal (E0 - Emin);
(* plane stress: sigma_xx = E/(1 - nu^2) (eps_xx + nu eps_yy) *)
sigXXVec = eeTest (epsXX + nu epsYY)/(1. - nu^2);

(* 3. Equilibrium check: bending moment of sigma_xx in every cross-section vs. P (L - x) *)
sigMat = Partition[sigXXVec, nely];                  (* nelx x nely, y index 1 = top *)
armVec = (nely + 1)/2. - Range[nely];                (* lever arm about mid-height *)
Mnum = Normal[sigMat.armVec];
Mref = -testForce (nelx - (Range[nelx] - 0.5));
selCols = Range[Ceiling[0.03 nelx], Floor[0.85 nelx]];
Print["Section equilibrium M_num/M_ref over x = 3%..85% of the span (should be ~1.00): mean = ",
  Mean[(Mnum/Mref)[[selCols]]], ", min = ", Min[(Mnum/Mref)[[selCols]]], ", max = ", Max[(Mnum/Mref)[[selCols]]]];

(* 4. Deformed geometry and quads *)
coords0 = Flatten[Table[{1. col, 1. (nely + 1 - row)}, {col, 0, nelx}, {row, 1, nely + 1}], 1];
coordsDef = coords0 + dispScale Partition[Utest, 2];
quadNodes = Transpose[{n1, n1 + 1, n2 + 1, n2}];     (* UL, LL, LR, UR node of every element *)
visible = Flatten[Position[UnitStep[x - rhoCutoff], 1]];
activeSig = sigXXVec[[visible]];
maxAbsSig = Max[10.^-12, Max[Abs[activeSig]]];

(* symmetric colour scale: -max = blue, 0 = white, +max = red *)
sigColorFunc = Blend[{{-maxAbsSig, Blue}, {0., White}, {maxAbsSig, Red}}, #] &;
elemPolysSigXX = Map[
   {EdgeForm[None], sigColorFunc[sigXXVec[[#]]],
     Polygon[coordsDef[[quadNodes[[#]]]]]} &, visible];

(* 5. Plot with legend *)
deformedSigXXPlot = Legended[
   Graphics[{{EdgeForm[{Dashed, Gray}], FaceForm[None], Rectangle[{0, 0}, {nelx, nely}]},
     elemPolysSigXX, {Black, PointSize[0.015], Point[coordsDef[[loadNodes]]]}},
    AspectRatio -> Automatic, Frame -> True, FrameLabel -> {"x", "y"},
       ImageSize -> 650],
   BarLegend[{sigColorFunc, {-maxAbsSig, maxAbsSig}}, ColorFunctionScaling -> False, LegendLabel -> Style[ToExpression["\\sigma_{xx}", TeXForm], 11]]];

Print["Min / Max sigma_xx (visible elements): ", {Min[activeSig], Max[activeSig]}];
deformedSigXXPlot


```

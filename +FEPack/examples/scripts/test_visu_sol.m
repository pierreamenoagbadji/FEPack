numbersNodes = [5, 10, 20, 30, 40];

set(0,'DefaultFigureWindowStyle','docked')

for idnodes = 1:length(numbersNodes)
  numNodes = numbersNodes(idnodes);

  load(['outputs/inputs_', int2str(numNodes), '.mat']);
  U2D = load(['outputs/U2D_', int2str(numNodes), '.mat']);
  Ur = load(['outputs/Ur_', int2str(numNodes), '.mat']);

  E.positive = U2D.positive - Ur.positive;
  E.negative = U2D.negative - Ur.negative;

  visu_sol(U2D, mesh2Dpos, mesh2Dneg, numCellsInfinite, numCellsSemiInfinite_pos, numCellsSemiInfinite_neg)
  % visu_sol(Ur, mesh2Dpos, mesh2Dneg, numCellsInfinite, numCellsSemiInfinite_pos, numCellsSemiInfinite_neg)
  % visu_sol(E, mesh2Dpos, mesh2Dneg, numCellsInfinite, numCellsSemiInfinite_pos, numCellsSemiInfinite_neg, [-1e-3, 1e-3])
end

%%
load('outputs/inputs_40.mat');
% load('outputs/inputs_copy_10.mat');
U2D_1 = load('outputs/U2D_40.mat');
% U2D_2 = load('outputs/U2D_copy_10_2.mat');
% U2D_3 = load('outputs/U2D_copy_10_3.mat');
visu_sol(U2D_1, mesh2Dpos, mesh2Dneg, numCellsInfinite, numCellsSemiInfinite_pos, numCellsSemiInfinite_neg, [-0.1, 1])
% visu_sol(U2D_2, mesh2Dpos, mesh2Dneg, numCellsInfinite, numCellsSemiInfinite_pos, numCellsSemiInfinite_neg, [-0.1, 1])
% visu_sol(U2D_3, mesh2Dpos, mesh2Dneg, numCellsInfinite, numCellsSemiInfinite_pos, numCellsSemiInfinite_neg, [-0.1, 1])

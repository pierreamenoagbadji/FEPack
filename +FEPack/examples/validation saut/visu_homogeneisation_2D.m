clear; clc;
import FEPack.*

% load output.mat
load output_1.mat
% load output_A_var_2D_homog_0_5.mat
% load output_A_var_2D_homog_0_25.mat
% load output_A_var_2D_homog_0_125.mat

plot_numCellsSemiInfinite_pos = 5;
plot_numCellsSemiInfinite_neg = 5;
plot_numCellsInfinite = 5;

%% Plot U
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

numCells = [1 1];
numCells(coSemiInf) = numCellsSemiInfinite_pos;
numCells(coInf) = 2*numCellsInfinite;
Icells = [1 1];
orientation = 1;
for idS = 1:plot_numCellsSemiInfinite_pos
  for idI = (-plot_numCellsInfinite):(plot_numCellsInfinite-1)
    Icells(coInf) = idI + numCellsInfinite + 1;
    Icells(coSemiInf) = idS;

    idcell = sub2ind(numCells, Icells(1), Icells(2));

    X = mesh_pos.points(:, 1) + (coSemiInf == 1) * (orientation * (idS - 1))*period + (coInf == 1) * idI * period;
    Y = mesh_pos.points(:, 2) + (coSemiInf == 2) * (orientation * (idS - 1))*period + (coInf == 2) * idI * period;

    trisurf(mesh_pos.triangles, X, Y, real(U.positive(:, idcell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    % set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
  idS
end
%
numCells = [1 1];
numCells(coSemiInf) = numCellsSemiInfinite_neg;
numCells(coInf) = 2*numCellsInfinite;
Icells = [1 1];
orientation = -1;
for idS = 1:plot_numCellsSemiInfinite_neg
  for idI = (-plot_numCellsInfinite):(plot_numCellsInfinite-1)
    Icells(coInf) = idI + numCellsInfinite + 1;
    Icells(coSemiInf) = idS;

    idcell = sub2ind(numCells, Icells(1), Icells(2));

    X = mesh_neg.points(:, 1) + (coSemiInf == 1) * (orientation * (idS - 1)) * period + (coInf == 1) * idI * period;
    Y = mesh_neg.points(:, 2) + (coSemiInf == 2) * (orientation * (idS - 1)) * period + (coInf == 2) * idI * period;

    trisurf(mesh_neg.triangles, X, Y, real(U.negative(:, idcell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    % set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
  idS
end


BBS = [-5, 5];
BBI = [-4, 4];
% BBS = [(-plot_numCellsSemiInfinite_neg + 1) * period, (plot_numCellsSemiInfinite_pos - 1) * period];
% BBI = [-plot_numCellsInfinite * period, plot_numCellsInfinite * period];
%
xlim((coInf == 1) * BBI + (coSemiInf == 1) * BBS);
ylim((coInf == 2) * BBI + (coSemiInf == 2) * BBS);
caxis([-1e-2, 1e-2]);


%% Plot Ueff
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

numCells = [1 1];
numCells(coSemiInf) = numCellsSemiInfinite_pos;
numCells(coInf) = 2*numCellsInfinite;
Icells = [1 1];
orientation = 1;
for idS = 1:plot_numCellsSemiInfinite_pos
  for idI = (-plot_numCellsInfinite):(plot_numCellsInfinite-1)
    Icells(coInf) = idI + numCellsInfinite + 1;
    Icells(coSemiInf) = idS;

    idcell = sub2ind(numCells, Icells(1), Icells(2));

    X = mesh_pos.points(:, 1) + (coSemiInf == 1) * (orientation * (idS - 1))*period + (coInf == 1) * idI * period;
    Y = mesh_pos.points(:, 2) + (coSemiInf == 2) * (orientation * (idS - 1))*period + (coInf == 2) * idI * period;

    trisurf(mesh_pos.triangles, X, Y, real(Ur.positive(:, idcell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    % set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
  idS
end

numCells = [1 1];
numCells(coSemiInf) = numCellsSemiInfinite_neg;
numCells(coInf) = 2*numCellsInfinite;
Icells = [1 1];
orientation = -1;
for idS = 1:plot_numCellsSemiInfinite_neg
  for idI = (-plot_numCellsInfinite):(plot_numCellsInfinite-1)
    Icells(coInf) = idI + numCellsInfinite + 1;
    Icells(coSemiInf) = idS;

    idcell = sub2ind(numCells, Icells(1), Icells(2));

    X = mesh_neg.points(:, 1) + (coSemiInf == 1) * (orientation * (idS - 1)) * period + (coInf == 1) * idI * period;
    Y = mesh_neg.points(:, 2) + (coSemiInf == 2) * (orientation * (idS - 1)) * period + (coInf == 2) * idI * period;

    trisurf(mesh_neg.triangles, X, Y, real(Ur.negative(:, idcell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    % set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
  idS
end


BBS = [-5, 5];
BBI = [-4, 4];
% BBS = [(-plot_numCellsSemiInfinite_neg + 1) * period, (plot_numCellsSemiInfinite_pos - 1) * period];
% BBI = [-plot_numCellsInfinite * period, plot_numCellsInfinite * period];
%%
xlim((coInf == 1) * BBI + (coSemiInf == 1) * BBS);
ylim((coInf == 2) * BBI + (coSemiInf == 2) * BBS);
caxis([-5e-3, 5e-3]);


%% Plot the error
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

numCells = [1 1];
numCells(coSemiInf) = numCellsSemiInfinite_pos;
numCells(coInf) = 2*numCellsInfinite;
Icells = [1 1];
orientation = 1;
for idS = 1:plot_numCellsSemiInfinite_pos
  for idI = (-plot_numCellsInfinite):(plot_numCellsInfinite-1)
    Icells(coInf) = idI + numCellsInfinite + 1;
    Icells(coSemiInf) = idS;

    idcell = sub2ind(numCells, Icells(1), Icells(2));

    X = mesh_pos.points(:, 1) + (coSemiInf == 1) * (orientation * (idS - 1))*period + (coInf == 1) * idI * period;
    Y = mesh_pos.points(:, 2) + (coSemiInf == 2) * (orientation * (idS - 1))*period + (coInf == 2) * idI * period;

    trisurf(mesh_pos.triangles, X, Y, real(E.positive(:, idcell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    % set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
  idS
end
%%
numCells = [1 1];
numCells(coSemiInf) = numCellsSemiInfinite_neg;
numCells(coInf) = 2*numCellsInfinite;
Icells = [1 1];
orientation = -1;
for idS = 1:plot_numCellsSemiInfinite_neg
  for idI = (-plot_numCellsInfinite):(plot_numCellsInfinite-1)
    Icells(coInf) = idI + numCellsInfinite + 1;
    Icells(coSemiInf) = idS;

    idcell = sub2ind(numCells, Icells(1), Icells(2));

    X = mesh_neg.points(:, 1) + (coSemiInf == 1) * (orientation * (idS - 1)) * period + (coInf == 1) * idI * period;
    Y = mesh_neg.points(:, 2) + (coSemiInf == 2) * (orientation * (idS - 1)) * period + (coInf == 2) * idI * period;

    trisurf(mesh_neg.triangles, X, Y, real(E.negative(:, idcell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    % set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
  idS
end


BBS = [-5, 5];
BBI = [-4, 4];
% BBS = [(-plot_numCellsSemiInfinite_neg + 1) * period, (plot_numCellsSemiInfinite_pos - 1) * period];
% BBI = [-plot_numCellsInfinite * period, plot_numCellsInfinite * period];
%%
xlim((coInf == 1) * BBI + (coSemiInf == 1) * BBS);
ylim((coInf == 2) * BBI + (coSemiInf == 2) * BBS);
caxis([-5e-3, 5e-3]);

clear; clc;
import FEPack.*

% load output_periodic_16.mat
% load output_periodic_32.mat
load output_periodic_64.mat

%% Plot U
figure;%(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% if (compareU)
%   subplot(1, 2, 1);
% end
for idS = 1:numCellsSemiInfinite_pos
  for idI = 1:2*numCellsInfinite
    Icell = sub2ind([numCellsSemiInfinite_pos, 2*numCellsInfinite], idS, idI);
    X = mesh2Dpos.points(:, 1) + (idS - 1);
    Y = mesh2Dpos.points(:, 2) + (idI - numCellsInfinite - 1);

    trisurf(mesh2Dpos.triangles, X, Y, real(U2D.positive(:, Icell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    % colormap jet
    % caxis([-0.02, 0.02]);
  end
end
%
for idS = 1:numCellsSemiInfinite_neg
  for idI = 1:2*numCellsInfinite
    Icell = sub2ind([numCellsSemiInfinite_neg, 2*numCellsInfinite], idS, idI);
    X = mesh2Dneg.points(:, 1) - (idS - 1);
    Y = mesh2Dneg.points(:, 2) + (idI - numCellsInfinite - 1);
    trisurf(mesh2Dneg.triangles, X, Y, real(U2D.negative(:, Icell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    % colormap jet
  end
end

xlim([-numCellsSemiInfinite_neg + 1, numCellsSemiInfinite_neg - 1]);
ylim([-numCellsInfinite, numCellsInfinite]);
caxis([-0.025, 0.025]);




% %% Plot Ur
% figure;
% for idS = 1:numCellsSemiInfinite_pos
%   for idI = 1:2*numCellsInfinite
%     Icell = sub2ind([numCellsSemiInfinite_pos, 2*numCellsInfinite], idS, idI);
%     X = mesh2Dpos.points(:, 1) + (idS - 1);
%     Y = mesh2Dpos.points(:, 2) + (idI - numCellsInfinite - 1);
% 
%     trisurf(mesh2Dpos.triangles, X, Y, real(Ur.positive(:, Icell)));
%     hold on;
%     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
%     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
%     % colormap jet
%   end
% end
% %
% for idS = 1:numCellsSemiInfinite_neg
%   for idI = 1:2*numCellsInfinite
%     Icell = sub2ind([numCellsSemiInfinite_neg, 2*numCellsInfinite], idS, idI);
%     X = mesh2Dneg.points(:, 1) - (idS - 1);
%     Y = mesh2Dneg.points(:, 2) + (idI - numCellsInfinite - 1);
%     trisurf(mesh2Dneg.triangles, X, Y, real(Ur.negative(:, Icell)));
%     hold on;
%     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
%     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
%     % colormap jet
%   end
% end
% 
% xlim([-numCellsSemiInfinite_neg + 1, numCellsSemiInfinite_pos - 1]);
% ylim([-numCellsInfinite, numCellsInfinite]);
% caxis([-0.025, 0.025]);




%% Plot error
figure;
for idS = 1:numCellsSemiInfinite_pos
  for idI = 1:2*numCellsInfinite
    Icell = sub2ind([numCellsSemiInfinite_pos, 2*numCellsInfinite], idS, idI);
    X = mesh2Dpos.points(:, 1) + (idS - 1);
    Y = mesh2Dpos.points(:, 2) + (idI - numCellsInfinite - 1);

    trisurf(mesh2Dpos.triangles, X, Y, real(E.positive(:, Icell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    % colormap jet
  end
end
%
for idS = 1:numCellsSemiInfinite_neg
  for idI = 1:2*numCellsInfinite
    Icell = sub2ind([numCellsSemiInfinite_neg, 2*numCellsInfinite], idS, idI);
    X = mesh2Dneg.points(:, 1) - (idS - 1);
    Y = mesh2Dneg.points(:, 2) + (idI - numCellsInfinite - 1);
    trisurf(mesh2Dneg.triangles, X, Y, real(E.negative(:, Icell)));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    % colormap jet
  end
end

xlim([-numCellsSemiInfinite_neg + 1, numCellsSemiInfinite_pos - 1]);
ylim([-numCellsInfinite, numCellsInfinite]);
caxis([-1e-2, 1e-2]);


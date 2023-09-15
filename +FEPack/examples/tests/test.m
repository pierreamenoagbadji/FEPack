
clear; clc;
load outputs/U_configB.mat
figure('Position', get(0, 'Screensize'), 'visible', 'off');

X = mesh2Dpos.points(:, 1); Y = mesh2Dpos.points(:, 2);
% trisurf(mesh2Dpos.triangles, X, Y, rho2Dpos([X, Y])); hold on;
% trisurf(mesh2Dpos.triangles, X+1, Y, rho2Dpos([X+1, Y]));
% trisurf(mesh2Dpos.triangles, X, Y-1, rho2Dpos([X, Y-1]));
% trisurf(mesh2Dpos.triangles, X+1, Y-1, rho2Dpos([X+1, Y-1]));
for idX = 0:2
  for idY=-2:1
    trisurf(mesh2Dpos.triangles, X+idX, Y-idY, mu2Dpos([X, Y-idY])); hold on;
  end
end

X = mesh2Dneg.points(:, 1); Y = mesh2Dneg.points(:, 2);
% trisurf(mesh2Dneg.triangles, X, Y, rho2Dneg([X, Y]));
% trisurf(mesh2Dneg.triangles, X-1, Y, rho2Dneg([X-1, Y]));
% trisurf(mesh2Dneg.triangles, X, Y-1, rho2Dneg([X, Y-1]));
% trisurf(mesh2Dneg.triangles, X-1, Y-1, rho2Dneg([X-1, Y-1]));
for idX = -2:0
  for idI=-2:2
    trisurf(mesh2Dneg.triangles, X+idX, Y+idI, mu2Dneg([X, Y+idI])); hold on;
  end
end
% axis([-2 2 -1 1])

% for idS = 1:numCellsSemiInfinite_pos
%   for idI = 1:2*numCellsInfinite_pos
%     Icell = sub2ind([numCellsSemiInfinite_pos, 2*numCellsInfinite_pos], idS, idI);
%     X = mesh2Dpos.points(:, 1) + (idS - 1);
%     Y = mesh2Dpos.points(:, 2) + (idI - numCellsInfinite_pos - 1) * period_pos;

%     trisurf(mesh2Dpos.triangles, X, Y, real(U.positive(:, Icell)));
%     hold on;
%     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
%     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
%   end
% end
% %
% for idS = 1:numCellsSemiInfinite_neg
%   for idI = 1:2*numCellsInfinite_neg
%     Icell = sub2ind([numCellsSemiInfinite_neg, 2*numCellsInfinite_neg], idS, idI);
%     X = mesh2Dneg.points(:, 1) - (idS - 1);
%     Y = mesh2Dneg.points(:, 2) + (idI - numCellsInfinite_neg - 1) * period_neg;
%     trisurf(mesh2Dneg.triangles, X, Y, real(U.negative(:, Icell)));
%     hold on;
%     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
%     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
%   end
% end

% for idS = 1:numCellsSemiInfinite_pos
%   for idI = 1:2*numCellsInfinite
%     Icell = sub2ind([numCellsSemiInfinite_pos, 2*numCellsInfinite], idS, idI);
%     X = mesh2Dpos.points(:, 1) + (idS - 1);
%     Y = mesh2Dpos.points(:, 2) + (idI - numCellsInfinite - 1);

%     trisurf(mesh2Dpos.triangles, X, Y, real(U2D.positive(:, Icell)));
%     hold on;
%     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
%     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
%     % colormap jet
%     % caxis([-0.02, 0.02]);
%   end
% end
% %%
% for idS = 1:numCellsSemiInfinite_neg
%   for idI = 1:2*numCellsInfinite
%     Icell = sub2ind([numCellsSemiInfinite_neg, 2*numCellsInfinite], idS, idI);
%     X = mesh2Dneg.points(:, 1) - (idS - 1);
%     Y = mesh2Dneg.points(:, 2) + (idI - numCellsInfinite - 1);
%     trisurf(mesh2Dneg.triangles, X, Y, real(U2D.negative(:, Icell)));
%     hold on;
%     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
%     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
%     % colormap jet
%   end
% end

%
axis([-3 3 -1 1])
caxis([0.5, 1.5]);
% colormap parula;
view(2);
shading interp;
colorbar off;
set(gca,'DataAspectRatio',[1 1 1]);
axis off;
%%
nomFichier = 'coeff_A_config_B';
print(nomFichier, '-dpng');
% Tronquer l'image construite
system(['convert ', nomFichier, '.png -trim ', nomFichier, '.png']);

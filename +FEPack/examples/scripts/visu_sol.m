function visu_sol(U, mesh2Dpos, mesh2Dneg, numCellsInfinite, numCellsSemiInfinite_pos, numCellsSemiInfinite_neg, clim)
  
  figure;
  set(groot,'defaultAxesTickLabelInterpreter','latex');
  set(groot,'defaulttextinterpreter','latex');
  set(groot,'defaultLegendInterpreter','latex');
  
  %% Positive side
  for idS = 1:numCellsSemiInfinite_pos
    for idI = 1:2*numCellsInfinite
      Icell = sub2ind([numCellsSemiInfinite_pos, 2*numCellsInfinite], idS, idI);
      X = mesh2Dpos.points(:, 1) + (idS - 1);
      Y = mesh2Dpos.points(:, 2) + (idI - numCellsInfinite - 1);

      trisurf(mesh2Dpos.triangles, X, Y, real(U.positive(:, Icell)));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    end
  end
  
  % %% Negative side
  % for idS = 1:numCellsSemiInfinite_neg
  %   for idI = 1:2*numCellsInfinite
  %     Icell = sub2ind([numCellsSemiInfinite_neg, 2*numCellsInfinite], idS, idI);
  %     X = mesh2Dneg.points(:, 1) - (idS - 1);
  %     Y = mesh2Dneg.points(:, 2) + (idI - numCellsInfinite - 1);
  %     trisurf(mesh2Dneg.triangles, X, Y, real(U.negative(:, Icell)));
  %     hold on;
  %     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  %     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  %   end
  % end

  % xlim([-numCellsSemiInfinite_neg + 1, numCellsSemiInfinite_pos - 1]);
  % ylim([-numCellsInfinite, numCellsInfinite]);

  if (nargin >= 7)
    caxis(clim);
  end
  
end
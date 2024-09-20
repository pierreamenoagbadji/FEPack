load([cheminDonnees, '/inputs']);
Usol.positive = cell(numCellsXpos, numCellsZ);
Usol.negative = cell(numCellsXneg, numCellsZ);

for idZ = 1:numCellsZsave
  idZ

  for idX = 1:numCellsXpos
    sol = load([cheminDonnees, '/guide_sol_pos_', int2str(idX), '_', int2str(idZ + (numCellsZ - numCellsZsave)/2)]);
    Usol.positive{idX, idZ} = sol.Usol;
  end

  for idX = 1:numCellsXneg
    sol = load([cheminDonnees, '/guide_sol_neg_', int2str(idX), '_', int2str(idZ + (numCellsZ - numCellsZsave)/2)]);
    Usol.negative{idX, idZ} = sol.Usol;
  end
end

save([cheminDonnees, '/inputs'], '-v7.3');

%% Plot Guide solution
load([cheminDonnees, '/inputs']);
mafig = figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

ZoriginSave = -sizeCellZ * 0.5 * numCellsZsave;

for idX = 1:numCellsXpos
  for idZ = 1:numCellsZsave
    fprintf('%d/%d et %d/%d\n', idX, numCellsXpos, idZ, numCellsZsave);
    cellZorigin = sizeCellZ * (idZ - 1) + ZoriginSave;

    trisurf(mesh_pos.triangles, mesh_pos.points(:, 1) + (idX - 1),...
                                mesh_pos.points(:, 2) + cellZorigin, real(Usol.positive{idX, idZ}));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
end

for idX = 1:numCellsXneg
  for idZ = 1:numCellsZsave
    fprintf('%d/%d et %d/%d\n', idX, numCellsXneg, idZ, numCellsZsave);
    cellZorigin = sizeCellZ * (idZ - 1) + ZoriginSave;

    trisurf(mesh_neg.triangles, mesh_neg.points(:, 1) - (idX - 1),...
                                mesh_neg.points(:, 2) + cellZorigin, real(Usol.negative{idX, idZ}));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
end

% Pour permettre à savefig de fonctionner malgré la taille des fichiers 
rootgroup = settings();
rootgroup.matlab.general.matfile.SaveFormat.PersonalValue = 'v7.3';
savefig(mafig, [cheminDonnees, '/solution'], 'compact');
% print([cheminDonnees, '/solution'], '-dpng');
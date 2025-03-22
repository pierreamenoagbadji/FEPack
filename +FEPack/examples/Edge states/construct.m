numKcell = 4;
numEcell = 4;

figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

hold on;

for idK = 1:numKcell
  for idE = 1:numEcell

    out = load(['outputs/disp_', num2str(idK), '_', num2str(idE)]);

    trisurf(out.dispmesh.triangles, out.dispmesh.points(:, 1), out.dispmesh.points(:, 2), log(abs(out.val))); hold on;

  end
end

% axis([kpars(1) kpars(end) Egies(1) Egies(end)]);
shading interp;
colormap jet;
view(2);
set(gca, 'FontSize', 16);
colorbar('TickLabelInterpreter', 'latex');
xlabel('$k_\parallel$');
ylabel('$E$');
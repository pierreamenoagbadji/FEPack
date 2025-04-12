figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
out = load('outputs/qt.mat');
qt = out.qt;
subplot(1, 2, 1);
qt.visualize_cache% (rootBB);
subplot(1, 2, 2);
qt.visualize% (rootBB);

%%
numKcell = 6;
numEcell = 4;

folder_name = '/';
% folder_name = 'out_1_0/';
% folder_name = 'out_1_1/';
% folder_name = 'out_2_1/';
% folder_name = 'out_3_1/';
% folder_name = 'out_4_1/';

figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

hold on;

for idK = 1:numKcell
  for idE = 1:numEcell
    try
      out = load(['outputs/' folder_name, 'disp_', num2str(idK), '_', num2str(idE)]);

      % val = abs(out.val);
      val = log(abs(out.val));
      val(isinf(val)) = -1;
      trisurf(out.dispmesh.triangles, out.dispmesh.points(:, 1), out.dispmesh.points(:, 2), val); hold on;

      % trisurf(out.dispmesh.triangles, out.dispmesh.points(:, 1), out.dispmesh.points(:, 2), log(abs(out.val))); hold on;
    catch ME
    end
  end
end

% axis([kpars(1) kpars(end) Egies(1) Egies(end)]);
shading interp;
colormap parula;
view(2);
set(gca, 'FontSize', 16);
colorbar('TickLabelInterpreter', 'latex');
xlabel('$k_\parallel$');
ylabel('$E$');
import FEPack.*
clear; clc;
 
% Transmission entre deux milieux constants
omega = 8 + 0.25i;
rho_pos = 1;
rho_neg = 1;

alphaG = 0.5;
G = @(z) 100 * FEPack.tools.cutoff(z, -alphaG, alphaG);
% alpha_G = 1;
% eps_G = 1e-8;
% supp_G = -log(eps_G) / alpha_G;
% Gint = @(x) exp(-alpha_G * x.^2) .* (abs(x) <= supp_G);
% Gint = @(x) FEPack.tools.cutoff(x(:, 1), -0.1, 0.1);
xibound = 30;
zbound = 4;

xi = linspace(-xibound, xibound, floor(100 * (2 * xibound))).';
Delta_xi = xi(2) - xi(1);

% Fourier transform of interface data
ypts = linspace(-alphaG, alphaG, ceil(10*xibound)).';
Delta_y = ypts(2) - ypts(1);
FT_G = Delta_y * exp(-1i * xi * ypts.') * G(ypts) / sqrt(2 * pi);

%% Solution
DtN_pos = sqrt(xi.^2 - rho_pos * omega^2);
DtN_neg = sqrt(xi.^2 - rho_neg * omega^2);
sol_int = FT_G ./ (DtN_pos + DtN_neg);

numNodes = 10;
struct_mesh = 1;
mesh_pos = meshes.MeshRectangle(struct_mesh, [0  1], [0 1], numNodes, numNodes);
mesh_neg = meshes.MeshRectangle(struct_mesh, [0 -1], [0 1], numNodes, numNodes);

exp_DtN_pos = exp(-DtN_pos * mesh_pos.points(:, 1).' + 1i * xi * mesh_pos.points(:, 2).');
exp_DtN_neg = exp( DtN_neg * mesh_neg.points(:, 1).' + 1i * xi * mesh_neg.points(:, 2).');

TFsol_pos = sol_int .* exp_DtN_pos / sqrt(2 * pi);
TFsol_neg = sol_int .* exp_DtN_neg / sqrt(2 * pi);

numCellsXpos = 4;
numCellsXneg = 4;
numCellsY = 4;
Usol.positive = cell(numCellsXpos, 2 * numCellsY);
Usol.negative = cell(numCellsXpos, 2 * numCellsY);

for idY = 1:2*numCellsY
  for idX = 1:numCellsXpos
    poids_pos = Delta_xi * exp(-DtN_pos *  (idX - 1) + 1i * xi * (idY - numCellsY - 1));

    Usol.positive{idX, idY} = poids_pos.' * TFsol_pos;
  end

  for idX = 1:numCellsXneg
    poids_neg = Delta_xi * exp( DtN_neg * -(idX - 1) + 1i * xi * (idY - numCellsY - 1));

    Usol.negative{idX, idY} = poids_neg.' * TFsol_neg;
  end
end

%% Plot solution
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%
for idY = 1:2*numCellsY
  for idX = 1:numCellsXpos
    X = mesh_pos.points(:, 1) + (idX - 1);
    Y = mesh_pos.points(:, 2) + idY - numCellsY - 1;

    trisurf(mesh_pos.triangles, X, Y, real(Usol.positive{idX, idY}));
    hold on;
  end
  
  for idX = 1:numCellsXneg
    X = mesh_neg.points(:, 1) - (idX - 1);
    Y = mesh_neg.points(:, 2) + idY - numCellsY - 1;

    trisurf(mesh_neg.triangles, X, Y, real(Usol.negative{idX, idY}));
    hold on;
  end
  
  view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
end

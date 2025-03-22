%% Import FEPack
clear; clc;
import FEPack.*

%% Problem-related parameters
HcObj = applications.HoneycombObject('none', 'none');

% Honeycomb lattice potentials
centers = [-1/sqrt(3), 1/sqrt(3); 0, 0];
ampsV   = [0;  0];
ampsW   = [10; 10];
rads    = [0.2; 0.2];
HcObj.V = @(x) FEPack.tools.atomicPotential(x, HcObj.vecPer1, HcObj.vecPer2, centers, ampsV, rads);
HcObj.W = @(x) 1 + FEPack.tools.atomicPotential(x, HcObj.vecPer1, HcObj.vecPer2, centers, ampsW, rads);
HcObj.V = @(x) 4.5 - (cos(x(:, 1:2) *  HcObj.dualVec1)...
                   +  cos(x(:, 1:2) *  HcObj.dualVec2)...
                   +  cos(x(:, 1:2) * (HcObj.dualVec1 + HcObj.dualVec2)));

HcObj.W = @(x) sin(x(:, 1:2) *  HcObj.dualVec1)...
             + sin(x(:, 1:2) *  HcObj.dualVec2)...
             - sin(x(:, 1:2) * (HcObj.dualVec1 + HcObj.dualVec2));

% Edge
HcObj.edge.a1 =  1;
HcObj.edge.b1 =  0;

% Domain wall 
HcObj.kappa = @(s) -1 + 2*FEPack.tools.domainwall(s+0.5);
pbinputs.delta = 1;

%% Energy/parallel quasi-momentum range
pbinputs.minKpar = -pi;
pbinputs.maxKpar = +pi;
pbinputs.numKpar = 8;
pbinputs.center_parallel_quasi_momentum_around_K = false;

% Energy range
pbinputs.minEgy = 0; % 17.5460-1e-1;
pbinputs.maxEgy = 25;% 17.5460+1e-1;
pbinputs.numEgy = 8;

%% Meshes and boundary conditions
% # nodes per X/Y edge for periodicity cells
pbinputs.numNodesXpos = 32;
pbinputs.numNodesYpos = 32; 
pbinputs.numNodesXneg = 32;
pbinputs.numNodesYneg = 32;
pbinputs.numNodesXint = 32;
pbinputs.numNodesYint = 32;

% Basis functions for boundary conditions ("Lagrange" or "Fourier")
pbinputs.bc_basis_functions_pos = "Lagrange";
pbinputs.bc_basis_functions_neg = "Lagrange";

%% Misc. 
% Should you use parallel computing?
opts.parallel_use = false;
opts.parallel_numcores = 4;

% plot the coefficients
opts.plot_coefficients = true;

% Plot dispersion functions, save them
opts.plot_disp = false;
opts.save_disp = true;

numKcell = 2;
numEcell = 2;

kpars = linspace(pbinputs.minKpar, pbinputs.maxKpar, numKcell+1);
Egies = linspace(pbinputs.minEgy,  pbinputs.maxEgy,  numEcell+1);

Kstep = kpars(2) - kpars(1);
Estep = Egies(2) - Egies(1);

if (opts.plot_disp)
  figure;
  set(groot,'defaultAxesTickLabelInterpreter','latex');
  set(groot,'defaulttextinterpreter','latex');
  set(groot,'defaultLegendInterpreter','latex');
  
  hold on;
end

for idK = 1:numKcell
  for idE = 1:numEcell

    clc;

    pbinputs.minKpar = kpars(idK);
    pbinputs.maxKpar = kpars(idK+1);
    pbinputs.minEgy  = Egies(idE);
    pbinputs.maxEgy  = Egies(idE+1);

    %% Computations
    [dispmesh, val] = launch_edge_states_rational(pbinputs, HcObj, opts);

    %% Plot the result
    if (opts.plot_disp)
      trisurf(dispmesh.triangles, dispmesh.points(:, 1), dispmesh.points(:, 2), log(abs(val))); hold on;
      % plot3( HcObj.highSymK'*HcObj.vecPer1*[1 1], [pbinputs.minEgy, pbinputs.maxEgy], [1e8 1e8], 'w--');
      % plot3(-HcObj.highSymK'*HcObj.vecPer1*[1 1], [pbinputs.minEgy, pbinputs.maxEgy], [1e8 1e8], 'w--');
    end

    %% Save the result
    if (opts.save_disp)
      save(['outputs/disp_', num2str(idK), '_', num2str(idE)], 'dispmesh', 'val', '-v7.3');
    end
    
  end
end

if (opts.plot_disp)
  axis([kpars(1) kpars(end) Egies(1) Egies(end)]);
  shading interp;
  colormap jet;
  view(2);
  set(gca, 'FontSize', 16);
  colorbar('TickLabelInterpreter', 'latex');
  xlabel('$k_\parallel$');
  ylabel('$E$');
end


%% Function
function [dispmesh, val] = launch_edge_states_rational(pbinputs, HcObj, opts)

  BBkpar = [pbinputs.minKpar, pbinputs.maxKpar];
  BBlambda = [pbinputs.minEgy, pbinputs.maxEgy];

  %% Operators
  % Make sure the edge coefficients are coprime integers
  a1 = HcObj.edge.a1;
  b1 = HcObj.edge.b1;
  [G, x, y] = gcd(a1, b1);  % a1*x + b1*y = G
  a2 = -y; b2 = x;

  if (G ~= 1)
    error('a1 et b1 doivent être premiers entre eux.');
  end

  % Edge vectors
  edge_vec1      =  a1 * HcObj.vecPer1  + b1 * HcObj.vecPer2;
  edge_vec2      =  a2 * HcObj.vecPer1  + b2 * HcObj.vecPer2;
  edge_dual_vec1 =  b2 * HcObj.dualVec1 - a2 * HcObj.dualVec2;
  edge_dual_vec2 = -b1 * HcObj.dualVec1 + a1 * HcObj.dualVec2;

  % Transformation matrices
  Rmat = [edge_vec1, edge_vec2];
  Tmat = [edge_dual_vec1'; edge_dual_vec2'] / (2*pi); % Inverse of Rmat
  Tvec = Tmat' * [1; 0];

  % Potentials
  Vpot = @(x) HcObj.V((Rmat * x(:, 1:2)')');
  Wpot = @(x) HcObj.W((Rmat * x(:, 1:2)')');

  % (+) half-guide
  op_pos.Tmat = Tmat;
  op_pos.vect = Tvec;
  op_pos.funQ = @(x) Vpot(x) + pbinputs.delta * Wpot(x);
  op_pos.funP = 1;
  op_pos.funR = 1;

  % (-) half-guide
  op_neg.Tmat = Tmat;
  op_neg.vect = Tvec;
  op_neg.funQ = @(x) Vpot(x) + pbinputs.delta * Wpot(x);
  op_neg.funP = 1;
  op_neg.funR = 1;

  % Interior domain
  yInt = ceil(1/(2*pi*pbinputs.delta));
  BBint = [-yInt, yInt];
  op_int.Tmat = Tmat;
  op_int.vect = Tvec;
  % op_int.funQ = @(x) Vpot(x) + pbinputs.delta * HcObj.kappa(2*pi*pbinputs.delta * x(:, 2)) .* Wpot(x);
  op_int.funQ = Wpot;
  op_int.funP = 1;
  op_int.funR = 1;

  % Launch it!
  [dispmesh, val] = edge_states_rational(...
    BBkpar, pbinputs.numKpar,...
    BBlambda, pbinputs.numEgy,...
    op_pos, pbinputs.numNodesXpos, pbinputs.numNodesYpos, pbinputs.bc_basis_functions_pos,...
    op_neg, pbinputs.numNodesXneg, pbinputs.numNodesYneg, pbinputs.bc_basis_functions_neg,...
    op_int, BBint, pbinputs.numNodesXint, pbinputs.numNodesYint,...
    opts.parallel_use, opts.parallel_numcores, opts.plot_coefficients);

end
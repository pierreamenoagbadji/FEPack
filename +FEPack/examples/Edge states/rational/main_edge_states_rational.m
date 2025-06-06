%% Import FEPack
clear; clc;
import FEPack.*

%% Problem-related parameters
HcObj = applications.HoneycombObject('none', 'none');

% Honeycomb lattice potentials
% centers = [-1/sqrt(3), 1/sqrt(3); 0, 0];
% ampsV   = [0;  0];
% ampsW   = [10; 10];
% rads    = [0.2; 0.2];
% HcObj.V = @(x) FEPack.tools.atomicPotential(x, HcObj.vecPer1, HcObj.vecPer2, centers, ampsV, rads);
% HcObj.W = @(x) 1 + FEPack.tools.atomicPotential(x, HcObj.vecPer1, HcObj.vecPer2, centers, ampsW, rads);

HcObj.V = @(x) 10 * cos(x(:, 1:2) *  HcObj.dualVec1) +...
               10 * cos(x(:, 1:2) *  HcObj.dualVec2) +...
               10 * cos(x(:, 1:2) * (HcObj.dualVec1  + HcObj.dualVec2));

HcObj.W = @(x)      sin(x(:, 1:2) *  HcObj.dualVec1) +...
                    sin(x(:, 1:2) *  HcObj.dualVec2) +...
                    sin(x(:, 1:2) * (HcObj.dualVec1  + HcObj.dualVec2));

% Edge
HcObj.edge.a1 =  1;
HcObj.edge.b1 =  2;

% Problem type ('interior' or 'interface')
pbinputs.problem_type = 'interface';

% Domain wall 
HcObj.kappa = @(s) -1 + 2*FEPack.tools.domainwall(s+0.5);
pbinputs.delta = 10;

%% Energy/parallel quasi-momentum range
pbinputs.minKpar = -pi;
pbinputs.maxKpar = +pi;
pbinputs.numKpar = 16;
pbinputs.center_parallel_quasi_momentum_around_K = false;

% Energy range
pbinputs.minEgy = 11.7341-10; % 17.5460-1e-1;
pbinputs.maxEgy = 11.7341+10;% 17.5460+1e-1;
pbinputs.numEgy = 16;

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
opts.parallel_numcores = 2;

% plot the coefficients
opts.plot_coefficients = false;

% Plot dispersion functions, save them
opts.plot_disp = true;
opts.save_disp = true;

%% Initialization
[mesh_pos, mesh_neg, mesh_int, op_pos, op_neg, op_int, BCstruct_pos, BCstruct_neg] = construct_mesh_and_FEmatrices(pbinputs, HcObj);

%% Compute indicator for edge state curves
rootBB = [pbinputs.minKpar, pbinputs.minEgy,...
          pbinputs.maxKpar, pbinputs.maxEgy];
maxDepth = 6;
force_depth = 2;
force_depth = min(force_depth, maxDepth);
threshold = 1e2;

on_band_border = @(fvals) (any(isnan(fvals), 2) & ~all(isnan(fvals), 2));
on_edge_border = @(fvals) (any(abs(fvals) > threshold, 2) & ~all(abs(fvals) > threshold, 2));
early_split    = @(depth) (depth <= force_depth);

% The rule: split leaves on the border of continuous spectrum XOR
% min depth hasn't been reached yet or the indicator function is 
% above a threshold.
rules = @(fvals, depth) ( on_band_border(fvals)) |...
                        (~on_band_border(fvals)  & (early_split(depth) | on_edge_border(fvals)));

qt = FEPack.meshes.Quadtree(maxDepth, rootBB, rules);

indicator_function = @(P) edge_states_rational(...
      P,...
      op_pos, mesh_pos, BCstruct_pos,...
      op_neg, mesh_neg, BCstruct_neg,...
      op_int, mesh_int, pbinputs.problem_type,...
      opts.parallel_use, opts.parallel_numcores, opts.plot_coefficients);

qt = qt.refine(indicator_function);

% Plot the result
if (opts.plot_disp)
  figure;
  set(groot,'defaultAxesTickLabelInterpreter','latex');
  set(groot,'defaulttextinterpreter','latex');
  set(groot,'defaultLegendInterpreter','latex');
  
  subplot(1, 2, 1);
  qt.visualize_cache(rootBB);
  subplot(1, 2, 2);
  qt.visualize(rootBB);
end

% Save the result
if (opts.save_disp)
  save('outputs/qt', 'qt', '-v7.3');
end


% %% Compute indicator for edge state curves
% numKcell = 2;
% numEcell = 2;

% kpars = linspace(pbinputs.minKpar, pbinputs.maxKpar, numKcell+1);
% Egies = linspace(pbinputs.minEgy,  pbinputs.maxEgy,  numEcell+1);

% Kstep = kpars(2) - kpars(1);
% Estep = Egies(2) - Egies(1);

% if (opts.plot_disp)
%   figure;
%   set(groot,'defaultAxesTickLabelInterpreter','latex');
%   set(groot,'defaulttextinterpreter','latex');
%   set(groot,'defaultLegendInterpreter','latex');
  
%   hold on;
% end

% for idK = 1:numKcell
%   for idE = 1:numEcell

%     clc;

%     pbinputs.minKpar = kpars(idK);
%     pbinputs.maxKpar = kpars(idK+1);
%     pbinputs.minEgy  = Egies(idE);
%     pbinputs.maxEgy  = Egies(idE+1);

%     BBkpar = [pbinputs.minKpar, pbinputs.maxKpar];
%     BBEgy  = [pbinputs.minEgy,  pbinputs.maxEgy];
%     dispmesh = FEPack.meshes.MeshRectangle(0, BBkpar, BBEgy, pbinputs.numKpar, pbinputs.numEgy);

%     %% Computations
%     val = edge_states_rational(...
%       dispmesh.points(:, 1:2),...
%       op_pos, mesh_pos, BCstruct_pos,...
%       op_neg, mesh_neg, BCstruct_neg,...
%       op_int, mesh_int, pbinputs.problem_type,...
%       opts.parallel_use, opts.parallel_numcores, opts.plot_coefficients);

%     %% Plot the result
%     if (opts.plot_disp)
%       trisurf(dispmesh.triangles, dispmesh.points(:, 1), dispmesh.points(:, 2), log(abs(val))); hold on;
%       % plot3( HcObj.highSymK'*HcObj.vecPer1*[1 1], [pbinputs.minEgy, pbinputs.maxEgy], [1e8 1e8], 'w--');
%       % plot3(-HcObj.highSymK'*HcObj.vecPer1*[1 1], [pbinputs.minEgy, pbinputs.maxEgy], [1e8 1e8], 'w--');
%     end

%     %% Save the result
%     if (opts.save_disp)
%       save(['outputs/disp_', num2str(idK), '_', num2str(idE)], 'dispmesh', 'val', '-v7.3');
%     end
    
%   end
% end

% if (opts.plot_disp)
%   axis([kpars(1) kpars(end) Egies(1) Egies(end)]);
%   shading interp;
%   colormap jet;
%   view(2);
%   set(gca, 'FontSize', 16);
%   colorbar('TickLabelInterpreter', 'latex');
%   xlabel('$k_\parallel$');
%   ylabel('$E$');
% end


%% Auxiliary functions
function [mesh_pos, mesh_neg, mesh_int, op_pos, op_neg, op_int, BCstruct_pos, BCstruct_neg] = construct_mesh_and_FEmatrices(pbinputs, HcObj)
  
  %% Coefficients
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
  op_neg.funQ = @(x) Vpot(x) - pbinputs.delta * Wpot(x);
  op_neg.funP = 1;
  op_neg.funR = 1;

  % Interior domain
  if strcmpi(pbinputs.problem_type, 'interior')
    yInt = ceil(1/(2*pi*pbinputs.delta));
    BBint = [-yInt, yInt];
    op_int.Tmat = Tmat;
    op_int.vect = Tvec;
    op_int.funQ = @(x) Vpot(x) + pbinputs.delta * HcObj.kappa(2*pi*pbinputs.delta * x(:, 2)) .* Wpot(x);
    op_int.funP = 1;
    op_int.funR = 1;
  else
    op_int = [];
    BBint = [];
  end

  %% Initialization
  [mesh_pos, mesh_neg, mesh_int, op_pos, op_neg, op_int, BCstruct_pos, BCstruct_neg] = initialize_edge_states_rational(...
    op_pos, pbinputs.numNodesXpos, pbinputs.numNodesYpos, pbinputs.bc_basis_functions_pos,...
    op_neg, pbinputs.numNodesXneg, pbinputs.numNodesYneg, pbinputs.bc_basis_functions_neg,...
    op_int, BBint, pbinputs.numNodesXint, pbinputs.numNodesYint,...
    pbinputs.problem_type...
  );

end

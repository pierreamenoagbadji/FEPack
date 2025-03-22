%% Import FEPack
clear; clc;
import FEPack.*

%% Should you use parallel computing?
parallel.use = false;
parallel.numcores = 4;

%% Problem-related parameters
HcObj = applications.HoneycombObject('none', 'none');

% Honeycomb lattice potentials
centers = [-1/sqrt(3), 1/sqrt(3); 0, 0];
ampsV   = [0;  0];
ampsW   = [1; -1];
rads    = [0.2; 0.2];
HcObj.V = @(x) FEPack.tools.atomicPotential(x, HcObj.vecPer1, HcObj.vecPer2, centers, ampsV, rads);
HcObj.W = @(x) FEPack.tools.atomicPotential(x, HcObj.vecPer1, HcObj.vecPer2, centers, ampsW, rads);

% Edge
HcObj.edge.a1 =  1;
HcObj.edge.b1 =  0;

% Domain wall 
HcObj.kappa = @(s) -1 + 2*FEPack.tools.domainwall(s+0.5);

% Small parameter
argspde.delta = 5;

% Parallel quasi-momentum range
argspde.minKpar = -0.01*argspde.delta;
argspde.maxKpar = +0.01*argspde.delta;
argspde.numKpar = 32;
argspde.center_parallel_quasi_momentum_around_K = true;

% Energy range
% If the endpoints of the interval are not specified,
% the code will chose ED + (-upsilon*delta, upsilon*delta)
argspde.minEgy = [];
argspde.maxEgy = [];
argspde.numEgy = 64;

%% Computation of Dirac points
% dirac.numNodesKx = 16;
% dirac.numNodesKy = 16;
dirac.numEigvals = 04; % Number of dispersion surfaces
dirac.numNodesX  = 64;
dirac.numNodesY  = 64;
dirac.numNodesK  = 512;

%% Meshes and boundary conditions
% # nodes per X/Y edge for periodicity cells
argspde.numNodesXpos = 64;
argspde.numNodesYpos = 64; 
argspde.numNodesXneg = 32;
argspde.numNodesYneg = 32;
argspde.numNodesXint = 32;
argspde.numNodesYint = 32;

% Basis functions for boundary conditions ("Lagrange" or "Fourier")
argspde.bc_basis_functions_pos = "Lagrange";
argspde.bc_basis_functions_neg = "Lagrange";

%% Tasks options
taskset.plot_coefficients = false;
taskset.spectrum_pert_bulk.compute = false;

%% Launch it!
[kpars, Egies, val] = edge_states_rational(parallel, HcObj, dirac, argspde, taskset);

%%
[X, Y] = meshgrid(kpars, Egies);
surf(X, Y, abs(val));
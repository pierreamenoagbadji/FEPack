load([cheminDonnees, '/inputs_', int2str(numNodes), '.mat']);

%% Positive half-space DtN
% ////////////////////////
TFBphiPosVec = zeros(mesh3Dpos.domain('xmin').numPoints, numBasis, numFloquetPoints);
for idI = 1:numBasis
  TFBphiPosVec(:, idI, :) = BlochTransform(mesh3Dpos.points(mesh3Dpos.domain('xmin').IdPoints, :),...
                                           FloquetPoints_pos, @(x) phis3Dpos(x, idI),...
                                           2, period, 1000);
end

% %%
% meshDom = mesh3Dpos.restrictToDomain(mesh3Dpos.domain('xmin'));
% for idI = 1:numBasis
%   trisurf(meshDom.triangles, meshDom.points(:, 2), meshDom.points(:, 3), real(TFBphiPosVec(:, idI, 2)));
%   view(2);
%   shading interp
%   set(gca, 'DataAspectRatio', [1 1 1])
%   pause;
% end
% %%

TFBphiPos = cell(numFloquetPoints, 1);
lambda_pos = zeros(numBasis);
wpos = (2*pi / period) / (numFloquetPoints - 1);

for idFB = 1:numFloquetPoints
  TFBphiPos{idFB} = BCstruct_pos.spB0.FE_to_spectral * TFBphiPosVec(:, :, idFB);
  
  TFBlambdaPos = load([cheminDonnees, '/TFBlambdaPos_', int2str(idFB), '.mat']);
  TFBlambdaPos = BCstruct_pos.spB0.massmat * TFBlambdaPos.TFBlambdaPos;

  lambda_pos = lambda_pos + wpos * TFBphiPos{idFB}' * TFBlambdaPos * TFBphiPos{idFB};
end

lambda_pos = lambda_pos / cutvecpos(1);

%% Negative half-space DtN
% ////////////////////////
TFBphiNegVec = zeros(mesh3Dneg.domain('xmin').numPoints, numBasis, numFloquetPoints);
for idI = 1:numBasis
  TFBphiNegVec(:, idI, :) = BlochTransform(mesh3Dneg.points(mesh3Dneg.domain('xmin').IdPoints, :),...
                                           FloquetPoints_neg, @(x) phis3Dneg(x, idI),...
                                           2, period, 1000);
end

TFBphiNeg = cell(numFloquetPoints, 1);
lambda_neg = zeros(numBasis);
wneg = (2*pi / period) / (numFloquetPoints - 1);

for idFB = 1:numFloquetPoints
  TFBphiNeg{idFB} = BCstruct_neg.spB0.FE_to_spectral * TFBphiNegVec(:, :, idFB);

  TFBlambdaNeg = load([cheminDonnees, '/TFBlambdaNeg_', int2str(idFB), '.mat']);
  TFBlambdaNeg = BCstruct_neg.spB0.massmat * TFBlambdaNeg.TFBlambdaNeg;

  lambda_neg = lambda_neg + wneg * TFBphiNeg{idFB}' * TFBlambdaNeg * TFBphiNeg{idFB};
end

lambda_neg = lambda_neg / cutvecneg(1);

%% Solve the integral equation on the interface
% /////////////////////////////////////////////
fprintf('Résolution équation interface\n');
numNodesEF = mesh2Dint.numPoints * mesh2Dtrs.numPoints;
[idIntEF, idTrsEF] = ind2sub([mesh2Dint.numPoints, mesh2Dtrs.numPoints], 1:numNodesEF);

projmat = spBint.projmat(idInt, idIntEF) .* spBtrs.projmat(idTrs, idTrsEF);
Gvec = Gint(mesh2Dint.points) * Gtrs(mesh2Dtrs.points).';
Gvec = Gvec(:);
mat_G_v = projmat * Gvec;

tic;
% -------------------------------------------------------------------- %
% The minus sign comes from the definition of the Lambda when they are %
% the DtN operators (see PeriodicHalfGuideBVP.m)                       %
% -------------------------------------------------------------------- %
trace_solution = -(lambda_pos + lambda_neg) \ mat_G_v;             %
% -------------------------------------------------------------------- %
toc;

%% Deduce the FB transform of the solution
% ////////////////////////////////////////
% FB transform of the solution's trace with period_pos
TFB_solution_pos = cell(numFloquetPoints, 1);
for idFB = 1:numFloquetPoints

  sol_pos_data = load([cheminDonnees, '/sol_pos_data_', int2str(idFB), '.mat']);

  TFB_solution_pos{idFB} = zeros(mesh3Dpos.numPoints, numCellsSemiInfinite_pos);
  R0Phi = TFBphiPos{idFB} * trace_solution;
  R1Phi = sol_pos_data.D * R0Phi;

  for idCell = 0:numCellsSemiInfinite_pos-1
    % Compute the solution in the current cell
    TFB_solution_pos{idFB}(:, idCell + 1) = sol_pos_data.E0 * R0Phi + sol_pos_data.E1 * R1Phi;

    % Update
    R0Phi = sol_pos_data.R * R0Phi;
    R1Phi = sol_pos_data.D * R0Phi;
  end

end

% FB transform of the solution's trace with period_neg
TFB_solution_neg = cell(numFloquetPoints, 1);
for idFB = 1:numFloquetPoints
  
  sol_neg_data = load([cheminDonnees, '/sol_neg_data_', int2str(idFB), '.mat']);

  TFB_solution_neg{idFB} = zeros(mesh3Dneg.numPoints, numCellsSemiInfinite_neg);
  R0Phi = TFBphiNeg{idFB} * trace_solution;
  R1Phi = sol_neg_data.D * R0Phi;

  for idCell = 0:numCellsSemiInfinite_neg-1
    % Compute the solution in the current cell
    TFB_solution_neg{idFB}(:, idCell + 1) = sol_neg_data.E0 * R0Phi + sol_neg_data.E1 * R1Phi;

    % Update
    R0Phi = sol_neg_data.R * R0Phi;
    R1Phi = sol_neg_data.D * R0Phi;
  end

end

%% Take the inverse Floquet transform: positive side
% //////////////////////////////////////////////////
numCells = [1 1 1];
numCells(infiniteDirection) = 2*numCellsInfinite;
numCells(semiInfiniteDirection) = numCellsSemiInfinite_pos;
Nu = prod(numCells);
[I1, I2, I3] = ind2sub(numCells, 1:Nu);
pointsIds = [I1; I2; I3]; % 3-by-Nu
tau = pointsIds(infiniteDirection, :) - numCellsInfinite' * ones(1, Nu) - 1; % Ni-by-Nu
W = prod((2*pi/period) ./ (numFloquetPoints - 1)); % 1-by-1
U3D.positive = zeros(mesh3Dpos.numPoints, Nu);

for idFB = 1:numFloquetPoints
  FloquetVar = FloquetPoints_pos(idFB);
  
  % The integral that defines the inverse Floquet-Bloch transform is computed
  % using a rectangular rule.
  exp_k_dot_x = exp(1i * mesh3Dpos.points(:, infiniteDirection) * FloquetVar); % N-by-1
  exp_k_dot_tau = exp(1i * FloquetVar * tau * period); % 1-by-Nu
  U_TFB = TFB_solution_pos{idFB}(:, pointsIds(semiInfiniteDirection, :)); % N-by-Nu

  U3D.positive = U3D.positive + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
end

U3D.positive = U3D.positive * sqrt(period / (2*pi));

%% Take the inverse Floquet transform: negative side
% //////////////////////////////////////////////////
numCells = [1 1 1];
numCells(infiniteDirection) = 2*numCellsInfinite;
numCells(semiInfiniteDirection) = numCellsSemiInfinite_neg;
Nu = prod(numCells);
[I1, I2, I3] = ind2sub(numCells, 1:Nu);
pointsIds = [I1; I2; I3]; % 3-by-Nu
tau = pointsIds(infiniteDirection, :) - numCellsInfinite' * ones(1, Nu) - 1; % Ni-by-Nu
W = prod((2*pi/period) ./ (numFloquetPoints - 1)); % 1-by-1
U3D.negative = zeros(mesh3Dneg.numPoints, Nu);

for idFB = 1:numFloquetPoints
  FloquetVar = FloquetPoints_neg(idFB);
  
  % The integral that defines the inverse Floquet-Bloch transform is computed
  % using a rectangular rule.
  exp_k_dot_x = exp(1i * mesh3Dneg.points(:, infiniteDirection) * FloquetVar); % N-by-1
  exp_k_dot_tau = exp(1i * FloquetVar * tau * period); % 1-by-Nu
  U_TFB = TFB_solution_neg{idFB}(:, pointsIds(semiInfiniteDirection, :)); % N-by-Nu

  U3D.negative = U3D.negative + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
end

U3D.negative = U3D.negative * sqrt(period / (2*pi));

%% Take the trace: Positive side
% //////////////////////////////
N2Dpos = mesh2Dpos.numPoints;
U2D.positive = zeros(N2Dpos, size(U3D.positive, 2));
dom = mesh3Dpos.domain('volumic');
for idI = 1:2*numCellsInfinite
  IcellY = (numCellsSemiInfinite_pos*(idI-1)+1):(numCellsSemiInfinite_pos*idI);
  X = mesh2Dpos.points(:, 1);
  Y = mesh2Dpos.points(:, 2); % ones(mesh2Dpos.numPoints, 1);
  Z = FEPack.tools.mymod(cutslopePos * (Y + idI - numCellsInfinite - 1));
  structLoc = dom.locateInDomain([X, Y, Z]);

  elts = dom.elements(structLoc.elements, :);
  elts = elts'; elts = elts(:);
  coos = structLoc.barycoos;
  coos = coos'; coos = coos(:);

  U2D.positive(:, IcellY) = reshape(sum(reshape(coos .* U3D.positive(elts, IcellY), dom.dimension+1, []), 1), N2Dpos, []);
end

%% Take the trace: Negative side
% //////////////////////////////
N2Dneg = mesh2Dneg.numPoints;
U2D.negative = zeros(N2Dneg, size(U3D.negative, 2));
dom = mesh3Dneg.domain('volumic');
for idI = 1:2*numCellsInfinite
  IcellY = (numCellsSemiInfinite_neg*(idI-1)+1):(numCellsSemiInfinite_neg*idI);
  X = mesh2Dneg.points(:, 1);
  Y = mesh2Dneg.points(:, 2); % ones(mesh2Dneg.numPoints, 1);
  Z = FEPack.tools.mymod(cutslopeNeg * (Y + idI - numCellsInfinite - 1));
  structLoc = dom.locateInDomain([X, Y, Z]);

  elts = dom.elements(structLoc.elements, :);
  elts = elts'; elts = elts(:);
  coos = structLoc.barycoos;
  coos = coos'; coos = coos(:);

  U2D.negative(:, IcellY) = reshape(sum(reshape(coos .* U3D.negative(elts, IcellY), dom.dimension+1, []), 1), N2Dneg, []);
end

% profile viewer
% profile OFF

%% Plot U
% ///////
mafig = figure;%(1);
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
%%
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

xlim([-numCellsSemiInfinite_neg, numCellsSemiInfinite_pos]);
ylim([-numCellsInfinite, numCellsInfinite]);


savefig(mafig, [cheminDonnees, '/solution_', int2str(numNodes)], 'compact');

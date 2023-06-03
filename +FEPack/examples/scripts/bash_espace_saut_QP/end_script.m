%

load('inputs.mat');

%% Get the Floquet-Bloch transform of the solution
TFBU = cell(numFloquetPoints, 1);

for idFB = 1:numFloquetPoints
  
  TFBU{idFB} = load(['outputs/TFBU_', int2str(idFB), '.mat']);

end

%% Take the inverse Floquet transform: positive side
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
  FloquetVar = FloquetPoints(idFB);
  
  % The integral that defines the inverse Floquet-Bloch transform is computed
  % using a rectangular rule.
  exp_k_dot_x = exp(1i * mesh3Dpos.points(:, infiniteDirection) * FloquetVar); % N-by-1
  exp_k_dot_tau = exp(1i * FloquetVar * tau * period); % 1-by-Nu
  U_TFB = TFBU{idFB}.positive(:, pointsIds(semiInfiniteDirection, :)); % N-by-Nu

  U3D.positive = U3D.positive + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
end

U3D.positive = U3D.positive * sqrt(period / (2*pi));

%% Take the inverse Floquet transform: negative side
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
  FloquetVar = FloquetPoints(idFB);
  
  % The integral that defines the inverse Floquet-Bloch transform is computed
  % using a rectangular rule.
  exp_k_dot_x = exp(1i * mesh3Dneg.points(:, infiniteDirection) * FloquetVar); % N-by-1
  exp_k_dot_tau = exp(1i * FloquetVar * tau * period); % 1-by-Nu
  U_TFB = TFBU{idFB}.negative(:, pointsIds(semiInfiniteDirection, :)); % N-by-Nu

  U3D.negative = U3D.negative + W * (exp_k_dot_x * exp_k_dot_tau) .* U_TFB;
end

U3D.negative = U3D.negative * sqrt(period / (2*pi));

%% Take the trace: Positive side
N2Dpos = mesh2Dpos.numPoints;
U2D.positive = zeros(N2Dpos, size(U3D.positive, 2));
dom = mesh3Dpos.domain('volumic');
for idI = 1:2*numCellsInfinite
  IcellY = (numCellsSemiInfinite_pos*(idI-1)+1):(numCellsSemiInfinite_pos*idI);
  X = mesh2Dpos.points(:, 1);
  Y = mesh2Dpos.points(:, 2); % ones(mesh2Dpos.numPoints, 1);
  Z = FEPack.tools.mymod(cutslope * (Y + idI - numCellsInfinite - 1));
  structLoc = dom.locateInDomain([X, Y, Z]);

  elts = dom.elements(structLoc.elements, :);
  elts = elts'; elts = elts(:);
  coos = structLoc.barycoos;
  coos = coos'; coos = coos(:);

  U2D.positive(:, IcellY) = reshape(sum(reshape(coos .* U3D.positive(elts, IcellY), dom.dimension+1, []), 1), N2Dpos, []);
end

%% Take the trace: Negative side
N2Dneg = mesh2Dneg.numPoints;
U2D.negative = zeros(N2Dneg, size(U3D.negative, 2));
dom = mesh3Dneg.domain('volumic');
for idI = 1:2*numCellsInfinite
  IcellY = (numCellsSemiInfinite_neg*(idI-1)+1):(numCellsSemiInfinite_neg*idI);
  X = mesh2Dneg.points(:, 1);
  Y = mesh2Dneg.points(:, 2); % ones(mesh2Dneg.numPoints, 1);
  Z = FEPack.tools.mymod(cutslope * (Y + idI - numCellsInfinite - 1));
  structLoc = dom.locateInDomain([X, Y, Z]);

  elts = dom.elements(structLoc.elements, :);
  elts = elts'; elts = elts(:);
  coos = structLoc.barycoos;
  coos = coos'; coos = coos(:);

  U2D.negative(:, IcellY) = reshape(sum(reshape(coos .* U3D.negative(elts, IcellY), dom.dimension+1, []), 1), N2Dneg, []);
end

profile viewer

%% Plot U
figure;%(1);
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

xlim([-numCellsSemiInfinite_neg + 1, numCellsSemiInfinite_pos - 1]);
ylim([-numCellsInfinite, numCellsInfinite]);
%caxis([-0.02, 0.02]);

%%

xlim([-numCellsSemiInfinite_neg + 1, numCellsSemiInfinite_neg - 1]);
ylim([-numCellsInfinite, numCellsInfinite]);
%caxis([-0.025, 0.025]);


% %% Compare U in the rational case
% if (compareU)
%   %
%   tic;
%   volBilinearIntg2D = @(muco, rhoco) (muco * grad2(u)) * grad2(v) - (opts.omega^2) * ((rhoco*id(u))*id(v));
% 
%   BCstruct2Dpos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh2Dpos.domain('xmin'));
%   BCstruct2Dpos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh2Dpos.domain('xmax'));
%   BCstruct2Dpos.BCdu = BCstruct_pos.BCdu;
%   BCstruct2Dpos.BCu = BCstruct_pos.BCu;
%   volBilinearIntg2Dpos = volBilinearIntg2D(mu2Dpos, rho2Dpos);
% 
%   %
%   BCstruct2Dneg.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh2Dneg.domain('xmin'));
%   BCstruct2Dneg.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh2Dneg.domain('xmax'));
%   BCstruct2Dneg.BCdu = BCstruct_neg.BCdu;
%   BCstruct2Dneg.BCu = BCstruct_neg.BCu;
%   volBilinearIntg2Dneg = volBilinearIntg2D(mu2Dneg, rho2Dneg);
% 
%   %
%   % jumpLinearIntg2D = G * id(v);
%   
%   %
%   Ur = PeriodicSpaceJumpBVP(1, 2, 1, ...
%                         volBilinearIntg2Dpos, mesh2Dpos, BCstruct2Dpos, numCellsSemiInfinite_pos,...
%                         volBilinearIntg2Dneg, mesh2Dneg, BCstruct2Dneg, numCellsSemiInfinite_neg,...
%                         G, numCellsInfinite, numFloquetPoints, opts);
%   toc;
%   %
%   E.positive = U2D.positive - Ur.positive;
%   E.negative = U2D.negative - Ur.negative;
% 
%   %%
%   % figure;
%   % for idS = 1:numCellsSemiInfinite_pos
%   %   for idI = 1:2*numCellsInfinite
%   %     Icell = sub2ind([numCellsSemiInfinite_pos, 2*numCellsInfinite], idS, idI);
%   %     X = mesh2Dpos.points(:, 1) + (idS - 1);
%   %     Y = mesh2Dpos.points(:, 2) + (idI - numCellsInfinite - 1);
% 
%   %     trisurf(mesh2Dpos.triangles, X, Y, real(Ur.positive(:, Icell)));
%   %     hold on;
%   %     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
%   %     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
%   %     % caxis([-0.02, 0.02]);
%   %     % colormap jet
%   %   end
%   % end
%   % %
%   % for idS = 1:numCellsSemiInfinite_neg
%   %   for idI = 1:2*numCellsInfinite
%   %     Icell = sub2ind([numCellsSemiInfinite_neg, 2*numCellsInfinite], idS, idI);
%   %     X = mesh2Dneg.points(:, 1) - (idS - 1);
%   %     Y = mesh2Dneg.points(:, 2) + (idI - numCellsInfinite - 1);
%   %     trisurf(mesh2Dneg.triangles, X, Y, real(Ur.negative(:, Icell)));
%   %     hold on;
%   %     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
%   %     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
%   %     % colormap jet
%   %     % caxis([-0.02, 0.02]);
%   %   end
%   % end
%   
%   % %
%   % xlim([-numCellsSemiInfinite_neg + 1, numCellsSemiInfinite_pos - 1]);
%   % ylim([-numCellsInfinite, numCellsInfinite]);
%   % caxis([-0.025, 0.025]);
% 
% end
% 
% %%
% % figure;
% % U2Dcst = load('U2D.mat');
% % E.positive = U2D.positive - U2Dcst.U2D.positive;
% % E.negative = U2D.negative - U2Dcst.U2D.negative;
% % for idS = 1:numCellsSemiInfinite_pos
% %   for idI = 1:2*numCellsInfinite
% %     Icell = sub2ind([numCellsSemiInfinite_pos, 2*numCellsInfinite], idS, idI);
% %     X = mesh2Dpos.points(:, 1) + (idS - 1);
% %     Y = mesh2Dpos.points(:, 2) + (idI - numCellsInfinite - 1);
% 
% %     trisurf(mesh2Dpos.triangles, X, Y, real(E.positive(:, Icell)));
% %     hold on;
% %     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
% %     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
% %     % caxis([-0.02, 0.02]);
% %     % colormap jet
% %   end
% % end
% % %
% % for idS = 1:numCellsSemiInfinite_neg
% %   for idI = 1:2*numCellsInfinite
% %     Icell = sub2ind([numCellsSemiInfinite_neg, 2*numCellsInfinite], idS, idI);
% %     X = mesh2Dneg.points(:, 1) - (idS - 1);
% %     Y = mesh2Dneg.points(:, 2) + (idI - numCellsInfinite - 1);
% %     trisurf(mesh2Dneg.triangles, X, Y, real(E.negative(:, Icell)));
% %     hold on;
% %     view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
% %     set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
% %     % colormap jet
% %     % caxis([-0.02, 0.02]);
% %   end
% % end
% 
% % %%
% % xlim([-numCellsSemiInfinite_neg + 1, numCellsSemiInfinite_pos - 1]);
% % ylim([-numCellsInfinite, numCellsInfinite]);
% % caxis([-0.0025, 0.00025]);
% 
% 
% save output_periodic.mat
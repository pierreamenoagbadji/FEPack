clear; clc;
import FEPack.*

% profile ON
opts.omega = 8 + 0.1i;
theta = pi/3;
opts.cutvec = [cos(theta), sin(theta)];
orientation = 1;
coInf = 2;
mu = @(x) 1 + 0.25*cos(2*pi*x(:, 1)) + 0.25*cos(2*pi*x(:, 2));
rho = @(x) 1 + 0.5*cos(2*pi*x(:, 1)).*sin(2*pi*x(:, 2));

N = 32;
BB = [0, 1; 0, 1]; BB(coInf, 2) = orientation;
mesh = meshes.MeshRectangle(1, BB(1, :), BB(2, :), N, N);
cell = mesh.domain('volumic');
Sigma0 = mesh.domains{2*coInf};
Sigma1 = mesh.domains{2*coInf-1};
%%
u = pdes.PDEObject; v = dual(u);
volIntg = (mu * gradDir(u, opts.cutvec))*gradDir(v, opts.cutvec) - (opts.omega^2) * ((rho*id(u))*id(v));

BoundaryStruct.basisfun0 = FEPack.spaces.PeriodicLagrangeBasis(Sigma0);
BoundaryStruct.basisfun1 = FEPack.spaces.PeriodicLagrangeBasis(Sigma1);
% FourierIds = [0 0]; FourierIds(3-coInf) = N/4;
% BoundaryStruct.basisfun0 = spaces.FourierBasis(Sigma0, FourierIds);
% BoundaryStruct.basisfun1 = spaces.FourierBasis(Sigma1, FourierIds);
% BoundaryStruct.basisfun = spaces.SpectralBasis(Sigma0, speye(Sigma0.numPoints), Sigma0.numPoints);
% BoundaryStruct.basisfun.computeBasisMatrices(0);

Nb = BoundaryStruct.basisfun0.numBasis;
BoundaryStruct.BCdu = 1.0;
BoundaryStruct.BCu = 1.0;
BoundaryStruct.representation = 'projection';
BoundaryStruct.phi = @(x) exp(2i*pi*x(:, 1));%(size(x, 1), 1);

% SS1 = FEPack.pdes.Form.intg_TU_V(Sigma1, 1*eye(Nb), BoundaryStruct.basisfun, BoundaryStruct.representation);
% SS0 = FEPack.pdes.Form.intg_TU_V(Sigma0, 1*eye(Nb), BoundaryStruct.basisfun, BoundaryStruct.representation);
% SS2 = FEPack.pdes.Form.intg_TU_V(Sigma0,     1.0);
%
% figure; spy(SS1); figure; spy(SS2)
% max(max(abs(SS1 - SS0)))

U = PeriodicHalfGuideBVP(mesh, orientation, coInf, volIntg, BoundaryStruct, 4, opts);
U = full(U);
% profile VIEWER

%%


% Plot U
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
for idI = 1:size(U, 2)
  % figure(3);
  subplot(-coInf+3, coInf, 1);
  trisurf(mesh.triangles, mesh.points(:, 1) + (coInf == 1) * orientation * (idI-1),...
                          mesh.points(:, 2) + (coInf == 2) * orientation * (idI-1), real(U(:, idI)));
  hold on;
  view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

  subplot(-coInf+3, coInf, 2);
  trisurf(mesh.triangles, mesh.points(:, 1) + (coInf == 1) * orientation * (idI-1),...
                          mesh.points(:, 2) + (coInf == 2) * orientation * (idI-1), imag(U(:, idI)));
  hold on;
  view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

  % %
  % Ur = BoundaryStruct.phi(mesh.points(:, 1) - (mesh.points(:, 2) + orientation*(idI-1))*opts.cutvec(1)/opts.cutvec(2)) .* exp(orientation*1i*opts.omega*(mesh.points(:, 2) + orientation*(idI-1))/opts.cutvec(2));
  % figure(4);
  % subplot(-coInf+3, coInf, 1);
  % trisurf(mesh.triangles, mesh.points(:, 1) + (coInf == 1) * orientation * (idI-1),...
  %                         mesh.points(:, 2) + (coInf == 2) * orientation * (idI-1), real(Ur));
  % hold on;
  % view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  % set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  %
  % subplot(-coInf+3, coInf, 2);
  % trisurf(mesh.triangles, mesh.points(:, 1) + (coInf == 1) * orientation * (idI-1),...
  %                         mesh.points(:, 2) + (coInf == 2) * orientation * (idI-1), imag(Ur));
  % hold on;
  % view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  % set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  % % pause;
end

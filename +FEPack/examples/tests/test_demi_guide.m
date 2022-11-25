clear; clc;
import FEPack.*

% profile ON
opts.omega = 8 + 0.1i;
theta = pi/3;
opts.cutvec = [cos(theta), sin(theta)];
orientation = 1;
coInf = 1;
mu = @(x) 1 + 0.25*cos(2*pi*x(:, 1)) + 0.25*cos(2*pi*x(:, 2));
rho = @(x) 1 + 0.5*cos(2*pi*x(:, 1)).*sin(2*pi*x(:, 2));

N = 32;
BB = [0, 1; 0, 1]; BB(coInf, 2) = orientation;
mesh = meshes.MeshRectangle(1, BB(1, :), BB(2, :), N, N);
cell = mesh.domain('volumic');
Sigma0 = mesh.domains{2*coInf};
Sigma1 = mesh.domains{2*coInf-1};
%
u = pdes.PDEObject; v = dual(u);
% volIntg = (mu * gradDir(u, opts.cutvec))*gradDir(v, opts.cutvec) - (opts.omega^2) * ((rho*id(u))*id(v));
muGrad = @(x) kron(eye(3), mu(x));
volIntg = (muGrad * grad(u))*grad(v) - (opts.omega^2) * ((rho*id(u))*id(v));

BCstruct.spB0 = FEPack.spaces.PeriodicLagrangeBasis(Sigma0);
BCstruct.spB1 = FEPack.spaces.PeriodicLagrangeBasis(Sigma1);

Nb = BCstruct.spB0.numBasis;
BCstruct.BCdu = 1.0;
BCstruct.BCu = 0.0;% @(x) 1 + 0.5*sin(2*pi*x(:, 1));% 1.0;
BCstruct.representation = '';
BCstruct.phi = @(x) exp(2i*pi*x(:, 1));%(size(x, 1), 1);

[U, BCstruct] = PeriodicHalfGuideBVP(mesh, orientation, coInf, volIntg, BCstruct, 4, opts);
U = full(U);
% profile VIEWER

% Plot U
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
for idI = 1:size(U, 2)
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
end

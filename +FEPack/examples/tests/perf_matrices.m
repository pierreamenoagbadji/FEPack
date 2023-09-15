clear; clc;
%%
import FEPack.*

%% Problem-related variables
opts.omega = 8 + 0.25i;
period = 1;
vecper = [-1, 1]; % [-sqrt(2), 1];
cutvec = [1/vecper(2); -vecper(1)/vecper(2)];
opts.cutmat = [[1; 0; 0], [0; cutvec]];
cutslope = cutvec(2) / cutvec(1);

%% Mesh
pregenerate_mesh = 1;
struct_mesh = 1;

%% Coefficients
% Cartesian coefficients
type_mu_pos = 2;  % 1: constant; 2: periodic; 
type_rho_pos = 2; % 1: constant; 2: periodic;  

switch (type_mu_pos)
case 1
  mu2DposCart = @(x) ones(size(x, 1), 1);
case 2
  mu2DposCart = @(x) 1 + 0.25*cos(2*pi*x(:, 1)) + 0.25*cos(2*pi*x(:, 2));
end

switch (type_rho_pos)
case 1
  rho2DposCart = @(x) ones(size(x, 1), 1);
case 2
  rho2DposCart = @(x) 1 + 0.5*cos(2*pi*x(:, 1)).*sin(2*pi*x(:, 2));
end

% 2D coefficients
mu2Dpos  = @(x) mu2DposCart( [x(:, 1) + cutvec(2)*x(:, 2), cutvec(1)*x(:, 2), zeros(size(x, 1), 1)]);
rho2Dpos = @(x) rho2DposCart([x(:, 1) + cutvec(2)*x(:, 2), cutvec(1)*x(:, 2), zeros(size(x, 1), 1)]);

% 3D coefficients
fun3D = @(fun2Dcart, x) fun2Dcart([x(:, 1) + x(:, 2), x(:, 3), zeros(size(x, 1), 1)]);
mu3Dpos = @(x) fun3D(mu2DposCart, x);
rho3Dpos = @(x) fun3D(rho2DposCart, x);

%% Bilinear and linear forms
u = pdes.PDEObject; v = dual(u);
gradu_gradv = @(muco) (muco * (opts.cutmat' * grad3(u))) * (opts.cutmat' * grad3(v));
gradu_vec1v = @(muco) (muco * (opts.cutmat' * grad3(u))) * (opts.cutmat' * [0; 1; 0] * v);
vec1u_gradv = @(muco) (muco * (opts.cutmat' * [0; 1; 0] * u)) * (opts.cutmat' * grad3(v));
vec1u_vec1v = @(muco) (muco * (opts.cutmat' * [0; 1; 0] * u)) * (opts.cutmat' * [0; 1; 0] * v);
u_v = @(rhoco) ((rhoco*u)*v);

%
nodes = [5, 10, 20, 30];%, 40, 50];
Nt = length(nodes);
temps1 = zeros(1, Nt);
temps2 = zeros(1, Nt);

FloquetVar = -pi + 2*pi*rand;

%%
for idI = 1:Nt
  idI
  
  numNodes3D = nodes(idI);

  if pregenerate_mesh

    if (struct_mesh)
      mesh_prefix = 'struct';
    else
      mesh_prefix = 'unstruct';
    end

    % Pick mesh from saved file
    m = load(['pregenMeshes/3D/', mesh_prefix, '_mesh_3D_', int2str(numNodes3D), '_positive.mat']);
    mesh3Dpos = m.mesh;

  else

    mesh3Dpos = meshes.MeshCuboid(struct_mesh, [0 1], [0 1], [0 1], numNodes3D, numNodes3D, numNodes3D);

  end

  %% Compute FE elementary matrices
  mat_gradu_gradv_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), gradu_gradv(mu3Dpos));
  mat_gradu_vec1v_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), gradu_vec1v(mu3Dpos));
  mat_vec1u_gradv_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), vec1u_gradv(mu3Dpos));
  mat_vec1u_vec1v_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), vec1u_vec1v(mu3Dpos));
  mat_u_v_pos         = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'),        u_v(rho3Dpos));

  %
  AA1 = mat_gradu_gradv_pos + 1i * FloquetVar * mat_vec1u_gradv_pos - 1i * FloquetVar * mat_gradu_vec1v_pos + FloquetVar * FloquetVar * mat_vec1u_vec1v_pos - (opts.omega^2) * mat_u_v_pos;

  AA2 = mat_gradu_gradv_pos + 1i * FloquetVar * mat_vec1u_gradv_pos - 1i * FloquetVar * mat_gradu_vec1v_pos + FloquetVar * FloquetVar * mat_vec1u_vec1v_pos + mat_u_v_pos;

  %
  LL = mat_u_v_pos * ones(mesh3Dpos.numPoints, 1);

  tic;
  UU1 = AA1 \ LL;
  temps1(idI) = toc;

  tic;
  UU2 = AA2 \ LL;
  temps2(idI) = toc;
end

%%
fid = fopen('outputs/temps.txt', 'w');
fprintf(fid, '%2d\t %0.5e\t %0.5e\n', [nodes; temps1; temps2]);
fclose(fid);
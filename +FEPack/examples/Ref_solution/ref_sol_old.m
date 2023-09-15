clear; clc;
%%
import FEPack.*
% profile ON

%% Problem-related variables
opts.omega = 8 + 0.25i;
period = 1;
opts.verbose = 0;
problem_setting = 'A'; % 'A' or 'B'
pregenerate_mesh = 0;
struct_mesh = 1;

Nx = 3; 
Nz = 3;
sizeCellX = 1.0;
sizeCellZ = 1.0;
Zorigin = -sizeCellZ * 0.5 * Nz;

% The number of cells in the x-direction must be even
Nx = 2 * ceil(Nx / 2);

if strcmpi(problem_setting, 'A')

  % 2D coefficients
  period_posFun = 1;
  period_negFun = 0.5 * sqrt(2);

  mu_pos  = @(x) ones(size(x, 1), 1); % 0.5 + perCutoffCircle(x, [1; 0], [0; period_posFun], [0.5, 0.5], [-0.2, 0.2]);
  rho_pos = @(x) ones(size(x, 1), 1); % 0.5 + perCutoffCuboid(x, [1; 0], [0; period_posFun], [0.5, 0.5], [-0.2, 0.2], [-0.2, 0.2], 0.5, 1);
  mu_neg  = @(x) ones(size(x, 1), 1); % 1 + 0.5 * cos(2*pi*x(:, 1)) .* cos(2*pi*x(:, 2)/period_negFun);
  rho_neg = @(x) ones(size(x, 1), 1); % 1 + 0.25 * sin(2*pi*x(:, 1)) + 0.25 * sin(2*pi*x(:, 2)/period_negFun);
  
else

  % 2D coefficients
  vecperFun = [-0.5*sqrt(2), 1]; % [-sqrt(2), 1];
  
  mu_pos  = @(x) 0.5 + perCutoffCircle(x, [1; 0], vecperFun, [0.5, 0.5], [-0.2, 0.2]);
  rho_pos = @(x) 0.5 + perCutoffCuboid(x, [1; 0], vecperFun, [0.5, 0.5], [-0.2, 0.2], [-0.2, 0.2], 0.5, 1);
  mu_neg  = @(x) ones(size(x, 1), 1);
  rho_neg = @(x) ones(size(x, 1), 1);

end

mu_glo  = @(x)  mu_pos(x) .* (x(:, 1) >= 0) +  mu_neg(x) .* (x(:, 1) < 0);
rho_glo = @(x) rho_pos(x) .* (x(:, 1) >= 0) + rho_neg(x) .* (x(:, 1) < 0);

% Jump data
% alpha_G = 3;
% eps_G = 1e-8;
% supp_G = -log(eps_G) / alpha_G;
% G = @(x) exp(-alpha_G * x(:, 2).^2) .* (abs(x(:, 2)) <= supp_G);
% G = @(x) FEPack.tools.cutoff(x(:, 2), -0.1, 0.1);
G = @(x) ones(size(x, 1), 1);

%% Mesh
numNodesX = 20;
numNodesZ = 20;

if pregenerate_mesh

  if (struct_mesh)
    mesh_prefix = 'struct';
  else
    mesh_prefix = 'unstruct';
  end

  % Pick mesh from saved file
  m = load(['pregenMeshes/2D/', mesh_prefix, '_mesh_2D_', int2str(numNodes), '_positive.mat']);
  meshcell = m.mesh;

else

  meshcell = meshes.MeshRectangle(struct_mesh, [0 sizeCellX], [0 sizeCellZ], numNodesX, numNodesZ);

end

% Domains
Sigma1 = meshcell.domain('xmin'); N1 = Sigma1.numPoints;
Sigma2 = meshcell.domain('xmax'); N2 = Sigma2.numPoints;
Sigma3 = meshcell.domain('ymin'); N3 = Sigma3.numPoints;
Sigma4 = meshcell.domain('ymax'); N4 = Sigma4.numPoints;
domCell = meshcell.domain('volumic');

%% Local cell problems
N = meshcell.numPoints;
u = pdes.PDEObject; v = dual(u);

% Dirichlet conditions on the boundary
ecs = ((u|Sigma1) == 0.0) & ((u|Sigma2) == 0.0) &...
      ((u|Sigma3) == 0.0) & ((u|Sigma4) == 0.0);

% Surfacic rhs
B1 = sparse(Sigma1.IdPoints, (1:N1), 1.0, N, N1);
B2 = sparse(Sigma2.IdPoints, (1:N2), 1.0, N, N2);
B3 = sparse(Sigma3.IdPoints, (1:N3), 1.0, N, N3);
B4 = sparse(Sigma4.IdPoints, (1:N4), 1.0, N, N4);

ecs.applyEcs;
ecs.b = [B1, B2, B3, B4];

% Solve and save local cell problems
fprintf('%3d et %3d\n', 0, 0);
for idX = 1:Nx
  for idZ = 1:Nz

    fprintf('\b\b\b\b\b\b\b\b\b\b\b%3d et %3d\n', idX, idZ);
    
    % Compute FE matrices
    cellXorigin = sizeCellX * (idX - 1 - 0.5 * Nx);
    cellZorigin = sizeCellZ * (idZ - 1) + Zorigin;

    mu_cell  = @(x)  mu_glo([x(:, 1) + cellXorigin, x(:, 2) + cellZorigin]);
    rho_cell = @(x) rho_glo([x(:, 1) + cellXorigin, x(:, 2) + cellZorigin]);

    MM = FEPack.pdes.Form.intg(domCell, rho_cell * u * v);
    KK = FEPack.pdes.Form.intg(domCell, (mu_cell * grad2(u)) * grad2(v));
    AA = KK - (opts.omega^2) * MM;

    % Elimination
    AA0 =  ecs.P * AA * ecs.P';
    LL0 = -ecs.P * AA * ecs.b;

    % Solve the linear system
    Ecell0 = ecs.b + ecs.P' * (AA0 \ LL0);
    
    % Local cell solutions
    Ecell = cell(4, 1);
    Ecell{1} = Ecell0(:,          1:N1);
    Ecell{2} = Ecell0(:,       N1+1:N1+N2);
    Ecell{3} = Ecell0(:,    N1+N2+1:N1+N2+N3);
    Ecell{4} = Ecell0(:, N1+N2+N3+1:N1+N2+N3+N4);

    % DtN operators
    DtNop = cell(4, 4);
    for idI = 1:4
      for idJ = 1:4
        DtNop{idI, idJ} = Ecell{idJ}' * (AA * Ecell{idI});
      end
    end
    
    % Save output
    save(['outputs/local_cell_sol_', int2str(idX), '_', int2str(idZ)], 'Ecell', 'DtNop');
  end
end


%% Plot some local cell solutions
plot_local_cell = 0;

if (plot_local_cell)
  figure;
  set(groot,'defaultAxesTickLabelInterpreter','latex');
  set(groot,'defaulttextinterpreter','latex');
  set(groot,'defaultLegendInterpreter','latex');

  LC1 = load(['outputs/local_cell_sol_', int2str(1), '_', int2str(1)]);

  for idI = 1:N1
    for idE = 1:4
      subplot(2, 2, idE);
      trisurf(meshcell.triangles, meshcell.points(:, 1), meshcell.points(:, 2), full(real(LC1.Ecell{idE}(:, idI))));
      hold off;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
    end
    pause;
  end
end

%% Solutions of local strip problems
% Npsi = Nz * (N1 - 1) + 1;
Npsi = Nz * N1;
Nphi = (Nz - 1) * N3;

psis0 = speye(Npsi, Npsi);
psis1 = speye(Npsi, Npsi);
% psis0 =  ones(Npsi, 1);
% psis1 = zeros(Npsi, 1);
fprintf('%3d\n', 0);

for idX = 1:Nx
  
  fprintf('\b\b\b\b%3d\n', idX);
  
  % Jacobi matrix 
  A = sparse(Nphi, Nphi);
  for idZ = 1:Nz-1 % Diagonal
    LC1 = load(['outputs/local_cell_sol_', int2str(idX), '_', int2str(idZ)]);
    LC2 = load(['outputs/local_cell_sol_', int2str(idX), '_', int2str(idZ+1)]);
    
    idPhi = 1+(idZ-1)*N3:idZ*N3;
    A(idPhi, idPhi) = LC2.DtNop{3, 3} + LC1.DtNop{4, 4}; % Diagonal
  end

  for idZ = 1:Nz-2 % Lower and upper diagonals
    LC = load(['outputs/local_cell_sol_', int2str(idX), '_', int2str(idZ+1)]);
    
    idPhi = 1+(idZ-1)*N3:idZ*N3;
    A(idPhi + N3, idPhi) = LC.DtNop{3, 4};
    A(idPhi, idPhi + N3) = LC.DtNop{4, 3};
  end

  % Rhs
  B0 = sparse(Nphi, size(psis0, 2));
  B1 = sparse(Nphi, size(psis1, 2));
  for idZ = 1:Nz-1
    LC1 = load(['outputs/local_cell_sol_', int2str(idX), '_', int2str(idZ)]);
    LC2 = load(['outputs/local_cell_sol_', int2str(idX), '_', int2str(idZ+1)]);

    idPhi = 1+(idZ-1)*N3:idZ*N3;
    % idPsi = 1+(idZ-1)*(N1-1):1+idZ*(N1-1);
    % B0(idPhi, :) = -LC1.DtNop{1, 4} * psis0(idPsi, :) - LC2.DtNop{1, 3} * psis0(idPsi+N1-1, :);
    % B1(idPhi, :) = -LC1.DtNop{2, 4} * psis1(idPsi, :) - LC2.DtNop{2, 3} * psis1(idPsi+N1-1, :);
    idPsi = 1+(idZ-1)*N1:idZ*N1;
    B0(idPhi, :) = -LC1.DtNop{1, 4} * psis0(idPsi, :) - LC2.DtNop{1, 3} * psis0(idPsi+N1, :);
    B1(idPhi, :) = -LC1.DtNop{2, 4} * psis1(idPsi, :) - LC2.DtNop{2, 3} * psis1(idPsi+N1, :);
  end

  % Solve system
  traceEstrip = A \ [B0, B1];
  traceE0strip = traceEstrip(:, 1:size(B0, 2));
  traceE1strip = traceEstrip(:, 1+size(B0, 2):end);

  % Add BC
  traceE0strip = [zeros(N3, size(B0, 2)); traceE0strip; zeros(N3, size(B0, 2))];
  traceE1strip = [zeros(N3, size(B1, 2)); traceE1strip; zeros(N3, size(B1, 2))];

  % Save output
  save(['outputs/local_strip_trace_', int2str(idX)], 'traceE0strip', 'traceE1strip');

  %
  % plot_local_strip = false;
  % if (plot_local_strip)
  %   figure;
  %   set(groot,'defaultAxesTickLabelInterpreter','latex');
  %   set(groot,'defaulttextinterpreter','latex');
  %   set(groot,'defaultLegendInterpreter','latex');
    
  %   E = cell(Nz, 1);
    
  %   for idZ = 1:Nz
  %     LC = load(['outputs/local_cell_sol_', int2str(idX), '_', int2str(idZ)]);
  %     idPhi = 1+(idZ-1)*N3:idZ*N3;
  %     % idPsi = 1+(idZ-1)*(N1-1):1+idZ*(N1-1);
  %     idPsi = 1+(idZ-1)*N1:idZ*N1;
  %     E{idZ} = LC.Ecell{1} * psis0(idPsi, :) + LC.Ecell{2} * psis1(idPsi, :) +...
  %              LC.Ecell{3} * (traceE0strip(idPhi, :) + traceE1strip(idPhi, :)) +...
  %              LC.Ecell{4} * (traceE0strip(idPhi+N3, :) + traceE1strip(idPhi+N3, :));
  %   end
    
  %   for idI = 1:size(psis0, 2)      
  %     for idZ = 1:Nz
  %       cellZorigin = sizeCellZ * (idZ - 1) + Zorigin;
  %       trisurf(meshcell.triangles, meshcell.points(:, 1), meshcell.points(:, 2) + cellZorigin, full(real(E{idZ}(:, idI))));
  %       hold on;
  %       view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  %       set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  %     end

      
  %     if (size(psis0, 2) > 1)
  %       caxis([0, 1]);
  %       hold off;
  %       pause;
  %     end
  %   end
  % end

  % validation = true;
  % if (validation)
  %   cellXorigin = sizeCellX * (idX - 1 - 0.5 * Nx);

  %   figure;
  %   set(groot,'defaultAxesTickLabelInterpreter','latex');
  %   set(groot,'defaulttextinterpreter','latex');
  %   set(groot,'defaultLegendInterpreter','latex');

  %   % Strip solution
  %   E = cell(Nz, 1);
  
  %   for idZ = 1:Nz
  %     LC = load(['outputs/local_cell_sol_', int2str(idX), '_', int2str(idZ)]);
  %     idPhi = 1+(idZ-1)*N3:idZ*N3;
  %     idPsi = 1+(idZ-1)*(N1-1):1+idZ*(N1-1);
  %     % idPsi = 1+(idZ-1)*N1:idZ*N1;
  %     E{idZ} = LC.Ecell{1} * psis0(idPsi, :) + LC.Ecell{2} * psis1(idPsi, :) +...
  %              LC.Ecell{3} * (traceE0strip(idPhi, :) + traceE1strip(idPhi, :)) +...
  %              LC.Ecell{4} * (traceE0strip(idPhi+N3, :) + traceE1strip(idPhi+N3, :));
  %   end
    
  %   subplot(1, 2, 1);
  %   for idI = 1:size(psis0, 2)      
  %     for idZ = 1:Nz
  %       cellZorigin = sizeCellZ * (idZ - 1) + Zorigin;
  %       trisurf(meshcell.triangles, meshcell.points(:, 1) + cellXorigin, meshcell.points(:, 2) + cellZorigin, full(real(E{idZ}(:, idI))));
  %       hold on;
  %       view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  %       set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  %     end
  %   end


  %   % Reference solution
  %   meshstrip = FEPack.meshes.MeshRectangle(struct_mesh, cellXorigin + [0,  sizeCellX], Zorigin + [0, Nz*sizeCellZ], numNodesX, Npsi);
    
  %   AAstrip = FEPack.pdes.Form.intg(meshstrip.domain('volumic'), (mu_glo * grad2(u)) * grad2(v))...
  %      - (opts.omega^2) * FEPack.pdes.Form.intg(meshstrip.domain('volumic'), rho_glo * u * v);
  %   ecs = ((u|meshstrip.domain('xmin')) == 0.0) & ((u|meshstrip.domain('xmax')) == 0.0) &...
  %         ((u|meshstrip.domain('ymin')) == 0.0) & ((u|meshstrip.domain('ymax')) == 0.0);
  %   Nstrip = meshstrip.domain('xmin').numPoints;
  %   B0strip = sparse(meshstrip.domain('xmin').IdPoints, (1:Nstrip), 1.0, meshstrip.numPoints, Nstrip);
  %   B1strip = sparse(meshstrip.domain('xmax').IdPoints, (1:Nstrip), 1.0, meshstrip.numPoints, Nstrip);
  %   ecs.applyEcs;
  %   ecs.b = [B0strip, B1strip];
  %   AA0 =  ecs.P * AAstrip * ecs.P';
  %   LL0 = -ecs.P * AAstrip * ecs.b;
  %   Estrip = ecs.b + ecs.P' * (AA0 \ LL0);
  %   E0strip = Estrip(:, 1:size(B0strip, 2));
  %   E1strip = Estrip(:, 1+size(B0strip, 2):end);

  %   % Plot solution
  %   subplot(1, 2, 2);
  %   trisurf(meshstrip.triangles, meshstrip.points(:, 1), meshstrip.points(:, 2), full(real(E0strip * psis0)));
  %   % hold on;
  %   % trisurf(meshneg.triangles, meshneg.points(:, 1), meshneg.points(:, 2), full(real(E0strip * psis0)));
  %   view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  %   set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  % end

end



%% Compute all traces
Ntrace = (Nx - 1) * Npsi;

% Jacobi matrix
A = sparse(Ntrace, Ntrace);

for idX = 1:Nx-1
  SS1 = load(['outputs/local_strip_trace_', int2str(idX)], 'traceE0strip', 'traceE1strip');
  SS2 = load(['outputs/local_strip_trace_', int2str(idX+1)], 'traceE0strip', 'traceE1strip');
  A0 = sparse(Npsi, Npsi);

  for idZ = 1:Nz
    LC1 = load(['outputs/local_cell_sol_', int2str(idX), '_', int2str(idZ)]);
    LC2 = load(['outputs/local_cell_sol_', int2str(idX+1), '_', int2str(idZ)]);

    idPhi = 1+(idZ-1)*N3:idZ*N3;
    idPsi = 1+(idZ-1)*N1:idZ*N1;

    A0(idPsi, idPsi) = LC1.DtNop{2, 2} + LC2.DtNop{1, 1};
    A0(idPsi, :) = A0(idPsi, :) +...
                   LC1.DtNop{3, 2} * SS1.traceE1strip(idPhi, :) + LC1.DtNop{4, 2} * SS1.traceE1strip(idPhi+N3, :) +...
                   LC2.DtNop{3, 1} * SS2.traceE0strip(idPhi, :) + LC2.DtNop{4, 1} * SS2.traceE0strip(idPhi+N3, :);
  end
  
  idI = 1+(idX-1)*Npsi:idX*Npsi;
  A(idI, idI) = A0;
end

for idX = 1:Nx-2 
  SS = load(['outputs/local_strip_trace_', int2str(idX+1)], 'traceE0strip', 'traceE1strip');
  A01 = sparse(Npsi, Npsi);
  A10 = sparse(Npsi, Npsi);
  
  for idZ = 1:Nz
    LC = load(['outputs/local_cell_sol_', int2str(idX+1), '_', int2str(idZ)]);

    idPhi = 1+(idZ-1)*N3:idZ*N3;
    idPsi = 1+(idZ-1)*N1:idZ*N1;

    A10(idPsi, idPsi) = LC.DtNop{2, 1};
    A10(idPsi, :) = A10(idPsi, :) +...
                    LC.DtNop{3, 1} * SS.traceE0strip(idPhi, :) +...
                    LC.DtNop{4, 1} * SS.traceE0strip(idPhi+N3, :);

    A01(idPsi, idPsi) = LC.DtNop{1, 2};
    A01(idPsi, :) = A01(idPsi, :) +...
                    LC.DtNop{3, 2} * SS.traceE1strip(idPhi, :) +...
                    LC.DtNop{4, 2} * SS.traceE1strip(idPhi+N3, :);
  end

  idI = 1+(idX-1)*Npsi:idX*Npsi;
  A(idI+Npsi, idI) = A01;
  A(idI, idI+Npsi) = A10;
end

% Right hand side
B = sparse(Ntrace, 1);
MMG = FEPack.pdes.Form.intg(Sigma1, u * v);
MMG = MMG(Sigma1.IdPoints, Sigma1.IdPoints); 
B0 = sparse(Npsi, 1);
for idZ = 1:Nz
  idPsi = 1+(idZ-1)*N1:idZ*N1;
  cellZorigin = sizeCellZ * (idZ - 1) + Zorigin;
  
  B0(idPsi) = MMG * G([meshcell.points(Sigma1.IdPoints, 1),...
                       meshcell.points(Sigma1.IdPoints, 2) + cellZorigin]);
end

idX = Nx/2;
idI = 1+(idX-1)*Npsi:idX*Npsi;
B(idI) = B0;

% Solve system and add BC
solPsi = A \ B;
solPsi = [zeros(Npsi, 1); solPsi; zeros(Npsi, 1)];

%% Save solution's traces
soltraceX = cell(Nx + 1, Nz);
soltraceZ = cell(Nx, Nz + 1);
solPsi = reshape(solPsi, Npsi, Nx + 1);

for idX = 1:Nx+1
  for idZ = 1:Nz
    idPsi = 1+(idZ-1)*N1:idZ*N1;
    soltraceX{idX, idZ} = solPsi(idPsi, idX);
  end
end

for idX = 1:Nx
  for idZ = 1:Nz+1
    SS = load(['outputs/local_strip_trace_', int2str(idX)], 'traceE0strip', 'traceE1strip');
    
    idPhi = 1+(idZ-1)*N3:idZ*N3;
    soltraceZ{idX, idZ} = SS.traceE0strip(idPhi, :) * solPsi(:, idX) +...
                          SS.traceE1strip(idPhi, :) * solPsi(:, idX+1);
  end
end

save('outputs/solution_trace', 'soltraceX', 'soltraceZ');

%% Compute solution
U = cell(Nx, Nz);

for idX = 1:Nx
  for idZ = 1:Nz
    LC = load(['outputs/local_cell_sol_', int2str(idX), '_', int2str(idZ)]);
    
    U{idX, idZ} = LC.Ecell{1} * soltraceX{idX, idZ} + LC.Ecell{2} * soltraceX{idX + 1, idZ} +...
                  LC.Ecell{3} * soltraceZ{idX, idZ} + LC.Ecell{4} * soltraceZ{idX, idZ + 1};
  end
end

% for idX = 1:Nx

  

%   for idZ = 1:Nz
%     LC = load(['outputs/local_cell_sol_', int2str(idX), '_', int2str(idZ)]);

%     % idPsi = 1+(idZ-1)*N1:idZ*N1;

%     U{idX, idZ} = LC.Ecell{1} * solPsi()
%   end
% end


%% Plot space solution
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

for idX = 1:Nx      
  for idZ = 1:Nz
    cellXorigin = sizeCellX * (idX - 1 - 0.5 * Nx);
    cellZorigin = sizeCellZ * (idZ - 1) + Zorigin;

    trisurf(meshcell.triangles, meshcell.points(:, 1) + cellXorigin, meshcell.points(:, 2) + cellZorigin, full(real(U{idX, idZ})));
    hold on;
    view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
  end
end

%% Validation
validation = true;
if (validation)
    
  % Positive side
  meshpos = FEPack.meshes.MeshRectangle(struct_mesh, [0,  0.5*Nx*sizeCellX], [Zorigin, Zorigin + Nz*sizeCellZ], numNodesX*Nx/2, numNodesZ*Nz);
  AApos = FEPack.pdes.Form.intg(meshpos.domain('volumic'), (mu_cell * grad2(u)) * grad2(v))...
        - (opts.omega^2) * FEPack.pdes.Form.intg(meshpos.domain('volumic'), u * v);
  ecs = ((u|meshpos.domain('xmin')) == 0.0) & ((u|meshpos.domain('xmax')) == 0.0) &...
        ((u|meshpos.domain('ymin')) == 0.0) & ((u|meshpos.domain('ymax')) == 0.0);
  Npos = meshpos.domain('xmin').numPoints;
  Bpos = sparse(meshpos.domain('xmin').IdPoints, (1:Npos), 1.0, meshpos.numPoints, Npos);
  ecs.applyEcs;
  ecs.b = Bpos;
  AA0 =  ecs.P * AApos * ecs.P';
  LL0 = -ecs.P * AApos * ecs.b;
  Upos = ecs.b + ecs.P' * (AA0 \ LL0);
  Lambdapos = Upos' * AApos * Upos;

  % Negative side
  meshneg = FEPack.meshes.MeshRectangle(struct_mesh, [0 -0.5*Nx*sizeCellX], [Zorigin, Zorigin + Nz*sizeCellZ], numNodesX*Nx/2, numNodesZ*Nz);
  AAneg = FEPack.pdes.Form.intg(meshneg.domain('volumic'), (mu_glo * grad2(u)) * grad2(v))...  
        - (opts.omega^2) * FEPack.pdes.Form.intg(meshneg.domain('volumic'), rho_glo * u * v);
  ecs = ((u|meshneg.domain('xmin')) == 0.0) & ((u|meshneg.domain('xmax')) == 0.0) &...
        ((u|meshneg.domain('ymin')) == 0.0) & ((u|meshneg.domain('ymax')) == 0.0);
  Nneg = meshneg.domain('xmin').numPoints;
  Bneg = sparse(meshneg.domain('xmin').IdPoints, (1:Nneg), 1.0, meshneg.numPoints, Nneg);
  ecs.applyEcs;
  ecs.b = Bneg;
  AA0 =  ecs.P * AAneg * ecs.P';
  LL0 = -ecs.P * AAneg * ecs.b;
  Uneg = ecs.b + ecs.P' * (AA0 \ LL0);
  Lambdaneg = Uneg' * AAneg * Uneg;

  % Interface problem
  MMG = FEPack.pdes.Form.intg(meshpos.domain('xmin'), u * v);
  MMG = MMG(meshpos.domain('xmin').IdPoints, meshpos.domain('xmin').IdPoints);
  vecG = MMG * G(meshpos.points(meshpos.domain('xmin').IdPoints, :));

  solint = (Lambdapos + Lambdaneg) \ vecG;

  % Plot solution
  figure;
  set(groot,'defaultAxesTickLabelInterpreter','latex');
  set(groot,'defaulttextinterpreter','latex');
  set(groot,'defaultLegendInterpreter','latex');
  trisurf(meshpos.triangles, meshpos.points(:, 1), meshpos.points(:, 2), full(real(Upos * solint)));
  hold on;
  trisurf(meshneg.triangles, meshneg.points(:, 1), meshneg.points(:, 2), full(real(Uneg * solint)));
  view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
  set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
end

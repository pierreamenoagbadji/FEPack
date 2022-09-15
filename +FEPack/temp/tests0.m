% % % % %%
% % % % clear; clc;
% % % % import FEPack.*
% % % % u = pdes.PDEObject;
% % % % rmesh = meshes.MeshRectangle;
% % % %
% % % %
% % % % xmin = rmesh.domain('xmin');
% % % % xmax = rmesh.domain('xmax');
% % % % ymin = rmesh.domain('ymin');
% % % % ymax = rmesh.domain('ymax');
% % % %
% % % % %
% % % % % ecs = assignEcs((u|ymin), sparse(ymin.numPoints, 1)) &...
% % % % %       assignEcs((u|ymax), sparse(ymax.numPoints, 1)) &...
% % % % %       assignEcs((u|xmax) - (u|xmin), sparse(xmin.numPoints, 1));
% % % % ecs = assignEcs((u|xmax) - (u|xmin), sparse(xmin.numPoints, 1)) &...
% % % %       assignEcs((u|ymax) - (u|ymin), sparse(ymin.numPoints, 1));
% % % %
% % % % ecs.applyEcs;
% % % %
% % % % % MM = FE.global_matrix('volumic', rmesh, 0, 0);
% % % % F = @(P) [cos(2*pi*P(:, 1)), sin(2*pi*P(:, 2)); zeros(size(P, 1), 1), exp(P(:, 1) + P(:, 2))];
% % % % MM = u.intg_DxU_V('volumic', rmesh, @(x) exp(x(:, 1)));
% % % % display(MM)
% % % % % % display(FEPack.pdes.PDEObject.U_V_elem(2, [0, 0; 1, 0; 0, 1]));
% % % % % Myx = FE.mat_elem(2, [0, 0, 0; 1, 0, 0; 0, 1, 0], 1, 3);
% % % % % Myy = FE.mat_elem(2, [0, 0, 0; 1, 0, 0; 0, 1, 0], 2, 2);
% % % % % display(Myx)
% % % %
% % % % % MM = FE.assembleFEmatrices(rmesh, @(P) FEPack.pdes.PDEObject.U_V_elem(2, P), 'volumic');
% % % % % display(MM);
% % % %
% % % % %%
% % % % clear; clc;
% % % % import FEPack.*
% % % % u = pdes.PDEObject;
% % % % cmesh = meshes.MeshCuboid(0);
% % % % % MM = u.intg_U_V('volumic', cmesh, @(x) exp(x(:, 1)));
% % % % % clc;
% % % % %%
% % % % N = 32;
% % % % points = rand(N, 2);
% % % % points(randi(N, 8, 1), 1) = 0.5;
% % % %
% % % % [~, cleX] = sort(points(:, 1));
% % % % [~, cleY] = sort(points(cleX, 2));
% % % % [~, cle] = sortrows(points, [1 2]);
% % % % sorted_points = points(cle, :);
% % % %
% % % % plot([0, 0], [1, 0], 'k'); hold on;
% % % % plot([1, 0], [1, 1], 'k');
% % % % plot([1, 1], [0, 1], 'k');
% % % % plot([0, 1], [0, 0], 'k');
% % % % plot(points(:, 1), points(:, 2), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
% % % %
% % % % for i = 1:N
% % % %   plot(sorted_points(i, 1), sorted_points(i, 2), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
% % % %   pause;
% % % % end
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % %
% % %
% % %
% % %
% % %%
% % clear; clc;
% % import FEPack.*
% % 
% % N = 8;
% % O = [0.38; 0.6]; r = 0.35;
% % f = @(x) FEPack.tools.cutoff(sqrt((x(:, 1)-O(1)).^2 + (x(:, 2)-O(2)).^2), -r, r, 0.5, 1e-8);% sin(2*pi*x(:, 1)) .* cos(2*pi*x(:, 2));
% % rmesh = meshes.MeshRectangle(1, [0, 1], [0, 1], N, N);
% % figure;
% % trisurf(rmesh.triangles, rmesh.points(:, 1), rmesh.points(:, 2), f(rmesh.points));
% % view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
% % % pause;
% % 
% % xmin = rmesh.domain('xmin');
% % xmax = rmesh.domain('xmax');
% % ymin = rmesh.domain('ymin');
% % ymax = rmesh.domain('ymax');
% % cell = rmesh.domain('volumic');
% % 
% % u = pdes.PDEObject; v = dual(u);
% % MM = pdes.Form.intg_U_V(cell);
% % 
% % matfun = @(P) [2+cos(P(:, 1)), 3+sin(P(:, 2)); exp(P(:, 1)), 2+zeros(size(P, 1), 1)];
% % 
% % KK = pdes.Form.intg_gradU_gradV(cell, matfun);
% % % [MMr, KKr] = matEF(rmesh);
% % AA = KK + MM;
% % LL = MM * f(rmesh.points);
% % 
% % % ecs = assignEcs((u|xmin), 1.0) & assignEcs((u|xmax), 1.0) & assignEcs((u|ymin), 1.0) & assignEcs((u|ymax), 1.0);
% % ecs = assignEcs((u|xmin) - (u|xmax), 0.0) & assignEcs((u|ymin) - (u|ymax), 0.0);
% % ecs.applyEcs;
% % 
% % % idInter = cell.points;
% % % idInter(unique([xmin.points; xmax.points; ymin.points; ymax.points])) = [];
% % % N0 = length(idInter);
% % % P0 = sparse((1:N0), idInter, 1, N0, rmesh.numPoints);
% % % ecs.P = P0;
% % 
% % LL = LL - AA * ecs.b;
% % 
% % AA0 = ecs.P * AA * ecs.P';
% % LL0 = ecs.P * LL;
% % U0 = AA0 \ LL0;
% % 
% % U = ecs.P' * U0 + ecs.b;
% % 
% % figure;
% % trisurf(rmesh.triangles, rmesh.points(:, 1), rmesh.points(:, 2), U);
% % view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
% % set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
% %%
% % % clc;
% % % import FEPack.*
% % % u = pdes.PDEObject;
% % % P = rand(2, 3);
% % % % P = [0, 0; 1, 0; 0, 1];
% % % [MelA, KelA] = mat_elem(P(:, 1), P(:, 2), P(:, 3));
% % % MelB = u.mat_elem(P, 2, 2, [1 0 0], [1 0 0]);
% % % KelB = u.mat_elem(P, 2, 2, [0 1 0], [0 1 0]) + ...
% % %        u.mat_elem(P, 2, 2, [0 0 1], [0 0 1]);
% % %
% % % [max(max(MelA - MelB)), max(max(KelA - KelB))]
% % % %%
% % % clc;
% % % import FEPack.*
% % % u = pdes.PDEObject;
% % % P = rand(2, 4);
% % % P = [P; rand(1, 4)];
% % % u.mat_elem(P, 3, 3, [0 1 1], [0 1 1])
% %
% %
% % % % %%
% % % clear; clc;
% % % import FEPack.*
% % %
% % % N = 8;
% % % rmesh = meshes.MeshRectangle(1, [0, 1], [0, 1], N, N);
% % % xmin = rmesh.domain('xmin');
% % % xmax = rmesh.domain('xmax');
% % % ymin = rmesh.domain('ymin');
% % % ymax = rmesh.domain('ymax');
% % % cell = rmesh.domain('volumic');
% % %
% % % sp = spaces.PeriodicLagrangeBasis(xmin);
% %
% % % u = pdes.PDEObject;
% % % % ecs = assignEcs((u|xmin), 1.0) & assignEcs((u|xmax), 1.0) & assignEcs((u|ymin), 1.0) & assignEcs((u|ymax), 1.0);
% % % ecs = assignEcs((u|xmin) - (u|xmax), 0.0) & assignEcs((u|ymin) - (u|ymax), 0.0);
% % % ecs.applyEcs;
% % % Q = ecs.P(:, xmin.IdPoints)
% %
% % %%
% % % x = linspace(0, 1, 8);
% % x = [0, sort(rand(1, 8-2)), 1];
% % plot([0, 0], [1 0]); hold on;
% % plot(x, 0, 'o', 'MarkerFaceColor', 'r');
% % ylim([-1 1]);
% %
% % y = mod(x + 0.25*sqrt(2), 1); y(end) = [];
% % for i = 1:length(y)
% %   plot(y(i), 0, 'o', 'MarkerFaceColor', 'b');
% %   pause;
% % end
% %
% % [~, I] = sort([x, y]);
% % % I = I(2:2:end)-length(x);
% % % for i = 1:length(y)
% % %   plot(y(I(i)), 0, 'o', 'MarkerFaceColor', 'b');
% % %   pause;
% % % end
% %
% % %%
% % clear; clc;
% % import FEPack.*
% %
% % N = 8;
% % rmesh = meshes.MeshRectangle(1, [0, 1], [0, 1], N, N);
% % xmin = rmesh.domain('xmin');
% % xmax = rmesh.domain('xmax');
% % ymin = rmesh.domain('ymin');
% % ymax = rmesh.domain('ymax');
% % cell = rmesh.domain('volumic');
% %
% % sp = spaces.FourierBasis(xmin, [0 1]);
% %
% % u = pdes.PDEObject;
% % fun = @(x) cos(2*pi*x);
% % sp.intg_shiftU_V([0 0.5])
% 
% clear; clc;
% import FEPack.*
% 
% N = 32;
% rmesh = meshes.MeshRectangle(1, [0, 1], [0, 1], N, N);
% xmin = rmesh.domain('xmin');
% xmax = rmesh.domain('xmax');
% ymin = rmesh.domain('ymin');
% ymax = rmesh.domain('ymax');
% cellule = rmesh.domain('volumic');
% u = pdes.PDEObject; v = dual(u);
% F = @(x) (1 + 8*pi*pi) * cos(2*pi*x(:, 1)) .* cos(2*pi*x(:, 2)); % cos(2*pi*x(:, 1)); 
% Fr = @(x, y) F([x, y]);
% % G = @(x, y)
% 
% AA = FEPack.pdes.Form.intg(cellule, grad(u)*grad(v) + id(u)*id(v));
% LL = FEPack.pdes.Form.intg(cellule, F * id(v));
% UU = AA \ LL;
% 
% UUr = cos(2*pi*rmesh.points(:, 1)) .* cos(2*pi*rmesh.points(:, 2));
% sqrt((UU - UUr)' * AA * (UU - UUr))
% % sp = FEPack.spaces.FourierBasis(xmin, [0, N/4], false);
% % % sp = FEPack.spaces.PeriodicLagrangeBasis(xmin);
% % % Nb = sp.numBasis;
% % % MM1 = FEPack.pdes.Form.intg_TU_V(xmin, eye(Nb), sp, 'weak evaluation');
% % % MM2 = FEPack.pdes.Form.intg_U_V(xmin);
% % MMphi1 = sp.intg_U_V(xmin);
% % MMphi2 = sp.intg_U_V(xmax);
% 
% % phis_mat = sp.phis(sp.domain.mesh.points(sp.domain.IdPoints, :), 1:sp.numBasis);
% % MM0 = matEF(rmesh);
% % MM0phi = phis_mat' * MM0 * phis_mat;
% 
% % max(max(abs(MMphi1 - MMphi2)))
% 
% 
% 
% % FEPack.pdes.Form.intg(cellule, grad(u)*grad(v))
% % theta = pi/3;
% % omega = 5 + 0.1i;
% % et = [cos(theta) sin(theta)];
% % 
% % AA = FEPack.pdes.Form.intg(cell, (F*(et * grad(u))) * (et * grad(v)) - omega*omega*id(u)*id(v));
% % Kt = FEPack.pdes.Form.intg(cell, (F*(et * grad(u))) * (et * grad(v)));
% % KK = FEPack.pdes.Form.intg(cell, ((F*dx(u)) * dx(v)) + ((F*dy(u)) * dy(v)));
% % [MMr, KKr, Ktr] = matEF(rmesh, Fr);
% % [MM0, KK0, Kt0] = matEF(rmesh);
% % max(max(abs(MM0 + KK0 - FEPack.pdes.Form.intg(cellule, grad(u)*grad(v) + id(u)*id(v)))))
% % 
% % AAr = Ktr - omega*omega*MMr;
% % AA0 = Kt0 - omega*omega*MM0;
% % 
% % display(max(max(abs(AA - AAr))));
% %
% % P = 3*[rand(2, 3); zeros(1, 3)];
% % 
% % % plot([P(1, 1), P(1, 2)], [P(2, 1), P(2, 2)]); hold on;
% % % plot([P(1, 2), P(1, 3)], [P(2, 2), P(2, 3)]);
% % % plot([P(1, 1), P(1, 3)], [P(2, 1), P(2, 3)]);
% % 
% % % Yquad = mapToRel.A * Xquad + mapToRel.B; plot(Yquad(1, :), Yquad(2, :), '*');
% % 
% % F = @(x) x(:, 1).^2 + x(:, 2).^2;% cos(2*pi*x(:, 1));
% % Fr = @(x, y) F([x, y]);
% % 
% % Kel = FEPack.pdes.Form.mat_elem(P, 2, 2, [0 1 0 0], [0 1 0 0], F) +...
% %       FEPack.pdes.Form.mat_elem(P, 2, 2, [0 0 1 0], [0 0 1 0], F);
% % 
% % [Mel, Kelr] = mat_elem(P(1:2, 1), P(1:2, 2), P(1:2, 3), Fr);
% % display(max(max(abs(Kel - Kelr))));
% 
% % F = @(x) x(:, 1).^2;% cos(2*pi*x(:, 1));
% % quadRule1 = FEPack.tools.QuadratureObject.symetrical_Gauss_triangle(6);
% % I1 = quadRule1.weights * F(quadRule1.points.')
% % 
% % quadRule2.points = [1/3 1/5 1/5 3/5; 1/3 1/5 3/5 1/5];
% % quadRule2.weights = [-9/32, 25/96, 25/96, 25/96];
% % I2 = quadRule2.weights * F(quadRule2.points.')



% clear; clc;
import FEPack.*

phi = @(x) exp(2i*pi*x(:, 2));
N = 32;
mesh = meshes.MeshRectangle(0, [0 1], [0 1], N, N);
dom = mesh.domain('xmin');
spB = spaces.FourierBasis(dom, [0 N/4]);

U1 = phi(mesh.points(dom.IdPoints, :));

U2 = spB.phis(mesh.points(dom.IdPoints, :), 1:spB.numBasis) * spB.FE_to_spectral * U1;

x = linspace(0, 1, size(U1, 1));
plot(x, imag(U1)); hold on;
plot(x, imag(U2));








%%
U1 = exp(2i*pi*mesh_int.points(mesh_int.domain('xmax').IdPoints, 2));
% U1 = U.interior(mesh_int.domain('xmax').IdPoints);
U2 = spBint_pos.phis(mesh_int.points(mesh_int.domain('xmax').IdPoints, :), 1:spBint_pos.numBasis) * spBint_pos.FE_to_spectral * U1;
x = linspace(0, 1, size(U1, 1));
plot(x, imag(U1)); hold on;
plot(x, imag(U2));

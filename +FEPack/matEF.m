function [M, K, Ktheta, S0, S1] = matEF(maillage)
% MATEF Assemblage des matrices elements finis.

II = zeros(maillage.numTriangles*9,1);
JJ = zeros(maillage.numTriangles*9,1);

VV_M  = zeros(maillage.numTriangles*9, 1);
VV_K  = zeros(maillage.numTriangles*9, 1);
VV_Kt = zeros(maillage.numTriangles*9, 1);

index_II = [1 2 3 1 2 3 1 2 3];
index_JJ = [1 1 1 2 2 2 3 3 3];


for t = 1:maillage.numTriangles

    S1 = maillage.points(maillage.triangles(t,1),1:2);
    S2 = maillage.points(maillage.triangles(t,2),1:2);
    S3 = maillage.points(maillage.triangles(t,3),1:2);

    % calcul des matrices elementaires du triangle t
    [Mloc, Kloc, KTloc] = mat_elem(S1, S2, S3);

    index = (9*(t-1)+1):(9*(t-1)+9);
    II(index) = maillage.triangles(t,index_II);
    JJ(index) = maillage.triangles(t,index_JJ);

    VV_M(index)  = Mloc(:);
    VV_K(index)  = Kloc(:);
    VV_Kt(index) = KTloc(:);


end

fprintf('\n');

M = sparse(II, JJ, VV_M, maillage.numPoints, maillage.numPoints);
K = sparse(II, JJ, VV_K, maillage.numPoints, maillage.numPoints);
Ktheta = sparse(II, JJ, VV_Kt, maillage.numPoints, maillage.numPoints);

% % Matrices de surface avec translation
% % S_{ij} = \int_0^1 w_j(x - tau) w_i(x) dx
% pts0 = maillage.points(maillage.Refneu == 1, 1); N0 = length(pts0);
% pts1 = maillage.points(maillage.Refneu == 3, 1); N1 = length(pts1);
% pts0 = sort([pts0; mod(pts0 - params.tau, 1)]);
% pts1 = sort([pts1; mod(pts1 - params.tau, 1)]);
%
% x0 = [pts0; 0.5*(pts0(1:end-1) + pts0(2:end))];
% x1 = [pts1; 0.5*(pts1(1:end-1) + pts1(2:end))];
% w0 = [1/6; ones(2*N0-2,1)/3; 1/6; 2*ones(2*N0-1, 1)/3];
% w1 = [1/6; ones(2*N1-2,1)/3; 1/6; 2*ones(2*N1-1, 1)/3];
%
% S0 = sparse(maillage.numPoints, maillage.numPoints);
% S1 = sparse(maillage.numPoints, maillage.numPoints);
%
% idS0 = find(maillage.Refneu == 1); [idS0, cle0] = sort(idS0);
% idS1 = find(maillage.Refneu == 3); [idS1, cle1] = sort(idS1);
%
% Q = @(x, idI, idJ) lagrangeP1(x - params.tau, idJ, maillage.points) .*...
%                    lagrangeP1(x,              idI, maillage.points);
%
% % VV_S0 = zeros(N0, N0);
% % VV_S1 = zeros(N1, N1);
% % for idI = 1:N0, for idJ = 1:N0, VV_S0(idI, idJ) = w0' * Q(x0, idI, idJ);end; end
% % for idI = 1:N1, for idJ = 1:N1, VV_S1(idI, idJ) = w1' * Q(x1, idI, idJ);end; end
% for idI = 1:N0
%   for idJ = 1:N0
%     S0(idS0(idI), idS0(idJ)) = w0' * Q(x0, (cle0(idI)), (cle0(idJ)));
%   end
% end
%
% for idI = 1:N1
%   for idJ = 1:N1
%     S1(idS1(idI), idS1(idJ)) = w1' * Q(x1, (cle1(idI)), (cle1(idJ)));
%   end
% end
%
%
% % S1(idI, idJ) = w1' * Q(x1, idI, idJ);
% % S0 = sparse()
end

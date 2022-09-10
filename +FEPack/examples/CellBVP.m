function U = CellBVP(mesh, bilinearIntg, linearIntg, ecs)
  % CellBVP(mesh, bilinearIntg, linearIntg, ecs) computes the solution of
  %
  %       Find u in V such that tau(u) = g, and
  %                A(u, v) = L(v) for any v in V such that tau(v) = 0,
  %
  % where tau is a 0-order linear operator.
  %
  % INPUTS: * mesh, Mesh object,
  %         * bilinearIntg, FE matrix associated to the bilinear form before
  %           the elimination process. bilinearIntg is an N-by-N matrix, where
  %           N is the number of degrees of freedom.
  %         * linearIntg, matrix associated to the linear forms before the
  %           elimination process. linearIntg is an N-by-mV matrix, where
  %           mV is the number of volumic right-hand-sides
  %         * ecs, EssentialConditions object. ecs might have mS surfacic
  %           right-hand sides
  %
  % OUTPUTS: * U, a N-by-mV*mS matrix containing the solution for each source
  %            and each surfacic data.
  %            For iS in {1,...,mS} and each iV in {1,...,mV}, the row index
  %            I = sub2ind(iS, iV) returns the solution U(:, I) computed for
  %            the iV-th source and for the iS-th surfacic data.

  % Default argument
  if (nargin < 4)
    ecs = FEPack.pdes.EssentialConditions;
  end
  if (nargin < 3)
    linearIntg = 0.0;
  end


  % FE matrix and linear form before elimination
  N = mesh.numPoints;
  AA = bilinearIntg;
  LL = linearIntg;

  if (LL == 0)
    LL = sparse(N, 1);
  end

  % Elimination
  if (size(ecs.C, 1) == 0)
    % There is no essential condition
    ecs.P = speye(N);
    ecs.b = sparse(N, 1);
  elseif (size(ecs.P, 1) == 0)
    % The projection matrix has not been computed yet
    ecs.applyEcs;
  end

  % Surfacic contribution to right-hand side
  mV = size(LL, 2);
  mS = size(ecs.b, 2);
  LL = kron(LL, ones(1, mS)) - kron(ones(1, mV), AA * ecs.b);

  % Eliminated FE matrix and right-hand side
  AA0 = ecs.P * AA * ecs.P';
  LL0 = ecs.P * LL;

  % Compute the solution
  U = ecs.P' * (AA0 \ LL0) + kron(ones(1, mV), ecs.b);

end

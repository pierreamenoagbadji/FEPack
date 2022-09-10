function [U, U0, dU0] = PeriodicHalfGuideBVP(mesh, orientation, infiniteDirection, volBilinearIntg, BoundaryStruct, numCells, opts)
  % PeriodicHalfGuideBVP(mesh, orientation, infiniteDirection, ...
  %                      volBilinearIntg, BoundaryStruct, numCells, opts) solves
  %
  %       Find u in V such that tau(u) = g, and
  %                A(u, v) + S(u, v) = Sf(v) for any v in V such that tau(v) = 0,
  %
  % where tau is a 0-order linear operator.
  %
  % INPUTS: * mesh, Mesh object.
  %
  %         * orientation, a number in {-1, 1} that represent whether the guide
  %           is infinite along a positive or a negative axis.
  %
  %         * infiniteDirection, direction along which the guide is infinite.
  %           infiniteDirection is between 1 and the mesh dimension.
  %
  %         * bilinearIntg, integrand associated to the volumic bilinear form.
  %            For instance, bilinearIntg = grad(u)*grad(v) - (omega^2)*id(u)*id(v).
  %
  %         * BoundaryStruct, a struct with multiple fields:
  %             - basisfun0, SpectralBasis object, set of basis functions
  %               associated to {x_i = 0} (if the guide is infinite along x_i).
  %             - basisfun1, same as basisfun0, but with {x_i = 1}. NOTE that
  %               basisfun0 and basisfun1 should correspond to the same
  %               basis on different domains.
  %             - BCdu, BCu and phi are such that the solution satisfies the
  %               boundary condition
  %                            BCdu * (du/dn) + BCu * u = phi.
  %               -> BCdu is a scalar.
  %               -> BCu can be either
  %                   ** a scalar,
  %                   ** a function_handle, or
  %                   ** an Nb-by-Nb matrix which represents an operator in the
  %                      spectral space spanned by the Nb basisfun0 functions.
  %                      NOTE: in this case, one has to specify if BCu is a
  %                      weak L2 representation or a matrix representation of
  %                      the operator, by means of BoundaryStruct.representation.
  %               -> phi is a function_handle.
  %             - representation: required if BCu is a matrix. It can be set
  %               either to 'weak evaluation' or 'projection' depending on
  %               how BCu is defined:
  %               -> if BCu_{i, j} = <Op(phi_j), phi_i>_L2, then 'weak-evaluation'
  %               -> If Op(phi_j) = BCu_{1, j} phi_1 + ... + BCu_{N, j} phi_N,
  %                  then 'projection'.
  %
  %         * numCells, Number of cells on which the solution is computed.
  %
  %         * opts, structure with fields
  %           -> omega, the frequency (needed to compute the flux);
  %           -> suffix for post-processing.
  %
  % OUTPUTS: * U, a N-by-numCells matrix containing the solution in each cell of periodicity.
  %
  %          * U0, the trace of U on the boundary {x_i = 0} (if the guide is
  %            infinite along x_i).
  %
  %          * dU0, the normal trace of U on the boundary {x_i = 0}.

  %% % ************************* %
  %  % Preliminary verifications %
  %  % ************************* %
  if (abs(orientation) ~= 1)
    % Produce an error if an inappropriate value is given for the
    % sign of the outward normal
    error('%d a ete donne comme orientation du demi-guide au lieu de -1 ou de 1.', orientation);
  end

  if min((1:mesh.dimension) ~= infiniteDirection)
    % The component along which the guide is infinite should be between
    % 1 and the mesh dimension
    error('Le demi-guide ne peut pas être infini dans la direction %d.', infiniteDirection);
  end

  %% % ************** %
  %  % Initialization %
  %  % ************** %
  N = mesh.numPoints;
  AA = FEPack.pdes.Form.intg(mesh.domain('volumic'), volBilinearIntg);

  % Find bounded directions and impose periodic condition on them
  boundedDirections = find((1:mesh.dimension) ~= infiniteDirection); % Bounded directions

  u = FEPack.pdes.PDEObject;
  v = dual(u);
  ecs = FEPack.pdes.EssentialConditions;
  for idI = 1:length(boundedDirections)
    ecs = ecs & assignEcs((u|mesh.domains{2*boundedDirections(idI)-1}) - (u|mesh.domains{2*boundedDirections(idI)}), 0.0);
  end

  % Domains
  Sigma0 = mesh.domains{2*infiniteDirection};
  Sigma1 = mesh.domains{2*infiniteDirection-1};

  % Boundary basis functions
  basisfun0 = BoundaryStruct.basisfun0;
  basisfun1 = BoundaryStruct.basisfun1;
  Nb = basisfun0.numBasis; % Number of basis functions

  %% % ***************************** %
  %  % Solve the local cell problems %
  %  % ***************************** %
  fprintf('1. Resolution des problemes de cellule locaux\n');
  % Periodic cell problems whose boundary conditions are similar to the one of
  % the half-guide problem

  % Surfacic right-hand sides
  B0 = sparse(N, Nb);
  if (basisfun0.is_interpolated)
    B0(Sigma0.IdPoints, :) = basisfun0.phis;
  else
    B0(Sigma0.IdPoints, :) = basisfun0.phis(mesh.points(Sigma0.IdPoints, :), 1:Nb);
  end

  B1 = sparse(N, Nb);
  if (basisfun1.is_interpolated)
    B1(Sigma1.IdPoints, :) = basisfun1.phis;
  else
    B1(Sigma1.IdPoints, :) = basisfun1.phis(mesh.points(Sigma1.IdPoints, :), 1:Nb);
  end

  if (abs(BoundaryStruct.BCdu) < eps)

    % Dirichlet boundary conditions %
    % ***************************** %
    fprintf('\tConditions de Dirichlet\n');

    % Homogeneous Dirichlet condition
    ecs = ecs & assignEcs(u|Sigma0, 0.0) & assignEcs(u|Sigma1, 0.0);
    ecs.applyEcs;

    % Add the surfacic contributions
    ecs.b = [B0, B1];

    % Solve the local cell problems
    Ecell = CellBVP(mesh, AA, 0.0, ecs);

    % Deduce the local cell solutions.
    E0 = Ecell(:, 1:Nb);
    E1 = Ecell; E1(:, 1:Nb) = [];

    fprintf('2. Traces et traces normales\n');
    % Traces and normal traces of the local cell solutions
    % Ekl is the trace of Ek on Sigmal
    E00 = speye(Nb);
    E10 = sparse(Nb, Nb);
    E01 = sparse(Nb, Nb);
    E11 = speye(Nb);

    % Fkl is the normal trace of Ek on Sigmal. It is evaluated weakly
    F00 = E0' * AA * E0;
    F10 = E0' * AA * E1;
    F01 = E1' * AA * E0;
    F11 = E1' * AA * E1;

  else

    % Robin (Neumann) boundary conditions %
    % *********************************** %
    fprintf('\tConditions de Robin-Neumann\n');

    % Surfacic contributions
    if isfield(BoundaryStruct, 'representation')
      representation = BoundaryStruct.representation;
    else
      representation = [];
    end
    SS0 = FEPack.pdes.Form.intg_TU_V(Sigma0, BoundaryStruct.BCu, basisfun0, representation);
    SS1 = FEPack.pdes.Form.intg_TU_V(Sigma1, BoundaryStruct.BCu, basisfun1, representation);
    AA = AA + SS0 + SS1;

    % Surfacic right-hand sides
    BB = [FEPack.pdes.Form.intg(Sigma0, id(u)*id(v)) * B0,...
          FEPack.pdes.Form.intg(Sigma1, id(u)*id(v)) * B1];

    % Solve the local cell problems
    Ecell = CellBVP(mesh, AA, BB, ecs);

    % Deduce the local cell solutions.
    E0 = Ecell(:, 1:Nb);
    E1 = Ecell; E1(:, 1:Nb) = [];

    fprintf('2. Traces et traces normales\n');
    % Traces and normal traces of the local cell solutions
    % Ekl is the trace of Ek on Sigmal
    E00 = basisfun0.FE_to_spectral * E0(Sigma0.IdPoints, :);
    E10 = basisfun1.FE_to_spectral * E1(Sigma0.IdPoints, :);
    E01 = basisfun0.FE_to_spectral * E0(Sigma1.IdPoints, :);
    E11 = basisfun1.FE_to_spectral * E1(Sigma1.IdPoints, :);

    % Fkl is the normal trace of Ek on Sigmal.
    invBCdu = 1.0 / BoundaryStruct.BCdu;
    F00 = -invBCdu * (BoundaryStruct.BCu * E00 - speye(Nb));
    F10 = -invBCdu *  BoundaryStruct.BCu * E10;
    F01 = -invBCdu *  BoundaryStruct.BCu * E01;
    F11 = -invBCdu * (BoundaryStruct.BCu * E11 - speye(Nb));

  end

  %% % ************************** %
  %  % Solve the Riccati equation %
  %  % ************************** %
  fprintf('3. Resolution du système de Riccati\n');
  % Solve the linearized eigenvalue problem associated to the Riccati equation
  flux = @(V) modesFlux(V, E00, E10, F00, F10, orientation, basisfun0, opts.omega);
  riccatiOpts.tol = 1.0e-2;
  if isfield(opts, 'suffix')
    riccatiOpts.suffix = opts.suffix;
  else
    riccatiOpts.suffix = '';
  end
  [R, D] = propagationOperators([E01, E11;  orientation*F01,  orientation*F11], ...
                                [E00, E10; -orientation*F00, -orientation*F10], flux, riccatiOpts);

  %% % ********************************* %
  %  % Compute the solution cell by cell %
  %  % ********************************* %
  fprintf('4. Construction de la solution\n');
  U = zeros(N, numCells);
  R0Phi = basisfun0.FE_to_spectral * BoundaryStruct.phi(mesh.points(Sigma0.IdPoints, :));
  R1Phi = D * R0Phi;

  for idCell = 0:numCells-1
    % Compute the solution in the current cell
    U(:, idCell + 1) = E0 * R0Phi + E1 * R1Phi;

    % Update the coefficients
    R0Phi = R * R0Phi;
    R1Phi = D * R0Phi;
  end

  %% % *********************************************** %
  %  % The transmission coefficient and the derivative %
  %  % *********************************************** %
  fprintf('5. Calcul de l''operateur de transmission global\n');
  % Compute the trace and the normal trace of the half-guide solution
  U0 = E00 + E10 * D;
  dU0 = F00 + F10 * D;

end

%% % ******** %
%  % Appendix %
%  % ******** %
% The flux function
function flux = modesFlux(V, E00, E10, F00, F10, orientation, basisfun, omega)

  Psi  = E00 + E10 * V;
  dPsi = F00 + F10 * V;

  flux = -orientation * imag((Psi' * basisfun.massmat * dPsi) / omega);

end

function U = PeriodicGuideJumpBVP(semiInfiniteDirection,...
                                  volBilinearIntg_pos, mesh_pos, BCstruct_pos, numCells_pos,...
                                  volBilinearIntg_neg, mesh_neg, BCstruct_neg, numCells_neg,...
                                  jumpData, opts)

  % PeriodicSpaceBVP
  % warning('Uniquement Dirichlet (saut Neumann)');

  % Preliminary set ups
  opts.computeSol = false;

  % Solve the problem in the positive half-guide
  % ////////////////////////////////////////////
  [~, E0pos, E1pos, Rpos, Dpos, BCstruct_pos, Lambda_pos] = PeriodicHalfGuideBVP(mesh_pos, +1, semiInfiniteDirection, volBilinearIntg_pos, BCstruct_pos, numCells_pos, opts);

  % Solve the problem in the negative half-guide
  % ////////////////////////////////////////////
  [~, E0neg, E1neg, Rneg, Dneg, BCstruct_neg, Lambda_neg] = PeriodicHalfGuideBVP(mesh_neg, -1, semiInfiniteDirection, volBilinearIntg_neg, BCstruct_neg, numCells_neg, opts); %#ok

  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % N = 32;
  % gamma = sqrt(4*pi*pi*(-N/4:N/4).^2 - opts.omega^2);
  % Lambda_pos = diag(gamma);
  % Lambda_neg = diag(gamma);
  % Nbpos = BCstruct_pos.spB0.numBasis;
  % Nbneg = BCstruct_neg.spB0.numBasis;
  %
  % for i = 1:numCells_pos
  %   Upos{i} = BCstruct_pos.spB0.phis(mesh_pos.points, 1:Nbpos) .* exp(-(mesh_pos.points(:, 1) + (i-1))*gamma);
  % end
  % for i = 1:numCells_neg
  %   Uneg{i} = BCstruct_neg.spB0.phis(mesh_neg.points, 1:Nbneg) .* exp((mesh_neg.points(:, 1) - (i-1))*gamma);
  % end
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Solve the problem at the interface
  % //////////////////////////////////
  Sigma0pos = mesh_pos.domains{2*semiInfiniteDirection};
  % Sigma0neg = mesh_neg.domains{2*semiInfiniteDirection};
  GG = BCstruct_pos.spB0.FE_to_spectral * jumpData(mesh_pos.points(Sigma0pos.IdPoints, :));
  
  % -------------------------------------------------------------------- %
  % The minus sign comes from the definition of the Lambda when they are %
  % the DtN operators (see PeriodicHalfGuideBVP.m)                       %
  % -------------------------------------------------------------------- %
  solphi = - (Lambda_pos + Lambda_neg) \ GG;  % see comment above        %
  % solphi = jumpData(mesh_pos.points(Sigma0pos.IdPoints, :));
  % solphi(end) = [];
  % -------------------------------------------------------------------- %
  % cond(Lambda_neg + Lambda_pos)
  
  % Construct the solution in the whole domain
  % //////////////////////////////////////////
  % Positive side
  U.positive = zeros(mesh_pos.numPoints, numCells_pos);
  R0Phi = solphi;
  R1Phi = Dpos * R0Phi;

  for idCell = 0:numCells_pos-1
    % Compute the solution in the current cell
    U.positive(:, idCell + 1) = E0pos * R0Phi + E1pos * R1Phi;

    % Update the coefficients
    R0Phi = Rpos * R0Phi;
    R1Phi = Dpos * R0Phi;
  end

  % for idCell = 1:numCells_pos
  %   U.positive(:, idCell) = Upos{idCell} * solphi;
  % end

  % Negative side
  U.negative = zeros(mesh_neg.numPoints, numCells_neg);
  R0Phi = solphi;
  R1Phi = Dneg * R0Phi;

  for idCell = 0:numCells_neg-1
    % Compute the solution in the current cell
    U.negative(:, idCell + 1) = E0neg * R0Phi + E1neg * R1Phi;

    % Update the coefficients
    R0Phi = Rneg * R0Phi;
    R1Phi = Dneg * R0Phi;
  end
  % for idCell = 1:numCells_neg
  %   U.negative(:, idCell) = Uneg{idCell} * solphi;
  % end

  % % Trace and normal trace of U
  % if (BCstruct_pos.spB0.is_interpolated)
  %   Uint = BCstruct_pos.spB0.phis * solphi;
  %   dUint.positive = BCstruct_pos.spB0.phis * Lambda_pos * solphi;
  %   dUint.negative = BCstruct_neg.spB0.phis * Lambda_neg * solphi;
  % else
  %   Sigma0neg = mesh_neg.domains{2*semiInfiniteDirection};

  %   Uint = BCstruct_pos.spB0.phis(mesh_pos.points(Sigma0pos.IdPoints, :), 1:BCstruct_pos.spB0.numBasis) * solphi;
  %   dUint.positive = BCstruct_pos.spB0.phis(mesh_pos.points(Sigma0pos.IdPoints, :), 1:BCstruct_pos.spB0.numBasis) * Lambda_pos * solphi;
  %   dUint.negative = BCstruct_neg.spB0.phis(mesh_pos.points(Sigma0neg.IdPoints, :), 1:BCstruct_neg.spB0.numBasis) * Lambda_neg * solphi;
  % end
end

function symBC = symNeumannBC(dimension, u, varargin)
    % Homogeneous Neumann conditions
    symBC = FEPack.applications.symBC.SymBoundaryCondition;

    for idI = 1:dimension

        symBC = symBC & ((dn(u)|varargin{2*idI-1}) == 0) & ((dn(u)|varargin{2*idI}) == 0);

    end
end
function symBC = symDirichletBC(dimension, u, varargin)
    % Homogeneous Dirichlet conditions
    symBC = FEPack.applications.symBC.SymBoundaryCondition;

    for idI = 1:dimension

        symBC = symBC & ((u|varargin{2*idI-1}) == 0) & ((u|varargin{2*idI}) == 0);

    end
end
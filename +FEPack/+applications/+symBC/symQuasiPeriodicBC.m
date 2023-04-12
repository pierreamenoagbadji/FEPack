function symBC = symQuasiPeriodicBC(dimension, u, FBvar, varargin)
    % Homogeneous quasiperiodic conditions
    symBC = FEPack.applications.symBC.SymBoundaryCondition;

    for idI = 1:dimension

        symBC = symBC & (((u|varargin{2*idI-1}) - exp(1i*FBvar)*(u|varargin{2*idI})) == 0) &...
                        (((dn(u)|varargin{2*idI-1}) + exp(1i*FBvar)*(dn(u)|varargin{2*idI})) == 0);

    end
end
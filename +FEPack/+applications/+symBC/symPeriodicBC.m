function symBC = symPeriodicBC(dimension, u, varargin)
    % Homogeneous periodic conditions
    
    symBC = FEPack.applications.symBC.symQuasiPeriodicBC(dimension, u, 0, varargin{:})

end
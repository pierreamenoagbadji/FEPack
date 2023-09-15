load([cheminDonnees, '/inputs.mat']);

fprintf('C. Transmission à l''interface\n');
solGuidePos = load([cheminDonnees, '/half_guide_solution_pos']);
solGuideNeg = load([cheminDonnees, '/half_guide_solution_neg']);

%% Transmission equation
% Right-hand side
numBasisX = size(solGuidePos.Lambda, 1);
Sigma0x = mesh_pos.domain('xmin'); N0x = Sigma0x.numPoints;
u = FEPack.pdes.PDEObject; v = dual(u);
GG = zeros(numBasisX, 1);

MM = FEPack.pdes.Form.intg(Sigma0x, u * v);
MM = MM(Sigma0x.IdPoints, Sigma0x.IdPoints);

%%
Sigma0xPoints = sort(mesh_pos.points(Sigma0x.IdPoints, 2));

for idZ = 1:numCellsZ
  idBasis = (1+(idZ-1)*(N0x-1):1+idZ*(N0x-1))';
  cellZorigin = sizeCellZ * (idZ - 1) + Zorigin;
  Gvec = G([mesh_pos.points(Sigma0x.IdPoints, 1), mesh_pos.points(Sigma0x.IdPoints, 2) + cellZorigin]);

  GG = GG + solGuidePos.B0(idBasis, :)' * MM * Gvec;
end

% Solve linear system
solphi = (solGuidePos.Lambda + solGuideNeg.Lambda) \ GG;

save([cheminDonnees, '/inputs.mat'], '-v7.3');

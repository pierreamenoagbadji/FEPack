import FEPack.*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Compute Dirac point
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HC = applications.HoneycombObject(true, true, 'atomic');

numNodesX = 50;
meshXY = meshes.MeshRectangle(1, [0 1], [0 1], numNodesX, numNodesX);
cellXY = meshXY.domain('volumic');
edge0x = meshXY.domain('xmin'); N0x = edge0x.numPoints;
edge1x = meshXY.domain('xmax'); N1x = edge1x.numPoints;
edge0y = meshXY.domain('ymin'); N0y = edge0y.numPoints;
edge1y = meshXY.domain('ymax'); N1y = edge1y.numPoints;

%% Represent coefficient
numX = 2;
numY = 2;
R = [v1, v2];
transQpot = @(x) HC.fun((R * x(:, 1:2)')');
figure;
for idX = -numX:numX-1
  for idY = -numY:numY-1
    X = meshXY.points(:, 1) + idX;
    Y = meshXY.points(:, 2) + idY;

    % trisurf(meshXY.triangles, X, Y, Qpot([X, Y]));
    trisurf(meshXY.triangles, X, Y, transQpot([X, Y])); 
    hold on;
    shading interp;
    view(2);
    set(gca, 'DataAspectRatio', [1 1 1]);
  end
end
error();
%%
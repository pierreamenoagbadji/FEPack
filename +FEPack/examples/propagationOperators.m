function [R, D] = propagationOperators(A, B, flux, opts)
% function [R, D] = PROPAGATIONOPERATORS(A, B, flux, opts)
% computes the propagation matrices by solving the problem
%
%   Find (R, D) such that A*[Id, D].' = B*[Id, D].'*R
%              + spectral radius / flux constraint
%
% using a spectral approach.
%
% One proceeds by finding the eigenpairs (V, E) such that
%
%   A * V = B * V * E
%
% The eigenpairs (E, V) of the matrix R are then choosen as
% the ones with a modulus strictly less than 1 or with a modulus
% equal to 1 and a negative corresponding flux.
%
% INPUT:  - A, B, 2N x 2N matrices, where N is the size of R,
%         - flux, a function which takes in argument a 2Nx1
%                 vector and returns a Nx1 vector containing the flux
%         - opts, a structure with field(s)
%                 + tol, a tolerance used to classify the eigenvalues
%                 + suffix, for exporting the eigenpairs
%
% OUTPUT: - R and D, N x N matrices

if (nargin < 4)

  % Default options
  opts.tol = 1.0e-2;
  opts.suffix = '';

end

% **************************** %
% Solve the eigenvalue problem %
% **************************** %
[V, E] = eig(full(A), full(B));

% Extract the eigenvalues and normalize the eigenvectors
E = diag(E);
V = V ./ (sqrt(diag(V'*V)).' * ones(size(V, 2), 1));

% ************************ %
% Classify the eigenvalues %
% ************************ %
idGoodEigs = find(abs(E) < 1.0 - opts.tol);
idBadEigs  = find(abs(E) > 1.0 + opts.tol);

isOnUnitCircle = abs(abs(E) - 1.0) < opts.tol;
idEigsImP = find(isOnUnitCircle & imag(E) > eps);
idEigsImN = find(isOnUnitCircle & imag(E) < eps);
idEigsIm0ReP = find(isOnUnitCircle & abs(imag(E)) < eps & real(E) > eps);
idEigsIm0ReN = find(isOnUnitCircle & abs(imag(E)) < eps & real(E) < eps);

% figure;
% % figRD = figure('Position', get(0, 'Screensize'), 'visible', 'off');
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% plot(exp(2i*pi*linspace(0, 1, 512))); hold on;
% plot(E(idGoodEigs), 'bo');
% % plot(E(idBadEigs),  'rx');
% plot(E(isOnUnitCircle),  'g*');
% xlim([-1.1, 1.1]); ylim([-1.1, 1.1]);
% title(['Spectre de $\mathcal{P}$; $N = ', num2str(floor(length(E)/2)), '$']);
% set(gca, 'DataAspectRatio',[1 1 1], 'FontSize', 20);
% % print(figRD, ['outputs/eigsR', num2str(floor(length(E)/2))], '-dpng');
% % pause(0.5);
% % close(figRD);

% Preliminary verifications
strReco = 'Revoir la tolerance ou le nombre de fonctions de base.';

if (length(idGoodEigs) ~= length(idBadEigs))

  % There has to be as many eigenvalues with a modulus < 1
  % as eigenvalues with a modulus > 1
  error(['[propagationOperators.m] ', ...
         'Les valeurs propres de module different de 1 ne peuvent ', ...
         'pas etre mises par paires :\nIl y en a %d de module < 1 ', ...
         'et %d de module > 1.\n', strReco],...
         length(idGoodEigs), length(idBadEigs));

end

if (length(idEigsImP) ~= length(idEigsImN))

  % The eigenvalues on the unit circle should be grouped by conjugate pairs
  error(['[propagationOperators.m] ', ...
         'Les valeurs propres de module 1 et de partie imaginaire ',...
         'non nulle ne peuvent pas etre classees par paires :\nIl ', ...
         'y en a %d a partie imaginaire > 0 et %d a partie ', ...
         'imaginaire < 0.\n', strReco],...
         length(idEigsImP), length(idEigsImN));

end

if (mod(length(idEigsIm0ReP), 2) || mod(length(idEigsIm0ReN), 2))

  % Eigenvalues 1 and -1 should have even multiplicity
  error(['[propagationOperators.m] ', ...
         'Les valeurs propres 1 et -1 ne sont pas de multiplicite ',...
         'paire : \n1 est de multiplicite %d et -1 est de ', ...
         'multiplicité %d.\n', strReco],...
         length(idEigsIm0ReP), length(idEigsIm0ReN));
end

N = floor(length(E)/2);

% *************************** %
% Choose the right eigenpairs %
% *************************** %
if (length(idGoodEigs) == N)

  % Case where the number of eigenvalues with a modulus < 1
  % coincides with the size of the propagation matrices
  IGood = idGoodEigs;

else

  % Some eigenvalues of the propagation operator are
  % on the unit circle.
  % They will be identified using the flux
  % Find the eigenpairs with the negative flux
  idGoodEigsOnUnitCircle = find(isOnUnitCircle & (sign(flux(V)) == -1));

  % Verify that we have the right number of good eigenvalues
  if (length(idGoodEigs) + length(idGoodEigsOnUnitCircle) ~= N)

    error(['[propagationOperators.m] ', ...
           '%d bonnes paires propres ont ete selectionnees. ', ...
           '%d bonnes paires propres etaient attendues'], ...
           length(idGoodEigs) + length(idGoodEigsOnUnitCircle), N);

  end

  IGood = [idGoodEigs;  idGoodEigsOnUnitCircle];

end

% Save good and bad eigenvalues
fid = fopen(['outputs/eigsP_QuasiD_', opts.suffix, '.txt'], 'w');
fprintf(fid, '%0.5e   %0.5e\n', [real(E(IGood)), imag(E(IGood))].');
fclose(fid);

fid = fopen(['outputs/badEigs_QuasiD_', opts.suffix, '.txt'], 'w');
Ebad = E; Ebad(IGood) = [];
fprintf(fid, '%0.5e   %0.5e\n', [real(Ebad), imag(Ebad)].');
fclose(fid);

% Construct the propagation operators
R = V(1:N, IGood) * diag(E(IGood)) * (V(1:N, IGood) \ eye(N));
D = V(N+1:end, IGood) * (V(1:N, IGood) \ eye(N));

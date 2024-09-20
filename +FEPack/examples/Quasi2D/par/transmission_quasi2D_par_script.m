% transmission_quasi2D_par('A', 1, 8, 'omega1', 'constant', 1, 'G1', true, false, true);
% transmission_quasi2D_par('A', 1,  10, 'omega1', 'constant', 1, 'G1', false, true, true);
% transmission_quasi2D_par('A', 1,  20, 'omega4', 'constant', 1, 'G1', false, true, true);

clear; clc;
!rm /UMA/tmp/amenoagbadji/solution_B_irrational_omega3_per11_G1/*
entrees = {'B', 1, 40, 'omega3', 'irrational', 1, 'G1', false, false, true, '1', 0}; 
thetas = linspace(0, 2*pi, 100);
% thetas = pi/4 + [-1e-2 0 1e-2];
% thetas = [0 pi/2];
for idT = 1:length(thetas)
  entrees{11} = int2str(idT);
  entrees{12} = thetas(idT);
  transmission_quasi2D_par_gen_medium
end
!tar -cvf solution.tar -C /UMA/tmp/amenoagbadji/solution_B_irrational_omega3_per11_G1 .

% clear; clc;
% entrees = {'A', 1, 20, 'omega1', 'irrational', 1, 'G1', false, false, true, ''}; transmission_quasi2D_par
% % entrees = {'A', 1, 10, 'omega1', 'irrational', [3, 2], 'G1', false, false, true, ''}; transmission_quasi2D_par
% entrees = {'A', 1, 20, 'omega1', 'irrational', 1, 'G2', false, false, true, ''}; transmission_quasi2D_par

% entrees = {'A', 1, 10, 'omega4', 'constant', 1, 'G1', false, true, true, ''}; transmission_quasi2D_par
% entrees = {'A', 1, 20, 'omega4', 'constant', 1, 'G1', false, true, true, ''}; transmission_quasi2D_par
% entrees = {'A', 1, 30, 'omega4', 'constant', 1, 'G1', false, true, true, ''}; transmission_quasi2D_par
% entrees = {'A', 1, 40, 'omega4', 'constant', 1, 'G1', false, true, true, ''}; transmission_quasi2D_par
% entrees = {'A', 1, 50, 'omega4', 'constant', 1, 'G1', false, true, true, ''}; transmission_quasi2D_par
% entrees = {'A', 1, 60, 'omega4', 'constant', 1, 'G1', false, true, true, ''}; transmission_quasi2D_par
% entrees = {'A', 1, 70, 'omega4', 'constant', 1, 'G1', false, true, true, ''}; transmission_quasi2D_par
% entrees = {'A', 1, 80, 'omega4', 'constant', 1, 'G1', false, true, true, ''}; transmission_quasi2D_par
% entrees = {'A', 1, 90, 'omega4', 'constant', 1, 'G1', false, true, true, ''}; transmission_quasi2D_par


% transmission_quasi2D_par('A', 1,  40, 'omega4', 'constant', 1, 'G1', false, true, true);
% transmission_quasi2D_par('A', 1,  50, 'omega4', 'constant', 1, 'G1', false, true, true);
% transmission_quasi2D_par('A', 1,  60, 'omega4', 'constant', 1, 'G1', false, true, true);
% transmission_quasi2D_par('A', 1,  70, 'omega1', 'constant', 1, 'G1', false, true, true);
% transmission_quasi2D_par('A', 1,  80, 'omega1', 'constant', 1, 'G1', false, true, true);
% transmission_quasi2D_par('A', 1,  90, 'omega1', 'constant', 1, 'G1', false, true, true);
% transmission_quasi2D_par('A', 1, 100, 'omega1', 'constant', 1, 'G1', false, true, true);
% transmission_quasi2D_par('B', 1, 10, 'omega1', 'constant', 1, 'G1', false, false, true);

% transmission_quasi2D_par('A', 1, 10, 'omega1', 'rational', 1, 'G1', false, false, true);
% transmission_quasi2D_par('B', 1, 10, 'omega1', 'rational', 1, 'G1', false, false, true);

% transmission_quasi2D_par('A', 1, 10, 'omega2', 'irrational', 1, 'G1', false, false, true);
% transmission_quasi2D_par('B', 1, 10, 'omega2', 'irrational', 1, 'G1', false, false, true);

% transmission_quasi2D_par('B', 1, 10, 'omega1', 'irrational', 1, 'G1', false, false, true);
% transmission_quasi2D_par('A', 0.075, 10, 'omega1', 'homogeneisation', 1, 'G1', false, false, true, '_eps_0_075');
% transmission_quasi2D_par('A', 0.25, 10, 'omega1', 'homogeneisation', 1, 'G1', false, false, true, '_eps_0_25');
% transmission_quasi2D_par('A', 0.25, 10, 'omega1', 'homogeneisation', 1, 'G1', false, false, true, '_eps_0_25');
% transmission_quasi2D_par('A', 0.1, 10, 'omega1', 'homogeneisation', 1, 'G1', false, false, true, '_eps_0_1');
% transmission_quasi2D_par('A', 0.05, 10, 'omega1', 'homogeneisation', 1, 'G1', false, false, true, '_eps_0_05');
% transmission_quasi2D_par('A', 0.01, 10, 'omega1', 'homogeneisation', 1, 'G1', false, false, true, '_eps_0_01');
% transmission_quasi2D_par('A', 0.01, 10, 'omega1', 'homogeneisation', 1, 'G1', false, false, true, '_eps_0_01');

% transmission_quasi2D_par('A', 1, 10, 'omega1', 'irrational', 1, 'G1', false, false, true);
% transmission_quasi2D_par('B', 1, 10, 'omega1', 'irrational', 1, 'G1', false, false, true);

% transmission_quasi2D_par('A', 1, 10, 'omega1', 'irrational', 2, 'G1', false, false, true);
% transmission_quasi2D_par('B', 1, 10, 'omega1', 'irrational', 2, 'G1', false, false, true);

% transmission_quasi2D_par('A', 1, 10, 'omega1', 'irrational', 1, 'G2', true, false, true);
% transmission_quasi2D_par('B', 1, 10, 'omega1', 'irrational', 1, 'G2', false, false, true);


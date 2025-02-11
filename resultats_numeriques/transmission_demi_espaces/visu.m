% for idI = 1:5
%   idI
%   H = openfig(['last_simus/solution_', int2str(idI), '.fig']);
%   set(H, 'Position', get(0, 'Screensize'), 'visible', 'off');
%   caxis([-2e-2, 2e-2]);
%   axis off;
%   colorbar off;
%   view(2);
%   print(H, ['last_simus/solution_', int2str(idI)], '-dpng');
%   system(['convert last_simus/solution_', int2str(idI), '.png -trim last_simus/solution_', int2str(idI), '.png']);
%   close(H);
% end
prefixe = '0_075';

H = openfig(['last_simus/solution_A_homogeneisation_omega1_per1_G1_eps_', prefixe, '/fig_solution_2D.fig']);
set(H, 'Position', get(0, 'Screensize'), 'visible', 'off');
shading interp;
caxis([-4e-2, 4e-2]);
axis off;
colorbar off;
view(2);
print(H, ['last_simus/solutions_homoge/fig_solution_2D_', prefixe], '-dpng');
system(['convert last_simus/solutions_homoge/fig_solution_2D_', prefixe, '.png -trim last_simus/solutions_homoge/fig_solution_2D_', prefixe, '.png']);
close(H);

H = openfig(['last_simus/solution_A_homogeneisation_omega1_per1_G1_eps_', prefixe, '/fig_erreur_2D.fig']);
set(H, 'Position', get(0, 'Screensize'), 'visible', 'off');
shading interp;
caxis([-3e-2, 3e-2]);
axis off;
colorbar off;
view(2);
print(H, ['last_simus/solutions_homoge/fig_erreur_2D_', prefixe], '-dpng');
system(['convert last_simus/solutions_homoge/fig_erreur_2D_', prefixe, '.png -trim last_simus/solutions_homoge/fig_erreur_2D_', prefixe, '.png']);
close(H);


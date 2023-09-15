function [coord, elem, val] = plotsolution(filenames, key, titre, xmin, xmax)
% ========================================================================
% READ AND PLOT DATA
% ========================================================================

% Valeurs par defaut
if (nargin < 3),    titre = '';   end
if (nargin < 4),    xmin  =  0;   end
if (nargin < 5),    xmax  =  1;   end


if (strcmp(key, '1Dvtu'))

    % Solutions 1D
    [coord, elem, val] = read_vtu(filenames);

%     if (coord(:,2) == 0)
%         xmin = min(coord(:, 1));
%         xmax = max(coord(:, 1));
%     end

    O1 = min(coord);
    O2 = max(coord);

    d  = sqrt( (O2(1) - O1(1))^2 + (O2(2) - O1(2))^2 );

    t = xmin + (xmax - xmin) * sqrt((coord(:, 1) - O1(1)).^2 + (coord(:, 2) - O1(2)).^2) ./ d;

    [t, I] = sort(t);
    val = val(I, :);

    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

    subplot(2, 1, 1); 
    plot(t, val(:, 1), 'color', 'r'); hold on;
    xlabel('$x$', 'Interpreter', 'Latex');
    title([titre, repmat(' -- ', ~strcmp(titre, '')), 'Partie reelle'], 'Interpreter', 'Latex');
    set(gca, 'FontSize', 16);

    subplot(2, 1, 2); 
    plot(t, val(:, 2), 'color', 'r'); hold on;
    xlabel('$x$', 'Interpreter', 'Latex');
    title([titre, repmat(' -- ', ~strcmp(titre, '')), 'Partie imaginaire'], 'Interpreter', 'Latex');
    set(gca, 'FontSize', 16);

end

if (strcmp(key, '1Dtxt'))

    % Solutions 1D enregistrées sous la forme [xi Re(u(xi)) Im(u(xi))]
    u = load(filenames);

    plot(u(:, 1), u(:, 2), 'color', rand(1, 3), 'LineWidth', 1.5);

end

if (strcmp(key, '2D'))

    % Solutions 2D
    [coord, elem, val] = read_vtu(filenames);

    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

    subplot(1, 2, 1); hold off;
    trisurf(elem + 1, coord(:, 1), coord(:, 2), val(:, 1));
    xlabel('$y_1$', 'Interpreter', 'Latex');
    ylabel('$y_2$', 'Interpreter', 'Latex');
    title([titre, repmat(' -- ', ~strcmp(titre, '')), 'Partie reelle'], 'Interpreter', 'Latex');
    % colormap jet;
    view(2);
    shading interp;
    colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

    subplot(1, 2, 2); hold off;
    trisurf(elem + 1, coord(:, 1), coord(:, 2), val(:, 2));
    xlabel('$y_1$', 'Interpreter', 'Latex');
    ylabel('$y_2$', 'Interpreter', 'Latex');
    title([titre, repmat(' -- ', ~strcmp(titre, '')), 'Partie imaginaire'], 'Interpreter', 'Latex');
    % colormap jet;
    view(2);
    shading interp;
    colorbar('TickLabelInterpreter', 'latex');
    set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

end

if (strcmp(key, '2Dguide-moins') || strcmp(key, '2Dguide-plus'))

    % Solution 2D construite sur plusieurs cellules
    if strcmp(key, '2Dguide-plus')
        nY =  1;
    else
        nY = -1;
    end

    for idxI = 1:length(filenames)

        [coord, elem, val{idxI}] = read_vtu(filenames{idxI});

        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');

        subplot(1, 2, 1);
        trisurf(elem + 1, coord(:, 1), coord(:, 2) + nY*(idxI-1), val{idxI}(:, 1));
        hold on;
        xlabel('$y_1$', 'Interpreter', 'Latex');
        ylabel('$y_2$', 'Interpreter', 'Latex');
        title([titre, repmat(' -- ', ~strcmp(titre, '')), 'Partie reelle'], 'Interpreter', 'Latex');
        % colormap jet;
        view(2);
        shading interp;
        colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

        subplot(1, 2, 2);
        trisurf(elem + 1, coord(:, 1), coord(:, 2) + nY*(idxI-1), val{idxI}(:, 2));
        hold on;
        xlabel('$y_1$', 'Interpreter', 'Latex');
        ylabel('$y_2$', 'Interpreter', 'Latex');
        title([titre, repmat(' -- ', ~strcmp(titre, '')), 'Partie imaginaire'], 'Interpreter', 'Latex');
        % colormap jet;
        view(2);
        shading interp;
        colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

    end

end

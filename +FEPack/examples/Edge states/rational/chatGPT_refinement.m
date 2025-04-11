% Quadtree Approximation with Reliable Refinement
blowup_function = @(x,y) 1 ./ ( ((x - 0.5).^2 + (y - 0.5).^2) - 0.2^2 + 1e-6 );
value_threshold = 1e2;  % Threshold for blow-up detection
max_depth = 8;         % Maximum depth
leaves = [0, 0, 1, 1, 0]; % Root node

figure;
set(gcf, 'Position', [100 100 600 600]);
axis equal tight;
box on;

for iter = 1:max_depth
    new_leaves = [];
    
    % Force more splits in early iterations
    force_split_depth = min(2, max_depth-1);
    
    for i = 1:size(leaves, 1)
        node = leaves(i, :);
        xmin = node(1); ymin = node(2);
        xmax = node(3); ymax = node(4);
        depth = node(5);
        
        % Only evaluate function if needed
        if depth <= force_split_depth || depth < max_depth
            corners = [xmin, ymin; xmax, ymin; xmax, ymax; xmin, ymax];
            f_corners = zeros(4, 1);
            for j = 1:4
                f_corners(j) = blowup_function(corners(j,1), corners(j,2));
            end
            has_blowup = any(abs(f_corners) > value_threshold);
        else
            has_blowup = false;
        end
        
        % Split conditions
        if (depth <= force_split_depth) || (has_blowup && depth < max_depth)
            xmid = (xmin + xmax)/2;
            ymid = (ymin + ymax)/2;
            new_depth = depth + 1;
            new_leaves = [new_leaves; 
                         xmin, ymin, xmid, ymid, new_depth;
                         xmid, ymin, xmax, ymid, new_depth;
                         xmin, ymid, xmid, ymax, new_depth;
                         xmid, ymid, xmax, ymax, new_depth];
        else
            new_leaves = [new_leaves; node];
        end
    end
    
    leaves = new_leaves;
    
    % Plot current quadtree
    clf;
    axis equal;
    axis([0 1 0 1]);
    hold on;
    
    % Plot all leaves
    for i = 1:size(leaves, 1)
        node = leaves(i, :);
        rectangle('Position', [node(1) node(2) (node(3)-node(1)) (node(4)-node(2))], ...
                 'EdgeColor', 'k');
    end
    
    % Highlight potential blow-up regions
    if iter == max_depth
        for i = 1:size(leaves, 1)
            node = leaves(i, :);
            if node(5) == max_depth
                corners = [node(1), node(2); node(3), node(2); node(3), node(4); node(1), node(4)];
                f_corners = arrayfun(@(x,y) blowup_function(x,y), corners(:,1), corners(:,2));
                if any(abs(f_corners) > value_threshold)
                    rectangle('Position', [node(1) node(2) (node(3)-node(1)) (node(4)-node(2))], ...
                             'FaceColor', 'r', 'EdgeColor', 'none');
                end
            end
        end
    end
    
    title(sprintf('Iteration %d - Depth %d - %d Nodes', iter, max_depth, size(leaves,1)));
    drawnow;
    pause;  % Clear pause between iterations
end
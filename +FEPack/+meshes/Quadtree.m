%> @file Quadtree.m
%> @brief Contains the meshes.Quadtree class.
% =========================================================================== %
%> @brief Quadtree structure used to locally refine rectangular
% grids based on splitting rules that may involve an indicator function
%>
%> Coded with the help of DeepSeek
% =========================================================================== %
classdef Quadtree < FEPack.FEPackObject
  % FEPack.meshes.Quadtree < FEPack.meshes.Mesh

  properties

    %> @brief Quadtree nodes [xmin, ymin, xmax, ymax, depth]
    leaves

    %> @brief Maximum refinement depth
    maxDepth 

    %> @brief Splitting rules
    splitting_rules

    %> @brief Nx2 matrix of [x,y] evaluation points
    cache_coords

    %> @brief Nx1 vector of function values
    cache_vals 

  end
  
  methods

    function obj = Quadtree(maxDepth, rootBB, rules)
      % function obj = Quadtree(maxDepth, rootBB)
      % Initialize quadtree with root node
      %
      % inputs: * maxDepth (integer), maximum depth
      %         * rootBB 4-sized vector containing 
      %                [xmin, ymin, xmax, ymax, depth]
      %         * rules (function handle) that takes 
      %                 as inputs 
      %               - a Nxp vector representing
      %               values (p = 4 is the number of corners)
      %               - a Nx1 vector containing leaves depths
      %               and returns a Nx1 boolean 
      if (nargin < 3)
        % Default rule: force early split, then split leaves where
        % the indicator function is above a certain threshold
        threshold = 1e2;
        force_depth = min(2, maxDepth);
        blow_up_rule = @(fvals) (any(abs(fvals) > threshold, 2) &...
                                ~all(abs(fvals) > threshold, 2));
        early_split_rule = @(depth) (depth <= force_depth);
        rules = @(fvals, depth) early_split_rule(depth) | blow_up_rule(fvals);
      end
      if (nargin < 2)
        rootBB = [0, 0, 1, 1];
      end

      obj.leaves = [rootBB, 0];
      obj.maxDepth = maxDepth;
      obj.cache_coords = zeros(0, 2);
      obj.cache_vals = zeros(0, 1);
      obj.set_splitting_rules(rules);
    end
    
    function obj = record_value(obj, P, val)
      % function obj = RECORD_VALUE(obj, P, val)
      % Records function evaluation in vector storage
      %
      % INPUTS: * obj, FEPack.meshes.Quadtree object.
      %         * P, a N-by-2 matrix containing the 
      %              coordinates of the points.
      %         * val, a N-by-1 vector containing 
      %                the corresponding values
      
      obj.cache_coords = [obj.cache_coords; P];
      obj.cache_vals   = [obj.cache_vals; val];
    end
    
    function [val, exists] = get_cached(obj, P, almostzero)
      % function [val, exists] = GET_CACHED(obj, P, almostzero)
      %
      % Check if coordinate is included in cache
      %
      % INPUTS: * obj, FEPack.meshes.Quadtree object.
      %         * P, a N-by-2 matrix containing the 
      %              coordinates of the points.
      %         * almostzero (optional) tolerance.
      %
      % OUPUTS: * val, a Ncached-by-1 vector containing
      %           the cached values. Ncached <= N is the
      %           number of points that are already in
      %           the cached coordinates.
      %
      %         * exists, a N-by-1 vector. exists(i) is one
      %           if P(i) is in the cache, and 0 otherwise.         
      
      if (nargin < 3)
        almostzero = 1e-8;
      end

      Pc = obj.cache_coords;
      N  = size(P,  1);
      Nc = size(Pc, 1);

      distance = (P(:, 1) * ones(1, Nc) - ones(N, 1) * Pc(:, 1)').^2 +...
                 (P(:, 2) * ones(1, Nc) - ones(N, 1) * Pc(:, 2)').^2;
      
      % Find points contained in the cache
      matches = (sqrt(distance) < almostzero);
      exists = any(matches, 2);
      
      % Modify value accordingly
      idExists = matches(exists, :)';
      Nex = size(idExists, 2);
      val = obj.cache_vals * ones(1, Nex);
      val = val(idExists);
    end
    
    function obj = refine(obj, indicator_function)
      % Main refinement loop with cached evaluations

      for iter = 1:obj.maxDepth
        newLeaves = [];
        nodes = obj.leaves;
        
        % Find leaves that should be split
        should_split = obj.check_node(nodes, indicator_function);

        % Nodes that are not split are added to the new leaves
        newLeaves = [newLeaves; nodes(~should_split, :)]; %#ok

        % Split leaves that should be
        newLeaves = [newLeaves; obj.split_node(nodes(should_split, :))]; %#ok

        obj.leaves = newLeaves;
        % obj.visualize();
        % pause;
      end
    end
    
    function set_splitting_rules(obj, rules)
      % obj = SET_SPLITTING_RULES(obj, rules)
      % Define the rules used for splitting
      %
      % inputs: * rules (function handle) that takes 
      %                 as inputs 
      %               - a Nxp vector representing
      %               values (p = 4 is the number of corners)
      %               - a Nx1 vector containing leaves depths
      %               and returns a Nx1 boolean 
      obj.splitting_rules = rules;
    end

    function should_split = check_node(obj, nodes, indicator_function)
      % Node evaluation with caching
      xmin  = nodes(:, 1); ymin = nodes(:, 2);
      xmax  = nodes(:, 3); ymax = nodes(:, 4);
      depth = nodes(:, 5);
      
      % Check corners
      corners{1} = [xmin, ymin];
      corners{2} = [xmax, ymin];
      corners{3} = [xmax, ymax];
      corners{4} = [xmin, ymax];

      N  = size(nodes, 1);
      fvals = zeros(N, 4);
      
      for idJ = 1:4
        % Check if function has been already evaluated at corners 
        [cached_val, exists] = obj.get_cached(corners{idJ});

        fvals( exists, idJ) = cached_val;
        fvals(~exists, idJ) = indicator_function(corners{idJ}(~exists, :));
        obj = obj.record_value(corners{idJ}(~exists, :), fvals(~exists, idJ));
      end
      
      should_split = obj.splitting_rules(fvals, depth);
    end
    
    function children = split_node(~, nodes)
      % Splits node into 4 children
      xmidpoints = (nodes(:, 1) + nodes(:, 3))/2;
      ymidpoints = (nodes(:, 2) + nodes(:, 4))/2;
      newDepth = nodes(:, 5) + 1;
      
      N = size(nodes, 1);
      children = zeros(4*N, 5);

      children(1:4:4*N, :) = [nodes(:,1), nodes(:,2), xmidpoints, ymidpoints, newDepth];
      children(2:4:4*N, :) = [xmidpoints, nodes(:,2), nodes(:,3), ymidpoints, newDepth];
      children(3:4:4*N, :) = [nodes(:,1), ymidpoints, xmidpoints, nodes(:,4), newDepth];
      children(4:4:4*N, :) = [xmidpoints, ymidpoints, nodes(:,3), nodes(:,4), newDepth];
    end
    
    function visualize(obj, axis_lim, max_visualization_depth)

      % Visualize quadtree up to a certain depth
      if (nargin < 3)
        max_visualization_depth = obj.maxDepth;
      end

      % figure;
      axis equal;
      % axis([0 1 0 1]);
      hold on;
      
      validId = (obj.leaves(:, 5) <= max_visualization_depth);
      valid_leaves = obj.leaves(validId, :);

      % Plot all nodes
      for idI = 1:size(valid_leaves, 1)
        node = valid_leaves(idI, :);
        rectangle('Position', [node(1) node(2) (node(3)-node(1)) (node(4)-node(2))], 'EdgeColor', 'k');
      end

      title(sprintf('Quadtree Approximation (Depth %d)', obj.maxDepth));
      drawnow;
      if (nargin >= 2) && ~isempty(axis_lim)
        xlim([axis_lim(1), axis_lim(3)]);
        ylim([axis_lim(2), axis_lim(4)]);
      end

    end

    function visualize_cache(obj, axis_lim)
      % Visualize cached values on the quadtree mesh
      % Uses same rendering approach as visualize_function but with cached values
      axis equal;
      hold on;
      
      % Create colormap
      cmap = parula(256);
      
      % Match cache points to leaves and compute mean values
      tic;
      leaf_vals = NaN(size(obj.leaves, 1), 1);
      for idI = 1:size(obj.leaves, 1)
        node = obj.leaves(idI, :);
        
        % Find all cache points within this leaf
        in_leaf = obj.cache_coords(:, 1) >= node(1) & ...
                  obj.cache_coords(:, 1) <= node(3) & ...
                  obj.cache_coords(:, 2) >= node(2) & ...
                  obj.cache_coords(:, 2) <= node(4);
        
        if any(in_leaf)
            leaf_vals(idI) = mean(obj.cache_vals(in_leaf));
        end
      end
      tps = toc;
      fprintf('Time 1: %f seconds.\n', tps);

      % Calculate value range using only valid leaves
      valid_vals = leaf_vals(~isnan(leaf_vals));
      if isempty(valid_vals)
        error('No cached values found in quadtree leaves');
      end
      
      max_val = prctile(valid_vals, 98); % 98th percentile to avoid outliers
      min_val = min(valid_vals);
      
      % Plot each cell with color based on cached values
      tic;
      for idI = 1:size(obj.leaves, 1)
        node = obj.leaves(idI, :);
        if ~isnan(leaf_vals(idI))
          % Normalize to colormap index
          cidx = round(255 * (leaf_vals(idI) - min_val) / (max_val - min_val)) + 1;
          cidx = max(1, min(256, cidx)); % Clamp to valid range
          
          % Draw rectangle
          rectangle('Position', [node(1),  node(2),...
                                 node(3) - node(1),...
                                 node(4) - node(2)],...
                    'FaceColor', cmap(cidx, :), ...
                    'EdgeColor', 'none');
        else
          % Draw gray rectangle
          rectangle('Position', [node(1),  node(2),...
                                 node(3) - node(1),...
                                 node(4) - node(2)],...
                    'FaceColor', [128 128 128]/255, ...
                    'EdgeColor', 'none');
        end
      end
      tps = toc;
      fprintf('Time 2: %f seconds.\n', tps);
      
      if (nargin >= 2)
        xlim([axis_lim(1), axis_lim(3)]);
        ylim([axis_lim(2), axis_lim(4)]);
      end

      % Add colorbar
      colormap(cmap);
      colorbar;
      
      title(sprintf('Cached Values Visualization (%d points)', size(obj.cache_coords, 1)));
      hold off;
    end
  end
end
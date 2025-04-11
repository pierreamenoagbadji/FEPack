%> @file Quadtree.m
%> @brief Contains the meshes.Quadtree class.
% =========================================================================== %
%> @brief Class for rectangular grids with possible local refinement
%>
%> Coded with DeepSeek
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

    %> @brief boolean decision for parfor
    parallel_use

    %> @brief integer number of cores
    parallel_numcores

  end
  
  methods

    function obj = Quadtree(maxDepth, rootBB, rules, parallel_use, parallel_numcores)
      % function obj = Quadtree(maxDepth, rootBB)
      % Initialize quadtree with root node
      %
      % inputs: maxDepth (integer), maximum depth
      %         rootBB 4-sized vector containing 
      %                [xmin, ymin, xmax, ymax, depth]
      %         rules (function handle) that takes 
      %               a Nx1 vector and returns a Nx1
      %               boolean 
      if (nargin < 4)
        parallel_use = false; 
        parallel_numcores = 0;
      end
      if (nargin == 4 && parallel_use)
        error('parallel_use has been set to true; you have to specify the number of cores through parallel_numcores. Use the command feature(''numcores'') to see the number of available cores.');
      end
      if (nargin < 3)
        threshold = 1e2;
        rules = @(fvals) (any(abs(fvals) > threshold) &&...
                         ~all(abs(fvals) > threshold));
      end
      if (nargin < 2)
        rootBB = [0, 0, 1, 1];
      end

      obj.leaves = [rootBB, 0];
      obj.maxDepth = maxDepth;
      obj.cache_coords = zeros(0, 2);
      obj.cache_vals = zeros(0, 1);
      obj.parallel_use = parallel_use;
      obj.parallel_numcores = parallel_numcores;
      obj.set_splitting_rules(rules);
    end
    
    function obj = record_value(obj, x, y, val)
      % function obj = RECORD_VALUE(obj, x, y, val)
      % Records function evaluation in vector storage

      obj.cache_coords(end+1, :) = [x, y];
      obj.cache_vals(end+1) = val;
    end
    
    function [val, exists] = get_cached(obj, x, y)
      % function [val, exists] = GET_CACHED(obj, x, y)
      % Checks cache using 1e-8 tolerance
      
      matches = all(abs(obj.cache_coords - [x, y]) < 1e-8, 2);
      exists = any(matches);
      if exists
        val = obj.cache_vals(matches);
      else
        val = NaN;
      end
    end
    
    function obj = refine(obj, indicator_function, minDepthForce)
      % Main refinement loop with cached evaluations

      if (nargin < 4)
        minDepthForce = 2;
      end

      if (obj.parallel_use)
        pool = gcp('nocreate'); % Get current parallel pool
        if isempty(pool)
          parpool('local', obj.parallel_numcores); % Start a parallel pool if none exists
        end
        fprintf('Using parfor (parallel execution)\n');
        
        for iter = 1:obj.maxDepth
          newLeaves = [];
          force_depth = min(minDepthForce, obj.maxDepth-1); % Force early splits
          
          parfor idI = 1:size(obj.leaves, 1)
            node = obj.leaves(idI,:); %#ok
            should_split = obj.check_node(node, indicator_function, force_depth);
              
            if should_split
              newLeaves = [newLeaves; obj.split_node(node)];
            else
              newLeaves = [newLeaves; node];
            end
          end
          
          obj.leaves = newLeaves;
        end
      else
        for iter = 1:obj.maxDepth
          newLeaves = [];
          force_depth = min(minDepthForce, obj.maxDepth-1); % Force early splits
          
          for idI = 1:size(obj.leaves, 1)
            node = obj.leaves(idI,:);
            should_split = obj.check_node(node, indicator_function, force_depth);
              
            if should_split
              newLeaves = [newLeaves; obj.split_node(node)]; %#ok
            else
              newLeaves = [newLeaves; node]; %#ok
            end
          end
          
          obj.leaves = newLeaves;
        end
      end
    end
    
    function set_splitting_rules(obj, rules)
      % obj = SET_SPLITTING_RULES(obj, rules)
      % Define the rules used for splitting
      %
      % inputs: rules (function handle) that
      % take a Nx1 vector and returns a Nx1
      % boolean 
      obj.splitting_rules = rules;
    end

    function should_split = check_node(obj, node, indicator_function, force_depth)
      % Node evaluation with caching
      xmin  = node(1); ymin = node(2);
      xmax  = node(3); ymax = node(4);
      depth = node(5);
      
      % Check corners
      corners = [xmin, ymin; xmax, ymin; xmax, ymax; xmin, ymax];
      fvals = zeros(4, 1);
      
      for idJ = 1:4
        % Check if function has been already evaluated at corners 
        [cached_val, exists] = obj.get_cached(corners(idJ, 1), corners(idJ, 2));

        if (exists)
          fvals(idJ) = cached_val;
        else
          fvals(idJ) = indicator_function(corners(idJ, 1), corners(idJ, 2));
          obj =  obj.record_value(corners(idJ, 1), corners(idJ, 2), fvals(idJ));
        end
      end
      
      should_split = (depth <= force_depth) || (obj.splitting_rules(fvals));
    end
    
    function children = split_node(~, node)
      % Splits node into 4 children
      xmid = (node(1) + node(3))/2;
      ymid = (node(2) + node(4))/2;
      newDepth = node(5) + 1;
      
      children = [
        node(1), node(2), xmid, ymid, newDepth
        xmid, node(2), node(3), ymid, newDepth
        node(1), ymid, xmid, node(4), newDepth
        xmid, ymid, node(3), node(4), newDepth
      ];
    end
    
    function visualize(obj, axis_lim)

      % Visualize quadtree
      % figure;
      axis equal;
      % axis([0 1 0 1]);
      hold on;
      
      % Plot all nodes
      for idI = 1:size(obj.leaves, 1)
        node = obj.leaves(idI, :);
        rectangle('Position', [node(1) node(2) (node(3)-node(1)) (node(4)-node(2))], 'EdgeColor', 'k');
      end

      title(sprintf('Quadtree Approximation (Depth %d)', obj.maxDepth));
      drawnow;
      if (nargin >= 2)
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
      
      % Calculate value range using only valid leaves
      valid_vals = leaf_vals(~isnan(leaf_vals));
      if isempty(valid_vals)
        error('No cached values found in quadtree leaves');
      end
      
      max_val = prctile(valid_vals, 95); % 95th percentile to avoid outliers
      min_val = min(valid_vals);
      
      % Plot each cell with color based on cached values
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
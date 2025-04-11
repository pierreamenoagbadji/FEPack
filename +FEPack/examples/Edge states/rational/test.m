qt = FEPack.meshes.BlowupQuadtree(1);
blowup_fun = @(x,y) 1./(abs(x.^2 + y.^2 - 0.25) + 1e-8);
qt = qt.refine(blowup_fun, 50);
figure;
qt.visualize_cache();

qt.visualize();

% % Configure quadtree
% max_depth = 10;  % Increase for finer resolution
% threshold = 10; % Lower for more sensitive detection
% 
% % Initialize figure
% figure('Position', [100 100 1200 600]);
% 
% % Example 1: Cardioid
% subplot(3, 3, 1);
% cardioid = @(x,y) 0.15*(1 + sin(atan2(y-0.5,x-0.5))) - sqrt((x-0.5).^2 + (y-0.5).^2);
% fimplicit(cardioid, [0 1 0 1], 'LineWidth', 2);
% title('Exact Cardioid Curve');
% axis equal;
% 
% subplot(3, 3, 4);
% qt = FEPack.meshes.BlowupQuadtree(max_depth);%, [-pi, pi, 0, 10]);
% qt = qt.refine(@(x,y) 1./(abs(cardioid(x,y)) + 1e-8), threshold);
% qt.visualize(@(x,y) 1./(abs(cardioid(x,y)) + 1e-8), threshold);
% title('Quadtree Approximation');
% 
% subplot(3, 3, 7);
% qt.visualize_cache;
% 
% % % Example 2: Double Spiral (corrected visualization)
% % subplot(3,3,2);
% % theta1 = @(x,y) atan2(y-0.3,x-0.3);
% % theta2 = @(x,y) atan2(y-0.7,x-0.7);
% % spiral1 = @(x,y) sqrt((x-0.3).^2 + (y-0.3).^2) - 0.1*(2+theta1(x,y)/pi);
% % spiral2 = @(x,y) sqrt((x-0.7).^2 + (y-0.7).^2) - 0.1*(2-theta2(x,y)/pi);
% % % Plot each spiral separately
% % fimplicit(spiral1, [0 1 0 1], 'LineWidth', 1.5, 'Color', 'b');
% % hold on;
% % fimplicit(spiral2, [0 1 0 1], 'LineWidth', 1.5, 'Color', 'r');
% % title('Exact Double Spiral');
% % axis equal;
% % hold off;
% 
% % subplot(3,3,5);
% % qt = FEPack.meshes.BlowupQuadtree(max_depth);
% % spiral_fn = @(x,y) min(abs(spiral1(x,y)), abs(spiral2(x,y)));
% % qt = qt.refine(@(x,y) 1./(spiral_fn(x,y) + 1e-8), threshold);
% % qt.visualize(@(x,y) 1./(spiral_fn(x,y) + 1e-8), threshold);
% % title('Quadtree Approximation');
% 
% % subplot(3, 3, 8);
% % qt.visualize_function(spiral_fn); % New blowup function visualization
% 
% % % 
% % % Example 3: Lemniscate (corrected visualization)
% % subplot(3, 3, 3);
% % lemniscate = @(x,y) ((x-0.5).^2 + (y-0.5).^2).^2 - ((x-0.5).^2 - (y-0.5).^2)/5;
% % fimplicit(lemniscate, [0 1 0 1], 'LineWidth', 2);
% % title('Exact Lemniscate');
% % axis equal;
% 
% % subplot(3, 3, 6);
% % qt = FEPack.meshes.BlowupQuadtree(max_depth);
% % qt = qt.refine(@(x,y) 1./(abs(lemniscate(x,y)) + 1e-8), threshold);
% % qt.visualize(@(x,y) 1./(abs(lemniscate(x,y)) + 1e-8), threshold);
% % title('Quadtree Approximation');
% 
% % subplot(3, 3, 9);
% % qt.visualize_function(lemniscate); % New blowup function visualization
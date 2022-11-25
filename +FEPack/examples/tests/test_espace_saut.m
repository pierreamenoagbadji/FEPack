clear; clc;
import FEPack.*

% profile ON
opts.omega = 8 + 0.1i;
coSemiInf = 1;
coInf = 2;

mu_pos = @(x) 1 + 0.25*cos(2*pi*x(:, 1)) + 0.25*cos(2*pi*x(:, 2));
rho_pos = @(x) 1 + 0.5*cos(2*pi*x(:, 1)).*sin(2*pi*x(:, 2));

mu_neg = @(x)  1 + 0.5*sin(2*pi*x(:, 1)).*sin(2*pi*x(:, 2));% ones(size(x, 1), 1);
rho_neg = @(x) 1 + 0.25*sin(2*pi*x(:, 1)) + 0.25*cos(2*pi*x(:, 2)); % ones(size(x, 1), 1);

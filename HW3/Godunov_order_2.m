% Traffic Flow Simulation with Second-Order Godunov Scheme
% This script implements the solution of the traffic flow equation
% using the second-order Godunov scheme with slope limiters.

clear all;
close all;
clc;

addpath Functions
%% Parameters
L = 1;          % Length of the domain [km]
T = 1;          % Total simulation time [s]
Nx = 100;       % Number of spatial cells
Nt = 500;       % Number of time steps

dx = L / Nx;    % Spatial step size
dt = T / Nt;    % Time step size

x = linspace(dx/2, L-dx/2, Nx)';  % Spatial grid (cell centers)
t = linspace(0, T+dt, Nt);           % Time grid

rho_max = 1;    % Maximum density (normalized)
u_max = 1;      % Maximum velocity (normalized)

% Check CFL condition
CFL = u_max * dt / dx;
if CFL > 0.5
    error('CFL condition for 2nd order method violated! CFL = %.4f > 0.5. Reduce dt or increase dx.', CFL);
else
    fprintf('CFL condition satisfied: %.4f < 0.5\n', CFL);
end

%% Define flux function and its derivatives
f = @(rho) u_max * (rho - rho.^2/rho_max);
df = @(rho) u_max * (1 - 2*rho/rho_max);
ddf = -2 * u_max / rho_max;

%% Select scenario
% scenario = 'Traffic Jam' / 'Green Light' / 'Traffic Flow' + BC and IC;
scenario = setScenario(rho_max,x, Nx);
fprintf('Simulating %s scenario with 2nd-order Godunov method\n', scenario.name);

% Initialize solution array
rho = zeros(Nx, Nt);
rho(:, 1) = scenario.rho_0;  % Set initial condition

% Choose slope limiter
% Options: 'minmod', 'superbee', 'vanLeer', 'MC', 'none'
limiter_type = 'MC';
fprintf('Using %s slope limiter\n', limiter_type);

[rho] = runOrder2(scenario, f, rho_max, Nx, Nt, dx, dt, limiter_type);

%% Visualize results
% Plot the final state
figure;
plot(x, scenario.rho_0, 'b--', 'LineWidth', 1.5);
hold on;
plot(x, rho(:, end), 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Position $x$ [km]');
ylabel('Density $\rho$ [vehicles/km]');
ylim([0 , rho_max*1.05])
title(sprintf('%s: 2nd-Order Godunov - Initial and Final Density', scenario.name));
legend('Initial Condition', 'Final Solution', Location='best');

%%  Create an animation of the solution evolution
figure;

fps = 24
t = floor(linspace(1,size(rho,2), fps));
% %%
% for n = 1:1:Nt
%     plot(x, rho(:, n), 'LineWidth', 1.5);
%     grid on;
%     xlabel('Position x [km]');
%     ylabel('Density \rho [vehicles/km]');
%     title(sprintf('%s: 2nd-Order Godunov - t = %.3f s', scenario, t(n)));
%     ylim([0, 1.1*rho_max]);
%     drawnow;
%     pause(0.05);
% end

%%
for n = t
    plot(x, rho(:, n), 'LineWidth', 1.5);
    grid on;
    xlabel('Position x [km]');
    ylabel('Density \rho [vehicles/km]');
    title(sprintf('%s: 2nd-Order Godunov - $t = %.3f$ s', scenario.name, (n-1)*T/Nt));
    ylim([0, 1.1*rho_max]);
    drawnow;
    pause(0.05);
end
%% Create a space-time plot
figure;
[X, T] = meshgrid(x, t);
surf(X, T, rho');  % Note the transpose of T and rho
colormap jet;
shading interp;
xlabel('Position x [km]');
ylabel('Time t [s]');
zlabel('Density \rho [vehicles/km]');
title(sprintf('%s: 2nd-Order Godunov - Density Evolution', scenario));
view(45, 30);
colorbar;

% Compare with first-order solution (if available)
% If you want to run first-order scheme for comparison, uncomment the following
% and add code to compute first-order solution

%% Runtime Performance Analysis
fprintf('Simulation completed.\n');
fprintf('CFL number used: %.4f\n', CFL);
fprintf('Grid resolution: %d cells\n', Nx);
fprintf('Time steps: %d\n', Nt);
fprintf('Limiter: %s\n', limiter_type);
% Traffic Flow Simulation with Godunov Schemes
% This script implements and compares 1st and 2nd order Godunov methods
% for the traffic flow equation.

clear all;
close all;
clc;
reset(groot)


% set(0,'DefaultFigureWindowStyle','docked');             % all figures as tabs in single window

set(groot, 'DefaultTextInterpreter', 'latex', ...           % interpreter Latex - text and annotations
           'DefaultAxesTickLabelInterpreter', 'latex', ...  % interpreter Latex - tick labels
           'DefaultLegendInterpreter', 'latex', ...         % interpreter Latex - legends
           'DefaultLineLineWidth', 1.5, ...                 % functions
           'DefaultAxesFontSize', 12, ...                   % axis and title
           'DefaultTextFontSize', 14, ...                   % sgtitle
           'DefaultAxesFontName', 'Times New Roman', ...    % axis and title
           'DefaultTextFontName', 'Times New Roman', ...    % sgtitle
           'DefaultAxesLineWidth', 1, ...                   % axis
           'DefaultConstantLineLineWidth', 1.2, ...         % xline and yline
           'DefaultAxesTitleFontSizeMultiplier', 1.2, ...   % title
           'DefaultFigureColor', 'w', ...                   % background color
           'DefaultLegendLocation', 'best' ...              % legend position
           );
           %'DefaultAxesBox', 'on', ...                    % plot box
 



%% Parameters
L = 1;          % Length of the domain [km]
T = 1;          % Total simulation time [s]
Nx = 100;       % Number of spatial cells
Nt = 500;       % Number of time steps

% dx = L / Nx;    % Spatial step size
% dt = T / Nt;    % Time step size
% 
% x = linspace(dx/2, L-dx/2, Nx)';  % Spatial grid (cell centers)
% t = linspace(0, T+dt, Nt+1);           % Time grid

rho_max = 1;                % Maximum density  (normalized)
rho_c = rho_max/2
u_max_km_h = 240 ;           % Maximum velocity in km/h
u_max = u_max_km_h/3600;    % Maximum velocity km/s ==
%% Define flux function and its derivatives
f = @(rho) u_max * (rho - rho.^2/rho_max);
df = @(rho) u_max * (1 - 2*rho/rho_max);
ddf = -2 * u_max / rho_max;

Mesh = createMesh(L, T, Nx, Nt);

% Check CFL condition for both schemes
CFL = u_max * Mesh.dt / Mesh.dx;
if CFL > 0.5
    error('CFL condition for 2nd order method violated! CFL = %.4f > 0.5. Reduce dt or increase dx.', CFL);
else
    fprintf('CFL condition satisfied: %.4f < 0.5\n', CFL);
end

%% Select scenario
% Uncomment the scenario you want to simulate
scenario = 'Traffic jam';
% scenario = 'Green light';
% scenario = 'Traffic flow';

% [scenario] = setScenario(rho_max, rho_c, Mesh.x, Mesh.Nx)

fprintf('Simulating %s scenario\n', scenario);

%% Set initial condition based on scenario
[rho_0, rho_L, rho_R] = init_condition(scenario, rho_max, Mesh.x, Mesh.Nx);

%% Select slope limiter for 2nd order scheme
% Options: 'minmod', 'superbee', 'vanLeer', 'MC', 'none'
limiter_type = 'MC';
fprintf('Using %s slope limiter for 2nd order scheme\n', limiter_type);

%% Solve using both schemes
% First-order solution
rho_1st = first_order_godunov(rho_0, rho_L, rho_R, scenario, Mesh, f);

% Second-order solution
rho_2nd = second_order_godunov(rho_0, rho_L, rho_R, scenario, Mesh, f, limiter_type);

%% Visualize results


% %% 3D Plots
[X, T] = meshgrid(Mesh.x, Mesh.t); % Create space-time grid
    
figure();
sgtitle(sprintf('%s solution for $\\rho(x,t)$, dx = %.3f km, dt = %.3f s', scenario, Mesh.dx, Mesh.dt));
 
subplot(1,2,1);
surf(X, T, rho_1st', 'EdgeColor', 'none'); % Transpose Rho to match dimensions
title('First order scheme');
xlabel('$x$ [km]');
ylabel('$t$ [s]');
zlabel('$\rho$ (vehicles/km)');
ylim([0, Mesh.T*1.01])
colorbar; % Add color scale
colormap(jet); % Use colormap for better visualization
clim([0,rho_max]); % Set common color limits
% Set common color limits
% zlim([-0.1, rho_max+0.1])
   
view(0,90); % Adjust view angle for better visualization
shading interp; % Smooth color transition
grid on;

subplot(1,2,2);
surf(X, T, rho_2nd', 'EdgeColor', 'none'); % Transpose Rho to match dimensions
title('Second order scheme');
colorbar; % Add color scale
colormap(jet); % Use colormap for better visualization
clim([0,rho_max]); % Set common color limits
xlabel('$x$ [km]');
ylabel('$t$ [s]');
zlabel('$\rho$ (vehicles/km)');
ylim([0, Mesh.T*1.01])
   
view(0,90);  % Adjust view angle for better visualization
shading interp; % Smooth color transition
grid on;



%%  Create an animation of the solution evolution
figure();

fps = 24;
t = floor(linspace(1,size(rho_1st,2), fps));

for n = t
    p1 = plot(Mesh.x, rho_1st(:, n), 'LineWidth', 1.5, 'Color','g');
    grid on;
    hold on;
    p2 = plot(Mesh.x, rho_2nd(:, n), 'LineWidth', 1.5, 'Color','b');


    xlabel('Position $x$ [km]');
    ylabel('Density $\rho$ [vehicles/km]');
    title(sprintf(['%s scenario animation ' ...
                 '$dx = %.3f$ km, $dt = %.3f$ s\n' ...
                 '$t = %.3f$ s'], scenario, Mesh.dx, Mesh.dt, (n-1)*Mesh.dt));

    % title(sprintf('1st order scheme - $t = %.3f$ s', (n-1)*Mesh.dt));
    ylim([0, 1.1*rho_max]);

    legend([p1, p2], {'1st order scheme', '2nd order scheme'}, 'Location', 'northeast');


    drawnow;

    pause(0.05);
    hold off;
end


% visualize_results(x, t, rho_0, rho_1st, rho_2nd, scenario, Nx, Nt, CFL, limiter_type);
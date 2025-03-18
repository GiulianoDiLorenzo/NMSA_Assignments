%% Reset and Set Up
clc
clear 
close all
reset(groot)


% set(0,'DefaultFigureWindowStyle','docked');             % all figures as tabs in single window
addpath Functions
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

L = 1000;                   % Length of the road
T = 1;                      % Total simulation time

rho_max = 200;              % Maximum density [car/km]
v_ref = 120;                 % Max speed [km/h]
u_max = v_ref;
% u_max = v_ref*1000/3600;    % Maximum speed   [m/s]


Nx = 80;                   % Number of spatial points
Nt = 10000;                   % Number of time steps


Mesh = createMesh(L, T, Nx, Nt);

CFL = u_max * Mesh.dt / Mesh.dx


if CFL > 1
    error('CFL condition violated!');
else 
    disp('CFL Condition respected')
end

%%
options = {'Traffic jam' , 'Green light', 'Traffic flow' };
        
% Create selection dialog
[selection, ok] = listdlg('PromptString', 'Select an option:', ...
                          'SelectionMode', 'single', ...
                          'ListString', options);

% Check if user made a selection
if ok
    selectedOption = options{selection};
    fprintf('You selected: %s\n', selectedOption);
    
    % Set variables based on selection
    switch selection
        case 1
            scenario = "Traffic jam";
        case 2
            scenario = "Green light";
        case 3
            scenario = "Traffic flow";
    end
    
else
    selectedOption = 'No selection';
    disp('No option selected. Exiting...');
end

%% First order scheme - Traffic Jam
% Setup & Run
scenario = 'Traffic jam';
[cond] = setConditions(scenario, rho_max, Mesh.L, Mesh.x);

[Rho_1_jam,~] = runOrder1Solution(scenario, cond, rho_max, u_max, Mesh);
[Rho_2_jam,~] = runOrder2Solution(scenario, cond, rho_max, u_max, Mesh);

%% First order scheme - Traffic Jam
% % 2D plot
draw2DSolution(scenario, Rho__jam, Mesh, rho_max, u_max);

%% First order scheme - Traffic Jam
% 3D plot
draw3DSolution(scenario, Rho_2_jam, Mesh);

%% First order scheme - Green light
% Setup & Run
scenario = 'Green light';
cond = setConditions(scenario, rho_max, Mesh.L, Mesh.x);

[Rho,t_c] = runOrder1Solution(scenario, cond, rho_max, u_max, Mesh);

%% First order scheme - Green light 
% 2D plot
draw2DSolution(scenario, Rho, Mesh, rho_max, u_max);

%% First order scheme - Green light
% 3D plot
draw3DSolution(scenario, Rho, Mesh);

%% First order scheme - Traffic flow
% Setup & Run
scenario = 'Traffic flow';
cond = setConditions(scenario, rho_max, Mesh.L, Mesh.x);

[Rho,t_c] = runSolution(scenario, cond, rho_max, u_max, Mesh);

%% First order scheme - Traffic flow
% % 2D plot
draw2DSolution(scenario, Rho, Mesh, rho_max, u_max);

%% First order scheme - Traffic flow
% 3D plot
draw3DSolution(scenario, Rho, Mesh);




%% 2nd Order scheme run

%% Param

L = 1;          % Length of the domain [km]
T = 1;          % Total simulation time [s]
Nx = 100;       % Number of spatial cells
Nt = 500;       % Number of time steps

dx = L / Nx;    % Spatial step size
dt = T / Nt;    % Time step size

x = linspace(dx/2, L-dx/2, Nx)';  % Spatial grid (cell centers)
t = linspace(0, T, Nt+1);           % Time grid

Mesh = createMesh(L, T, Nx, Nt);

rho_max = 1;    % Maximum density (normalized)
u_max = 1;      % Maximum velocity (normalized)

% Check CFL condition
CFL = u_max * dt / dx;
if CFL > 0.5
    error('CFL condition for 2nd order method violated! CFL = %.4f > 0.5. Reduce dt or increase dx.', CFL);
else
    fprintf('CFL condition satisfied: %.4f < 0.5\n', CFL);
end

%% Run

% scenario = 'Traffic Jam' / 'Green Light' / 'Traffic Flow' + BC and IC;
scenario = setScenario(rho_max,x, Nx);
fprintf('Simulating %s scenario with 2nd-order Godunov method\n', scenario.name);

% Choose slope limiter
% Options: 'minmod', 'superbee', 'vanLeer', 'MC', 'none'
limiter_type = 'MC';
fprintf('Using %s slope limiter\n', limiter_type);

[rho] = runOrder2(scenario, u_max, rho_max, Nx, Nt, dx, dt, limiter_type);

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

fps = 24;
t = floor(linspace(1,size(rho,2), fps));

for n = t
    plot(x, rho(:, n), 'LineWidth', 1.5);
    grid on;
    xlabel('Position x [km]');
    ylabel('Density \rho [vehicles/km]');
    title(sprintf('%s: 2nd-Order Godunov - $t = %.3f$ s', scenario.name, (n-1)*T/Nt));
    ylim([0, 1.1*rho_max]);
    drawnow;
    pause(1/24);
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
title(sprintf('%s: 2nd-Order Godunov - Density Evolution', scenario.name));
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


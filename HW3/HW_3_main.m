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

L = 1;          % Length of the domain [km]
T = 1;          % Total simulation time [s]
Nx = 100;       % Number of spatial cells
Nt = 500;       % Number of time steps

Mesh = createMesh(L, T, Nx, Nt);

rho_max = 1;    % Maximum density (normalized)
u_max = 1;      % Maximum velocity (normalized)

% Check CFL condition
CFL = u_max * Mesh.dt / Mesh.dx;
if CFL > 0.5
    error('CFL condition for 2nd order method violated! CFL = %.4f > 0.5. Reduce dt or increase dx.', CFL);
else
    fprintf('CFL condition satisfied: %.4f < 0.5\n', CFL);
end


%% Run

% scenario = 'Traffic Jam' / 'Green Light' / 'Traffic Flow' + BC and IC;
scenario = setScenario(rho_max, Mesh.x, Mesh.Nx);

[Rho_1,~] = runOrder1Solution(scenario, rho_max, u_max, Mesh);

% Choose slope limiter
% Options: 'minmod', 'superbee', 'vanLeer', 'MC', 'none'
limiter_type = 'MC';
fprintf('Using %s slope limiter\n', limiter_type);

[Rho_2] = runOrder2Solution(scenario, u_max, rho_max, Mesh, limiter_type);

%% Plot the initial and final states
figure;
sgtitle(sprintf(['%s scenario, Initial and Final density, ' ...
                 '$dx = %.3f$ km, $dt = %.3f$ s'], scenario.name, Mesh.dx, Mesh.dt));

subplot(2,1,1)
plot(Mesh.x, scenario.rho_0, 'b--', 'LineWidth', 1.5);
hold on;
plot(Mesh.x, Rho_1(:, end), 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Position $x$ [km]');
ylabel('Density $\rho$ [vehicles/km]');
ylim([0 , rho_max*1.05])
title('1st-Order Godunov');
legend('Initial Condition', 'Final Solution', Location='northeast');

subplot(2,1,2)
plot(Mesh.x, scenario.rho_0, 'b--', 'LineWidth', 1.5);
hold on;
plot(Mesh.x, Rho_2(:, end), 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Position $x$ [km]');
ylabel('Density $\rho$ [vehicles/km]');
ylim([0 , rho_max*1.05])
title(sprintf('2nd-Order Godunov '));
% legend('Initial Condition', 'Final Solution', Location='bestoutside');


%% 3D Plots
[X, T] = meshgrid(Mesh.x, Mesh.t); % Create space-time grid
    
figure();
sgtitle(sprintf('%s solution, dx = %.3f km, dt = %.3f s', scenario.name, Mesh.dx, Mesh.dt));
 
subplot(1,2,1);
surf(X, T, Rho_1', 'EdgeColor', 'none'); % Transpose Rho to match dimensions
title('First order scheme');
xlabel('$x$ [km]');
ylabel('$t$ [s]');
zlabel('$\rho$ (vehicles/km)');
% zlim([-0.1, rho_max+0.1])
   
view(0,90); % Adjust view angle for better visualization
shading interp; % Smooth color transition
grid on;

subplot(1,2,2);
surf(X, T, Rho_2', 'EdgeColor', 'none'); % Transpose Rho to match dimensions
title('Second order scheme');
colorbar; % Add color scale
colormap(jet); % Use colormap for better visualization
xlabel('$x$ [km]');
ylabel('$t$ [s]');
zlabel('$\rho$ (vehicles/km)');
% zlim([-0.1, rho_max+0.1])
   
view(0,90);  % Adjust view angle for better visualization
shading interp; % Smooth color transition
grid on;



%%  Create an animation of the solution evolution
figure();
sgtitle(sprintf(['%s scenario animation' ...
                 '$dx = %.3f$ km, $dt = %.3f$ s'], scenario.name, Mesh.dx, Mesh.dt));


fps = 24;
t = floor(linspace(1,size(Rho_1,2), fps));

for n = t
    subplot(2,1,1);
    plot(Mesh.x, Rho_1(:, n), 'LineWidth', 1.5);
    grid on;
    xlabel('Position x [km]');
    ylabel('Density \rho [vehicles/km]');
    title(sprintf('1st order scheme - $t = %.3f$ s', (n-1)*Mesh.dt));
    ylim([0, 1.1*rho_max]);


    subplot(2,1,2);
    plot(Mesh.x, Rho_2(:, n), 'LineWidth', 1.5);
    grid on;
    xlabel('Position x [km]');
    ylabel('Density \rho [vehicles/km]');
    title(sprintf('2nd order scheme - $t = %.3f$ s', (n-1)*Mesh.dt));
    ylim([0, 1.1*rho_max]);
    drawnow;

    
    pause(0.2);
end


%%  NEW Create an animation of the solution evolution
figure();

fps = 24;
t = floor(linspace(1,size(Rho_1,2), fps));

for n = t
    p1 = plot(Mesh.x, Rho_1(:, n), 'LineWidth', 1.5, 'Color','g');
    grid on;
    hold on;
    p2 = plot(Mesh.x, Rho_2(:, n), 'LineWidth', 1.5, 'Color','b');


    xlabel('Position $x$ [km]');
    ylabel('Density $\rho$ [vehicles/km]');
    title(sprintf(['%s scenario animation ' ...
                 '$dx = %.3f$ km, $dt = %.3f$ s\n' ...
                 '$t = %.3f$ s'], scenario.name, Mesh.dx, Mesh.dt, (n-1)*Mesh.dt));

    % title(sprintf('1st order scheme - $t = %.3f$ s', (n-1)*Mesh.dt));
    ylim([0, 1.1*rho_max]);

    legend([p1, p2], {'1st order scheme', '2nd order scheme'}, 'Location', 'northeast');


    drawnow;

    pause(0.1);
    hold off;
end



% %% Param
% 
% L = 1;          % Length of the domain [km]
% T = 1;          % Total simulation time [s]
% Nx = 100;       % Number of spatial cells
% Nt = 500;       % Number of time steps
% 
% dx = L / Nx;    % Spatial step size
% dt = T / Nt;    % Time step size
% 
% x = linspace(dx/2, L-dx/2, Nx)';  % Spatial grid (cell centers)
% t = linspace(0, T, Nt+1);           % Time grid
% 
% Mesh = createMesh(L, T, Nx, Nt);
% 
% rho_max = 1;    % Maximum density (normalized)
% u_max = 1;      % Maximum velocity (normalized)
% 
% % Check CFL condition
% CFL = u_max * dt / dx;
% if CFL > 0.5
%     error('CFL condition for 2nd order method violated! CFL = %.4f > 0.5. Reduce dt or increase dx.', CFL);
% else
%     fprintf('CFL condition satisfied: %.4f < 0.5\n', CFL);
% end
% 
% %% Run
% 
% % scenario = 'Traffic Jam' / 'Green Light' / 'Traffic Flow' + BC and IC;
% scenario = setScenario(rho_max,x, Nx);
% fprintf('Simulating %s scenario with 2nd-order Godunov method\n', scenario.name);
% 
% % Choose slope limiter
% % Options: 'minmod', 'superbee', 'vanLeer', 'MC', 'none'
% limiter_type = 'MC';
% fprintf('Using %s slope limiter\n', limiter_type);
% 
% [rho] = runOrder2(scenario, u_max, rho_max, Nx, Nt, dx, dt, limiter_type);
% 
% %% Visualize results
% % Plot the final state
% figure;
% plot(x, scenario.rho_0, 'b--', 'LineWidth', 1.5);
% hold on;
% plot(x, rho(:, end), 'r-', 'LineWidth', 1.5);
% grid on;
% xlabel('Position $x$ [km]');
% ylabel('Density $\rho$ [vehicles/km]');
% ylim([0 , rho_max*1.05])
% title(sprintf('%s: 2nd-Order Godunov - Initial and Final Density', scenario.name));
% legend('Initial Condition', 'Final Solution', Location='best');
% 
% %%  Create an animation of the solution evolution
% figure;
% 
% fps = 24;
% t = floor(linspace(1,size(rho,2), fps));
% 
% for n = t
%     plot(x, rho(:, n), 'LineWidth', 1.5);
%     grid on;
%     xlabel('Position x [km]');
%     ylabel('Density \rho [vehicles/km]');
%     title(sprintf('%s: 2nd-Order Godunov - $t = %.3f$ s', scenario.name, (n-1)*T/Nt));
%     ylim([0, 1.1*rho_max]);
%     drawnow;
%     pause(1/24);
% end
% %% Create a space-time plot
% figure;
% [X, T] = meshgrid(x, t);
% surf(X, T, rho');  % Note the transpose of T and rho
% colormap jet;
% shading interp;
% xlabel('Position x [km]');
% ylabel('Time t [s]');
% zlabel('Density \rho [vehicles/km]');
% title(sprintf('%s: 2nd-Order Godunov - Density Evolution', scenario.name));
% view(45, 30);
% colorbar;
% 
% % Compare with first-order solution (if available)
% % If you want to run first-order scheme for comparison, uncomment the following
% % and add code to compute first-order solution
% 
% %% Runtime Performance Analysis
% fprintf('Simulation completed.\n');
% fprintf('CFL number used: %.4f\n', CFL);
% fprintf('Grid resolution: %d cells\n', Nx);
% fprintf('Time steps: %d\n', Nt);
% fprintf('Limiter: %s\n', limiter_type);
% 

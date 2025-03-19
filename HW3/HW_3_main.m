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
 
%% Parameters and set up

L = 1;          % Length of the domain [km]
T = 1;          % Total simulation time [s]
Nx = 100;       % Number of spatial cells
Nt = 400;       % Number of time steps

rho_max = 1;    % Maximum density (normalized)
u_max = 1;      % Maximum velocity (normalized)

f   =  @(rho) u_max * (rho - rho.^2/rho_max); 
df  =  @(rho) u_max * (1 - 2*rho/rho_max);
ddf =  - 2 * u_max/rho_max;

rho_c = rho_max/2;  % critical rho, for which df=0


Mesh = createMesh(L, T, Nx, Nt);
scenario = setScenario(rho_max, rho_c, Mesh.x, Mesh.Nx); % scenario = 'Traffic Jam' / 'Green Light' / 'Traffic Flow' + BC + IC;


%% Numerical Run

% Choose slope limiter
% Options: 'minmod', 'superbee', 'vanLeer', 'MC', 'none'
limiter_type = 'MC';

Rho_1 = runOrder1Solution(scenario, Mesh, f, rho_c);
[Rho_2] = runOrder2Solution(scenario, f, rho_max, rho_c, u_max, Mesh, limiter_type);
fprintf('Using %s slope limiter\n', limiter_type);

% %%
% draw2DSolution(scenario.name, Rho_2, Mesh, rho_max, u_max);
% draw2DSolution(scenario.name, Rho_1, Mesh, rho_max, u_max);

    

%% 3D Plots
[X, T] = meshgrid(Mesh.x, Mesh.t); % Create space-time grid
    
figure();
sgtitle(sprintf('%s solution for $\\rho(x,t)$, dx = %.3f km, dt = %.3f s', scenario.name, Mesh.dx, Mesh.dt));
 
subplot(1,2,1);
surf(X, T, Rho_1', 'EdgeColor', 'none'); % Transpose Rho to match dimensions
title('First order scheme');
xlabel('$x$ [km]');
ylabel('$t$ [s]');
zlabel('$\rho$ (vehicles/km)');
colorbar; % Add color scale
colormap(jet); % Use colormap for better visualization
clim([min(scenario.rho_L, scenario.rho_R), rho_max]); % Set common color limits
% zlim([-0.1, rho_max+0.1])
   
% view(0,90); % Adjust view angle for better visualization
shading interp; % Smooth color transition
grid on;

subplot(1,2,2);
surf(X, T, Rho_2', 'EdgeColor', 'none'); % Transpose Rho to match dimensions
title('Second order scheme');
colorbar; % Add color scale
colormap(jet); % Use colormap for better visualization
clim([min(scenario.rho_L, scenario.rho_R), rho_max]); % Set common color limits
xlabel('$x$ [km]');
ylabel('$t$ [s]');
zlabel('$\rho$ (vehicles/km)');
% zlim([-0.1, rho_max+0.1])
   
% view(0,90);  % Adjust view angle for better visualization
shading interp; % Smooth color transition
grid on;


%%  Create an animation of the solution evolution
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
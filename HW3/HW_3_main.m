%% Reset and Set Up
clc
clear 
close all
reset(groot)


set(0,'DefaultFigureWindowStyle','docked');             % all figures as tabs in single window

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
           'DefaultAxesBox', 'on', ...                      % plot box
           'DefaultLegendLocation', 'best' ...              % legend position
           );

%% Parameters

L = 1;                  % Length of the road
T = 1;                  % Total simulation time

Nx = 10;               % Number of spatial points
Nt = 20;               % Number of time steps

dx = L / Nx;            % Spatial step
dt = T / Nt;            % Time step (CFL condition should be checked)
x = linspace(0, L, Nx); % Spatial grid

rho_max = 1;            % Maximum density (normalized)
u_max = 1;              % Maximum speed   (normalized)


%% Boundary and Initial Conditions 
% =======================================
% ============= TRAFFIC JAM =============
% =======================================
[rho_L, rho_R, rho_x_0] = setConditions('jam', rho_max, L, x)




%% Boundary and Initial Conditions 
% =======================================
% ============= GREEN LIGHT =============
% =======================================
[rho_L, rho_R, rho_x_0] = setConditions('green', rho_max, L, x)



%% Boundary and Initial Conditions 
% =======================================
% ============= TRAFFIC FLOW ============
% =======================================
[rho_L, rho_R, rho_x_0] = setConditions('flow', rho_max, L, x)



%% First order scheme - Traffic Jam
[rho_L, rho_R, rho_x_0] = setConditions('jam', rho_max, L, x)

rho = rho_x_0;      % rho(x,t) at t = 0
f = zeros(Nx,1);    % Initialize f(x,t)

for n = 2:Nt        % rho contains the value for t=0 and is already computed
    for i = 2:Nx-1  % rho contains the value for x = {0 , L} and they stay constant

        % Update rule for the cell on the right
        if rho(i-1,n-1) <= rho(i+1,n-1)
            fR = min( f(i-1,n-1) , f(i,n) );
        else
            fR = max( f(i-1,n-1) , f(i,n) );
        end


        % Update rule for the cell on the left
        if rho(i-1,n-1) <= rho(i+1,n-1)
            fR = min( f(i-1,n-1) , f(i,n) );
        else
            fR = max( f(i-1,n-1) , f(i,n) );
        end
        
        fL ;'

        rho(i,n) = 

    end
end
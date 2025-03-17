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

L = 1000;                      % Length of the road
T = 1;                      % Total simulation time

Nx = 200;                   % Number of spatial points
Nt = 200;                    % Number of time steps

Mesh = createMesh(L, T, Nx, Nt);

rho_max = 200;                % Maximum density (normalized)
u_max = 150;                  % Maximum speed   (normalized)

CFL = u_max * Mesh.dt / Mesh.dx


if CFL > 1
    error('CFL condition violated!');
else 
    disp('CFL Condition respected')
end



%% First order scheme - Traffic Jam
% Setup & Run
scenario = 'Traffic jam';
[cond] = setConditions(scenario, rho_max, Mesh.L, Mesh.x);

[Rho_1_jam,~] = runOrder1Solution(scenario, cond, rho_max, u_max, Mesh);
[Rho_2_jam,~] = runOrder2Solution(scenario, cond, rho_max, u_max, Mesh);

%% First order scheme - Traffic Jam
% % 2D plot
draw2DSolution(scenario, Rho_1_jam, Mesh, rho_max);

%% First order scheme - Traffic Jam
% 3D plot
draw3DSolution(scenario, Rho_1_jam, Mesh);

%% First order scheme - Green light
% Setup & Run
scenario = 'Green light';
cond = setConditions(scenario, rho_max, Mesh.L, Mesh.x);

[Rho,t_c] = runOrder1Solution(scenario, cond, rho_max, u_max, Mesh);

%% First order scheme - Green light 
% 2D plot
draw2DSolution(scenario, Rho, Mesh, rho_max);

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
draw2DSolution(scenario, Rho, Mesh, rho_max);

%% First order scheme - Traffic flow
% 3D plot
draw3DSolution(scenario, Rho, Mesh);



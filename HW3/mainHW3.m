% Traffic Flow Simulation with Godunov Schemes
% This script implements and compares 1st and 2nd order Godunov methods
% for the traffic flow equation.

clear all;
close all;
clc;
reset(groot)

addpath Functions/'Full Study'
addpath Functions/'Show Results'
addpath Functions/'Update Values'

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

%% Parameters
L = 1;          % Length of the domain [km]
T = 5;          % Total simulation time [s]
Nx = 100;       % Number of spatial cells
Nt = 1000*T;       % Number of time steps

rho_max = 1;                % Maximum density  (normalized)
rho_c = rho_max/2;
u_max_km_h = 200 ;           % Maximum velocity in km/h
u_max = u_max_km_h/3600;    % Maximum velocity km/s ==
%% Define flux function and its derivatives
% [f1, f2, rho_c_f1, rho_c_f2] = setFluxes(rho_max, u_max, false);
plotFlux = false;
fluxes = setFluxes(rho_max, u_max, plotFlux);
% saveas(gcf, 'Pictures/Fluxes_graph.png');

% Select scenario - Uncomment the scenario you want to simulate
sce = 'Traffic jam';
% sce = 'Green light';
% sce = 'Traffic flow';



[results] = runFullStudy(L,T,Nx,Nt,sce,fluxes);

%% 3D Plots - Visualize results

drawDensity(results);
saveas(gcf, sprintf('Pictures/%s density comp %ds.png', results.scenario.name, results.Mesh.T));


drawFlux(results, fluxes);
saveas(gcf, sprintf('Pictures/%s flux comp %ds.png', results.scenario.name, results.Mesh.T));

%% 2D Animation

animateDensity(results, rho_max);
animateFlux(results, rho_max, u_max, fluxes);



%% Create and save animation as GIF
makeDensityGIF(results, rho_max);
makeFluxGIF(results, fluxes);
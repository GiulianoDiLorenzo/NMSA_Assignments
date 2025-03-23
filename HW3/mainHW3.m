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
T = 1;          % Total simulation time [s]
Nx = 100;       % Number of spatial cells
Nt = 500*T;       % Number of time steps

rho_max = 1;                % Maximum density  (normalized)
rho_c = rho_max/2;
u_max_km_h = 400 ;           % Maximum velocity in km/h
u_max = u_max_km_h/3600;    % Maximum velocity km/s 


%% Define flux function and its derivatives
plotFlux = false;
fluxes = setFluxes(rho_max, u_max, plotFlux);
% saveas(gcf, 'Pictures/Fluxes_graph.png');

% Select scenario - Uncomment the scenario you want to simulate
% sce = 'Traffic jam';
sce = 'Green light';
% sce = 'Traffic flow';

[results] = runFullStudy(L,T,Nx,Nt,sce,fluxes);

%% 2D plots - Initial & Final states
saveMe = false;

drawStartAndStop(results,fluxes, 'density', saveMe);
drawStartAndStop(results,fluxes, 'flux', saveMe);


%% 3D Plots - Visualize results
saveMe = false;

drawResults(results, fluxes, 'density', saveMe);
drawResults(results, fluxes, 'flux', saveMe);

%% 2D Animation
animateResults(results, fluxes, 'density');
animateResults(results, fluxes, 'flux');

%% Create and save animation as GIF
makeDensityGIF(results, rho_max);
makeFluxGIF(results, fluxes);
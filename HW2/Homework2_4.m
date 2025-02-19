% Solution of HOMEWORK 2
%     S utt = gamma^2 (S ux)x + f
%     S(x) variable cross-section
%     f(x) force
%     u(x,t) solution

clc
clear 
close all
reset(groot)

addpath('functions')

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

%% Point 4, parameters
% velocity, length and gamma
c = 2;
L = 1;
gamma = c/L;

% space domain (0,1) and time domain (0,T)
T = 5;

% number of space steps
NX = 1000;
% space step / spatial mesh size
dx = L/NX;

% number of time steps (minimum 2*T*NX for stability)
NT = 2*T*NX;
% time step / temporal mesh size
dt = T/NT;

% lambda
lambda = gamma*dt/dx;

% space domain definition: 0 to L in NX elements
x = linspace(0, L, NX).';

% time domain definition: 0 to T in NT elements
t = linspace(0, T, NT);

% % initial conditions
% u_0 = zeros(NX, 1);
% u_1 = zeros(NX, 1);

ux_0 = -1/2 * (sin((pi*(t-0.1))/0.05) + abs(sin((pi*(t-0.1))/0.05)));
% extending ux_0
ux_0 = [-1/2 * (sin((pi*(-dt-0.1))/0.05) + abs(sin((pi*(-dt-0.1))/0.05))), ux_0];


% % variable section definitions
% S = ...;
% Sx = ...;

%% Point 4, constant profile
% this block computes the numerical solution for S(x)=1
% prints the solution surface

tic

% initializing the solution;
sol = zeros(NX+1, NT+1);

% applying the boundary conditions
sol(1,:) = sol(2,:) - dx * ux_0;        % 1st condition
sol(:, 2) = 0;                      % 3rd condition
sol(:, 1) = sol(:, 2);          % 4th condition
sol(end,:) = 0;                         % 2nd condition

for n = 2:NT      % evaluating n+1, up to NT to consider virtual line

    for k = 2:NX
        sol(k,n+1) = 2*sol(k,n) - sol(k,n-1) + lambda^2 * (sol(k+1,n) - 2*sol(k,n) + sol(k-1,n));
    end

    sol(1,n+1) = sol(2,n+1) - dx * ux_0(n+1);                % 1st condition
end

% discarding virtual lines (1,:) and (:,1)
sol = sol(2:end, 2:end);

toc

% printing solution as surface
figure()
[X_axis, T_axis] = ndgrid(x, t);
surf(X_axis, T_axis, sol);
shading interp
title("Numerical solution with constant cross-section S(x) = 1 m ($N_X$ = " + NX + ", $N_T$ = "+ NT + ")");
xlabel('$x$ [m]');
ylabel('$t$ [s]');
zlabel('$u_h(x,t)$');
colorbar

% %% EXTRA - Plotting solution comparison with constant cross-section
% figure()
% for n = 1:NT
%     plot(x, sol(:,n));
%     title("Solution evolution (point 4) with constant cross-section S(x) = 1 m ($N_X$ = " + NX + ", $N_T$ = "+ NT + ") at t = "+ (n-1)*dt+ " s");
%     xlabel('$x$ [m]');
%     ylabel('$u_h(x,t)$');
%     ylim([0,0.7]);
%     grid on;
%     pause(1e-5);
%     hold off;
% end

%% Point 4, variable profile S_1
% this block computes the numerical solution for S(x) given in the file Sections.m
% prints the solution surface

% Sections.m
S1 = [0 1;0.09 0.4;0.11 2.4;0.24 2.4;0.26 3.2;0.29 3.2;0.32 4.2;...
0.41 4.2;0.47 3.2;0.59 1.8;0.65 1.6;0.71 1.6;0.74 1;0.76 0.8;...
0.82 0.8;0.88 2;0.91 2;0.94 3.2;1 3.2];

% interpolating sections with new x-axis
S_1 = interp1(S1(:,1), S1(:,2), x.', "linear").';

% elongating S_1 to take into account virtual line
S_1 = [S_1(1) + (S_1(1) - S_1(2)); S_1];

tic

% initializing the solution;
sol_var_1 = zeros(NX+1, NT+1);

% applying the boundary conditions
sol_var_1(1,:) = sol_var_1(2,:) - dx * ux_0;    % 1st condition
sol_var_1(:, 2) = 0;                          % 3rd condition
sol_var_1(:, 1) = sol_var_1(:, 2);              % 4th condition
sol_var_1(end,:) = 0;                         % 2nd condition

for n = 2:NT      % evaluating n+1, up to NT to consider virtual line

    for k = 2:NX
        sol_var_1(k,n+1) = 2*sol_var_1(k,n) - sol_var_1(k,n-1) + lambda^2/S_1(k) * (((S_1(k+1)-S_1(k-1))/4 + S_1(k)) * sol_var_1(k+1,n) - 2*S_1(k) * sol_var_1(k,n) + ((S_1(k-1)-S_1(k+1))/4 + S_1(k)) * sol_var_1(k-1,n));
    end

    sol_var_1(1,n+1) = sol_var_1(2,n+1) - dx * ux_0(n+1);                % 1st condition
end

% discarding virtual lines (1,:) and (:,1)
sol_var_1 = sol_var_1(2:end, 2:end);

toc

% printing solution as surface
figure()
[X_axis, T_axis] = ndgrid(x, t);
surf(X_axis, T_axis, sol_var_1);
shading interp
title("Numerical solution with variable cross-section $S_1$ ($N_X$ = " + NX + ", $N_T$ = "+ NT + ")");
xlabel('$x$ [m]');
ylabel('$t$ [s]');
zlabel('$u_h(x,t)$');
colorbar

%% Point 4, variable profile S_2
% this block computes the numerical solution for S(x) given in the file Sections.m
% prints the solution surface

% Sections.m
S2 = [0 1;0.03 0.60;0.09 0.4;0.12 1.6;0.18 0.6;0.29 0.2;0.35 0.4;...
    0.41 0.8;0.47 1;0.50 0.6;0.59 2;0.65 3.2;0.85 3.2;0.94 2;1 2];

% interpolating sections with new x-axis
S_2 = interp1(S2(:,1), S2(:,2), x.', "linear").';

% elongating S_2 to take into account virtual line
S_2 = [S_2(1) + (S_2(1) - S_2(2)); S_2];

tic

% initializing the solution;
sol_var_2 = zeros(NX+1, NT+1);

% applying the boundary conditions
sol_var_2(1,:) = sol_var_2(2,:) - dx * ux_0;    % 1st condition
sol_var_2(:, 2) = 0;                          % 3rd condition
sol_var_2(:, 1) = sol_var_2(:, 2);              % 4th condition
sol_var_2(end,:) = 0;                         % 2nd condition

for n = 2:NT      % evaluating n+1, up to NT to consider virtual line

    for k = 2:NX
        sol_var_2(k,n+1) = 2*sol_var_2(k,n) - sol_var_2(k,n-1) + lambda^2/S_2(k) * (((S_2(k+1)-S_2(k-1))/4 + S_2(k)) * sol_var_2(k+1,n) - 2*S_2(k) * sol_var_2(k,n) + ((S_2(k-1)-S_2(k+1))/4 + S_2(k)) * sol_var_2(k-1,n));
    end

    sol_var_2(1,n+1) = sol_var_2(2,n+1) - dx * ux_0(n+1);                % 1st condition
end

% discarding virtual lines (1,:) and (:,1)
sol_var_2 = sol_var_2(2:end, 2:end);

toc

% printing solution as surface
figure()
[X_axis, T_axis] = ndgrid(x, t);
surf(X_axis, T_axis, sol_var_2);
shading interp
title("Numerical solution with variable cross-section $S_2$ ($N_X$ = " + NX + ", $N_T$ = "+ NT + ")");
xlabel('$x$ [m]');
ylabel('$t$ [s]');
zlabel('$u_h(x,t)$');
colorbar

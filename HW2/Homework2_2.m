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

%% Point 2, parameters
% velocity, length and gamma
c = 1;
L = 1;
gamma = c/L;

% space domain (0,1) and time domain (0,T)
T = 5;

% number of space steps
NX = 800;
% space step / spatial mesh size
dx = L/NX;

% number of time steps (minimum T*NX for stability)
NT = T*NX;
% time step / temporal mesh size
dt = T/NT;

% lambda
lambda = gamma*dt/dx;

% space domain definition: 0 to L in NX elements
x = linspace(0, L, NX).';

% time domain definition: 0 to T in NT elements
t = linspace(0, T, NT);

% section definitions
S_const = 1;
S = (1+2*x).^2;
Sx = 4 + 8*x;   % derivative in x of S

%% Point 2, exact solution
% exact solution for computing the error
uex     = - cos(pi/2*x) * cos(3*pi*t);

% 1st derivatives
ux      = pi/2 * sin(pi/2*x) * cos(3*pi*t);
ut      = 3*pi * cos(pi/2*x) * sin(3*pi*t);

% 2nd derivatives
uxx     = (pi/2)^2 * cos(pi/2*x) * cos(3*pi*t);
utt     = (3*pi)^2 * cos(pi/2*x) * cos(3*pi*t);

% initial conditions
u_0 = uex(:,1);             % u_0 = u(x,0)
u_1 = ut(:,1);              % u_1 = du/dt(x,0)

% excitation forces
f  = S_const * (utt - gamma^2 * uxx);
f_var = S .* utt - gamma^2 .* (Sx .* ux + S .* uxx);

%% Point 2, constant profile
% this block computes the numerical solution for S(x)=1
% prints the solution surface
% computes the L2-error
% prints the approximation error surface

sol = computeSolutionConstant(NX, NT, dt, u_0, u_1, lambda, f);

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

% approximation difference (in space and time)
err = sol - uex;
e_L2 = sqrt(sum(sum(err.^2)) * dx * dt / (L * T));

% plotting error (in abs) as surface
figure()
[X_axis, T_axis] = ndgrid(x, t);
surf(X_axis, T_axis, abs(err));
shading interp
title("Numerical error (in abs) with constant cross-section S(x) = 1 m ($N_X$ = " + NX + ", $N_T$ = "+ NT + ")");
xlabel('$x$ [m]');
ylabel('$t$ [s]');
zlabel('$|u_h(x,t) - u_{ex}(x,t)|$');
colorbar

% %% EXTRA - Plotting solution comparison with constant cross-section
% figure()
% for n = 1:NT
%     plot(x, uex(:,n), "LineStyle", ":");
%     hold on;
%     plot(x, sol(:,n));
%     legend('$u_{ex}$', '$u$');
%     title("Solution comparison with constant cross-section S(x) = " + S_const + " m ($N_X$ = " + NX + ", $N_T$ = "+ NT + ") at t = "+ (n-1)*dt+ " s");
%     xlabel('$x$ [m]');
%     ylabel('$u(x,t), \quad u_{ex}(x,t)$');
%     %ylim([-1,1]);
%     grid on;
%     pause(1e-5);
%     hold off;
% end
% 
% %% EXTRA - Plotting error with constant cross-section
% figure()
% for n = 1:NT
%     plot(x, err(:,n));
%     title("Error evaluation with constant cross-section S(x) = " + S_const + " m ($N_X$ = " + NX + ", $N_T$ = "+ NT + ") at t = "+ (n-1)*dt+ " s");
%     xlabel('$x$ [m]');
%     ylabel('$u(x,t) - u_{ex}(x,t)$');
%     %ylim([-1,1]);
%     grid on;
%     pause(1e-5);
%     hold off;
% end

%% Point 2, instability example (if you run this, take care of running the previous blocks again)
NX = 800;
NT = 2000;

dx = L/NX;
dt = T/NT;
lambda = gamma*dt/dx;

x = linspace(0, L, NX).';
t = linspace(0, T, NT);

uex     = - cos(pi/2*x) * cos(3*pi*t);
ux      = pi/2 * sin(pi/2*x) * cos(3*pi*t);
ut      = 3*pi * cos(pi/2*x) * sin(3*pi*t);
uxx     = (pi/2)^2 * cos(pi/2*x) * cos(3*pi*t);
utt     = (3*pi)^2 * cos(pi/2*x) * cos(3*pi*t);

u_0 = uex(:,1);
u_1 = ut(:,1);
f  = S_const * (utt - gamma^2 * uxx);

sol = computeSolutionConstant(NX, NT, dt, u_0, u_1, lambda, f);

% printing instable solution
figure()
m = 6;
plot(x, sol(:, m));
grid on
title("Numerical solution with constant cross-section S(x) = 1 m ($N_X$ = " + NX + ", $N_T$ = "+ NT + ") at time t = " + m*dt + " s");
xlabel("$x$ [m]");
ylabel("$u_h(x,t)$");

%% Point 3, variable profile
% this block computes the numerical solution for S(x)=(1+2x)^2
% prints the solution surface
% computes the L2-error
% prints the approximation error surface

sol_var = computeSolutionVariable(NX, NT, dx, dt, S, u_0, u_1, lambda, f_var);

% printing solution as surface
figure()
[X_axis, T_axis] = ndgrid(x, t);
surf(X_axis, T_axis, sol);
shading interp
title("Numerical solution with variable cross-section S(x) = $(1+2x)^2$ ($N_X$ = " + NX + ", $N_T$ = "+ NT + ")");
xlabel('$x$ [m]');
ylabel('$t$ [s]');
zlabel('$u_h(x,t)$');
colorbar

% approximation difference (in space and time)
err_var = sol_var - uex;
e_L2_var = sqrt(sum(sum(err.^2)) * dx * dt / (L * T));

% plotting error (in abs) as surface
figure()
[X_axis, T_axis] = ndgrid(x, t);
surf(X_axis, T_axis, abs(err_var));
shading interp
title("Numerical error (in abs) with variable cross-section S(x) = $(1+2x)^2$ ($N_X$ = " + NX + ", $N_T$ = "+ NT + ")");
xlabel('$x$ [m]');
ylabel('$t$ [s]');
zlabel('$|u_h(x,t) - u_{ex}(x,t)|$');
colorbar

% %% EXTRA - Plotting solution comparison with variable cross-section
% figure()
% for n = 1:NT
%     plot(x, uex(:,n), "LineStyle", ":");
%     hold on;
%     plot(x, sol_var(:,n));
%     legend('$u_{ex}$', '$u$');
%     title("Solution comparison with variable cross-section S(x) = $(1+2x)^2$ ($N_X$ = " + NX + ", $N_T$ = "+ NT + ") at t = "+ (n-1)*dt+ " s");
%     xlabel('$x$ [m]');
%     ylabel('$u(x,t), \quad u_{ex}(x,t)$');
%     %ylim([-1,1]);
%     grid on;
%     pause(1e-5);
%     hold off;
% end
% 
% %% EXTRA - Plotting error with variable cross-section
% figure()
% for n = 1:NT
%     plot(x, err_var(:,n));
%     title("Error evaluation with variable cross-section S(x) = $(1+2x)^2$ ($N_X$ = " + NX + ", $N_T$ = "+ NT + ") at t = "+ (n-1)*dt+ " s");
%     xlabel('$x$ [m]');
%     ylabel('$u(x,t) - u_{ex}(x,t)$');
%     %ylim([-1,1]);
%     grid on;
%     pause(1e-5);
%     hold off;
% end

%% Looping with different NX and NT (using NT = k*NX)
% this block loops the computation for both constant and variable profile S(x)
% takes different NX and computes NT = k*NX for each case
% where k=mesh_const is a refinement factor for the time domain

mesh_const = 5;     % refinement factor (tried 5, 10 and 15)
NX = 100*linspace(1,20,20);
NT = mesh_const*NX;

e_L2 = zeros(1, length(NX));
e_L2_var = zeros(1, length(NX));

for i = 1:length(NX)

    dx = L/NX(i);
    dt = T/NT(i);
    lambda = gamma*dt/dx;
    
    x = linspace(0, L, NX(i)).';
    t = linspace(0, T, NT(i));
    
    S = (1+2*x).^2;
    Sx = 4 + 8*x;
    
    uex     = - cos(pi/2*x) * cos(3*pi*t);
    
    ux      = pi/2 * sin(pi/2*x) * cos(3*pi*t);
    ut      = 3*pi * cos(pi/2*x) * sin(3*pi*t);
    uxx     = (pi/2)^2 * cos(pi/2*x) * cos(3*pi*t);
    utt     = (3*pi)^2 * cos(pi/2*x) * cos(3*pi*t);
    
    u_0 = uex(:,1);
    u_1 = ut(:,1);
    
    f  = S_const * (utt - gamma^2 * uxx);
    f_var = S .* utt - gamma^2 .* (Sx .* ux + S .* uxx);
    
    sol = computeSolutionConstant(NX(i), NT(i), dt, u_0, u_1, lambda, f);
    e_L2(i) = sqrt(sum(sum((sol - uex).^2)) * dx * dt / (L * T));
    
    sol_var = computeSolutionVariable(NX(i), NT(i), dx, dt, S, u_0, u_1, lambda, f_var);
    e_L2_var(i) = sqrt(sum(sum((sol_var - uex).^2)) * dx * dt / (L * T));

end

% printing norm errors in function of dx/L
figure()
% plot(flip(L./NX), flip(e_L2), '-o');
plot(flip(1./NX)*100, flip(e_L2), '-o');
title("Norm error as function of relative space step, $N_T$ = " + mesh_const + " $N_X$");
% xlabel("$\Delta x$ [m]");
xlabel("$\Delta x/L$ [\%]");
ylabel("$||u_h-u_{ex}||_{L^2}$");
grid on
hold on
% plot(flip(L./NX), flip(e_L2_var), '-o');
plot(flip(1./NX)*100, flip(e_L2_var), '-o');
legend("Constant profile", "Variable profile");

% printing norm errors in function of dt/T
figure()
% plot(flip(T./NT), flip(e_L2), '-o');
plot(flip(1./NT)*100, flip(e_L2), '-o');
title("Norm error as function of relative time step, $N_T$ = " + mesh_const + " $N_X$");
% xlabel("$\Delta t$ [m]");
xlabel("$\Delta t/T$ [\%]");
ylabel("$||u_h-u_{ex}||_{L^2}$");
grid on
hold on
% plot(flip(T./NT), flip(e_L2_var), '-o');
plot(flip(1./NT)*100, flip(e_L2_var), '-o');
legend("Constant profile", "Variable profile");

%% Looping with different NT, keeping NX
% this block loops the computation for both constant and variable profile S(x)
% takes a single NX = space_steps and computes NT = k*NX for each case
% where k increases each time and is a refinement factor for the time domain

space_steps = 100;       % (100)=1% of space domain
NX = space_steps * ones(1, 20);
NT = NX .* round(linspace(5, 25, 20));

e_L2 = zeros(1, length(NX));
e_L2_var = zeros(1, length(NX));

for i = 1:length(NX)

    dx = L/NX(i);
    dt = T/NT(i);
    lambda = gamma*dt/dx;
    
    x = linspace(0, L, NX(i)).';
    t = linspace(0, T, NT(i));
    
    S = (1+2*x).^2;
    Sx = 4 + 8*x;
    
    uex     = - cos(pi/2*x) * cos(3*pi*t);
    
    ux      = pi/2 * sin(pi/2*x) * cos(3*pi*t);
    ut      = 3*pi * cos(pi/2*x) * sin(3*pi*t);
    uxx     = (pi/2)^2 * cos(pi/2*x) * cos(3*pi*t);
    utt     = (3*pi)^2 * cos(pi/2*x) * cos(3*pi*t);
    
    u_0 = uex(:,1);
    u_1 = ut(:,1);
    
    f  = S_const * (utt - gamma^2 * uxx);
    f_var = S .* utt - gamma^2 .* (Sx .* ux + S .* uxx);
    
    sol = computeSolutionConstant(NX(i), NT(i), dt, u_0, u_1, lambda, f);
    e_L2(i) = sqrt(sum(sum((sol - uex).^2)) * dx * dt / (L * T));
    
    sol_var = computeSolutionVariable(NX(i), NT(i), dx, dt, S, u_0, u_1, lambda, f_var);
    e_L2_var(i) = sqrt(sum(sum((sol_var - uex).^2)) * dx * dt / (L * T));

end

% printing norm errors in function of dt/T
figure()
% plot(flip(T./NT), flip(e_L2), '-o');
plot(flip(1./NT)*100, flip(e_L2), '-o');
title("Norm error as function of relative time step, $N_X$ = " + space_steps + " (" + L*100/space_steps + " \%)");
% xlabel("$\Delta t$ [m]");
xlabel("$\Delta t/T$ [\%]");
ylabel("$||u_h-u_{ex}||_{L^2}$");
grid on
hold on
% plot(flip(T./NT), flip(e_L2_var), '-o');
plot(flip(1./NT)*100, flip(e_L2_var), '-o');
legend("Constant profile", "Variable profile");

%% Looping with different NX, keeping NT
% this block loops the computation for both constant and variable profile S(x)
% takes a single NT = time_steps and computes NX = NT/k for each case
% where k decreases each time and is a refinement factor for the space domain

time_steps = 500;       % (500)=1% of space domain
NT = time_steps * ones(1, 20);
NX = round(NT .* (linspace(1, 20, 20)/100));

e_L2 = zeros(1, length(NX));
e_L2_var = zeros(1, length(NX));

for i = 1:length(NX)

    dx = L/NX(i);
    dt = T/NT(i);
    lambda = gamma*dt/dx;
    
    x = linspace(0, L, NX(i)).';
    t = linspace(0, T, NT(i));
    
    S = (1+2*x).^2;
    Sx = 4 + 8*x;
    
    uex     = - cos(pi/2*x) * cos(3*pi*t);
    
    ux      = pi/2 * sin(pi/2*x) * cos(3*pi*t);
    ut      = 3*pi * cos(pi/2*x) * sin(3*pi*t);
    uxx     = (pi/2)^2 * cos(pi/2*x) * cos(3*pi*t);
    utt     = (3*pi)^2 * cos(pi/2*x) * cos(3*pi*t);
    
    u_0 = uex(:,1);
    u_1 = ut(:,1);
    
    f  = S_const * (utt - gamma^2 * uxx);
    f_var = S .* utt - gamma^2 .* (Sx .* ux + S .* uxx);
    
    sol = computeSolutionConstant(NX(i), NT(i), dt, u_0, u_1, lambda, f);
    e_L2(i) = sqrt(sum(sum((sol - uex).^2)) * dx * dt / (L * T));
    
    sol_var = computeSolutionVariable(NX(i), NT(i), dx, dt, S, u_0, u_1, lambda, f_var);
    e_L2_var(i) = sqrt(sum(sum((sol_var - uex).^2)) * dx * dt / (L * T));

end

% printing norm errors in function of dx/L
figure()
plot(flip(1./NX)*100, flip(e_L2), '-o');
title("Norm error as function of relative time step, $N_T$ = " + time_steps + " (" + T*100/time_steps + " \%)");
xlabel("$\Delta x/L$ [\%]");
ylabel("$||u_h-u_{ex}||_{L^2}$");
grid on
hold on
plot(flip(1./NX)*100, flip(e_L2_var), '-o');
legend("Constant profile", "Variable profile");

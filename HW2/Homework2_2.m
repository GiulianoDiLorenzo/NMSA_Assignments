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
NX = 1000;
% space step / spatial mesh size
dx = L/(NX-1);

% number of time steps (minimum T*NX for stability)
NT = T*NX;
% time step / temporal mesh size
dt = T/(NT-1);

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
title("Numerical solution with constant cross-section $ S(x) = 1 $ m ($N_X$ = " + NX + ", $N_T$ = "+ NT + ")");
xlabel('$x$ [m]');
ylabel('$t$ [s]');
zlabel('$u_h(x,t)$');
colorbar

% approximation difference (in space and time)
err = sol - uex;

% plotting error (in abs) as surface
figure()
[X_axis, T_axis] = ndgrid(x, t);
surf(X_axis, T_axis, abs(err));
shading interp
title("Numerical residual (in abs) with constant cross-section $ S(x) = 1 $ m ($N_X$ = " + NX + ", $N_T$ = "+ NT + ")");
xlabel('$x$ [m]');
ylabel('$t$ [s]');
zlabel('$|u_h(x,t) - u_{ex}(x,t)|$');
colorbar

%% Point 2, instability example (if you run this, take care of running the previous blocks again)
NX = 800;
NT = 2000;

dx = L/(NX-1);
dt = T/(NT-1);
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
title("Numerical solution with variable cross-section $ S(x) = (1+2x)^2 $ ($N_X$ = " + NX + ", $N_T$ = "+ NT + ")");
xlabel('$x$ [m]');
ylabel('$t$ [s]');
zlabel('$u_h(x,t)$');
colorbar

% approximation difference (in space and time)
err_var = sol_var - uex;

% plotting error (in abs) as surface
figure()
[X_axis, T_axis] = ndgrid(x, t);
surf(X_axis, T_axis, abs(err_var));
shading interp
title("Numerical error (in abs) with variable cross-section $ S(x) = (1+2x)^2 $ ($N_X$ = " + NX + ", $N_T$ = "+ NT + ")");
xlabel('$x$ [m]');
ylabel('$t$ [s]');
zlabel('$|u_h(x,t) - u_{ex}(x,t)|$');
colorbar

%% Looping with different NX and NT
% this block loops the computation for both constant and variable profile S(x)
% takes different NX and NT and evaluates the L2-error

NX = [5, 10, 100, 200, 1000];
NT = [5, 10, 15, 20, 25]*1e3;

e_L2 = zeros(length(NX));
e_L2_var = zeros(length(NX));

for i = 1:length(NX)
    for j = 1:length(NT)
        dx = L/(NX(i)-1);
        dt = T/(NT(j)-1);
        lambda = gamma*dt/dx;
        
        x = linspace(0, L, NX(i)).';
        t = linspace(0, T, NT(j));
        
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
        
        sol = computeSolutionConstant(NX(i), NT(j), dt, u_0, u_1, lambda, f);
        e_L2(i, j) = sqrt(sum(sum((sol - uex).^2)) * dx * dt);
        
        sol_var = computeSolutionVariable(NX(i), NT(j), dx, dt, S, u_0, u_1, lambda, f_var);
        e_L2_var(i, j) = sqrt(sum(sum((sol_var - uex).^2)) * dx * dt);
    end

end

%% Looping with different NT, keeping NX
% this block loops the computation for both constant and variable profile S(x)
% takes different NT, keeps NX and evaluates the L2-error

NX = 100;       % (100)=1% of space domain
NT = NX * round(linspace(5, 25, 20));

dx = L/(NX-1);
x = linspace(0, L, NX).';
Dt = T./(NT-1);

e_L2 = zeros(1, length(NT));
e_L2_var = zeros(1, length(NT));

dim = 1;    % L2-norm in space

for i = 1:length(NT)

    dt = Dt(i);
    lambda = gamma*dt/dx;
    
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
    
    sol = computeSolutionConstant(NX, NT(i), dt, u_0, u_1, lambda, f);
    e_L2(i) = mean(sqrt(sum((sol - uex).^2, dim) * dx));
  
    sol_var = computeSolutionVariable(NX, NT(i), dx, dt, S, u_0, u_1, lambda, f_var);
    e_L2_var(i) = mean(sqrt(sum((sol_var - uex).^2, dim) * dx));

end

% scheme convergence
p_time = log(e_L2(1:end-1)./e_L2(2:end)) ./ log(Dt(1:end-1)./Dt(2:end));
p_time_var = log(e_L2_var(1:end-1)./e_L2_var(2:end)) ./ log(Dt(1:end-1)./Dt(2:end));
p_time_avg = sum(p_time)/length(p_time);
p_time_var_avg = sum(p_time_var)/length(p_time_var);

% printing norm errors in function of dt
figure()
plot(log(Dt), log(e_L2), '-o');
title("Error convergence with time ($N_X$ = " + NX + ")");
xlabel("$\log(\Delta t)$");
ylabel("$\log(e)$");
grid on
hold on
plot(log(Dt), log(e_L2_var), '-o');
legend("Constant profile, p = " + round(p_time_avg), "Variable profile, p = " + round(p_time_var_avg));

%% Looping with different NX, keeping NT
% this block loops the computation for both constant and variable profile S(x)
% takes different NX, keeps NT and evaluates the L2-error

NT = 500;       % (500)=1% of time domain
NX = round(NT * (linspace(1, 20, 20)/100));

dt = T/(NT-1);
t = linspace(0, T, NT);
Dx = L./(NX-1);

e_L2 = zeros(1, length(NX));
e_L2_var = zeros(1, length(NX));

dim = 2;    % L2-norm in time

for i = 1:length(NX)

    dx = Dx(i);
    lambda = gamma*dt/dx;
    
    x = linspace(0, L, NX(i)).';
    
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
    
    sol = computeSolutionConstant(NX(i), NT, dt, u_0, u_1, lambda, f);
    e_L2(i) = mean(sqrt(sum((sol - uex).^2, dim) * dt));
    
    sol_var = computeSolutionVariable(NX(i), NT, dx, dt, S, u_0, u_1, lambda, f_var);
    e_L2_var(i) = mean(sqrt(sum((sol_var - uex).^2, dim) * dt));

end

% scheme convergence
p_space = log(e_L2(1:end-1)./e_L2(2:end)) ./ log(Dx(1:end-1)./Dx(2:end));
p_space_var = log(e_L2_var(1:end-1)./e_L2_var(2:end)) ./ log(Dx(1:end-1)./Dx(2:end));
p_space_avg = sum(p_space)/length(p_space);
p_space_var_avg = sum(p_space_var)/length(p_space_var);

% printing norm errors in function of dx
figure()
plot(log(Dx), log(e_L2), '-o');
title("Error convergence with space ($N_T$ = " + NT + ")");
xlabel("$\log(\Delta x)$");
ylabel("$\log(e)$");
grid on
hold on
plot(log(Dx), log(e_L2_var), '-o');
legend("Constant profile", "Variable profile");
legend("Constant profile, p = " + round(p_space_avg), "Variable profile, p = " + round(p_space_var_avg));

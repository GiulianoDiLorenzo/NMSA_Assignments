% Solution of HOMEWORK 2
%     S utt = gamma^2 (S ux)x + f
%     S(x) variable cross-section
%     f(x) force
%     u(x,t) solution

% Neumann conditions
%     ux(0) = ux(1) = 0
% 
% Initial profile and derivative
%     u(x,0), u_t(x,0)

% Finite difference method (SIMPLIFIED VERSION)
%     u(k,n+1) - 2*u(k,n) + u(k,n-1) = 
%     = lambda^2 *(u(k+1,n) - 2*u(x,n) + u(k-1,n)) + f(k,n)

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
% velocity, length and gamma
c = 1;
L = 1;
gamma = c/L;

% space domain (0,1) and time domain (0,T)
T = 5;
I = [0 L];

% number of time steps (minimum 5*NX for stability)
NT = 4000;
% time step / temporal mesh size
dt = T/NT;

% number of space steps
NX = 800;
% space step / spatial mesh size
dx = (I(2) - I(1))/NX;

% lambda
lambda = gamma*dt/dx;

% space domain definition: 0 to L in NX elements
x = linspace(I(1), I(2), NX).';

% time domain definition: 0 to T in NT elements
t = linspace(0, T, NT);

% section definition
S_const = 1;
S = (1+2*x).^2;
Sx = 4 + 8*x;   % derivative in x of S_var

%% Point 2, exact solution
% exact solution for computing the error
% uex   = cos(pi*(x/2+1)) * cos(3*pi*t);
uex     = - cos(pi/2*x) * cos(3*pi*t);

% 1st derivatives
ux      = pi/2 * sin(pi/2*x) * cos(3*pi*t);
ut      = 3*pi * cos(pi/2*x) * sin(3*pi*t);

% 2nd derivatives
% utt = -(3*pi)^2 * cos(pi*(x/2+1)) * cos(3*pi*t);
% uxx = -(pi/2)^2 * cos(pi*(x/2+1)) * cos(3*pi*t);
uxx     = (pi/2)^2 * cos(pi/2*x) * cos(3*pi*t);
utt     = (3*pi)^2 * cos(pi/2*x) * cos(3*pi*t);

% initial conditions
u_0 = uex(:,1);             % u_0 = u(x,0)
u_1 = ut(:,1);              % u_1 = du/dt(x,0)
u_t_1 = zeros(size(t));     % u_t_1 = du/dt(1,t) 
u_x_0 = zeros(size(t));     % u_x_0 = du/dx(0,t)

% excitation forces
f  = S_const * (utt - gamma^2 * uxx);
f_var = S .* utt - gamma^2 .* (Sx .* ux + S .* uxx);

%% Point 2, constant cross-section, leap-frog scheme
tic

% initializing the solution;
sol = zeros(NX+1, NT+1);

% applying the boundary conditions
sol(2:end, 2) = u_0;                        % 3rd condition
sol(2:end, 1) = sol(2:end, 2) - dt * u_1;   % 4th condition
sol(end,:) = u_0(end);                      % 2nd condition
sol(1,:) = sol(2,:);                        % 1st condition

for n = 2:NT      % evaluating n+1, up to NT to consider virtual line

    for k = 2:NX-1
        sol(k,n+1) = 2*sol(k,n) - sol(k,n-1) + lambda^2 * (sol(k+1,n) - 2*sol(k,n) + sol(k-1,n)) + dt^2 * f(k-1,n);
    end

    sol(1,n+1) = sol(2,n+1);                % 1st condition
end

% discarding virtual lines (1,:) and (:,1)
sol = sol(2:end, 2:end);

toc

% error evaluation
err = sol - uex;    % distance between numerical and exact solution
[max_err, ~] = max(abs(err(:)));        % maximum error (in abs) for numerical scheme

%% Plotting solution comparison with constant cross-section
figure()
for n = 1:NT
    plot(x, uex(:,n), "LineStyle", ":");
    hold on;
    plot(x, sol(:,n));
    legend('$u_{ex}$', '$u$');
    title("Solution comparison with constant cross-section S(x) = " + S_const + " m ($N_X$ = " + NX + ", $N_T$ = "+ NT + ") at t = "+ (n-1)*dt+ " s");
    xlabel('$x$ [m]');
    ylabel('$u(x,t), \quad u_{ex}(x,t)$');
    %ylim([-1,1]);
    grid on;
    pause(1e-5);
    hold off;
end

%% Plotting solution as surface
figure()
[X, T] = ndgrid(x, t);
surf(X, T, sol);
shading interp
title("Numerical solution with constant cross-section S(x) = " + S_const + " m ($N_X$ = " + NX + ", $N_T$ = "+ NT + ")");
xlabel('$x$ [m]');
ylabel('$t$ [s]');
zlabel('$u(x,t)$');
colorbar

%% Plotting error with constant cross-section
figure()
for n = 1:NT
    plot(x, err(:,n));
    title("Error evaluation with constant cross-section S(x) = " + S_const + " m ($N_X$ = " + NX + ", $N_T$ = "+ NT + ") at t = "+ (n-1)*dt+ " s");
    xlabel('$x$ [m]');
    ylabel('$u(x,t) - u_{ex}(x,t)$');
    %ylim([-1,1]);
    grid on;
    pause(1e-5);
    hold off;
end

%% Plotting error (in abs) as surface
figure()
[X, T] = ndgrid(x, t);
surf(X, T, abs(err));
shading interp
title("Numerical error (in abs) with constant cross-section S(x) = " + S_const + " m ($N_X$ = " + NX + ", $N_T$ = "+ NT + ")");
xlabel('$x$ [m]');
ylabel('$t$ [s]');
zlabel('$|u(x,t) - u_{ex}(x,t)|$');
colorbar

%% Point 2, variable cross-section, leap-frog scheme
tic

% initializing the solution;
sol_var = zeros(NX+1, NT+1);

% applying the boundary conditions
sol_var(2:end, 2) = u_0;                        % 3rd condition
sol_var(2:end, 1) = sol_var(2:end, 2) - dt * u_1;   % 4th condition
sol_var(end,:) = u_0(end);                      % 2nd condition
sol_var(1,:) = sol_var(2,:);                        % 1st condition

% extending S(x) to be coherent with the virtual lines
S = [(1-2*dx)^2; S];

for n = 2:NT      % evaluating n+1, up to NT to consider virtual line

    for k = 2:NX-1
        sol_var(k,n+1) = 2*sol_var(k,n) - sol_var(k,n-1) + lambda^2/S(k) * (((S(k+1)-S(k-1))/4 + S(k)) * sol_var(k+1,n) - 2*S(k) * sol_var(k,n) + ((S(k-1)-S(k+1))/4 + S(k)) * sol_var(k-1,n)) + dt^2/S(k) * f_var(k-1,n);
    end

    sol_var(1,n+1) = sol_var(2,n+1);                % 1st condition
end

% discarding virtual lines
sol_var = sol_var(2:end, 2:end);
S = S(2:end);

toc

% error evaluation
err_var = sol_var - uex;    % distance between numerical and exact solution
[max_err, ~] = max(abs(err_var(:)));        % maximum error (in abs) for numerical scheme

%% Plotting solution comparison with variable cross-section
figure()
for n = 1:NT
    plot(x, uex(:,n), "LineStyle", ":");
    hold on;
    plot(x, sol_var(:,n));
    legend('$u_{ex}$', '$u$');
    title("Solution comparison with variable cross-section S(x) = $(1+2x)^2$ ($N_X$ = " + NX + ", $N_T$ = "+ NT + ") at t = "+ (n-1)*dt+ " s");
    xlabel('$x$ [m]');
    ylabel('$u(x,t), \quad u_{ex}(x,t)$');
    %ylim([-1,1]);
    grid on;
    pause(1e-5);
    hold off;
end

%% Plotting solution as surface
figure()
[X, T] = ndgrid(x, t);
surf(X, T, sol_var);
shading interp
title("Numerical solution with variable cross-section S(x) = $(1+2x)^2$ ($N_X$ = " + NX + ", $N_T$ = "+ NT + ")");
xlabel('$x$ [m]');
ylabel('$t$ [s]');
zlabel('$u(x,t)$');
colorbar

%% Plotting error with variable cross-section
figure()
for n = 1:NT
    plot(x, err_var(:,n));
    title("Error evaluation with variable cross-section S(x) = $(1+2x)^2$ ($N_X$ = " + NX + ", $N_T$ = "+ NT + ") at t = "+ (n-1)*dt+ " s");
    xlabel('$x$ [m]');
    ylabel('$u(x,t) - u_{ex}(x,t)$');
    %ylim([-1,1]);
    grid on;
    pause(1e-5);
    hold off;
end

%% Plotting error (in abs) as surface
figure()
[X, T] = ndgrid(x, t);
surf(X, T, abs(err_var));
shading interp
title("Numerical error (in abs) with variable cross-section S(x) = $(1+2x)^2$ ($N_X$ = " + NX + ", $N_T$ = "+ NT + ")");
xlabel('$x$ [m]');
ylabel('$t$ [s]');
zlabel('$|u(x,t) - u_{ex}(x,t)|$');
colorbar




%%






a1 =  c;
a2 = -c;

% Initial data
% u0  = @(x) ...;     % u(x,0)
% ut0 = @(x) ...;     % u_t(x,0)

% % Initial conditions for system u_t = A u_x  (u_t + A u_x = 0)
% v0  = @(x) 2*pi*sin(2*pi*x); % u_t(t=0)
% up0 = @(x) 0.*x;             % u_x(t=0)

% % Initial conditions for w_t = D w_x  (w_t + D w_x = 0)
% % by using the transformation w(x,0) = T^(-1) u(x,0)
% 
% w0 = @(x) sqrt(1+c^2)/(2*c).*( 2*pi*sin(2*pi*x));
% w1 = @(x) sqrt(1+c^2)/(2*c).*(-2*pi*sin(2*pi*x));

% Neumann boundary conditions
gp1 = @(t) 0.*t;
gp2 = @(t) 0.*t;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% menu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flag = menu('Please select a method',...
%     'Forward-Euler/Central',...
%     'Lax-Friedrichs',...
%     'Lax-Wendroff',...
%     'Upwind');
% 
% disp('Please wait ...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Initial conditions
% w = [w_1, w_2]'; 
% u = [u_1, u_2]';
% SOL_w(1:NX+1,:) -> w_1; % SOL_w(NX+2:end,:) -> w_2;
% the same for SOL_u


SOL_w = zeros(2*(NX+1),NT+1);
SOL_u = zeros(2*(NX+1),NT+1);

% initial values
for j = 1 : NX+1
    SOL_w(j,1)      = w0(I(1) + (j-1)*dx);
    SOL_w(NX+1+j,1) = w1(I(1) + (j-1)*dx);
    SOL_u(j,1)      = v0(I(1) + (j-1)*dx);
    SOL_u(NX+1+j,1) = up0(I(1) + (j-1)*dx);

end


for n = 1:NT
    for k = 2:NX
        
        if flag == 1 % FCE scheme
            SOL_w(j,n+1)      = SOL_w(j,n)      - 0.5*lambda*a1*(SOL_w(j+1,n)-SOL_w(j-1,n)); % w_1
            SOL_w(NX+1+j,n+1) = SOL_w(NX+1+j,n) - 0.5*lambda*a2*(SOL_w(NX+1+j+1,n)-SOL_w(NX+1+j-1,n)); % w_2
        elseif flag == 2 % LF scheme
            SOL_w(j,n+1)      = 0.5*(SOL_w(j+1,n)+SOL_w(j-1,n))           - 0.5*lambda*a1*(SOL_w(j+1,n)-SOL_w(j-1,n)); % w_1
            SOL_w(NX+1+j,n+1) = 0.5*(SOL_w(NX+1+j+1,n)+SOL_w(NX+1+j-1,n)) - 0.5*lambda*a2*(SOL_w(NX+1+j+1,n)-SOL_w(NX+1+j-1,n)); %w_2
        elseif flag == 3 % LW scheme
            SOL_w(j,n+1) = SOL_w(j,n) - 0.5*lambda*a1*(SOL_w(j+1,n)-SOL_w(j-1,n)) ...
                              + 0.5*lambda^2*a1^2*(SOL_w(j+1,n)-2*SOL_w(j,n)+SOL_w(j-1,n));
            SOL_w(NX+1+j,n+1) = SOL_w(NX+1+j,n) - 0.5*lambda*a2*(SOL_w(NX+1+j+1,n)-SOL_w(NX+1+j-1,n)) ...
                              + 0.5*lambda^2*a2^2*(SOL_w(NX+1+j+1,n)-2*SOL_w(NX+1+j,n)+SOL_w(NX+1+j-1,n));
        elseif flag == 4 % UPWIND
            SOL_w(j,n+1) = SOL_w(j,n) - 0.5*lambda*a1*(SOL_w(j+1,n)-SOL_w(j-1,n)) ...
                              + 0.5*lambda*abs(a1)*(SOL_w(j+1,n)-2*SOL_w(j,n)+SOL_w(j-1,n));
            SOL_w(NX+1+j,n+1) = SOL_w(NX+1+j,n) - 0.5*lambda*a2*(SOL_w(NX+1+j+1,n)-SOL_w(NX+1+j-1,n)) ...
                              + 0.5*lambda*abs(a2)*(SOL_w(NX+1+j+1,n)-2*SOL_w(NX+1+j,n)+SOL_w(NX+1+j-1,n));

        end
    end
    
    % Boundary conditions
    % can be improved by adding external end points --> SOL_w(1:2*(NX+1)+4,:)
    SOL_w(NX+2,n+1) = SOL_w(NX+3,n+1); %constant extrapolation w_2(0,t) = w_2(h,t)
    % update w_1(0,t)
    SOL_w(1,n+1) = sqrt(1+c^2)/c*gp1((n+1)*dt) + SOL_w(NX+2,n+1);


    SOL_w(NX+1,n+1) = SOL_w(NX,n+1);  %constant extrapolation w_1(L,t) = w_1(L-h,t)
    % update w_2(L,t)
    SOL_w(end,n+1)  = SOL_w(NX+1,n+1)-sqrt(1+c^2)/c*gp2((n+1)*dt);

    
    % compute u by using the transformation u = T w
    % w = [w_1, w_2] = [u_t, u_x];
    SOL_u(1:NX+1,n+1)   = c/sqrt(1+c^2)*(SOL_w(1:NX+1,n+1)-SOL_w(NX+2:end,n+1));
    SOL_u(NX+2:end,n+1) = 1/sqrt(1+c^2)*(SOL_w(1:NX+1,n+1)+SOL_w(NX+2:end,n+1));
end

% Computation of the L2-error for t=T
SOLu_EX  = zeros(NX+1,1);
SOLut_EX = zeros(NX+1,1);
SOLux_EX = zeros(NX+1,1);

for j = 1:NX+1
    SOLu_EX(j,1)  = u(I(1) + (j-1)*dx,NT*dt);
    SOLut_EX(j,1) = ut(I(1) + (j-1)*dx,NT*dt);
    SOLux_EX(j,1) = ux(I(1) + (j-1)*dx,NT*dt);

end
L2ut_ERR   = norm(SOLut_EX-SOL_u(1:NX+1,n+1),2)*dx^0.5;
L1ut_ERR   = norm(SOLut_EX-SOL_u(1:NX+1,n+1),1)*dx;
LINFut_ERR = norm(SOLut_EX-SOL_u(1:NX+1,n+1),Inf);

L2ux_ERR   = norm(SOLux_EX-SOL_u(NX+2:end,n+1),2)*dx^0.5;
L1ux_ERR   = norm(SOLux_EX-SOL_u(NX+2:end,n+1),1)*dx;
LINFux_ERR = norm(SOLux_EX-SOL_u(NX+2:end,n+1),Inf);

Err_t = [L2ut_ERR, L1ut_ERR, LINFut_ERR]
Err_x = [L2ux_ERR, L1ux_ERR, LINFux_ERR]


% Visualization
time = ones(NX+1,1)*linspace(0,T,NT+1);
space = linspace(I(2),I(1),NX+1)'*ones(1,NT+1);
surf(time,space,SOL_u(1:NX+1,:),'EdgeColor','none');
title('u_t','FontSize',16);
xlabel('time t','FontSize',16);
ylabel('space x','FontSize',16);
view(2);

figure;
surf(time,space,SOL_u(NX+2:end,:),'EdgeColor','none');
title('u_x','FontSize',16);
xlabel('time t','FontSize',16);
ylabel('space x','FontSize',16);
view(2);

%charachteristic variables
figure;
subplot(1,2,1);
surf(time,space,SOL_w(1:NX+1,:),'EdgeColor','none');
title('w_1','FontSize',16);
xlabel('time t','FontSize',16);
ylabel('space x','FontSize',16);
view(2);
subplot(1,2,2);
surf(time,space,SOL_w(NX+2:end,:),'EdgeColor','none');
title('w_2','FontSize',16);
xlabel('time t','FontSize',16);
ylabel('space x','FontSize',16);
view(2);





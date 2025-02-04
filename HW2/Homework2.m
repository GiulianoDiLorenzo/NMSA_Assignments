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

clear; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% exact solution for computing the error
u  = @(x,t) ...;
ut = @(x,t) ...;
ux = @(x,t) ...;


% velocity
c = 1;
a1 =  c;
a2 = -c;

% Initial data
u0  = @(x) ...;     % u(x,0)
ut0 = @(x) ...;     % u_t(x,0)

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

% Space domain (0,1) and time domain (0,T)
T = 5;
I = [0 1];

% Number of time steps
NT = 8000;
% Number of space steps
NX = 400;


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

dt = T/NT;              % time step (Delta t)
dx = (I(2)-I(1))/NX;    % mesh size (h)
lambda = dt/dx;

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
        
        % if flag == 1 % FCE scheme
        %     SOL_w(j,n+1)      = SOL_w(j,n)      - 0.5*lambda*a1*(SOL_w(j+1,n)-SOL_w(j-1,n)); % w_1
        %     SOL_w(NX+1+j,n+1) = SOL_w(NX+1+j,n) - 0.5*lambda*a2*(SOL_w(NX+1+j+1,n)-SOL_w(NX+1+j-1,n)); % w_2
        % elseif flag == 2 % LF scheme
        %     SOL_w(j,n+1)      = 0.5*(SOL_w(j+1,n)+SOL_w(j-1,n))           - 0.5*lambda*a1*(SOL_w(j+1,n)-SOL_w(j-1,n)); % w_1
        %     SOL_w(NX+1+j,n+1) = 0.5*(SOL_w(NX+1+j+1,n)+SOL_w(NX+1+j-1,n)) - 0.5*lambda*a2*(SOL_w(NX+1+j+1,n)-SOL_w(NX+1+j-1,n)); %w_2
        % elseif flag == 3 % LW scheme
        %     SOL_w(j,n+1) = SOL_w(j,n) - 0.5*lambda*a1*(SOL_w(j+1,n)-SOL_w(j-1,n)) ...
        %                       + 0.5*lambda^2*a1^2*(SOL_w(j+1,n)-2*SOL_w(j,n)+SOL_w(j-1,n));
        %     SOL_w(NX+1+j,n+1) = SOL_w(NX+1+j,n) - 0.5*lambda*a2*(SOL_w(NX+1+j+1,n)-SOL_w(NX+1+j-1,n)) ...
        %                       + 0.5*lambda^2*a2^2*(SOL_w(NX+1+j+1,n)-2*SOL_w(NX+1+j,n)+SOL_w(NX+1+j-1,n));
        % elseif flag == 4 % UPWIND
        %     SOL_w(j,n+1) = SOL_w(j,n) - 0.5*lambda*a1*(SOL_w(j+1,n)-SOL_w(j-1,n)) ...
        %                       + 0.5*lambda*abs(a1)*(SOL_w(j+1,n)-2*SOL_w(j,n)+SOL_w(j-1,n));
        %     SOL_w(NX+1+j,n+1) = SOL_w(NX+1+j,n) - 0.5*lambda*a2*(SOL_w(NX+1+j+1,n)-SOL_w(NX+1+j-1,n)) ...
        %                       + 0.5*lambda*abs(a2)*(SOL_w(NX+1+j+1,n)-2*SOL_w(NX+1+j,n)+SOL_w(NX+1+j-1,n));
        % 
        % end

        SOL_u(k,n+1) = lambda^2 *(SOL_u(k+1,n) - 2*SOL_u(k,n) + SOL_u(k-1,n)) + F(k,n) + 2*SOL_u(k,n) - SOL_u(k,n-1);
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





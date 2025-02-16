%% Data for the following PDE
%   rho(x)* w * u(x)  + d(mu*du/dx)/dx = f in (0,L)
%   Dirichlet condition: u(0) = g_D
%   Robin condition : u(L) - j/(rho(x)*w) * du(x=L)/dx = g_A

function [Data] = NewDataTest(TestName, numEl, omega)

% if strcmp(TestName,'Test1')
%     %% Input Data Test1
%     % mu*u''(x) + u(x) = (4*pi^2+1)*cos(2*pi*x)  x in (0,1)
%     %     u'(0) = u'(1) = 0;
% 
%     Data.name = TestName;
%     Data.domain = [0,1];
%     Data.boundary = 'NN';
% 
%     % Parameters and external forces 
%     Data.mu = 1;
%     Data.sigma = 1;
%     Data.force = @(x) (4*pi^2+1)*cos(2*pi*x);
%     Data.gN = @(x) 0.*x;
%     Data.gD = @(x) 0.*x;
% 
%     % Exact solution for error analysis
%     Data.uex = @(x) cos(2*pi*x);
%     Data.uex_x = @(x) -2*pi*sin(2*pi*x);
% 
% 
% elseif strcmp(TestName,'Test2')
%     %% Input Data Test2
%     % mu * u''(x) + rho*omega*u= 0  x in (0,L)
%     % u(0) = u(L) = 0;
%     L = 1; % Domain end
% 
%     Data.name = TestName;
%     Data.domain = [0,L];
%     Data.L = L;
%     Data.boundary = 'DD';
% 
%     % Parameters and external forces 
%     Data.mu = 1;
%     Data.rho = 1;
%     Data.omega = Data.mu / Data.rho * (pi/Data.L)^2;
%     Data.force = @(x) 0.*x;
%     Data.gN = @(x) 0.*x; % Dirichlet condition on the left part of the domain
%     Data.gD = @(x) 0.*x; % Dirichlet condition on the right part of the domain
% 
%     % Exact solution for error analysis
%     Data.uex = @(x) sin(sqrt(Data.rho*Data.omega/Data.L)*x);
%     Data.uex_x = @(x) sqrt(Data.rho*Data.omega/Data.L)*cos(sqrt(Data.rho*Data.omega/Data.L)*x);


if strcmp(TestName,'TestHW1_1')
    % INPUT DATA HERE
    %   rho(x)* w * u(x)  + d(mu*du/dx)/dx = f in (0,L)
    %   Dirichlet condition: u(0) = g_D
    %   Robin condition : u(L) - j/(rho(x)*w) * du(x=L)/dx = g_A

    L = 1; % Domain end

    Data.name = TestName;
    Data.domain = [0,L];
    Data.L = L;
    Data.boundary = 'DR';

    % Exact solution for error analysis
    Data.uex = @(x) sin(2*pi*x);
    Data.uex_x = @(x) 2*pi*cos(2*pi*x);
    Data.uex_xx = @(x) -(2*pi)^2*uex(x);
    
    % Parameters and external forces, simple case with constant terms rho,mu
    Data.mu = @(x) 2.*ones(size(x));
    Data.rho = @(x) 1.*ones(size(x));
    Data.omega = omega;

    x = linspace(0,L, numEl);
    Data.x = x;

    % Compute mu(x) values
    mu_values = Data.mu(x);
    % Compute numerical derivative using forward finite difference dmu/dx
    dmu_dx = diff(mu_values) ./ diff(x);

    % Putting the last point equal its predecessor
    Data.mu_x =  [dmu_dx, dmu_dx(end)];  % NaN for the last value (since diff reduces size)
    
    % d(mu * d(u_ex/dx))/dx
    d_mu_u_ex_x = diff(Data.mu(x) .* Data.uex_x(x) ) ./ diff(x);

    Data.d_mu_u_ex_x  =  [d_mu_u_ex_x , d_mu_u_ex_x(end)];  % NaN for the last value (since diff reduces size)

    

    Data.alpha = Data.mu(L) .* Data.rho(L) * Data.omega ; % We're interested only in alpha(end), i.e. at x=L
    
    Data.force = @(x) (Data.rho(x)*Data.omega - Data.mu(x)*(2*pi)^2) .* Data.uex(x);


    %  % Non constant profiles
    % Data.force = @(x) Data.rho(x) .* Data.omega .* Data.uex(x) + Data.d_mu_u_ex
    % Data.d_mu_u_ex()x
    % (Data.rho(x)*Data.omega - Data.mu(x)*(2*pi)^2) .* Data.uex;

    % Dirichlet condition on the left part of the domain

    Data.gD = 0 ;

    % Robin on the right part of the domain
    Data.gR = Data.uex(L) - 1i * Data.uex_x(L) ./ (Data.rho(L) .* Data.omega); 

end



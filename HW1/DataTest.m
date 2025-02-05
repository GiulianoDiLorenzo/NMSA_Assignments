%% Data for the forllowing PDE
%  mu*u'' + rho*omega*u = f in (a,b)
%   Dirichlet or Neumann (homo) conditions  

function [Data] = DataTest(TestName)

if strcmp(TestName,'Test1')
    %% Input Data Test1
    % mu*u''(x) + u(x) = (4*pi^2+1)*cos(2*pi*x)  x in (0,1)
    %     u'(0) = u'(1) = 0;
    
    Data.name = TestName;
    Data.domain = [0,1];
    Data.boundary = 'NN';
    
    % Parameters and external forces 
    Data.mu = 1;
    Data.sigma = 1;
    Data.force = @(x) (4*pi^2+1)*cos(2*pi*x);
    Data.gN = @(x) 0.*x;
    Data.gD = @(x) 0.*x;
    
    % Exact solution for error analysis
    Data.uex = @(x) cos(2*pi*x);
    Data.graduex = @(x) -2*pi*sin(2*pi*x);
    
    
elseif strcmp(TestName,'Test2')
    %% Input Data Test2
    % mu * u''(x) + rho*omega*u= 0  x in (0,L)
    % u(0) = u(L) = 0;
    L = 1; % Domain end

    Data.name = TestName;
    Data.domain = [0,L];
    Data.L = L;
    Data.boundary = 'DD';
    
    % Parameters and external forces 
    Data.mu = 1;
    Data.rho = 1;
    Data.omega = Data.mu / Data.rho * (pi/Data.L)^2;
    Data.force = @(x) 0.*x;
    Data.gN = @(x) 0.*x; % Dirichlet condition on the left part of the domain
    Data.gD = @(x) 0.*x; % Dirichlet condition on the right part of the domain
    
    % Exact solution for error analysis
    Data.uex = @(x) sin(sqrt(Data.rho*Data.omega/Data.L)*x);
    Data.graduex = @(x) sqrt(Data.rho*Data.omega/Data.L)*cos(sqrt(Data.rho*Data.omega/Data.L)*x);


elseif strcmp(TestName,'TestHW1_1')
    %% Input Data Test2
    % mu * u''(x) + rho*omega*u= 0  x in (0,L)
    % u(0) = 0;
    % g_A = du/dx(x=L) = 0
    L = 1; % Domain end

    Data.name = TestName;
    Data.domain = [0,L];
    Data.L = L;
    Data.boundary = 'DR';
    
    % Parameters and external forces, simple case with constant terms rho,mu
    Data.mu = @(x) ones(size(x));
    Data.rho = @(x) ones(size(x));
    Data.omega = 1;
    Data.alpha = @(x) Data.mu(x) .* Data.rho(x) * Data.omega ./ 1i;

    Data.force = @(x) (Data.rho(x) * Data.omega - Data.mu(x)*(2*pi)^2) .* sin(2*pi*x);
    Data.gD = @(x) 0.*x;                % Dirichlet condition on the left part of the domain
    Data.gR = @(x) [zeros(size(x)-1), sin(2*pi*Data.L) - 1i*2*pi/(Data.rho(L).*Data.omega)*cos(2*pi*L)]  % Dirichlet condition on the right part of the domain
    
    % Exact solution for error analysis
    Data.uex = @(x) sin(2*pi*x)
    Data.graduex = @(x) 2*pi*cos(2*pi*x);

elseif strcmp(TestName,'TestHW1_2')
    %% INCOMPLETE ASK MAZZIERI
    %% Input Data Test2
    % mu * u''(x) + rho*omega*u= 0  x in (0,L)
    % u(0) = 0;
    % g_A = du/dx(x=L) = 0
    L = 1; % Domain end

    Data.name = TestName;
    Data.domain = [0,L];
    Data.L = L;
    Data.boundary = 'DR';
    
    % Parameters and external forces, simple case with constant terms rho,mu
    mu1 = 1;
    mu2 = 1;
    rho1 = 1;
    rho2 = 1;
    Data.mu = @(x) mu1*(x <= L/2) + mu2 * (x > L/2);
    Data.rho = @(x) rho1*(x <= L/2) + rho2 * (x > L/2);
    Data.omega = 1;
    Data.alpha = @(x) Data.mu(x) .* Data.rho(x) * Data.omega ./ 1i;

    Data.force1 = @(x) ((rho1*Data.omega)* sin(2*pi*x) + mu1*2*pi*cos(2*pi*x)) .* (x <= L/2);
    Data.force2 = @(x) ((rho2*Data.omega)* sin(2*pi*x) + mu2*2*pi*cos(2*pi*x)) .* (x > L/2);
    Data.force = @(x) Data.force1(x) + Data.force2(x);

    Data.gD = @(x) 0.*x;                % Dirichlet condition on the left part of the domain
    Data.gR = @(x) [zeros(size(x)-1), sin(2*pi*Data.L) - 1i*2*pi/(rho2*Data.omega)*cos(2*pi*L)]  % Dirichlet condition on the right part of the domain
    
    % Exact solution for error analysis
    Data.uex = @(x) sin(2*pi*x)
    Data.graduex = @(x) 2*pi*cos(2*pi*x);


elseif strcmp(TestName,'TestHW1_3a')
    %% Input Data Test2
    % mu * u''(x) + rho*omega*u= 0  x in (0,L)
    % u(0) = 0;
    % g_A = du/dx(x=L) = 0
    L = 5; % Domain end
    Data.f_0 = 3;
    Data.name = TestName;
    Data.domain = [0,L];
    Data.L = L;
    Data.boundary = 'DR';
    
    % Parameters and external forces, simple case with constant terms rho,mu
    Data.mu = @(x) ones(size(x));
    Data.c = @(x) ones(size(x));

    Data.rho = @(x) Data.mu(x) / (Data.c(x).^2);
    Data.omega = 1;
    Data.alpha = @(x) Data.mu(x) .* Data.rho(x) * Data.omega ./ 1i;


    Data.force = @(x) 0.*x;

    Data.gD = @(x) 2*Data.omega^2/(sqrt(pi)*Data.f_0^2)*exp(-Data.omega^2/Data.f_0^2) .*ones(size(x));                % Dirichlet condition on the left part of the domain
    Data.gR = @(x) 0.*x;  % Dirichlet condition on the right part of the domain
    
    % Exact solution for error analysis
    %Data.uex = @(x) sin(2*pi*x)
    %Data.graduex = @(x) 2*pi*cos(2*pi*x);
end

%% Data for the forllowing PDE
% - mu*u'' + sigma*u = f in (a,b)
%   Dirichlet or Neumann (homo) conditions  

function [Data] = DataTest(TestName)

if strcmp(TestName,'Test1')
    %% Input Data Test1
    % - u''(x) + u(x) = (4*pi^2+1)*cos(2*pi*x)  x in (0,1)
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
    %- u''(x) = 16*cos(4*x)  x in (0,2*pi)
    %    u(0) = u(2*pi) = 0;

    Data.name = TestName;
    Data.domain = [0,2*pi];
    Data.boundary = 'DD';
    
    % Parameters and external forces 
    Data.mu = 1;
    Data.sigma = 0;
    Data.force = @(x) 16*cos(4*x);
    Data.gN = @(x) 0.*x;
    Data.gD = @(x) 1+0.*x;
    
    % Exact solution for error analysis
    Data.uex = @(x) cos(4*x);
    Data.graduex = @(x) -4*sin(4*x);
    
end

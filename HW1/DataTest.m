%% Data for the following PDE

function [Data] = DataTest(TestName)

if strcmp(TestName,'Test1')
    %% Input Data Test1
    
    Data.name = TestName;
    Data.domain = [0,1];
    % D = Dirichlet; N = Neumann; R = Robin; I = Impedance; P = periodic; 
    Data.boundary = 'NN';
    
    % Parameters and external forces 
    Data.ro    = 1.225;
    Data.omega = 170*2*pi;      % 100*2*pi and 500*2*pi
    Data.vel   = 340; 
    
    Data.force = @(x) 0.*x;
    Data.gN1   = @(omega,ro,vel) (-1.e-3*omega.^2.*ro) + 0.*vel; % Neumann at x = 0
    Data.gN2  =  @(omega,ro,vel) 0.*omega.*ro + 0.*vel;          % Neumann at x = 1
    Data.gD1  =  @(omega,ro,vel) 0.*ro.*vel.*omega;  %x = 0
    Data.gD2  =  @(omega,ro,vel) 0.*ro.*vel.*omega;  %x = 1
    

    % Exact solution for error analysis
    Data.uex     = @(x,omega,ro,vel) 1.e-3.*ro.*vel.*omega.*(sin(omega.*x./vel)+ 1/tan(omega*1./vel).*cos(omega.*x./vel));
    Data.graduex = @(x,omega,ro,vel) 1.e-3.*ro.*vel.*omega.*(omega./vel.*cos(omega.*x./vel) - 1/tan(omega.*1./vel).*omega./vel.*sin(omega.*x./vel));
    
elseif strcmp(TestName,'Test2')
    %% Input Data Test1

    % p''(x) + omega^2 p(x) = f;
    % p(0) = p(3) = 0;
           
    % solution is p = sin(x*2*pi/3)*1/(omega^2-(2*pi/3)^2)
    % f = sin(2*pi/3 * x)
    
    Data.name = TestName;
    Data.domain = [0,3];
    % D = Dirichlet; N = Neumann; R = Robin; I = Impedance; P = periodic; 
    Data.boundary = 'DD';
    
    % Parameters and external forces 
    Data.ro    = 1;
    Data.omega = 2*pi/3;
    Data.vel   = 1; 
    
    Data.force = @(x) sin(x*2*pi/3);
    Data.gN1   = @(omega,ro,vel) 2*pi/3*1/(omega^2-(2*pi/3)^2) + 0.*ro.*vel; %x = 0 ;
    Data.gN2  =  @(omega,ro,vel) 2*pi/3*1/(omega^2-(2*pi/3)^2) + 0.*ro.*vel;  %x = 1
    Data.gD1  =  @(omega,ro,vel) 0.*ro.*vel.*omega;  %x = 0
    Data.gD2  =  @(omega,ro,vel) 0.*ro.*vel.*omega;  %x = 1

    
    % Exact solution for error analysis
    Data.uex     = @(x,omega,ro,vel) sin(x*2*pi/3).*1./(omega.^2-(2*pi/3)^2) + 0.*ro.*vel;
    Data.graduex = @(x,omega,ro,vel) 2*pi/3.*cos(x*2*pi/3).*1./(omega.^2-(2*pi/3)^2) + 0.*ro.*vel;

end

           




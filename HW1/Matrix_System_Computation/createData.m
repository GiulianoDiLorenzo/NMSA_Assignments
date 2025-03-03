function [Data] = createData(TestName, L, N_pts, omega, mu_vals, rho_vals )
% ========================================================================
%   OUTPUT : Structure of the Galerkin formulation with problem's functions
%            discretization

%   INPUTS : 
%       - TestName  --> Identification for questions of HW1
%       - L         --> Length of the domin [m]
%       - N_pts     --> Number of nodes in the mesh
%       - omega     --> Scalar value for omega parameter
%       - mu_vals   --> array of (1,2) for the mu values on each half of
%                       the domain
%       - rho_vals  --> array of (1,2) for the rho values on each half of
%                       the domain
% ========================================================================

Data.name = TestName;
Data.L = L;
Data.domain = [0,L];
Data.boundary = 'DR';
Data.omega = omega;


Data.mu_vals = mu_vals;
mu1 = mu_vals(1);    % Value of mu in [0, L/2[
mu2 = mu_vals(2);    % Value of mu in [L/2, L]

% Adding the mu values to Data
Data.mu1 = mu1;  
Data.mu2 = mu2;    


Data.rho_vals = rho_vals;
rho1 = rho_vals(1);    % Value of rho in [0, L/2[
rho2 = rho_vals(2);    % Value of rho in [L/2, L]

% Adding the rho values to Data
Data.rho1 = rho1; 
Data.rho2 = rho2;    

Data.alpha = mu2 * rho2 * omega ;

x = linspace(0,L, N_pts);
Data.x = x.';


if strcmp(TestName,'TestHW1_2a')
    
    % Exact solution for error analysis
    Data.uex = @(x)     sin(2*pi*x);
    Data.uex_x = @(x)   2*pi .* cos(2*pi*x);
    Data.uex_xx = @(x)  - (2*pi)^2 .* Data.uex(x);
    
    % Function handler for mu
    Data.mu = @(x) mu_vals(1).*ones(size(x));

    % Function handler for rho
    Data.rho = @(x) rho_vals(1).*ones(size(x));

    % Force vector
    Data.force = @(x) Data.rho(x) * Data.omega .* Data.uex(x) + Data.mu(x) .* Data.uex_xx(x);

    % Dirichlet condition on the left part of the domain
    Data.gD = Data.uex(0) ;

    % Robin on the right part of the domain
    Data.gR = Data.uex(L) - 1i * Data.uex_x(L) ./ (Data.rho(L) .* Data.omega); 

elseif strcmp(TestName,'TestHW1_2b')
    
    % [Data.uex , Data.uex_x, Data.uex_xx] = setExactSolution(Data);
    Data.uex = @(x)     sin(2*pi*x);
    Data.uex_x = @(x)   2*pi .* cos(2*pi*x);
    Data.uex_xx = @(x)  - (2*pi)^2 .* Data.uex(x);

    % Function handler for mu
    Data.mu = @(x) (x >= 0  & x < L/2) .* mu1 ...
                 + (x >= L/2 & x <= L) .* mu2; %...
                 %+ (x == L/2)   .* (mu1 + mu2)/2;
    
    % Function handler for rho
    Data.rho = @(x) (x >= 0   & x <L/2) .* rho1 ...
                  + (x >= L/2 & x <= L) .* rho2; %...
                  %+ (x == L/2)   .* (rho1 + rho2)/2;

    % Compute numerical derivative using forward finite difference dmu/dx
    dmu_dx = diff(Data.mu(Data.x)) ./ diff(Data.x);

    % Duplicating the last derivative for the last point (since diff reduces size)
    dmu_dx =  [dmu_dx; dmu_dx(end)];
    
    % d(mu * d(u_ex/dx))/dx
    double_deriv = diff(Data.mu(Data.x) .* Data.uex_x(Data.x) ) ./ diff(Data.x);

    % Duplicating the last derivative for the last point (since diff reduces size)
    double_deriv = [double_deriv ; double_deriv(end)]; 

    % Force vector
    Data.force = @(x) Data.rho(x).* Data.omega .* Data.uex(x) + double_deriv;

    % Dirichlet condition on the left part of the domain
    Data.gD = Data.uex(0) ;

    % Robin on the right part of the domain
    Data.gR = Data.uex(L) - 1i * Data.uex_x(L) ./ (Data.rho(L) .* Data.omega); 

elseif strcmp(TestName,'TestHW1_3a')
    f0 = 3;
    Data.f0 = 3;

    % Wave speed
    Data.c = @(x) ones(size(x));

    % Density
    Data.rho = @(x) rho1* ones(size(x));

    % Tension (??)
    Data.mu = @(x)  Data.rho(x);

    Data.mu1 = rho1;
    Data.mu2 = rho1;

    % We're interested only in alpha(end), i.e. at x=L
    Data.alpha = Data.mu(L) .* Data.rho(L) ;

    % Force vector
    Data.force = @(x) zeros(size(x));

    % Dirichlet condition on the left part of the domain
    Data.gD = 2 * omega^2 / (sqrt(pi) * f0^2) * exp( -(omega/f0)^2 ) ;

    % Robin on the right part of the domain
    Data.gR = 0; 

elseif strcmp(TestName,'TestHW1_3b')

    f0 = 3;
    Data.f0 = 3;

    Data.c = @(x) (x >= 0  & x < L/2) .* 0.1 ...
                + (x >= L/2 & x <= L)   .* 2;

    Data.rho = @(x) rho1 * ones(size(x));

    Data.mu = @(x) Data.rho(x) .* (Data.c(x)).^2;

    mu_vals = unique(Data.mu(x), 'stable');

    Data.mu1 = mu_vals(1);
    Data.mu2 = mu_vals(2);

    % We're interested only in alpha(end), i.e. at x=L
    Data.alpha = Data.mu(L) .* Data.rho(L) * Data.omega ;

    % Force vector
    Data.force = @(x) zeros(size(x));

    % Dirichlet condition on the left part of the domain
    Data.gD = 2 * omega^2 / (sqrt(pi) * f0^2) * exp( -(omega/f0)^2 ) ;

    % Robin on the right part of the domain
    Data.gR = 0; 
end



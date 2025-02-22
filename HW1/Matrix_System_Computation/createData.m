function [Data] = createData(TestName, N_pts, omega, mu_vals, rho_vals )
% ========================================================================
%   OUTPUT : Structure of the Galerkin formulation with problem's functions
%            discretization

%   INPUTS : 
%       - TestName  --> Identification for questions of HW1
%       - N_pts     --> Number of nodes in the mesh
%       - omega     --> Scalar value for omega parameter
%       - mu_vals   --> array of (1,2) for the mu values on each half of
%                       the domain
%       - rho_vals  --> array of (1,2) for the rho values on each half of
%                       the domain
% ========================================================================

if strcmp(TestName,'TestHW1_2a')
    L = 1; % Domain end

    Data.name = TestName;
    Data.domain = [0,L];
    Data.L = L;
    Data.boundary = 'DR';

    % Exact solution for error analysis
    Data.uex = @(x)     sin(2*pi*x);
    Data.uex_x = @(x)   2*pi .* cos(2*pi*x);
    Data.uex_xx = @(x)  - (2*pi)^2 .* Data.uex(x);

    Data.omega = omega;

    Data.mu_vals = mu_vals;
    Data.rho_vals = rho_vals;

    
    % Parameters and external forces, simple case with constant terms rho,mu
    Data.mu = @(x) mu_vals(1).*ones(size(x));
    Data.rho = @(x) rho_vals(1).*ones(size(x));
    
    x = linspace(0,L, N_pts);
    Data.x = x.';

    % Compute mu(x) values
    mu_values = Data.mu(x);

    % Compute numerical derivative using forward finite difference dmu/dx
    dmu_dx = diff(mu_values) ./ diff(x);

    % Duplicating the last derivative for the last point (since diff reduces size)
    Data.mu_x =  [dmu_dx, dmu_dx(end)];  % NaN for the last value (since diff reduces size)
    
    % d(mu * d(u_ex/dx))/dx
    d_mu_u_ex_x = @(x) diff(Data.mu(x) .* Data.uex_x(x) ) ./ diff(x);

    % Duplicating the last derivative for the last point (since diff reduces size)
    % Data.d_mu_u_ex_x  =  [d_mu_u_ex_x(x) , d_mu_u_ex_x(end)];
   
    % We're interested only in alpha(end), i.e. at x=L
    Data.alpha = Data.mu(L) .* Data.rho(L) * Data.omega ; 

    % Force vector
    Data.force = @(x) Data.rho(x) * Data.omega .* Data.uex(x) + Data.mu(x) .* Data.uex_xx(x);

    % Dirichlet condition on the left part of the domain
    Data.gD = Data.uex(0) ;

    % Robin on the right part of the domain
    Data.gR = Data.uex(L) - 1i * Data.uex_x(L) ./ (Data.rho(L) .* Data.omega); 

elseif strcmp(TestName,'TestHW1_2b')
    L = 1; % Domain end

    Data.name = TestName;
    Data.domain = [0,L];
    Data.L = L;Data.omega = omega;
    Data.boundary = 'DR';

    [Data.uex , Data.uex_x, Data.uex_xx] = setExactSolution(Data);

    Data.omega = omega;

    % Parameters and external forces, simple case with piece-wise constant terms rho,mu
    mu1 = mu_vals(1);    % Value of mu in [0, L/2]
    mu2 = mu_vals(2);    % Value of mu in [L/2, L]
    Data.mu1 = mu1;
    Data.mu2 = mu2;

    % Define the function handle
    Data.mu = @(x) (x >= 0  & x < L/2) .* mu1 ...
                 + (x >= L/2 & x <= L)   .* mu2;

    rho1 = rho_vals(1);       % Value of rho in [0, L/2]
    rho2 = rho_vals(2);     % Value of rho in [L/2, L]
    Data.rho1 = rho1;
    Data.rho2 = rho2;

    % Define the function handle
    Data.rho = @(x) (x >= 0  & x <= L/2) .* rho1 ...
                  + (x > L/2 & x <= L)   .* rho2;


    x = linspace(0,L, N_pts).';
    Data.x = x;

    % Compute mu(x) values
    mu_values = Data.mu(x);
    
    % Compute numerical derivative using forward finite difference dmu/dx
    dmu_dx = diff(mu_values) ./ diff(x);

    % Duplicating the last derivative for the last point (since diff reduces size)
    Data.mu_x =  [dmu_dx; dmu_dx(end)];
    
    % d(mu * d(u_ex/dx))/dx
    d_mu_u_ex_x = diff(Data.mu(x) .* Data.uex_x(x) ) ./ diff(x);

    % Duplicating the last derivative for the last point (since diff reduces size)
    Data.d_mu_u_ex_x  =  [d_mu_u_ex_x ; d_mu_u_ex_x(end)]; 

    % We're interested only in alpha(end), i.e. at x=L
    Data.alpha = Data.mu(L) .* Data.rho(L) * Data.omega ; 

    % Force vector
    Data.force = @(x) Data.rho(x).* Data.omega .* Data.uex(x) + Data.d_mu_u_ex_x;

    % Dirichlet condition on the left part of the domain
    Data.gD = Data.uex(0) ;

    % Robin on the right part of the domain
    Data.gR = Data.uex(L) - 1i * Data.uex_x(L) ./ (Data.rho(L) .* Data.omega); 

elseif strcmp(TestName,'TestHW1_3a')

    L = 5; % Domain end
    f0 = 3;
    Data.f0 = 3;

    Data.name = TestName;
    Data.domain = [0,L];
    Data.L = L;
    Data.boundary = 'DR';

    x = linspace(0,L, N_pts).';
    Data.x = x;

    Data.omega = omega;

    rho1 = rho_vals(1);       % Value of rho in [0, L]
 
    % Wave speed
    Data.c = @(x) ones(size(x));
    Data.rho = @(x) rho1* ones(size(x));

    Data.mu = @(x)  Data.rho(x);

    % Compute mu(x) values
    mu_values = Data.mu(x);

    % Compute numerical derivative using forward finite difference dmu/dx
    dmu_dx = diff(mu_values) ./ diff(x);

    % Putting the last point equal its predecessor
    Data.mu_x =  [dmu_dx; dmu_dx(end)];  % NaN for the last value (since diff reduces size)
    
    % d(mu * d(u_ex/dx))/dx
    % d_mu_u_ex_x = diff(Data.mu(x) .* Data.uex_x(x) ) ./ diff(x);
    % 
    % Data.d_mu_u_ex_x  =  [d_mu_u_ex_x , d_mu_u_ex_x(end)];  % NaN for the last value (since diff reduces size)

    % We're interested only in alpha(end), i.e. at x=L
    Data.alpha = Data.mu(L) .* Data.rho(L) * Data.omega ;

    % Force vector
    Data.force = @(x) zeros(size(x));


    % Dirichlet condition on the left part of the domain
    Data.gD = 2 * omega^2 / (sqrt(pi) * f0^2) * exp( -(omega/f0)^2 ) ;

    % Robin on the right part of the domain
    Data.gR = 0; 

elseif strcmp(TestName,'TestHW1_3b')
    L = 5; % Domain end
    f0 = 3;
    Data.f0 = 3;

    Data.name = TestName;
    Data.domain = [0,L];
    Data.L = L;
    Data.boundary = 'DR';

    x = linspace(0,L, N_pts).';
    Data.x = x;

    Data.omega = omega;
    rho1 = rho_vals(1);       % Value of rho in [0, L]
    Data.c = @(x) (x >= 0 & x <= L/2) .* 0.1 + (x > L/2 & x <= L) .* 2;
    Data.rho = @(x) rho1 * ones(size(x));

    Data.mu = @(x)  (Data.c(x)).^2 .*Data.rho(x);

    % Compute mu(x) values
    mu_values = Data.mu(x);

    % Compute numerical derivative using forward finite difference dmu/dx
    dmu_dx = diff(mu_values) ./ diff(x);

    % Putting the last point equal its predecessor
    Data.mu_x =  [dmu_dx; dmu_dx(end)];

    % d(mu * d(u_ex/dx))/dx
    % d_mu_u_ex_x = diff(Data.mu(x) .* Data.uex_x(x) ) ./ diff(x);
    % 
    % Data.d_mu_u_ex_x  =  [d_mu_u_ex_x , d_mu_u_ex_x(end)];  % NaN for the last value (since diff reduces size)

    % We're interested only in alpha(end), i.e. at x=L
    Data.alpha = Data.mu(L) .* Data.rho(L) * Data.omega ;

    % Force vector
    Data.force = @(x) zeros(size(x));

    % Dirichlet condition on the left part of the domain
    Data.gD = 2 * omega^2 / (sqrt(pi) * f0^2) * exp( -(omega/f0)^2 ) ;

    % Robin on the right part of the domain
    Data.gR = 0; 
end



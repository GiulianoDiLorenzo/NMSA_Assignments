%% Data for the following PDE
%   rho(x)* w * u(x)  + d(mu*du/dx)/dx = f in (0,L)
%   Dirichlet condition: u(0) = g_D
%   Robin condition : u(L) - j/(rho(x)*w) * du(x=L)/dx = g_A

function [Data] = createData(TestName, numEl, omega)

if strcmp(TestName,'TestHW1_2a')
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
    Data.mu = @(x) 5.*ones(size(x));
    Data.rho = @(x) 1.*ones(size(x));
    Data.omega = omega;

    x = linspace(0,L, numEl);
    Data.x = x.';

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


    % Non constant profiles
    % Data.force = @(x) Data.rho(x) .* Data.omega .* Data.uex(x) + Data.d_mu_u_ex
    % Data.d_mu_u_ex()x
    % (Data.rho(x)*Data.omega - Data.mu(x)*(2*pi)^2) .* Data.uex;

    % Dirichlet condition on the left part of the domain

    Data.gD = 0 ;

    % Robin on the right part of the domain
    Data.gR = Data.uex(L) - 1i * Data.uex_x(L) ./ (Data.rho(L) .* Data.omega); 

elseif strcmp(TestName,'TestHW1_2b')
    % INPUT DATA HERE
    %   rho(x)* w * u(x)  + d(mu*du/dx)/dx = f in (0,L)
    %   Dirichlet condition: u(0) = g_D
    %   Robin condition : u(L) - j/(rho(x)*w) * du(x=L)/dx = g_A

    L = 1; % Domain end

    Data.name = TestName;
    Data.domain = [0,L];
    Data.L = L;
    Data.boundary = 'DR';

    [Data.uex , Data.uex_x, Data.uex_xx] = setExactSolution(Data);

    % Exact solution for error analysis
    % Data.uex = @(x)  (x >= 0 & x <= L/2) .*  sin(2*pi*x) + (x > L/2 & x <= L) .*  sin(4*pi*x);
    % Data.uex_x = @(x) (x >= 0 & x <= L/2) .*  2*pi .* cos(2*pi*x) + (x > L/2 & x <= L) .*  4*pi .* cos(4*pi*x);
    % Data.uex_xx = @(x) (x >= 0 & x <= L/2) .*  (-(2*pi)^2) .* sin(2*pi*x) + (x > L/2 & x <= L) .* (-(4*pi)^2) .* sin(4*pi*x);

    % Data.uex = @(x) sin(2*pi*x);
    % Data.uex_x = @(x) 2*pi*cos(2*pi*x);
    % Data.uex_xx = @(x) -(2*pi)^2*uex(x);
    

    % Parameters and external forces, simple case with piece-wise constant terms rho,mu
    mu1 = 1;    % Value of mu in [0, L/2]
    mu2 = 10;    % Value of mu in [L/2, L]
    Data.mu1 = mu1;
    Data.mu2 = mu2;
    % Define the function handle
    Data.mu = @(x) (x >= 0 & x <= L/2) .* mu1 + (x > L/2 & x <= L) .* mu2;

    rho1 = 1;    % Value of mu in [0, L/2]
    rho2 = 150;    % Value of mu in [L/2, L]
    Data.rho1 = rho1;
    Data.rho2 = rho2;
    % Define the function handle
    Data.rho = @(x) (x >= 0 & x <= L/2) .* rho1 + (x > L/2 & x <= L) .* rho2;

    Data.omega = omega;

    x = linspace(0,L, numEl).';
    Data.x = x;

    % Compute mu(x) values
    mu_values = Data.mu(x);
    % Compute numerical derivative using forward finite difference dmu/dx
    dmu_dx = diff(mu_values) ./ diff(x);

    % Putting the last point equal its predecessor
    Data.mu_x =  [dmu_dx; dmu_dx(end)];  % NaN for the last value (since diff reduces size)
    
    % d(mu * d(u_ex/dx))/dx
    d_mu_u_ex_x = diff(Data.mu(x) .* Data.uex_x(x) ) ./ diff(x);
    Data.d_mu_u_ex_x  =  [d_mu_u_ex_x ; d_mu_u_ex_x(end)];  % NaN for the last value (since diff reduces size)

    Data.alpha = Data.mu(L) .* Data.rho(L) * Data.omega ; % We're interested only in alpha(end), i.e. at x=L
    
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

    Data.c = @(x) ones(size(x));
    Data.rho = @(x) 5* ones(size(x));

    x = linspace(0,L, numEl).';
    Data.x = x;

    Data.mu =@(x)  Data.rho(x);

    Data.omega = omega;

   
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


    Data.alpha = Data.mu(L) .* Data.rho(L) * Data.omega ; % We're interested only in alpha(end), i.e. at x=L
    
    Data.force = @(x) zeros(size(x));

    %  % Non constant profiles
    % Data.force = @(x) Data.rho(x) .* Data.omega .* Data.uex(x) + Data.d_mu_u_ex
    % Data.d_mu_u_ex()x
    % (Data.rho(x)*Data.omega - Data.mu(x)*(2*pi)^2) .* Data.uex;

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

    Data.c = @(x) (x >= 0 & x <= L/2) .* 0.1 + (x > L/2 & x <= L) .* 2;
    Data.rho = @(x) 5* ones(size(x));

    x = linspace(0,L, numEl).';
    Data.x = x;

    Data.mu =@(x)  (Data.c(x)).^2 .*Data.rho(x);

    Data.omega = omega;

   
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

    Data.alpha = Data.mu(L) .* Data.rho(L) * Data.omega ; % We're interested only in alpha(end), i.e. at x=L
    
    Data.force = @(x) zeros(size(x));

    %  % Non constant profiles
    % Data.force = @(x) Data.rho(x) .* Data.omega .* Data.uex(x) + Data.d_mu_u_ex
    % Data.d_mu_u_ex()x
    % (Data.rho(x)*Data.omega - Data.mu(x)*(2*pi)^2) .* Data.uex;

    % Dirichlet condition on the left part of the domain
    Data.gD = 2 * omega^2 / (sqrt(pi) * f0^2) * exp( -(omega/f0)^2 ) ;

    % Robin on the right part of the domain
    Data.gR = 0; 
end



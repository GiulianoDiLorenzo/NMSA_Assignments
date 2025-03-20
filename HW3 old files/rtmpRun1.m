function [Rho] = rtmpRun1(scenario, Mesh, f, rho_c)

% FIRST_ORDER_GODUNOV Solves the traffic flow equation using first-order Godunov scheme
%
% Inputs:
%       scenario: Structure containing initial and boundary conditions
%       density distribution
%       Mesh: structure with the details and characteristics of [0,L]x[0,T]
%   f handler for the flux function
%   rho_c critical density

% Output:
%       Rho: matrix with the space and time domain of the numerical
%       solution

    % Initialize solution
    Rho = zeros(Mesh.Nx+1, Mesh.Nt+1);
    Rho(:,1) = scenario.rho_0;        %Initial condition
   
    % Bouindary condition
    Rho(1,:) = Rho(1,1);
    
    fprintf('Running order 1.....\n');
    for n = 1:Mesh.Nt
        %Apply boundary conditions for the given time-step
        Rho = apply_BC(Rho, n, scenario, f);
        
        % Update interior point
        for i = 2:Mesh.Nx-1  % rho contains the value for x = {0 , L} and they stay constant
           
            % Calculating fluxes at cell interfaces
            fR = godunov_flux(Rho(i,n) , Rho(i+1,n), rho_c, f);
            fL = godunov_flux(Rho(i-1,n) , Rho(i,n), rho_c, f);

            % Update solution
            Rho(i,n+1) = Rho(i,n) - Mesh.dt/Mesh.dx * (fR - fL);
        end
    
        % Apply boundary conditions  based on scenario
        % Rho = apply_BC(Rho, n+1, scenario, f);

        if strcmp(scenario.name, 'Traffic jam')
            % Already computed
            Rho(end, n+1) = scenario.rho_R;

        elseif strcmp(scenario.name, 'Green light')
           % Free outflow boundary condition
           % Rho(end, n+1) = Rho(end-1, n+1);
           Rho(end, n+1) = Rho(end, n) - dt/dx * (fR - fL);

        elseif strcmp(scenario.name, 'Traffic flow')
            % Free outflow boundary condition
            Rho(end, n+1) = Rho(end-1, n+1);
        else
            error('Unrecognized scenario')


    end
end


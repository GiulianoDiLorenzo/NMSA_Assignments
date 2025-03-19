function rho = first_order_godunov(rho_0, rho_L, rho_R, rho_c, scenario, Nx, Nt, dx, dt, f)
% FIRST_ORDER_GODUNOV Solves the traffic flow equation using first-order Godunov scheme
%
% Inputs:
%   rho_0    - Initial density distribution
%   rho_L    - Left boundary density value
%   rho_R    - Right boundary density value
%   scenario - String specifying the scenario
%   Nx       - Number of grid points
%   Nt       - Number of time steps
%   dx       - Grid spacing
%   dt       - Time step size
%   f        - Flux function handle
%
% Output:
%   rho      - Solution array with density at all grid points and time steps

    % Initialize solution array
    rho = zeros(Nx, Nt);
    rho(:, 1) = rho_0;  % Set initial condition
    
    % Main time-stepping loop
    for n = 1:Nt-1
        % Apply boundary conditions for the current time step
        % rho = apply_boundary_conditions(rho, n, Nx, scenario, rho_L, rho_R);
        
        % Update interior points using first-order Godunov scheme
        for i = 2:Nx-1
            % Calculate fluxes at cell interfaces
            fR = godunov_flux(rho(i, n), rho(i+1, n), rho_c, f);
            fL = godunov_flux(rho(i-1, n), rho(i, n), rho_c, f);
            
            % Update solution
            rho(i, n+1) = rho(i, n) - dt/dx * (fR - fL);
        end
        
        % Apply boundary conditions to the new time step
        switch scenario
            case 'Traffic Jam'
                % Fixed boundary conditions
                rho(1, n+1) = rho_L;
                rho(Nx, n+1) = rho_R;
                
            case 'Green Light'
                % Fixed left boundary
                rho(1, n+1) = rho_L;
                
                % One-sided update for right boundary (outflow)
                fL = godunov_flux(rho(Nx-2, n), rho(Nx-1, n), f);
                fR = godunov_flux(rho(Nx-1, n), rho(Nx, n), f);
                rho(Nx, n+1) = rho(Nx, n) - dt/dx * (fR - fL);
                
            case 'Traffic Flow'
                % Fixed left boundary
                rho(1, n+1) = rho_L;
                
                % One-sided update for right boundary (outflow)
                fL = godunov_flux(rho(Nx-2, n), rho(Nx-1, n), f);
                fR = godunov_flux(rho(Nx-1, n), rho(Nx, n), f);
                rho(Nx, n+1) = rho(Nx, n) - dt/dx * (fR - fL);
        end
        
        % Ensure physical bounds
        rho(:, n+1) = max(0, min(1, rho(:, n+1)));
    end
end
function rho = new_apply_BC(rho, n, Nx, scenario, f)
% APPLY_BOUNDARY_CONDITIONS Sets the boundary values based on the scenario
%
% Inputs:
%   rho      - Density array at all grid points and time steps
%   n        - Current time step index
%   Nx       - Number of grid points
%   scenario - String specifying the scenario
%   rho_L    - Left boundary density value
%   rho_R    - Right boundary density value
%
% Output:
%   rho      - Updated density array with boundary conditions applied

    switch scenario.name
        case 'Traffic jam'
            % Fixed boundary conditions
            rho(1, n) = scenario.rho_L;
            rho(Nx, n) = scenario.rho_R;
            
        case 'Green light'
            % Fixed left boundary, free outflow right boundary
            rho(1, n) = scenario.rho_L;
            rho(Nx, n) = rho(Nx-1, n);
             
        case 'Traffic flow'
            % Fixed left boundary, free outflow right boundary
            rho(1, n) = scenario.rho_L;
            rho(Nx, n) = rho(Nx-1, n);  % Zero-gradient (free outflow)
            % fL = godunov_flux_1(rho(Nx-2, n), rho(Nx-1, n), f);
            % fR = godunov_flux_1(rho(Nx-1, n), rho(Nx, n), f);
            % rho(Nx, n+1) = rho(Nx, n) - dt/dx * (fR - fL);
    end
end
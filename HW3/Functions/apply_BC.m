function Rho = apply_BC(Rho, n, scenario, f)
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
        case 'Traffic Jam'
            % Fixed boundary conditions
            Rho(1, n) = scenario.rho_L;
            Rho(Nx, n) = scenario.rho_R;
            
        case 'Green Light'
            % Fixed left boundary
            Rho(1, n) = scenario.rho_L;

            % rho(end, n) = rho(end-1, n);  % Zero-gradient (free outflow)
            % 
            % % One-sided update for right boundary (outflow)
            % fL = godunov_flux(Rho(end-2, n), Rho(end-1, n), rho_c, f);
            % fR = godunov_flux(Rho(end-1, n), Rho(end, n), rho_c, f);
            Rho(end, n+1) = Rho(end, n) - dt/dx * (fR - fL);
            
        case 'Traffic Flow'
            % Fixed left boundary, free outflow right boundary
            Rho(1, n) = scenario.rho_L;
            Rho(end, n) = Rho(end-1, n);  
    end
end
function rho = apply_boundary_conditions(rho, n, Nx, scenario, rho_L, rho_R)
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

    switch scenario
        case 'Traffic Jam'
            % Fixed boundary conditions
            rho(1, n) = rho_L;
            rho(Nx, n) = rho_R;
            
        case 'Green Light'
            % Fixed left boundary, free outflow right boundary
            rho(1, n) = rho_L;
            rho(Nx, n) = rho(Nx-1, n);  % Zero-gradient (free outflow)
            
        case 'Traffic Flow'
            % Fixed left boundary, free outflow right boundary
            rho(1, n) = rho_L;
            rho(Nx, n) = rho(Nx-1, n);  % Zero-gradient (free outflow)
    end
end
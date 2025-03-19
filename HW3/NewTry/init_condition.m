function [rho_0, rho_L, rho_R] = init_condition(scenario, rho_max, x, Nx)
% INIT_CONDITION Sets up initial and boundary conditions based on scenario
%
% Inputs:
%   scenario - String specifying the scenario ('Traffic Jam', 'Green Light', or 'Traffic Flow')
%   rho_max  - Maximum density value
%   x        - Grid points array
%   Nx       - Number of grid points
%
% Outputs:
%   rho_0    - Initial density distribution
%   rho_L    - Left boundary density value
%   rho_R    - Right boundary density value

    switch scenario
        case 'Traffic jam'
            rho_L = 0.8 * rho_max;
            rho_R = rho_max;
            rho_0 = zeros(Nx, 1);
            rho_0(x < 0.5) = rho_L;
            rho_0(x >= 0.5) = rho_R;
            
        case 'Green light'
            rho_L = 0.3 * rho_max;
            rho_R = 0.6 * rho_max;
            rho_0 = zeros(Nx, 1);
            rho_0(x < 0.5) = rho_L;
            rho_0(x >= 0.5) = rho_R;
            
        case 'Traffic flow'
            rho_L = 0.8 * rho_max;
            rho_R = 0.5 * rho_max;
            rho_0 = zeros(Nx, 1);
            rho_0(x < 0.5) = rho_L;
            rho_0(x >= 0.5) = rho_R;
            
        otherwise
            error('Unknown scenario. Choose "Traffic Jam", "Green Light", or "Traffic Flow".');
    end
end
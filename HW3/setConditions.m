function [rho_L, rho_R, rho_x_0] = setConditions(scenario, rho_max, L, x)

if strcmp(scenario, 'jam')
    % ============= TRAFFIC JAM =============
    rho_L = 0.2* rho_max;
    rho_R = rho_max;
    rho_x_0 = (x <= L/2) .* rho_L + (x > L/2) .* rho_R;

elseif strcmp(scenario, 'green')
    % ============= GREEN LIGHT =============
    rho_L = 0.7* rho_max;
    rho_R = 0.2* rho_max;
    rho_x_0 = (x <= L/2) .* rho_L + (x > L/2) .* rho_R;

elseif strcmp(scenario, 'flow')
    % ============= TRAFFIC FLOW ============
    rho_L = rho_max;
    rho_R = 0.5* rho_max;
    rho_x_0 = (x <= L/2) .* rho_L + (x > L/2) .* rho_R;

else 
    error('Unrecognized scenario, possibilities: jam / green / flow');
end
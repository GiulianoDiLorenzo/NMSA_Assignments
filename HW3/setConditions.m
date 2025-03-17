function [cond] = setConditions(scenario, rho_max, L, x)
% ========================================================================
%   OUTPUT : 
%       - rho_L     -->  value of rho on the left part of the domain
%       - rho_R     --> value of rho on the right part of the domain
%       - rho_x_0   -->  array of size (Nx,1) containing the initial 
%                         condition: rho(x,0)
%
%   INPUTS : 
%       - scenario  --> str indicating the considered scenario, either
%                       ('jam', 'green' , 'flow')
%       - rho_max   --> Maximum value of rho(x,t)
%       - L         --> scalar, length of the domain
%       - x         --> array of size (Nx,1) spatial mesh
% ========================================================================

cond = struct();
    if strcmp(scenario, 'Traffic jam')
        % ============= TRAFFIC JAM =============
        cond.rho_L = 0.2* rho_max;
        cond.rho_R = rho_max;
        % rho_x_0 = (x <= L/2) .* rho_L + (x > L/2) .* rho_R;
    
    elseif strcmp(scenario, 'Green light')
        % ============= GREEN LIGHT =============
        cond.rho_L = 0.7* rho_max;
        cond.rho_R = 0.2* rho_max;
        % rho_x_0 = (x <= L/2) .* rho_L + (x > L/2) .* rho_R;
    
    elseif strcmp(scenario, 'Traffic flow')
        % ============= TRAFFIC FLOW ============
        cond.rho_L = rho_max;
        cond.rho_R = 0.5* rho_max;
        % rho_x_0 = (x <= L/2) .* rho_L + (x > L/2) .* rho_R;
    
    else 
        error('Unrecognized scenario, possibilities: Traffic jam / Green light / Traffic flow');
    
    
    end

    cond.rho_x_0 = (x <= L/2) .* cond.rho_L + (x > L/2) .* cond.rho_R;

end



%% Define Godunov flux function
function F = godunov_flux(rho_L, rho_R, f)
    % Compute the Godunov flux at the interface between rho_L and rho_R
    
    % Find value at which df(rho) = 0 (critical density)
    rho_c = 0.5;  % For the flux function f(rho) = u_max * (rho - rho^2/rho_max), rho_c = rho_max/2
    
    if rho_L <= rho_R
        % Rarefaction wave
        if rho_L <= rho_c && rho_c <= rho_R
            % Critical density is between rho_L and rho_R
            F = f(rho_c);
        else
            % Take minimum of the two fluxes
            F = min(f(rho_L), f(rho_R));
        end
    else
        % Shock wave
        F = max(f(rho_L), f(rho_R));
    end
end

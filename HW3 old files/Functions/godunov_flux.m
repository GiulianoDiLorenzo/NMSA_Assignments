%% Define Godunov flux function
function F = godunov_flux(rho_L, rho_R, rho_c, f)
    % Compute the Godunov flux at the interface between rho_L and rho_R
    
    
    
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

function F = godunov_flux(rho_L, rho_R, f, rho_c)
% GODUNOV_FLUX Computes the Godunov flux at the interface between two states
%
% Inputs:
%   rho_L - Density value on the left side of the interface
%   rho_R - Density value on the right side of the interface
%   f     - Flux function handle
%
% Output:
%   F     - Computed Godunov flux

    % Find value at which df(rho) = 0 (critical density)
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
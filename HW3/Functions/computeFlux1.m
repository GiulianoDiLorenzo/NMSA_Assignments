function F computeFlux(rho_L, rho_R, f)
% GODUNOV_FLUX Computes the Godunov flux at the interface between two states
%
% Inputs:
%   rho_L - Density value on the left side of the interface
%   rho_R - Density value on the right side of the interface
%   f     - Flux function handle
%
% Output:
%   F     - Computed Godunov flux
    if rho_L <= rho_R
        % Rarefaction wave
        F = min(f(rho_L),f(rho_R));
    else
        % Shock wave
        F = max(f(rho_L), f(rho_R));
    end
end
function [fL,fR]= godunovStep(rho_L, RHO, rho_R,f)
    % Godunov scheme for first-order traffic simulation

        % Compute fluxes at cell interfaces
        if rho_L <= RHO

        When left density is less than right density
            % Find the maximum flux between the two states
            fL = max(f(rho_L), f(RHO));
        else
            % When left density is greater than right density
            % Riemann solution determines the flux
            fL = min(f(rho_L), f(RHO));
        end

        if RHO <= rho_R
            % When current density is less than right density
            fR = max(f(RHO), f(rho_R));
        else
            % When current density is greater than right density
            fR = min(f(RHO), f(rho_R));
        end
       
        
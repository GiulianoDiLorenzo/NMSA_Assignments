function Rho = runOrder1Solution(scenario, Mesh, f, rho_c)

    % Initialize solution
    Rho = zeros(length(Mesh.x), length(Mesh.Nt));
    Rho(:,1) = scenario.rho_0;
    
    % Boundary condition
    
    fprintf('Running order 1.....\n');
    for n = 1:Mesh.Nt
        %Apply boundary conditions for the given time-step
        % Rho = apply_BC(Rho, n, scenario, f);
        
        % Update interior point
        for i = 2:Mesh.Nx  % rho contains the value for x = {0 , L} and they stay constant
           
            % Calculating fluxes at cell interfaces
            fR = godunov_flux(Rho(i,n) , Rho(i+1,n), rho_c, f);
            fL = godunov_flux(Rho(i-1,n) , Rho(i,n), rho_c, f);

            % Update solution
            Rho(i,n+1) = Rho(i,n) - Mesh.dt/Mesh.dx * (fR - fL);
        end
    
        % Apply boundary conditions  based on scenario
        % Rho = apply_BC(Rho, n+1, scenario, f);

        if strcmp(scenario.name, 'Traffic jam')
            % Already computed
            Rho(end, n+1) = scenario.rho_R;

        elseif strcmp(scenario.name, 'Green light')
           % Free outflow boundary condition
           Rho(end, n+1) = Rho(end-1, n+1);

           % One-sided update for right boundary (outflow)
            fL = godunov_flux(Rho(Nx-2, n), Rho(Nx-1, n), f);
            fR = godunov_flux(Rho(Nx-1, n), Rho(Nx, n), f);
            Rho(Nx, n+1) = Rho(Nx, n) - dt/dx * (fR - fL);

        elseif strcmp(scenario.name, 'Traffic flow')
            % Free outflow boundary condition
            Rho(end, n+1) = Rho(end-1, n+1);

            % One-sided update for right boundary (outflow)
                fL = godunov_flux(Rho(Nx-2, n), Rho(Nx-1, n), f);
                fR = godunov_flux(Rho(Nx-1, n), Rho(Nx, n), f);
                Rho(Nx, n+1) = Rho(Nx, n) - dt/dx * (fR - fL);
        else
            error('Unrecognized scenario')
        end

    end
end


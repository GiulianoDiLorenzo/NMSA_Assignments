function rho = GodunovOrder1(scenario, Mesh, f, rho_c)
% FIRST_ORDER_GODUNOV Solves the traffic flow equation using first-order Godunov scheme
%
% Inputs:
%   rho_0    - Initial density distribution
%   rho_L    - Left boundary density value
%   rho_R    - Right boundary density value
%   scenario - String specifying the scenario
%   Nx       - Number of grid points
%   Nt       - Number of time steps
%   dx       - Grid spacing
%   dt       - Time step size
%   f        - Flux function handle
%
% Output:
%   rho      - Solution array with density at all grid points and time steps


    Nx = Mesh.Nx;
    Nt = Mesh.Nt;
    dx = Mesh.dx;
    dt = Mesh.dt;

    rho_0 = scenario.rho_0;
    rho_L = scenario.rho_L;
    rho_R = scenario.rho_R;

    % Initialize solution array
    rho = zeros(Nx, Nt+1);
    rho(:, 1) = rho_0;  % Set initial condition


    rho_ext = zeros(Nx+2, Nt+1);
    rho_ext(1, :) = rho(1,:); 
    rho_ext(2:end-1,:) = rho;
    rho_ext(end, :) = rho(end,:);  

    
    % Main time-stepping loop
    for n = 1:Nt
        % Apply boundary conditions for the current time step
        % rho = new_apply_BC(rho, n, Nx, scenario, rho_L, rho_R, f);
        
        % Update interior points using first-order Godunov scheme
        for i = 2:Nx+1
            % Calculate fluxes at cell interfaces
            % fR = computeFlux1(rho(i, n), rho(i+1, n), f);
            % fL = computeFlux1(rho(i-1, n), rho(i, n), f);
            % fR = computeFlux1(rho_ext(i, n), rho_ext(i+1, n), f);
            % fL = computeFlux1(rho_ext(i-1, n), rho_ext(i, n), f);
            rL = rho_ext(i-1,n);
            rH = rho_ext(i,n);
            rR = rho_ext(i+1,n);

            fL = godunov_flux(rL, rH, f, rho_c);
            fR = godunov_flux(rH,rR, f, rho_c);
            % Update solution
            % rho(i, n+1) = rho(i, n) - dt/dx * (fR - fL);
            rho_ext(i, n+1) = rho_ext(i, n) - dt/dx * (fR - fL);
        end
        
        % rho = applyBC(rho, n+1, Mesh, scenario, f);
        rho_ext = applyBC(rho_ext, n+1, Mesh, scenario, f);
        
        rho_ext(:, n+1) = max(0, min(1, rho_ext(:, n+1)));
    end

    rho = rho_ext(2:end-1,:);
end
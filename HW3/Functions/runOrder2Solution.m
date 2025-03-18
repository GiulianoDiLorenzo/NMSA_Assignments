function [rho] = runOrder2Solution(scenario, u_max, rho_max, Mesh, limiter_type)
    
    Nx = Mesh.Nx;
    Nt = Mesh.Nt;
    dx = Mesh.dx;
    dt = Mesh.dt;

    % Initialize solution array
    rho = zeros(Nx, Nt+1);
    rho(:, 1) = scenario.rho_0;  % Set initial condition

    % Define flux function and its derivatives
    f = @(rho) u_max * (rho - rho.^2/rho_max);
    df = @(rho) u_max * (1 - 2*rho/rho_max);
    ddf = -2 * u_max / rho_max;

    fprintf('Running order 2.....\n');
    % Main time-stepping loop with second-order Godunov update
    for n = 1:Nt
        % Apply boundary conditions appropriate for the scenario
        switch scenario.name
            case 'Traffic Jam'
                % Fixed boundary conditions
                rho(1, n) = scenario.rho_L;
                rho(Nx, n) = scenario.rho_R;
                
            case 'Green Light'
                % Fixed left boundary, free outflow right boundary
                rho(1, n) = scenario.rho_L;
                
            case 'Traffic Flow'
                % Fixed left boundary, free outflow right boundary
                rho(1, n) = scenario.rho_L;
        end
        
        % Step 1: Compute limited slopes for each cell
        slopes = zeros(Nx, 1);
        
        % Extended state vector with ghost cells for proper slope computation at boundaries
        rho_ext = [rho(1, n); rho(:, n); rho(Nx, n)];
        
        for i = 2:Nx+1  % Loop over interior cells (including ghost cells)
            % Compute forward and backward differences
            forward_diff = (rho_ext(i+1) - rho_ext(i)) / dx;
            backward_diff = (rho_ext(i) - rho_ext(i-1)) / dx;
            
            % Compute ratio of slopes for limiting
            if abs(backward_diff) < 1e-10
                if abs(forward_diff) < 1e-10
                    r = 1;  % Both slopes are effectively zero
                else
                    r = 100 * sign(forward_diff);  % Large r with sign of forward_diff
                end
            else
                r = forward_diff / backward_diff;
            end
            
            % Apply selected slope limiter
            phi = slope_limiter(r, limiter_type);
            
            % Compute limited slope
            slopes(i-1) = phi * backward_diff;
        end
        
        % Step 2: Data reconstruction at cell interfaces
        rho_interfaces_L = zeros(Nx+1, 1);
        rho_interfaces_R = zeros(Nx+1, 1);
        
        % Compute left and right states at all interfaces
        for i = 1:Nx-1
            % Right boundary of cell i (left state at interface i+1/2)
            rho_interfaces_L(i+1) = rho(i, n) + slopes(i) * dx/2;
            
            % Left boundary of cell i+1 (right state at interface i+1/2)
            rho_interfaces_R(i+1) = rho(i+1, n) - slopes(i+1) * dx/2;
        end
        
        % Handle domain boundaries
        % Left boundary
        rho_interfaces_L(1) = rho(1, n) - slopes(1) * dx/2;
        rho_interfaces_R(1) = rho(1, n);  % Fixed boundary condition
        
        % Right boundary
        rho_interfaces_L(Nx+1) = rho(Nx, n);  % Fixed boundary condition
        rho_interfaces_R(Nx+1) = rho(Nx, n) + slopes(Nx) * dx/2;
        
        % Enforce physically meaningful values (density between 0 and rho_max)
        rho_interfaces_L = max(0, min(rho_max, rho_interfaces_L));
        rho_interfaces_R = max(0, min(rho_max, rho_interfaces_R));
        
        % Step 3: Compute fluxes at all interfaces using Godunov solver
        fluxes = zeros(Nx+1, 1);
        for i = 1:Nx+1
            fluxes(i) = godunov_flux(rho_interfaces_L(i), rho_interfaces_R(i), f);
        end
        
        % Step 4: Update solution using conservative formula
        for i = 1:Nx
            rho(i, n+1) = rho(i, n) - dt/dx * (fluxes(i+1) - fluxes(i));
        end
        
        % Apply boundary conditions to the new time step
        switch scenario.name
            case 'Traffic Jam'
                % Already set at the beginning of the loop
                
            case 'Green Light'
                % Free outflow boundary condition
                rho(Mesh.Nx, n+1) = rho(Nx-1, n+1);
                
            case 'Traffic Flow'
                % Free outflow boundary condition
                rho(Mesh.Nx, n+1) = rho(Nx-1, n+1);
        end
        
        % Apply left boundary condition to the new time step
        rho(1, n+1) = scenario.rho_L;
        
        % Ensure physical bounds are maintained (due to numerical errors)
        rho(:, n+1) = max(0, min(rho_max, rho(:, n+1)));
    end


end
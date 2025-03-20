function rho = second_order_godunov(scenario, Mesh, f, limiter_type)
% SECOND_ORDER_GODUNOV Solves the traffic flow equation using second-order Godunov scheme
%
% Inputs:
%   rho_0       - Initial density distribution
%   rho_L       - Left boundary density value
%   rho_R       - Right boundary density value
%   scenario    - String specifying the scenario
%   Mesh        - Structure of the Mesh 
%   Nx          - Number of grid points
%   Nt          - Number of time steps
%   dx          - Grid spacing
%   dt          - Time step size
%   f           - Flux function handle
%   limiter_type - String specifying the slope limiter type
%
% Output:
%   rho         - Solution array with density at all grid points and time steps
    
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
    
    % Main time-stepping loop
    for n = 1:Nt
        % Apply boundary conditions for the current time step
        % rho = new_apply_BC(rho, n, Nx, scenario, rho_L, rho_R, f);
        
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
        
        % Enforce physically meaningful values
        rho_interfaces_L = max(0, min(1, rho_interfaces_L));
        rho_interfaces_R = max(0, min(1, rho_interfaces_R));
        
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
            case 'Traffic_jam'
                % Fixed boundary conditions
                rho(1, n+1) = rho_L;
                rho(Nx, n+1) = rho_R;
                
            case 'Green_light'
                % Fixed left boundary
                rho(1, n+1) = rho_L;
                
                % Update right boundary with the calculated value (outflow)
                % We keep the calculated value from the flux update
                
            case 'Traffic_flow'
                % Fixed left boundary
                rho(1, n+1) = rho_L;
                
                % Update right boundary with the calculated value (outflow)
                % We keep the calculated value from the flux update
        end
        
        % Ensure physical bounds
        rho(:, n+1) = max(0, min(1, rho(:, n+1)));
    end
end
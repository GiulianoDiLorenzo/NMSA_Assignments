%% Function to solve the traffic flow with exact flux expression
function [rho_godunov, rho_lw, cfl] = run_exact_flux(Mesh, rho_max, u_max, scenario)
    % Grid setup

    h = Mesh.dx;         % Spatial step size
    dt = Mesh.dt;        % Time step size
    x = Mesh.x;  % Spatial grid
    t = Mesh.t;  % Time grid
    
    % Initialize density
    rho = zeros(Mesh.Nx, Mesh.Nt+1);
    rho_godunov = first_order_godunov(rho_0, rho_L, rho_R, scenario, Mesh, f)
    % 
    % % Apply IC
    % rho(:,1) = scenario.rho_0;
    % % for i = 1:Mesh.Nx
    % %     if x(i) >= Mesh.L/2
    % %         rho(i, 1) = rho_R;
    % %     else
    % %         rho(i, 1) = rho_L;
    % %     end
    % % end
    % 
    % rho_R = scenario.rho_R;
    % rho_L = scenario.rho_L;
    % % Define flux function f(ρ) = ρ log(ρmax/ρ)
    % flux = @(rho) rho .* log(rho_max ./ rho);
    % 
    % % Define flux derivative (needed for characteristic speed)
    % flux_derivative = @(rho) log(rho_max ./ rho) - 1;
    % 
    % % Compute the solution using first-order Godunov method
    % rho_godunov = rho;
    % for n = 1:Mesh.Nt
    %     % Apply boundary conditions
    %     rho_godunov(1, n) = rho_L;    % Left boundary
    %     rho_godunov(Mesh.Nx, n) = rho_R;  % Right boundary
    % 
    %     % Calculate numerical fluxes
    %     F = zeros(Mesh.Nx, 1);
    %     for i = 1:Mesh.Nx-1
    %         % Calculate characteristic speed at interface
    %         if rho_godunov(i, n) == rho_godunov(i+1, n)
    %             a = flux_derivative(rho_godunov(i, n));
    %         else
    %             a = (flux(rho_godunov(i+1, n)) - flux(rho_godunov(i, n))) / (rho_godunov(i+1, n) - rho_godunov(i, n));
    %         end
    % 
    %         % Godunov flux
    %         if a >= 0
    %             F(i) = flux(rho_godunov(i, n));
    %         else
    %             F(i) = flux(rho_godunov(i+1, n));
    %         end
    %     end
    % 
    %     % Update solution
    %     for i = 2:Mesh.Nx
    %         rho_godunov(i, n+1) = rho_godunov(i, n) - dt/h * (F(i) - F(i-1));
    %     end
    % end
    % 
    % % Compute the solution using Lax-Wendroff scheme (second-order)
    % rho_lw = rho;
    % for n = 1:Mesh.Nt
    %     % Apply boundary conditions
    %     rho_lw(1, n) = rho_L;    % Left boundary
    %     rho_lw(Mesh.Nx, n) = rho_R;  % Right boundary
    % 
    %     % Calculate numerical fluxes
    %     F = zeros(Mesh.Nx, 1);
    %     for i = 1:Mesh.Nx-1
    %         % Linear reconstruction
    %         rho_minus = rho_lw(i, n) + 0.5 * (rho_lw(i+1, n) - rho_lw(i, n));
    %         rho_plus = rho_lw(i+1, n) - 0.5 * (rho_lw(i+1, n) - rho_lw(i, n));
    % 
    %         % Calculate characteristic speed at interface
    %         if rho_lw(i, n) == rho_lw(i+1, n)
    %             a = flux_derivative(rho_lw(i, n));
    %         else
    %             a = (flux(rho_lw(i+1, n)) - flux(rho_lw(i, n))) / (rho_lw(i+1, n) - rho_lw(i, n));
    %         end
    % 
    %         % Lax-Wendroff flux
    %         if a >= 0
    %             F(i) = flux(rho_minus);
    %         else
    %             F(i) = flux(rho_plus);
    %         end
    %     end
    % 
    %     % Update solution
    %     for i = 2:Mesh.Nx
    %         rho_lw(i, n+1) = rho_lw(i, n) - dt/h * (F(i) - F(i-1));
    %     end
    % end
    % 
    % % Check CFL condition
    % max_speed = max(abs(flux_derivative(rho(rho > 0))));
    % cfl = max_speed * dt / h;
    % fprintf('CFL number: %.4f\n', cfl);
    % if cfl > 1
    %     warning('CFL condition violated! Reduce time step or increase spatial step.');
    % end
end
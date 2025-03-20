function Claude_trial
    % Parameters
    L = 1;              % Domain length [km]
    T = 1;              % Simulation time [hours]
    M = 100;            % Number of spatial cells
    N = 500;            % Number of time steps
    
    % Grid setup
    h = L/M;            % Spatial step size
    dt = T/N;           % Time step size
    x = linspace(0, L, M+1);  % Cell boundaries
    x_center = 0.5*(x(1:end-1) + x(2:end));  % Cell centers
    
    % Traffic parameters
    rho_max = 1.0;      % Maximum density [vehicles/km]
    u_max = 1.0;        % Maximum velocity [km/hour]
    
    % Select case (1: Traffic Jam, 2: Green Light, 3: Traffic flow)
    caseNum = 3;
    
    % Initialize density based on selected case
    rho = initialize_density(x_center, rho_max, caseNum);
    
    % Time integration
    for n = 1:N
        % Apply boundary conditions
        rho = apply_boundary_conditions(rho, caseNum, rho_max);
        
        % Update density using second-order Lax-Wendroff
        rho = update_density_LW(rho, dt, h, u_max, rho_max);
        
        % Visualize every 10 steps
        if mod(n, 10) == 0
            visualize_solution(x_center, rho, n, dt, caseNum);
        end
    end
end

function rho = initialize_density(x, rho_max, caseNum)
    % Initialize density based on case number
    switch caseNum
        case 1  % Traffic Jam
            rho_L = 0.3 * rho_max;
            rho_R = rho_max;
            rho = rho_L * ones(size(x));
            rho(x >= 0.5) = rho_R;
            
        case 2  % Green Light
            rho_L = 0.7 * rho_max;
            rho_R = 0.3 * rho_max;
            rho = rho_L * ones(size(x));
            rho(x >= 0.5) = rho_R;
            
        case 3  % Traffic flow
            rho_L = rho_max;
            rho_R = 0.5 * rho_max;
            rho = rho_L * ones(size(x));
            rho(x >= 0.5) = rho_R;
            
        otherwise
            error('Invalid case number');
    end
end

function rho = apply_boundary_conditions(rho, caseNum, rho_max)
    % Apply appropriate boundary conditions based on case
    switch caseNum
        case 1  % Traffic Jam: outflow boundary conditions
            rho(1) = rho(2);
            rho(end) = rho(end-1);
            
        case 2  % Green Light: outflow boundary conditions
            rho(1) = rho(2);
            rho(end) = rho(end-1);
            
        case 3  % Traffic flow: periodic boundary conditions
            rho(1) = rho(end-1);
            rho(end) = rho(2);
            
        otherwise
            error('Invalid case number');
    end
end

function f = flux(rho, u_max, rho_max)
    % Compute the flux function: f(rho) = rho * u(rho)
    f = u_max * rho .* (1 - rho/rho_max);
end

function a = characteristic_speed(rho, u_max, rho_max)
    % Compute the characteristic speed: f'(rho)
    a = u_max * (1 - 2*rho/rho_max);
end

function rho_new = update_density_LW(rho, dt, h, u_max, rho_max)
    % Second-order Lax-Wendroff scheme with linear reconstruction
    
    % Number of cells
    M = length(rho);
    
    % Create extended arrays for boundary handling
    rho_ext = [rho(1), rho, rho(end)];
    
    % Initialize numerical fluxes
    F = zeros(1,M+1);
    
    % Compute numerical fluxes at cell interfaces
    for i = 1:M+1
        % Linear reconstruction to get left and right states
        rho_minus = rho_ext(i) + 0.5 * (rho_ext(i+1) - rho_ext(i));
        rho_plus = rho_ext(i+1) - 0.5 * (rho_ext(i+1) - rho_ext(i));
        
        % Ensure densities stay within physical bounds
        rho_minus = max(0, min(rho_max, rho_minus));
        rho_plus = max(0, min(rho_max, rho_plus));
        
        % Solve the Riemann problem (Godunov flux)
        F(i) = godunov_flux(rho_minus, rho_plus, u_max, rho_max);
    end
    
    % Update cell averages
    rho_new = rho - (dt/h) .* (F(2:end) - F(1:end-1));
    
    % Ensure solutions stay within physical bounds
    rho_new = max(0.*rho, min(rho_max, rho_new));
end

function F = godunov_flux(rho_l, rho_r, u_max, rho_max)
    % Compute the Godunov flux for the Riemann problem with left and right states
    
    % Critical density at which flux is maximum
    rho_critical = rho_max / 2;
    
    % Compute characteristic speeds
    a_l = characteristic_speed(rho_l, u_max, rho_max);
    a_r = characteristic_speed(rho_r, u_max, rho_max);
    
    % Godunov flux computation
    if rho_l <= rho_r
        % Rarefaction wave
        if rho_l <= rho_critical && rho_r >= rho_critical
            % Sonic point included, flux is maximum
            F = flux(rho_critical, u_max, rho_max);
        elseif rho_r <= rho_critical
            % Both states left of critical point, use right flux
            F = flux(rho_r, u_max, rho_max);
        else
            % Both states right of critical point, use left flux
            F = flux(rho_l, u_max, rho_max);
        end
    else
        % Shock wave
        a_shock = (flux(rho_r, u_max, rho_max) - flux(rho_l, u_max, rho_max)) / (rho_r - rho_l);
        
        if a_shock >= 0
            F = flux(rho_l, u_max, rho_max);
        else
            F = flux(rho_r, u_max, rho_max);
        end
    end
end

function visualize_solution(x, rho, n, dt, caseNum)
    % Create descriptive title based on case
    time = n * dt;
    
    switch caseNum
        case 1
            title_str = 'Traffic Jam';
        case 2
            title_str = 'Green Light';
        case 3
            title_str = 'Traffic Flow';
        otherwise
            title_str = 'Unknown Case';
    end
    
    % Plot the solution
    figure(1);
    plot(x, rho, 'b-', 'LineWidth', 2);
    grid on;
    xlabel('Position (km)');
    ylabel('Density (vehicles/km)');
    title(sprintf('%s: t = %.3f', title_str, time));
    ylim([0, 1.1]);
    
    % Add pause to visualize the evolution
    pause(0.05);
end
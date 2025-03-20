function [] = draw3DSolution(scenario_name, Rho, Mesh)


    figure();
    [X, T] = meshgrid(Mesh.x, Mesh.t); % Create space-time grid
    
    surf(X, T, Rho', 'EdgeColor', 'none'); % Transpose Rho to match dimensions
    
    colorbar; % Add color scale
    colormap(jet); % Use colormap for better visualization
    
    xlabel('$x$ [km]');
    ylabel('$t$ [s]');
    zlabel('$\rho$ (vehicles/km)');
    % zlim([-0.1, rho_max+0.1])
    title(sprintf('%s solution, dx = %.4f m, dt = %.4f s', scenario_name, Mesh.dx, Mesh.dt));
    
    % view(135, 30); % Adjust view angle for better visualization
    shading interp; % Smooth color transition
    grid on;

end
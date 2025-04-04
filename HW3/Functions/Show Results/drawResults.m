function [] = drawResults(results, fluxes, target, saveMe)
    Mesh        = results.Mesh;
    scenario    = results.scenario;

    [X, T] = meshgrid(Mesh.x, Mesh.t); % Create space-time grid

    %% Density plots
    if strcmp(target, 'density')

        rho_1st     = results.rho_1st;
        rho_2nd     = results.rho_2nd;
        rho_1st_ex  = results.rho_1st_ex;
        rho_2nd_ex  = results.rho_2nd_ex;
    
        rho_L = results.scenario.rho_L;
        rho_R = results.scenario.rho_R;
            
        figure();
        sgtitle(sprintf('%s solution for $\\rho(x,t)$\nWith $dx = %.2f$ km, $dt = %.3f$ s, $T = %d$ s', scenario.name, Mesh.dx, Mesh.dt, Mesh.T));
         
        subplot(2,2,1);
            surf(X, T, rho_1st.', 'EdgeColor', 'none'); % Transpose Rho to match dimensions
            title('1st order scheme with $f(\rho)= u_{\max} \left( \rho-\rho^2 / \rho _{max} \right)$');
            
            xlabel('$x$ [km]');
            ylabel('$t$ [s]');
            zlabel('$\rho$ (vehicles/km)');
            ylim([0, Mesh.T*1.01])
            colorbar; % Add color scale
            colormap(jet); % Use colormap for better visualization
            clim([min(rho_L, rho_R), max(rho_L, rho_R)]); % Set common color limits
            % Set common color limits
            % zlim([-0.1, rho_max+0.1])
               
            view(0,90); % Adjust view angle for better visualization
            shading interp; % Smooth color transition
            grid on;
        
        
        subplot(2,2,2);
            surf(X, T, rho_2nd.', 'EdgeColor', 'none'); % Transpose Rho to match dimensions
            title('2nd order scheme with $f(\rho)= u_{\max} \left( \rho-\rho^2 / \rho _{max} \right)$');
            colorbar; % Add color scale
            colormap(jet); % Use colormap for better visualization
            clim([min(rho_L, rho_R),max(rho_L, rho_R)]); % Set common color limits
            xlabel('$x$ [km]');
            ylabel('$t$ [s]');
            zlabel('$\rho$ (vehicles/km)');
            ylim([0, Mesh.T*1.01])
               
            view(0,90);  % Adjust view angle for better visualization
            shading interp; % Smooth color transition
            grid on;
        
        subplot(2,2,3);
            surf(X, T, rho_1st_ex.', 'EdgeColor', 'none'); % Transpose Rho to match dimensions
            title('1st order scheme with $f(\rho) = \rho \log(\rho_{max} / \rho)$');
            colorbar; % Add color scale
            colormap(jet); % Use colormap for better visualization
            clim([min(rho_L, rho_R),max(rho_L, rho_R)]); % Set common color limits
            xlabel('$x$ [km]');
            ylabel('$t$ [s]');
            zlabel('$\rho$ (vehicles/km)');
            ylim([0, Mesh.T*1.01])
               
            view(0,90);  % Adjust view angle for better visualization
            shading interp; % Smooth color transition
            grid on;
        
        subplot(2,2,4);
            surf(X, T, rho_2nd_ex.', 'EdgeColor', 'none'); % Transpose Rho to match dimensions
            title('2nd order scheme with $f(\rho) = \rho \log(\rho_{max} / \rho)$');
            colorbar; % Add color scale
            colormap(jet); % Use colormap for better visualization
            clim([min(rho_L, rho_R),max(rho_L, rho_R)]); % Set common color limits
            xlabel('$x$ [km]');
            ylabel('$t$ [s]');
            zlabel('$\rho$ (vehicles/km)');
            ylim([0, Mesh.T*1.01])
           
            view(0,90);  % Adjust view angle for better visualization
            shading interp; % Smooth color transition
            grid on;

    %% Flux plots
    elseif strcmp(target, 'flux')
        flux_1st    = results.flux_1st;
        flux_2nd    = results.flux_2nd;
        flux_1st_ex = results.flux_1st_ex;
        flux_2nd_ex = results.flux_2nd_ex;
    
        f1_max = fluxes.f1_max;
        f2_max = fluxes.f2_max;
    
        [X, T] = meshgrid(Mesh.x, Mesh.t); % Create space-time grid
            
        figure();
        sgtitle(sprintf('%s solution for $f(\\rho)$\nWith $dx = %.2f$ km, $dt = %.3f$ s, $T = %d$ s', scenario.name, Mesh.dx, Mesh.dt, Mesh.T));
        
        subplot(2,2,1);
            surf(X, T, flux_1st.', 'EdgeColor', 'none'); % Transpose Rho to match dimensions
            title('1st order scheme with $f(\rho)= u_{\max} \left( \rho-\rho^2 / \rho _{max} \right)$');
            
            xlabel('$x$ [km]');
            ylabel('$t$ [s]');
            zlabel('$\rho$ (vehicles/km)');
            ylim([0, Mesh.T*1.01])
            colorbar; % Add color scale
            colormap(jet); % Use colormap for better visualization
            clim([0,f1_max]); % Set common color limits
               
            view(0,90); % Adjust view angle for better visualization
            shading interp; % Smooth color transition
            grid on;
        
        subplot(2,2,2);
            surf(X, T, flux_2nd.', 'EdgeColor', 'none'); % Transpose Rho to match dimensions
            title('2nd order scheme with $f(\rho)= u_{\max} \left( \rho-\rho^2 / \rho _{max} \right)$');
            colorbar; % Add color scale
            colormap(jet); % Use colormap for better visualization
            clim([0,f1_max]); % Set common color limits
            xlabel('$x$ [km]');
            ylabel('$t$ [s]');
            zlabel('$\rho$ (vehicles/km)');
            ylim([0, Mesh.T*1.01])
               
            view(0,90);  % Adjust view angle for better visualization
            shading interp; % Smooth color transition
            grid on;
        
        subplot(2,2,3);
            surf(X, T, flux_1st_ex.', 'EdgeColor', 'none'); % Transpose Rho to match dimensions
            title('1st order scheme with $f(\rho) = \rho \log(\rho_{max} / \rho)$');
            colorbar; % Add color scale
            colormap(jet); % Use colormap for better visualization
            clim([0,f2_max]); % Set common color limits
            xlabel('$x$ [km]');
            ylabel('$t$ [s]');
            zlabel('$\rho$ (vehicles/km)');
            ylim([0, Mesh.T*1.01])
           
            view(0,90);  % Adjust view angle for better visualization
            shading interp; % Smooth color transition
            grid on;
        
        subplot(2,2,4);
            surf(X, T, flux_2nd_ex.', 'EdgeColor', 'none'); % Transpose Rho to match dimensions
            title('2nd order scheme with $f(\rho) = \rho \log(\rho_{max} / \rho)$');
            colorbar; % Add color scale
            colormap(jet); % Use colormap for better visualization
            clim([0,f2_max]); % Set common color limits
            xlabel('$x$ [km]');
            ylabel('$t$ [s]');
            zlabel('$\rho$ (vehicles/km)');
            ylim([0, Mesh.T*1.01])
               
            view(0,90);  % Adjust view angle for better visualization
            shading interp; % Smooth color transition
            grid on;

    else
        error('Unrecognized target, select density or flux');
    end

    if saveMe
        saveas(gcf, sprintf('Pictures/%s %s comp %ds.png', results.scenario.name, target, results.Mesh.T));
    end
end
function [fluxes] = setFluxes(rho_max, u_max, plotMe)

    fluxes = struct();

    rho = linspace(1e-2, rho_max, 1e3);

    f1 = @(rho) u_max * (rho - rho.^2/rho_max);
    df1 = @(rho) u_max * (1 - 2*rho/rho_max);
    % ddf1 = -2 * u_max / rho_max;

    rho_c_f1 = rho_max/2;
  
    
    f2 = @(rho) rho .* log10(rho_max./rho);
    df2 = @(rho) log10(rho_max./rho) - 1;
    % ddf_2 = @(rho) 1/rho;
    rho_c_f2 = rho_max/10;

    fluxes.u_max = u_max;
    fluxes.rho_max = rho_max;
    fluxes.f1 = f1;
    fluxes.f1_max = max(f1(rho));
    fluxes.df1 = df1;
    fluxes.f2 = f2;
    fluxes.f2_max = max(f2(rho));
    fluxes.df2 = df2;
    fluxes.rho_c_f1 = rho_c_f1;
    fluxes.rho_c_f2 = rho_c_f2;
    

    

    if plotMe
        figure();
        % Plots of the fluxes
        subplot(2,1,1);
            plot(rho, f1(rho), "Color",'b')
            hold on;
            plot(rho, f2(rho))
            legend('$f_1(\rho)$','$f_2(\rho)$', Location='best')
            xlabel('$\rho$ [car/km]');
            ylabel('$f$ [car/s]');
            title('Graph of the different fluxes')
            grid on;
        
            % Modify x-axis labels
            xt = xticks(); % Get current x-tick positions
            xt(end) = rho_max; % Ensure last tick is at rho_max
            xticklabels([string(xt(1:end-1)), '$\rho_{\max}$']); % Replace last tick with LaTeX
        
            yt = yticks(); % Get current x-tick positions
            yt(end) = max(f2(rho)); % Ensure last tick is at rho_max
            xticklabels([string(xt(1:end-1)), '$\rho_{\max}$']); % Replace last tick with LaTeX
        
            yticklabels([string(yt(1:end-1)), '$f_{ex}^{\max}$']); % Replace last tick with LaTeX
            
            set(gca, 'TickLabelInterpreter', 'latex'); % Ensure LaTeX rendering
    
        % Plots of the characteristic speed
        subplot(2,1,2)
            plot(rho, df1(rho), "Color",'b')
            hold on;
            plot(rho, df2(rho))
            
            scatter(rho_c_f1, df1(rho_c_f1), 50, 'r', 'filled') % Red marker for first flux
            scatter(rho_c_f2, df2(rho_c_f2), 50, 'g', 'filled') % Green marker for second flux
             % Add text annotations
            text(rho_c_f1, df1(rho_c_f1), '$\rho_{c1}$', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 12);
            text(rho_c_f2, df2(rho_c_f2), '$\rho_{c2}$', 'Interpreter', 'latex', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontSize', 12);
    
            legend('$df_1 / d\rho$','$df_2 / d\rho$', Location='best')
            xlabel('$\rho$ [car/km]');
            ylabel('$fc(\rho)$ [km/s]');
            title('Graph of the different characteristic speeds')
            grid on;
        
            % Modify x-axis labels
            xt = xticks(); % Get current x-tick positions
            xt(end) = rho_max; % Ensure last tick is at rho_max
            xticklabels([string(xt(1:end-1)), '$\rho_{\max}$']); % Replace last tick with LaTeX
            
            set(gca, 'TickLabelInterpreter', 'latex'); % Ensure LaTeX rendering
    end
    
end
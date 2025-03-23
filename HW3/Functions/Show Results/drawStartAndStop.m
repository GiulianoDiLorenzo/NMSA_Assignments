function [] = drawStartAndStop(results, fluxes, target, saveMe)
    Mesh        = results.Mesh;
    scenario    = results.scenario;

    rho_max = fluxes.rho_max;

    fps = 24;
    fps_scaled = round(fps* (-Mesh.T/18 + 19/18));
    num_frames = fps_scaled * Mesh.T;
    t = floor(linspace(1, results.Mesh.Nt+1, num_frames));
    
    %% Plotting Density
    if strcmp(target, 'density')
        rho_1st     = results.rho_1st;
        rho_2nd     = results.rho_2nd;
        rho_1st_ex  = results.rho_1st_ex;
        rho_2nd_ex  = results.rho_2nd_ex;

        figure();
        sgtitle(sprintf(['%s solution - graph of $\\rho(x,t)$\nWith $dx = %.2f$ km, ' ...
                         '$dt = %.3f$ s, $T = %d$ s'], scenario.name, ...
                         Mesh.dx, Mesh.dt, Mesh.T));
            subplot(2,1,1)
            p1 = plot(Mesh.x, rho_1st(:,1), 'LineWidth', 1.5, 'Color','b');
            grid on;
            hold on;
            p2 = plot(Mesh.x, rho_1st(:, round(end/2)), 'LineWidth', 1.5, 'Color','b', 'LineStyle','--');
            p3 = plot(Mesh.x, rho_1st(:, end), 'LineWidth', 1.5, 'Color','b', 'LineStyle',':');
    
            p4 = plot(Mesh.x, rho_2nd(:, 1), 'LineWidth', 1.5, 'Color','r');
            p5 = plot(Mesh.x, rho_2nd(:, round(end/2)), 'LineWidth', 1.5, 'Color','r', 'LineStyle','--');
            p6 = plot(Mesh.x, rho_2nd(:, end), 'LineWidth', 1.5, 'Color','r', 'LineStyle',':');
    
            xlabel('Position $x$ [km]');
            ylabel('Density $\rho$ [vehicles/km]');
            title('$f(\rho)= u_{\max} \left( \rho-\rho^2 / \rho _{max} \right)$');
    
            ylim([0, 1.1*rho_max]);
    
            legend([p1, p2, p3, p4, p5, p6], ...
                   {'1st order scheme, $t=0s$', ...
                    '1st order scheme, $t=T/2$ s', ...
                    '1st order scheme, $t=T$ s', ...
                    '2nd order scheme, $t=0s$', ...
                    '2nd order scheme, $t=T/2$ s', ...
                    '2nd order scheme, $t=T$ s'}, ...
                    'Location', 'best');

            subplot(2,1,2)
            p1 = plot(Mesh.x, rho_1st_ex(:,1), 'LineWidth', 1, 'Color','b');
            grid on;
            hold on;
            p2 = plot(Mesh.x, rho_1st_ex(:, round(end/2)), 'LineWidth', 1, 'Color','b', 'LineStyle','--');
            p3 = plot(Mesh.x, rho_1st_ex(:, end), 'LineWidth', 1, 'Color','b', 'LineStyle',':');
    
            p4 = plot(Mesh.x, rho_2nd_ex(:, 1), 'LineWidth', 1.5, 'Color','r');
            p5 = plot(Mesh.x, rho_2nd_ex(:, round(end/2)), 'LineWidth', 1.5, 'Color','r', 'LineStyle','--');
            p6 = plot(Mesh.x, rho_2nd_ex(:, end), 'LineWidth', 1.5, 'Color','r','LineStyle',':');
    
            xlabel('Position $x$ [km]');
            ylabel('Density $\rho$ [vehicles/km]');
            title ('$f(\rho) = \rho \log(\rho_{max} / \rho)$');
        
            ylim([0, 1.1*rho_max]);
    

     %% Plotting Flux        
    elseif strcmp(target,'flux')
        flux_1st    = results.flux_1st;
        flux_2nd    = results.flux_2nd;
        flux_1st_ex = results.flux_1st_ex;
        flux_2nd_ex = results.flux_2nd_ex;

        figure();
        sgtitle(sprintf(['%s solution - graph of $f(\\rho(x,t))$\nWith $dx = %.2f$ km, ' ...
                         '$dt = %.3f$ s, $T = %d$ s'], scenario.name, ...
                         Mesh.dx, Mesh.dt, Mesh.T));
            subplot(2,1,1)
            p1 = plot(Mesh.x, flux_1st(:,1), 'LineWidth', 1.5, 'Color','b');
            grid on;
            hold on;
            p2 = plot(Mesh.x, flux_1st(:, round(end/2)), 'LineWidth', 1.5, 'Color','b', 'LineStyle','--');
            p3 = plot(Mesh.x, flux_1st(:, end), 'LineWidth', 1.5, 'Color','b', 'LineStyle',':');
    
            p4 = plot(Mesh.x, flux_2nd(:, 1), 'LineWidth', 1.5, 'Color','r');
            p5 = plot(Mesh.x, flux_2nd(:, round(end/2)), 'LineWidth', 1.5, 'Color','r','LineStyle','--');
            p6 = plot(Mesh.x, flux_2nd(:, end), 'LineWidth', 1.5, 'Color','r', 'LineStyle',':');
    
            xlabel('Position $x$ [km]');
            ylabel('Flux $f(\rho(x,t))$ [vehicles/s]');
            title('$f(\rho)= u_{\max} \left( \rho-\rho^2 / \rho _{max} \right)$');
    
            ylim([0, 1.01*fluxes.f1_max]);

            % Add horizontal lines at f1_max
            yline(fluxes.f1_max, '--k', '$f1_{max}$', 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5, Interpreter='latex');
        
    
            legend([p1, p2, p3, p4, p5, p6], ...
                   {'1st order scheme, $t=0s$', ...
                    '1st order scheme, $t=T/2$ s', ...
                    '1st order scheme, $t=T$ s', ...
                    '2nd order scheme, $t=0s$', ...
                    '2nd order scheme, $t=T/2$ s', ...
                    '2nd order scheme, $t=T$ s'}, ...
                    'Location', 'best');

            subplot(2,1,2)
            p1 = plot(Mesh.x, flux_1st_ex(:,1), 'LineWidth', 1, 'Color','b');
            grid on;
            hold on;
            p2 = plot(Mesh.x, flux_1st_ex(:, round(end/2)), 'LineWidth', 1, 'Color','b', 'LineStyle','--');
            p3 = plot(Mesh.x, flux_1st_ex(:, end), 'LineWidth', 1, 'Color','b', 'LineStyle',':');
    
            p4 = plot(Mesh.x, flux_2nd_ex(:, 1), 'LineWidth', 1.5, 'Color','r');
            p5 = plot(Mesh.x, flux_2nd_ex(:, round(end/2)), 'LineWidth', 1.5, 'Color','r','LineStyle','--');
            p6 = plot(Mesh.x, flux_2nd_ex(:, end), 'LineWidth', 1.5, 'Color','r','LineStyle',':');
    
            xlabel('Position $x$ [km]');
            ylabel('Flux $f(\rho(x,t))$ [vehicles/s]');
            title ('$f(\rho) = \rho \log(\rho_{max} / \rho)$');

            % Add horizontal lines at f1_max
            yline(fluxes.f2_max, '--k', '$f2_{max}$', 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5, Interpreter='latex');
        
            ylim([0, 1.01*fluxes.f2_max]);
    
            % legend([p1, p2, p3, p4, p5, p6], ...
            %        {'1st order scheme, $t=0s$', ...
            %         '1st order scheme, $t=T/2$ s', ...
            %         '1st order scheme, $t=T$ s', ...
            %         '2nd order scheme, $t=0s$', ...
            %         '2nd order scheme, $t=T/2$ s', ...
            %         '2nd order scheme, $t=T$ s'}, ...
            %         'Location', 'best');
    
    else 
        error('Unrecognized target, select density or flux');
    end

    if saveMe
        saveas(gcf, sprintf('Pictures/%s %s start stop %ds.png', results.scenario.name, target, results.Mesh.T));
    end
    
    
    
end
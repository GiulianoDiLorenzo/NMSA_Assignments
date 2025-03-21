function [] = animateDensity(results, rho_max)

    Mesh        = results.Mesh;
    scenario    = results.scenario;
    rho_1st     = results.rho_1st;
    rho_2nd     = results.rho_2nd;
    rho_1st_ex  = results.rho_1st_ex;
    rho_2nd_ex  = results.rho_2nd_ex;

    fps = 24;
    fps_scaled = round(fps* (-Mesh.T/18 + 19/18));
    num_frames = fps_scaled * Mesh.T;
    t = floor(linspace(1, size(rho_1st, 2), num_frames));
    
    figure();
    for n = t
        p1 = plot(Mesh.x, rho_1st(:, n), 'LineWidth', 1.5, 'Color',[0.5, 0, 0.5]);
        grid on;
        hold on;
        p2 = plot(Mesh.x, rho_2nd(:, n), 'LineWidth', 1.5, 'Color','b');
        p3 = plot(Mesh.x, rho_1st_ex(:,n), 'LineWidth', 1.5, 'Color','r');
        p4 = plot(Mesh.x, rho_2nd_ex(:,n), 'LineWidth', 1.5, 'Color', [1, 0.5, 0]); % Orange color
    
    
        xlabel('Position $x$ [km]');
        ylabel('Density $\rho$ [vehicles/km]');
        title(sprintf(['%s scenario animation\n' ...
                     '$dx = %.3f$ km, $dt = %.3f$ s, $T = %d$ s\n' ...
                     '$t = %.3f$ s'], scenario.name, Mesh.dx, Mesh.dt, Mesh.T, (n-1)*Mesh.dt));
    
        % title(sprintf('1st order scheme - $t = %.3f$ s', (n-1)*Mesh.dt));
        ylim([0.99*min(scenario.rho_L, scenario.rho_R), 1.01*rho_max]);
    
        legend([p1, p2, p3, p4], {'1st order scheme with $f_1(\rho)$', ...
                                  '2nd order scheme with $f_1(\rho) \ $', ...
                                  '1st order scheme with $f_2(\rho)$', ...
                                  '2nd order scheme with $f_2(\rho) \ $'}, 'Location', 'best');
    
    
        drawnow;
    
        pause(1/24);
        hold off;
    end

end
% Helper function to create GIF animation
function create_gif_animation(rho_1st, rho_2nd, rho_1st_ex, rho_2nd_ex, ...
                             Mesh, scenario, rho_max, L, fps, duration, filename)
    % Set up figure for GIF
    fig = figure('Position', [100, 100, 800, 600]);
    
    % Calculate frame indices
    total_frames = fps * duration;
    frame_indices = floor(linspace(1, size(rho_1st, 2), total_frames));
    
    % Create GIF animation
    for i = 1:length(frame_indices)
        n = frame_indices(i);
        
        % Plot the different solutions
        p1 = plot(Mesh.x, rho_1st(:, n), 'LineWidth', 2, 'Color', [0.0, 0.6, 0.0]);
        grid on;
        hold on;
        p2 = plot(Mesh.x, rho_2nd(:, n), 'LineWidth', 2, 'Color', [0.0, 0.0, 0.8]);
        p3 = plot(Mesh.x, rho_1st_ex(:, n), 'LineWidth', 2, 'Color', [0.8, 0.0, 0.0]);
        p4 = plot(Mesh.x, rho_2nd_ex(:, n), 'LineWidth', 2, 'Color', [0.8, 0.6, 0.0]);
        
        % Customize plot appearance
        xlabel('Position $x$ [km]', 'Interpreter', 'latex', 'FontSize', 14);
        ylabel('Density $\rho$ [vehicles/km]', 'Interpreter', 'latex', 'FontSize', 14);
        title(sprintf('%s Scenario - $t = %.3f$ s', scenario.name, (n-1)*Mesh.dt), ...
              'Interpreter', 'latex', 'FontSize', 16);
        
        ylim([0, 1.1*rho_max]);
        xlim([0, L]);
        
        legend([p1, p2, p3, p4], ...
               {'1st order scheme', '2nd order scheme', ...
                '1st order with exact flux', '2nd order with exact flux'}, ...
               'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
        
        % Add parameters
        annotation('textbox', [0.02, 0.02, 0.4, 0.05], 'String', ...
                   sprintf('$dx = %.3f$ km, $dt = %.3f$ s', Mesh.dx, Mesh.dt), ...
                   'Interpreter', 'latex', 'FitBoxToText', 'on', 'LineStyle', 'none', ...
                   'FontSize', 12, 'BackgroundColor', [1 1 1 0.7]);
        
        % Improve overall appearance
        set(gca, 'TickLabelInterpreter', 'latex');
        set(gca, 'FontSize', 12);
        box on;
        
        % Capture frame for GIF
        frame = getframe(fig);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        
        % Write to GIF file
        if i == 1
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 1/fps);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1/fps);
        end
        
        hold off;
        
        % Display progress less frequently for GIF creation
        if mod(i, floor(total_frames/5)) == 0
            fprintf('GIF progress: %.0f%%\n', 100*i/total_frames);
        end
    end
    fprintf('GIF animation saved as "%s"\n', filename);
end
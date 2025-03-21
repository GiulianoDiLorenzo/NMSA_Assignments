function [] = makeDensityGIF(results, rho_max)

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
    
    % Create Pictures subfolder if it doesn't exist
    picturesFolder = fullfile('GIFs');
    if ~exist(picturesFolder, 'dir')
        mkdir(picturesFolder);
    end
    
    % Set file path and name
    filename = fullfile(picturesFolder, sprintf('%s density animation %ds.gif', scenario.name, Mesh.T));
    
    
    filename = sprintf('GIFs/%s density animation %ds.gif', scenario.name, Mesh.T);
    figure('Position', [100 100 800 600]);  % Set figure size for better quality
    
    % First iteration to set up the gif file
    n = t(1);
    
    
    p1 = plot(Mesh.x, rho_1st(:, n), 'LineWidth', 1.5, 'Color',[0.5, 0, 0.5]);
    grid on;
    hold on;
    p2 = plot(Mesh.x, rho_2nd(:, n), 'LineWidth', 1.5, 'Color','b');
    p3 = plot(Mesh.x, rho_1st_ex(:,n), 'LineWidth', 1.5, 'Color','r');
    p4 = plot(Mesh.x, rho_2nd_ex(:,n), 'LineWidth', 1.5, 'Color', [1, 0.5, 0]); 
    
    xlabel('Position $x$ [km]', 'Interpreter', 'latex');
    ylabel('Density $\rho$ [vehicles/km]', 'Interpreter', 'latex');
    title(sprintf(['%s scenario density animation\n' ...
                   '$dx = %.3f$ km, $dt = %.3f$ s, $T = %d$ s\n' ...
                   '$t = %.3f$ s'], scenario.name, Mesh.dx, Mesh.dt, Mesh.T, (n-1)*Mesh.dt), 'Interpreter', 'latex');
    
    ylim([0.99*min(scenario.rho_L, scenario.rho_R), 1.01*rho_max]);
    legend([p1, p2, p3, p4], {'1st order scheme with $f_1(\rho)$', ...
                              '2nd order scheme with $f_1(\rho) \ $', ...
                              '1st order scheme with $f_2(\rho)$', ...
                              '2nd order scheme with $f_2(\rho) \ $'}, 'Location', 'best');
    drawnow;
    
    
    % Capture the frame and save it to the gif file
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 1/fps_scaled);
    
    
    fprintf('Density GIF Progress : 0 %%');
    % Subsequent frames
    for i = 2:length(t)
        n = t(i);
        % fprintf('Animation progress %d  %', round(i/length(t))*100);
    
    
        % Display progress less frequently for GIF creation
        fprintf('\b\b\b%2d%%', round(100*i/length(t)));  % Overwrites the same line
        
        % Create the plot
        p1 = plot(Mesh.x, rho_1st(:, n), 'LineWidth', 1.5, 'Color',[0.5, 0, 0.5]);
        grid on;
        hold on;
        p2 = plot(Mesh.x, rho_2nd(:, n), 'LineWidth', 1.5, 'Color','b');
        p3 = plot(Mesh.x, rho_1st_ex(:,n), 'LineWidth', 1.5, 'Color','r');
        p4 = plot(Mesh.x, rho_2nd_ex(:,n), 'LineWidth', 1.5, 'Color', [1, 0.5, 0]); 
        xlabel('Position $x$ [km]', 'Interpreter', 'latex');
        ylabel('Density $\rho$ [vehicles/km]', 'Interpreter', 'latex');
        title(sprintf(['%s scenario animation\n' ...
                   '$dx = %.3f$ km, $dt = %.3f$ s, $T = %d$ s\n' ...
                   '$t = %.3f$ s'], scenario.name, Mesh.dx, Mesh.dt, Mesh.T, (n-1)*Mesh.dt), 'Interpreter', 'latex');
        ylim([0.99*min(scenario.rho_L, scenario.rho_R), 1.01*rho_max]);
        legend([p1, p2, p3, p4], {'1st order scheme with $f_1(\rho)$', ...
                                  '2nd order scheme with $f_1(\rho) \ $', ...
                                  '1st order scheme with $f_2(\rho)$', ...
                                  '2nd order scheme with $f_2(\rho) \ $'}, 'Location', 'best');
        drawnow;
        
        % Capture the frame and append it to the gif file
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1/fps_scaled);
        
        hold off;
    end
    fprintf('\n');
    fprintf('Animation saved as %s\n\n', filename);

end
function [] = drawFluxes(rho_max, u_max)

    f = @(rho) u_max * (rho - rho.^2/rho_max);
    % df = @(rho) u_max * (1 - 2*rho/rho_max);
    % ddf = -2 * u_max / rho_max;
    
    f_ex = @(rho) rho .* log10(rho_max./rho);
    % df_ex = @(rho) log10(rho/rho_max) + 1;
    % ddf_ex = @(rho) 1/rho;

    rho = linspace(1e-2, rho_max, 1e3);

    figure();
    plot(rho, f(rho), "Color",'b')
    hold on;
    plot(rho, f_ex(rho))
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
    yt(end) = max(f_ex(rho)); % Ensure last tick is at rho_max
    xticklabels([string(xt(1:end-1)), '$\rho_{\max}$']); % Replace last tick with LaTeX

    yticklabels([string(yt(1:end-1)), '$f_{ex}^{\max}$']); % Replace last tick with LaTeX
    
    set(gca, 'TickLabelInterpreter', 'latex'); % Ensure LaTeX rendering
    
end
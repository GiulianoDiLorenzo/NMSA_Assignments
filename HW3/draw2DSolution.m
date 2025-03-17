function [] = draw2DSolution(scenario, Rho, Mesh, rho_max)
    figure();
    grid on;
    xlabel('$x$ [m]');
    ylabel('$\rho (x)$ [car/km]');
    ylim([-0.1 1.2*rho_max]);
    
    for n = 1 : Mesh.Nt
        plot(Mesh.x , Rho(:,n));
        hold on;
    
        grid on;
        xlabel('$x$ [m]');
        ylabel('$\rho (x)$ [car/km]');
        ylim([-0.1, 1.2*rho_max]);
        title(sprintf('%s solution, dx = %.4f m, dt = %.4f s, t = %.4f s', scenario, Mesh.dx, Mesh.dt, (n-1)*Mesh.dt));
        pause(0.05);
        hold off;
    end
end
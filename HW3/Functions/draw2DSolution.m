function [] = draw2DSolution(scenario, Rho, Mesh, rho_max)
    figure();
    grid on;
    xlabel('$x$ [m]');
    ylabel('$\rho (x)$ [car/km]');
    legend('$\rho(x,0)$','$\rho(x,t)$');
    ylim([-0.1 1.2*rho_max]);

    t = round(linspace (1,Mesh.Nt,100));
    
    for n = 1:5:Mesh.Nt+1
        k=n;

        if k >1
            k = n-1;
        end

        plot(Mesh.x , Rho(:,1));

        hold on;
        % plot(Mesh.x, Rho(:,1), Color='r', LineStyle='-.')
        plot(Mesh.x , Rho(:,k));
        grid on;
        xlabel('$x$ [m]');
        ylabel('$\rho (x)$ [car/km]');
        ylim([-0.1, 1.2*rho_max]);
        title(sprintf('%s solution, dx = %.4f m, dt = %.4f s, t = %.4f s', scenario, Mesh.dx, Mesh.dt, k*Mesh.dt));
        pause(0.05);
        legend('$\rho(x,0)$', '$\rho(x,t)$');
        hold off;
        
    end
end
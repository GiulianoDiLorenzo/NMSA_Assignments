function [] = draw2DSolution(scenario, Rho, Mesh, rho_max, u_max)
    dt = Mesh.dt;
    figure();
    grid on;
    xlabel('$x$ [m]');
    ylabel('$\rho (x)$ [car/km]');
    legend('$\rho(x,0)$','$\rho(x,t)$');
    ylim([-0.1 1.2*rho_max]);
    

    t = round(linspace (1,Mesh.Nt,Mesh.Nt/10));
    
    for n = 1:length(t)
        p1 = plot(Mesh.x, Rho(:,1), Color='r', LineStyle='-.');
        hold on;
        grid on;
        p2 = plot(Mesh.x , Rho(:,t(n)));
        
        xlabel('$x$ [km]');
        ylabel('$\rho (x)$ [car/km]');
        ylim([-0.1, 1.2*rho_max]);
        % titleText = ['%s solution, $\rho_{max}= %.f[car/km]$ , $u_{max}= %.f[m/s]$, \n ' ...
        %                '$t =$ %.4f s', scenario.name, rho_max, u_max, t(n)*dt];
        % title(sprintf('%s solution, $\rho_{max}=%.1f$ [car/km], $u_{max}=%.1f$ s \n' ...
        %        '$t = %.4f$ s', scenario.name, rho_max, u_max, t(n)*Mesh.dt));
        % 

        legend([p1, p2],'$\rho(x,0)$', '$\rho(x,t)$','Location', 'northeast');

        drawnow;
        pause(0.1);
        hold off;
        
    end
end
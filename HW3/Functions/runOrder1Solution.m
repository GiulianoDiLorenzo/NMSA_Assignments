function [Rho,t_c] = runOrder1Solution(scenario, cond, rho_max, u_max, Mesh)

    Rho = zeros(Mesh.Nx, Mesh.Nt);

    rho_x_0 = cond.rho_x_0;

    Rho(:,1) = rho_x_0;
    
    f   =  @(rho) u_max * (rho - rho.^2/rho_max); 
    rho = linspace(0,rho_max,100);
    plot(rho)
    df  =  @(rho) u_max * (1 - 2*rho/rho_max);
    ddf =  - 2 * u_max/rho_max;
    
    d_rho = ( rho_x_0(2:end) - rho_x_0(1:end-1) ) / Mesh.dx;
    t_c = min(abs(d_rho.*ddf));
    fprintf('Solution valid up to time t_c = %.3f s \n', t_c);
    
    for n = 1:Mesh.Nt-1       % rho contains the value for t=0 and is already computed
        
        Rho(1, n+1) = Rho(1, n);
        

        for i = 2:Mesh.Nx-1  % rho contains the value for x = {0 , L} and they stay constant
           
            % Update rule for the right cell
            if Rho(i,n) <= Rho(i+1,n)
                % Rarefaction wave
                fR = min( f(Rho(i,n)) ,  f(Rho(i+1,n)) );
            else
                % Shock wave
                fR = max( f(Rho(i,n)) ,  f(Rho(i+1,n)) );
            end
    
            % Update rule for the left cell
            if Rho(i-1,n) <= Rho(i,n)
                % Rarefaction wave
                fL = min( f(Rho(i,n)) ,  f(Rho(i-1,n)) );
            else 
                % Shock wave
                fL = max( f(Rho(i,n)) ,  f(Rho(i-1,n)) );
            end
    
            Rho(i,n+1) = Rho(i,n) - Mesh.dt/Mesh.dx * (fR - fL);
        end
    
        % Apply boundary conditions  based on scenario
        if strcmp(scenario, 'Traffic jam')
            % Already computed
            Rho(end, n+1) = cond.rho_R;
        elseif strcmp(scenario, 'Green light')
           % Free outflow boundary condition
           Rho(end, n+1) = Rho(end-1, n+1);
        elseif strcmp(scenario, 'Traffic flow')
            % Free outflow boundary condition
            Rho(end, n+1) = Rho(end-1, n+1);
        else
            error('Unrecognized scenario');
        end


    end
end
function [Rho,t_c] = runOrder2Solution(scenario, cond, rho_max, u_max, Mesh)

    Rho = zeros(Mesh.Nx, Mesh.Nt);

    rho_x_0 = cond.rho_x_0;

    Rho(:,1) = rho_x_0;
    
    f   =  @(rho) u_max * (rho - rho.^2/rho_max); 
    df  =  @(rho) u_max * (1 - 2*rho/rho_max);
    ddf =  - 2 * u_max/rho_max;
    
    d_rho = ( rho_x_0(2:end) - rho_x_0(1:end-1) ) / Mesh.dx;
    t_c = min(abs(d_rho.*ddf));
    fprintf('Solution valid up to time t_c = %.3f s \n', t_c);
    
    for n = 1:Mesh.Nt-1       % rho contains the value for t=0 and is already computed
        
        Rho(1, n+1) = Rho(1, n);
        

        for i = 2:Mesh.Nx-1  % rho contains the value for x = {0 , L} and they stay constant
            

            rho_Rm = Rho(i,n) + 0.5 * (Rho(i+1,n) - Rho(i,n));
            rho_Rp = Rho(i,n) - 0.5 * (Rho(i+1,n) - Rho(i,n));


            rho_Lm = Rho(i,n) + 0.5 * (Rho(i,n) - Rho(i-1,n));
            rho_Lp = Rho(i,n) - 0.5 * (Rho(i,n) - Rho(i-1,n));

            rho_minus_ip = Rho(i,n) + 0.5 * (Rho(i+1,n) - Rho(i,n));  % ρ^n-_{i+1/2}
            rho_plus_ip = Rho(i+1,n) - 0.5 * (Rho(i+1,n) - Rho(i,n)); % ρ^n+_{i+1/2}
        
            rho_minus_im = Rho(i-1,n) + 0.5 * (Rho(i,n) - Rho(i-1,n)); % ρ^n-_{i-1/2}
            rho_plus_im = Rho(i,n) - 0.5 * (Rho(i,n) - Rho(i-1,n));    % ρ^n+_{i-1/2}
        
           
            % Update rule for the right cell
             if rho_Rm <= rho_Rp
                fR = min(f(rho_Rm), f(rho_Rp));
            else
                fR = max(f(rho_Rm), f(rho_Rp));
            end
    
            % Update rule for the left cell
            if rho_Lm <= rho_Lp
                fL = min(f(rho_Lm), f(rho_Lp));
            else
                fL = max(f(rho_Lm), f(rho_Lp));
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
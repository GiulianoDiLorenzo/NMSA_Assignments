function [scenario] = setScenario(sce, rho_max, rho_c, x, Nx)

    scenario = struct();

    if strcmp(sce, 'Traffic jam')
        scenario.name = "Traffic jam";
        scenario.rho_L = 0.2 * rho_max;
        scenario.rho_R = rho_max;        

    elseif strcmp(sce, 'Green light')
        scenario.name = "Green light";
        scenario.rho_L = 0.9 * rho_max;
        scenario.rho_R = 0.2 * rho_max;

    elseif strcmp(sce, 'Traffic flow')
        scenario.name = "Traffic flow";
        scenario.rho_L = rho_max;
        scenario.rho_R = rho_max/2;
    end
               
    scenario.rho_c= rho_c;

    scenario.rho_0 = zeros(Nx, 1);
    scenario.rho_0(x < 0.5) = scenario.rho_L;
    scenario.rho_0(x >= 0.5) = scenario.rho_R;

    if(scenario.rho_L < rho_c)
        fprintf('Information travels Right (-->) on the left part\n');
    elseif  (scenario.rho_L == rho_c)
        fprintf('Information stands still on the left part\n');
    else  
        fprintf('Information travels Left (<--) on the left part\n');
    end

    if (scenario.rho_R < rho_c)
        fprintf('Information travels Right (-->) on the right part\n\n')
    elseif  (scenario.rho_R == rho_c)
        fprintf('Information stands still on the right part\n\n');
    else  
        fprintf('Information travels Left (<--) on the right part\n\n');
    end

end
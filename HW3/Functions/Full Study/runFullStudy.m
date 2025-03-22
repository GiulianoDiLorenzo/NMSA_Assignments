function [results] = runFullStudy(L,T,Nx,Nt,sce, fluxes)

    u_max = fluxes.u_max;
    rho_max = fluxes.rho_max;
    rho_c_f1 = fluxes.rho_c_f1;
    rho_c_f2 = fluxes.rho_c_f2;
    f1 = fluxes.f1;
    f2 = fluxes.f2;
    

    %% Creating Mesh, and Scenario
    Mesh = createMesh(L, T, Nx, Nt);

    %% Check CFL condition for both schemes
    CFL = u_max * Mesh.dt / Mesh.dx;
    if CFL > 0.5
        error('CFL condition for 2nd order method violated! CFL = %.4f > 0.5. Reduce dt or increase dx.\n', CFL);
    else
        fprintf('CFL condition satisfied: %.4f < 0.5\n\n', CFL);
    end

     % Options: 'minmod', 'superbee', 'vanLeer', 'MC', 'none'
    limiter_type = 'MC';
    fprintf('Using %s slope limiter for 2nd order scheme\n\n', limiter_type);
    
    %% Run the simulation with f1
    fprintf('Simulating %s scenario with f1\n', sce);
    [scenario] = setScenario(sce, rho_max, rho_c_f1, Mesh.x, Mesh.Nx);
    rho_1st = GodunovOrder1(scenario, Mesh, f1,rho_c_f1);
    rho_2nd = GodunovOrder2(scenario, Mesh, f1, rho_c_f1, rho_max, limiter_type);
    
    %% Run the simulation with f2
    fprintf('Simulating %s scenario with f2\n', sce);
    [scenario] = setScenario(sce, rho_max, rho_c_f2, Mesh.x, Mesh.Nx);
    rho_1st_ex = GodunovOrder1(scenario, Mesh, f2, rho_c_f2);
    rho_2nd_ex = GodunovOrder2(scenario, Mesh, f2, rho_c_f2, rho_max,limiter_type);

    rho_1st = real(rho_1st);
    rho_2nd = real(rho_2nd);
    rho_1st_ex = real(rho_1st_ex);
    rho_2nd_ex = real(rho_2nd_ex);


    %% Extracting only the real part and storing in results
    results.Mesh        = Mesh;
    results.scenario    = scenario;
    results.rho_1st     = rho_1st;
    results.rho_2nd     = rho_2nd;
    results.rho_1st_ex  = rho_1st_ex;
    results.rho_2nd_ex  = rho_2nd_ex;

    results.flux_1st     = u_max .* rho_1st .* (1 - rho_1st/rho_max);
    results.flux_2nd     = u_max .* rho_2nd .* (1 - rho_2nd/rho_max);
    results.flux_1st_ex  = rho_1st_ex .* log10(rho_max./rho_1st_ex);
    results.flux_2nd_ex  = rho_2nd_ex .* log10(rho_max./rho_2nd_ex);
end
function [result] = runNumericalSolution(TestName, omega, N_pts, computeError)
    result = struct();
    % Generate data
    [Data] = createData(TestName, N_pts, omega);
    
    % Compute results
    [sol, mesh, data] = systemAssembly(Data, N_pts);

    % Compute error
    if computeError
        error = Data.uex(Data.x).' - sol;
        L2_err = compute_L2_Error(error, mesh.h);
        result.err = L2_err;  % Store error
    end
   
    
    % Store results in the datas structure
    result.n_pts = N_pts;  % Store nEl value
    result.sol = sol;  % Store solution
    result.mesh = mesh; % Store mesh structure
    result.data = data; % Store data structure
end
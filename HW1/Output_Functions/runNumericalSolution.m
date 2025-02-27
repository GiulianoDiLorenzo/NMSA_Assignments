function [result] = runNumericalSolution(TestName, L, omega, N_pts, compute_Err, mu_vals, rho_vals)
% ========================================================================
%   OUTPUT : Structure with all data : problem formulation and solution

%   INPUTS : 
%       - TestName  --> Identification for questions of HW1
%       - N_pts     --> Number of nodes in the mesh
%       - omega     --> Scalar value for omega parameter
%       - compute_Err   --> Boolean for L2 error computation (==True)
%       - mu_vals   --> array of (1,2) for the mu values on each half of
%                       the domain
%       - rho_vals  --> array of (1,2) for the rho values on each half of
%                       the domain
% ========================================================================
    
    result = struct();

    %==========================================================================
    % GENERATE DATA
    %==========================================================================
    % [Data] = createData(TestName, N_pts, omega, mu_vals, rho_vals);
    [Data] = createData(TestName, L, N_pts, omega, mu_vals, rho_vals);
   
    %==========================================================================
    % MESH GENERATION
    %==========================================================================
    [Mesh] = createMesh(Data,N_pts);

    %==========================================================================
    % BUILD BASIS FUNCTIONS
    %==========================================================================
    plotPhi = false;
    [Phi] = buildPhi(Mesh, plotPhi);
    
    %==========================================================================
    % BUILD MATRICES AND SOLVE SYSTEM
    %==========================================================================
    [M_star,F] = systemAssembly(Data, Mesh, Phi, N_pts);

    %==========================================================================
    % INVERTING THE SYSTEM TO GET THE SOLUTION
    %==========================================================================
    fprintf('Computing numerical solution U = inv(M*) * F* \n');
    sol = M_star\F;

    sol = real(sol); %Extracting only the real part of the solution
    

    %==========================================================================
    % STORE RELEVANT DATA IN THE OUTPUT STRUCTURE
    %==========================================================================
    result.n_pts = N_pts;  % Store nEl value
    result.sol = sol;  % Store solution
    result.mesh = Mesh; % Store mesh structure
    result.data = Data; % Store data structure

    %==========================================================================
    % COMPUTE ERROR IF NEEDED
    %==========================================================================
    if compute_Err
        error = Data.uex(Data.x) - sol;
        L2_err = compute_L2_Error(error, Mesh.h);
        result.err = L2_err;  % Store error
    end 

end
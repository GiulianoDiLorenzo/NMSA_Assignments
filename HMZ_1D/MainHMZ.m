function [Errors,Solutions,Femregion,Data] = MainHMZ(Data, nEl)
%%
%    INPUT:
%          Data    : (struct) Data struct
%          nEl     : (int)    Number of mesh elements
%
%    OUTPUT:
%          errors      : (struct) contains the computed errors
%          solutions   : (sparse) nodal values of the computed and exact
%                        solution
%          femregion   : (struct) finite element space
%
%          Data        : (struct)  Data struct
%

fprintf('============================================================\n')
fprintf(['Solving test ', Data.name, ' with ',num2str(nEl),' elements \n']);

%==========================================================================
% MESH GENERATION
%==========================================================================

[Region] = CreateMesh(Data,nEl);

%==========================================================================
% FINITE ELEMENT REGION
%==========================================================================

[Femregion] = CreateFemregion(Data,Region);

%==========================================================================
% BUILD FINITE ELEMENT MATRICES and RIGHT-HAND SIDE
%==========================================================================

[A_nbc,M_nbc] = Matrix1D(Data,Femregion);


%==========================================================================
% BUILD FINITE ELEMENTS RHS
%==========================================================================

[b_nbc] = Rhs1D(Data,Femregion);


%==========================================================================
% SOLVE LINEAR SYSTEM - COMPUTE BOUNDARY CONDITIONS
%==========================================================================

% Solve the reduced system
K_nbc =  - A_nbc + Data.omega^2/Data.vel^2 * M_nbc;

if (strcmp(Data.boundary,'NN'))
    p = K_nbc\b_nbc;
else
    [K,b,u_g] = BoundaryConditions(K_nbc,b_nbc,Femregion,Data);
    p = K\b;
    p = p + u_g;
end


[p] = Snapshot(Femregion, p, 'pressure');

% Computation of the velocity 
%n = length(p);
%vel(2:n-1,1) = (p(3:n,1)-p(1:n-2,1))/(2*Femregion.h);
%vel(1) = Data.gN1(Data.omega,Data.ro,Data.vel);
%vel(n) = Data.gN2(Data.omega,Data.ro,Data.vel);
% dis = 1/(Data.ro*Data.omega).^2 * vel;

% Plot of displacement
%[vel] = Snapshot(Femregion, vel, 'velocity');


% Impedance in x=0
% Z_imp = p(1)/Data.gN1(Data.omega,Data.ro,Data.vel);

%==========================================================================
% POST-PROCESSING OF THE SOLUTION
%==========================================================================

[Solutions] = PostProcessing(Data,Femregion,p);


%==========================================================================
% ERROR ANALYSIS
%==========================================================================
Errors = [];
if (Data.calc_errors)
    [Errors] = ComputeErrors(Data,Femregion,Solutions);
end



%==========================================================================
% EIGENVALUE ANALYSIS
%==========================================================================

if Data.eig_analysis
    % eigenmodes(resonant frequencies) and eingevectors
    if (strcmp(Data.boundary,'NN'))
        M = M_nbc; A = A_nbc;
    else
        [M,~, ~] = BoundaryConditions(M_nbc, b_nbc, Femregion, Data);
        [A,~, ~] = BoundaryConditions(A_nbc, b_nbc, Femregion, Data);
    end
    
    
    % eigenvalues D -- eigenvector V
    % A x = k^2 M x
    [V,D] = eigs(A,M,6,0.1);
    
    k2 = diag(D);
    omega_h = sqrt(k2*Data.vel^2);
    modes = omega_h/(2*pi);
    kex = [0:5]'*pi/Data.domain(2);
    modes_ex = sqrt(kex.^2*Data.vel^2)/(2*pi);
    
    %% plot of 6 eigenvectors
    if Data.plot_eigvct
        for i = 1 : 6
            v = V(:,i);
            
            strTitle = ['eigenvector n ', num2str(i)];
            [v] = Snapshot(Femregion, v, strTitle);
        end
    end
    
    disp([modes modes_ex]);
    
end


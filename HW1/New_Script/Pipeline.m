function [Errors,Solutions,FemRegion,Data] = Pipeline(Data, nEl)
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

[FemRegion] = CreateFemregion(Data,Region);

%==========================================================================
% BUILD FINITE ELEMENT MATRICES and RIGHT-HAND SIDE
%==========================================================================

% [Phi] = buildPhi(Region, true);
% % Mass matrix assembly;
% M = buildM(Data, Femregion, Phi);
% 
% % Stiffness matrix assembly;
% A = buildA(Data, Femregion, Phi);
% 
% M_star = M - A;
% 
% % Computing load vector F = F - M_star * R + D
% [F,R] = buildF(Data, Femregion, M_star);

plotPhi = false;
[Phi] = buildPhi(Region, plotPhi);

% % Mass matrix assembly;
fprintf('Building Matrices and Vectors ...\n');
[M] = buildM(Data, FemRegion, Phi);

% Stiffness matrix assembly;
% fprintf('Building A ...');
[A] = buildA(Data, FemRegion, Phi);

% Boundary conditions matrix assembly
% fprintf('Building R* ...\n');
[R_star] = buildR_star(Data,FemRegion);

% fprintf('Building M* ...\n');
M_star = Data.omega * M - A - R_star;

% Excitation force vector
% fprintf('Building F...\n');
F = Data.force(Data.x)';
% Boundary conditions vector 
[R_bc] = buildR_bc(Data,FemRegion);
F = F + R_bc;

%==========================================================================
% Invert system for solution and extract u;
%==========================================================================
fprintf('Computing numerical solution...\n\n');
u_num = (M_star)\F;
u_num = u_num ./ ( max(abs(real(u_num))) );

%==========================================================================
% Gathering the exact and numerical solutions in 1 structure
%==========================================================================

uex_comp = Data.uex(Data.x)';

[Solutions] = struct('u_ex',uex_comp,'u_num',real(u_num));



figure();
sol_num = real(Solutions.u_num);

% sol_num = abs(sol_num).* exp(1i * angle(sol_num));


plot(Data.x, Solutions.u_ex);
hold on 
plot(Data.x, sol_num );
grid on;
xlabel('x [m]');
ylabel('Amplitude');
title(['Exact and Numerical and Solutions, with h = ', num2str(FemRegion.h)]);
legend('u_{ex}(x)', 'u_{num}(x)');

%==========================================================================
% ERROR ANALYSIS
%==========================================================================
[Errors] = ComputeErrors(Data,FemRegion,Solutions);






end
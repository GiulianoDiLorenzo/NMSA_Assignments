function [sol,Mesh,Data] = systemAssembly(Data, N_pts)
% ========================================================================
%   OUTPUT : 
%       - sol    --> Numerical approximation of the solution u_num of size
%                    (N_pts, 1)
%       - Mesh  -->  Structure of the Mesh of the system 
%       - Data   --> Structure of the Galerkin formulation with problem's
% 
%   INPUTS : 
%       - Data   --> Structure of the Galerkin formulation with problem's 
%       - N_pts  --> Number of nodes in the domain
% ========================================================================

fprintf('\n============================================================\n')
fprintf(['Solving test ', Data.name, ' with ',num2str(N_pts),' elements \n']);

%==========================================================================
% MESH GENERATION
%==========================================================================

[Mesh] = CreateMesh(Data,N_pts);

%==========================================================================
% BUILD FINITE ELEMENT MATRICES and RIGHT-HAND SIDE
%==========================================================================

plotPhi = false;
[Phi] = buildPhi(Mesh, plotPhi);

% % Mass matrix assembly;
fprintf('Building M ...\n');
[M] = buildM(Data, Mesh, Phi);

% Stiffness matrix assembly;
fprintf('Building A ...\n');
[A] = buildA(Data, Mesh);

% Boundary conditions matrix assembly
fprintf('Building R* ...\n');
[R_star] = buildR_star(Data,Mesh);

fprintf('Building M* = ( w * M - A - R* ) ...\n');
M_star = Data.omega * M - A - R_star;

% Excitation force vector
fprintf('Building F* = F + R_bc ...\n');
F = Mesh.h * Data.force(Data.x);

% Boundary conditions vector 
[R_bc] = buildR_bc(Data,Mesh);

F = F + R_bc;

%==========================================================================
% Inverting the system for the numerical solution 
%==========================================================================
fprintf('Computing numerical solution U = inv(M*) * F* \n');
sol = (M_star)\F;

sol = real(sol);

end
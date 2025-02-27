function [M_star, F] = systemAssembly(Data, Mesh, Phi, N_pts)
% ========================================================================
%   OUTPUT : 
%       - M_star --> Matrix of the left hand side of the system,
%                    size (N_pts, N_pts)
%       - F      --> Vector of the right hand side of the matrix system, 
%                    size (N_pts,1))   
% 
%   INPUTS : 
%       - Data   --> Structure of the Galerkin formulation with problem's
%       - Mesh   --> Structure of the Mesh of the system 
%       - Phi    --> Structure with handlers for the N_pts basis functions
%       - N_pts  --> Number of nodes in the domain
% ========================================================================

fprintf('\n============================================================\n')
fprintf(['Solving test ', Data.name, ' with ',num2str(N_pts),' elements \n']);


% Mass matrix assembly;
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

% %==========================================================================
% % Inverting the system for the numerical solution 
% %==========================================================================
% fprintf('Computing numerical solution U = inv(M*) * F* \n');
% sol = (M_star)\F;
% 
% sol = real(sol);

end
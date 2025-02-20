function [sol,Mesh,Data] = systemAssembly(Data, nEl)

fprintf('\n============================================================\n')
fprintf(['Solving test ', Data.name, ' with ',num2str(nEl),' elements \n']);

%==========================================================================
% MESH GENERATION
%==========================================================================

[Mesh] = CreateMesh(Data,nEl);

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
% Invert system for solution and extract u;
%==========================================================================
fprintf('Computing numerical solution U = inv(M*) * F* \n');
sol = (M_star)\F;

sol = real(sol);

end
function [L2_err,sol,Mesh,Data] = getResults(Data, nEl)

fprintf('============================================================\n')
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
fprintf('Building Matrices and Vectors ...\n');
[M] = buildM(Data, Mesh, Phi);

% Stiffness matrix assembly;
fprintf('Building A ...\n');
[A] = buildA(Data, Mesh);

% Boundary conditions matrix assembly
fprintf('Building R* ...\n');
[R_star] = buildR_star(Data,Mesh);

fprintf('Building M* ...\n');
M_star = Data.omega * M - A - R_star;

% Excitation force vector
fprintf('Building F...\n');
F = Mesh.h * Data.force(Data.x)';

% Boundary conditions vector 
[R_bc] = buildR_bc(Data,Mesh);

F = F + R_bc;

%==========================================================================
% Invert system for solution and extract u;
%==========================================================================
fprintf('Computing numerical solution...\n\n');
sol = (M_star)\F;

% sol = real(sol);
% sol = sol ./ ( max(sol)) );
% sol = abs(sol).* exp(1i * angle(sol));

%==========================================================================
% ERROR COMPUTATION
%==========================================================================
error = Data.uex(Data.x).' - sol;
L2_err = compute_L2_Error(error, Mesh.h);

end
function [A,b,u_g] = BoundaryConditions(A,b,femregion,Data)
%% [A,b,u_g] = BoundaryConditions(A,b,femregion,Data)
%==========================================================================
% Assign Dirchlet boundary conditions, after assembling the matrices
%==========================================================================
%    called in Main.m
%
%    INPUT:
%          A           : (sparse(ndof,ndof) real) stiffness matrix
%          b           : (sparse(ndof,1) real) rhs vector
%          femregion   : (struct)  see CreateFemregion.m
%          Data        : (struct)  see DataTest.m

%
%    OUTPUT:
%          A           : (sparse(ndof,ndof) real) stiffness matrix
%          b           : (sparse(ndof,1) real) rhs vector
%          u_g         : (sparse(ndof,1) real) evaluation of Dirichlet conditions 
%

fprintf('Assign Dirichlet boundary conditions ... \n');


ndof = length(b);
u_g = sparse(ndof,1);

if strcmp(Data.boundary,'DD')   
    boundary_points = femregion.boundary_points;
elseif strcmp(Data.boundary,'DN') 
    boundary_points = femregion.boundary_points(1);
elseif strcmp(Data.boundary,'ND')
    boundary_points = femregion.boundary_points(end);  
elseif strcmp(Data.boundary, 'DR')
    boundary_points = femregion.boundary_points(1); 
    D = zeros(ndof,1);
    D(end) = Data.mu(femregion.coord(end))*Data.gR(femregion.coord(end))*Data.alpha(femregion.coord(end));
end

x = femregion.dof(boundary_points,1);
u_g(boundary_points) = Data.gD(x); % Compute the lifting operator ug

A_0 = A;
b_0 = b-A*u_g; % modify the load vector --> F(v) = F(v) - a(ug,v)


% Reduce the system A in order to solve the pb with 
% homogeneous Dirichlet conditions 
for k = 1:length(boundary_points)
    A_0(boundary_points(k),:) = 0;
    A_0(:,boundary_points(k)) = 0;
    A_0(boundary_points(k),boundary_points(k)) = 1;
    b_0(boundary_points(k)) = 0;
end


b = b_0 + D;
A = A_0;

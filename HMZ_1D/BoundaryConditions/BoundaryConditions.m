function [A,b,u_g] = BoundaryConditions(A,b,femregion,Data)
%% [A,b,u_g] = BoundaryConditions(A,b,femregion,Data,t)
%==========================================================================
% Assign Dirchlet boundary conditions
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


ndof = length(b);
u_g = sparse(ndof,1);

if strcmp(Data.boundary,'DD')   
    boundary_points = femregion.boundary_points;
    x = femregion.dof(boundary_points,1);
    u_g(boundary_points(1)) = Data.gD1(Data.omega,Data.ro,Data.vel); 
    u_g(boundary_points(2)) = Data.gD2(Data.omega,Data.ro,Data.vel);
elseif strcmp(Data.boundary(1),'D') 
    boundary_points = femregion.boundary_points(1);
    x = femregion.dof(boundary_points,1);
    u_g(boundary_points(1)) = Data.gD1(Data.omega,Data.ro,Data.vel); 
elseif strcmp(Data.boundary(2),'D')
    boundary_points = femregion.boundary_points(end);  
    x = femregion.dof(boundary_points,1);
    u_g(boundary_points(1)) = Data.gD2(Data.omega,Data.ro,Data.vel); 
end


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

b = b_0;
A = A_0;

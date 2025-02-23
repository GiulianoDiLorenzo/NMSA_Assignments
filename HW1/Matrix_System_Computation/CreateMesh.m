function [Mesh] = CreateMesh(Data,N_pts)
% ========================================================================
%   OUTPUT : Structure of the Mesh of the system 

%   INPUTS : 
%       - Data      --> Structure of the Galerkin formulation with problem's 
%                       functions discretization
%       - N_pts     --> Number of nodes in the mesh
%       - omega     --> Scalar value for omega parameter
%       - mu_vals   --> array of (1,2) for the mu values on each half of
%                       the domain
%       - rho_vals  --> array of (1,2) for the rho values on each half of
%                       the domain
% ========================================================================


x0 = Data.domain(1);
xL = Data.domain(2);

p = linspace(x0,xL,N_pts);         % points of the mesh
 MeshSize = (xL-x0)./(N_pts-1);     % mesh size of the domain
%================================================

% Mesh data structure
Mesh = struct('dim',1,...                   % Dimension of the domain
               'domain',Data.domain,...     % Interval of the domain
               'h',MeshSize,...             % Mesh size of the systerm
               'n_pts',N_pts,...            % Number of nodes in the mesh
               'coord',p',...               % Coordinates of the mesh points
               'boundary_points',[x0,xL]);  % Boundaries of the domain
           
           
end
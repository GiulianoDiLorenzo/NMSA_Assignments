function [R_bc] = buildR_bc(Data,Mesh)
% ========================================================================
%   OUTPUT : Boundary Conditions "Force" vector Array of size (num_pts,1)

%   INPUTS : 
%       - Data --> Structure of the Galerkin formulation of the pb
%       - Mesh --> Structure of the Mesh representation 
% ========================================================================


num_pts = Mesh.n_pts;
h =  Mesh.h;
beta = Data.mu(0) / h;

% M of size (N_h+1)x(N_h+1)
R_bc = zeros( num_pts , 1);

R_bc(1) = -beta * Data.gD;
R_bc(end) = - 1i* Data.alpha * Data.gR;

end

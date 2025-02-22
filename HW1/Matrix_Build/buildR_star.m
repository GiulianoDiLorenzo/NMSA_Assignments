function [R_star] = buildR_star(Data,Mesh)
% ========================================================================
%   OUTPUT : Extra residual Matrix of size (num_pts,num_pts)

%   INPUTS : 
%       - Data --> Structure of the Galerkin formulation of the pb
%       - Mesh --> Structure of the Mesh representation 
% ========================================================================

num_pts = Mesh.n_pts;

h =  Mesh.h;
beta = Data.mu(0) / h;


% R_star of size (N_h+1)x(N_h+1)
R_star=zeros( num_pts , num_pts);

R_star(1,2) =  beta;
R_star(end,end) = 1i * Data.alpha;
end

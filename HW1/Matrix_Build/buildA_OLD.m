function [A] = buildA_OLD(Data,Mesh)
% ========================================================================
%   OUTPUT : Stiffness matrix of dimension (num_pts , num_pts)

%   INPUTS : 
%       - Data --> Structure of Galerkin formulation of the pb
%       - Mesh --> Structure of the Mesh representation 
% ========================================================================
    

num_pts = Mesh.n_pts;
h =  Mesh.h;
mu = Data.mu;
x = Mesh.coord;

% A of size (N_h+1)x(N_h+1)
A = zeros( num_pts , num_pts );

% Computing the edges of A
A(1,1) =  ( 2*mu(x(1)) + mu(x(2)) ) / (2*h);
A(num_pts,num_pts) = ( 2*mu(x(end-1)) + mu(x(end)) ) / (2*h);

%Computing the superior triangular zone of A
for i = 2:num_pts-1  

    % Storing points of interests
    x_prev = x(i-1);
    x_now  = x(i);
    x_next = x(i+1);


    % ===================================================================
    % ====================== ELEMENTS A(i,i-1) ==========================
    % ===================================================================
    % Computing m_{i,i-1} using Trapezoidal rule
    A(i,i-1) = - (mu(x_prev) + mu(x_now)) / (2*h); 


    % ===================================================================
    % ====================== ELEMENTS A(i,i) ===========================
    % ===================================================================
    % Computing m_{i,i} using Trapezoidal rule
    A(i,i) = (mu(x_prev) + 2*mu(x_now) + mu(x_next) ) / (2*h) ;

    % ===================================================================
    % ====================== ELEMENTS A(i,i+1) ==========================
    % ===================================================================
    % Computing m_{i,i+1} using Trapezoidal rule
    A(i,i+1) = - (mu(x_now) + mu(x_next)) / (2*h); 
end


% ===================================================================
% ===================== MAKING A SYMMETRIC===========================
% ===================================================================
A(1,2) = A(2,1);
A(end,end-1) = A(end-1, end);

end

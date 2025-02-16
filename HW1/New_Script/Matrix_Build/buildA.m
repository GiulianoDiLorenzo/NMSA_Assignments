function [A] = buildA(Data,Femregion, Phi)
    %==========================================================================
% Build the local mass matrix for the term (uv)
%==========================================================================
%    called in C_matrix1D.m
%
%    INPUT:
%          dphiq       : (array real) evaluation of the basis function on
%                        quadrature nodes
%          w_1D        : (array real) quadrature weights
%          nln         : (integer) number of local unknowns
%          BJ          : Jacobian of the map 
%
%    OUTPUT:
%          M_loc       :  (array real) Local mass matrix

num_pts = length(Femregion.dof);
h =  Femregion.h;
mu = Data.mu;
x = Femregion.coord;

% A of size (N_h+1)x(N_h+1)
A=zeros( num_pts , num_pts );

% Computing the boundaries of A
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
    % ====================== ELEMENTS A(i,i1) ===========================
    % ===================================================================

    % Computing m_{i,i} using Trapezoidal rule
    A(i,i) = (mu(x_prev) + 2*mu(x_now) + mu(x_next) ) / (2*h) ;

    % ===================================================================
    % ====================== ELEMENTS A(i,i+1) ==========================
    % ===================================================================

    % Computing m_{i,i+1} using Trapezoidal rule
    A(i,i+1) = - (mu(x_now) + mu(x_next)) / (2*h); 
end

%Making A symmetric
A(1,2) = A(2,1);
for i = 1 : num_pts-1
    A(i+1,i) = A(i,i+1);
end

% % % New considerations
% A(1,:) = zeros(1, num_pts);
% A(1,1) = 1;
% 
% A(end,:) = zeros(1, num_pts);
% A(end,end) = 1;
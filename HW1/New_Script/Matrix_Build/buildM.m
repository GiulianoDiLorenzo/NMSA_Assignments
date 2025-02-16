function [M] = buildM(Data,Femregion, Phi)
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
% w = Data.omega;
rho = Data.rho;
x = Femregion.coord;


% M of size (N_h+1)x(N_h+1)
M=zeros( num_pts , num_pts );


phi_1 = Phi{1};
% phi_1_val = phi_1(x);


M(1,1) =  h * rho(0);


phi_end = Phi{end};


M(num_pts,num_pts) = h * rho(1) * phi_end(Data.L)^2  / 2 ;

%Computing the superior triangular zone of M
for i = 2:num_pts-1
    
    phi_i_now = Phi{i};
  
    % Storing points of interests
    x_prev = x(i-1);
    x_now  = x(i);
    x_next = x(i+1);

    % ===================================================================
    % ====================== ELEMENTS M(i,i-1) ==========================
    % ===================================================================

    % Points of reference for [x_prev,x_now] = int_prev
    x_j_p_prev = ( x_prev + x_now )/2 + h/(2*sqrt(3));
    x_j_n_prev = ( x_prev + x_now )/2 - h/(2*sqrt(3));


    % Computing phi_{i-1}
    phi_i_prev = Phi{i-1}; 

    prod1= rho(x_j_p_prev) * phi_i_prev(x_j_p_prev) * phi_i_now(x_j_p_prev);
    prod2= rho(x_j_n_prev) * phi_i_prev(x_j_n_prev) * phi_i_now(x_j_n_prev);

    sum1 = prod1 + prod2;

    M(i,i-1) = h/2 * sum1;

    % ===================================================================
    % ====================== ELEMENTS M(i,i) ===========================+
    % ===================================================================
    % int_prev = [x_prev ; x_now];
    % int_next = [x_now ; x_prev]; 
    % 
    % % Applying trapezoidal rule for diagonal elements
    % sum_prev = sum(rho(int_prev) .* (phi_i_now(int_prev)).^2)
    % sum_next = sum(rho(int_next) .* (phi_i_now(int_next)).^2)

    % sum_prev = rho(x_prev)*(phi_i_now(x_prev)).^2 + rho(x_now)*(phi_i_now(x_now)).^2;
    % sum_next = rho(x_next)*(phi_i_now(x_next))^2 + rho(x_now)*(phi_i_now(x_now))^2;
    M(i,i)   = h * rho(x(i)) ;  

    % ===================================================================
    % ====================== ELEMENTS M(i,i+1) ==========================
    % ===================================================================
    % Points of rference for [x_now,x_next] = int_next
    x_j_p_next = ( x_next + x_now )/2 + h/(2*sqrt(3));
    x_j_n_next = ( x_next + x_now )/2 - h/(2*sqrt(3));

    % x_j_next = [x_j_n_next; x_j_p_next];

    % Computing phi_{i+1}
    phi_i_next = Phi{i+1}; 

    prod3 = rho(x_j_p_next) * phi_i_next(x_j_p_next) * phi_i_now(x_j_p_next);
    prod4 = rho(x_j_n_next) * phi_i_next(x_j_n_next) * phi_i_now(x_j_n_next);
    sum2 = prod3 + prod4;

    % Computing m_{i,i+1} using Gauss-Legendre Polynomials
    % sum1 = sum( rho(x_j_next) .* phi_i_next(x_j_next) .* phi_i_now(x_j_next));
    
    M(i,i+1) = h/2 * sum2;
   
end
    

%Making M symmetric
M(1,2) = M(2,1);
for i = 1 : num_pts-1
    M(i+1,i) = M(i,i+1);
end
% 
% % % New considerations
% M(1,:) = zeros(1, num_pts);
% M(1,1) = 1;
% 
% M(end,:) = zeros(1, num_pts);
% M(end,end) = 1;


end
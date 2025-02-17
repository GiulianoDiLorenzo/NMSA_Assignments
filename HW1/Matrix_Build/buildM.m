function [M] = buildM(Data,Mesh, Phi)

num_pts = Mesh.n_pts;

h =  Mesh.h;
% w = Data.omega;
rho = Data.rho;
x = Mesh.coord;


% M of size (N_h+1)x(N_h+1)
M=zeros( num_pts , num_pts );

M(1,1) =  h * rho(0);
M(num_pts,num_pts) = h * rho(1)  / 2 ;


%Computing the superior triangular zone of M
for i = 2:num_pts-1

    % Storing points of interests
    x_prev = x(i-1);
    x_now  = x(i);
    x_next = x(i+1);

    % Set of basis functions to consider
    phi_i_now = Phi{i};
    phi_i_prev = Phi{i-1};
    phi_i_next = Phi{i+1};

    % ===================================================================
    % ====================== ELEMENTS M(i,i-1) ==========================
    % ===================================================================

    % Points of reference for [x_prev,x_now] = int_prev
    x_j_p_prev = ( x_prev + x_now )/2 + h/(2*sqrt(3));
    x_j_n_prev = ( x_prev + x_now )/2 - h/(2*sqrt(3));

    prod1= rho(x_j_p_prev) * phi_i_prev(x_j_p_prev) * phi_i_now(x_j_p_prev);
    prod2= rho(x_j_n_prev) * phi_i_prev(x_j_n_prev) * phi_i_now(x_j_n_prev);
    sum1 = prod1 + prod2;

    M(i,i-1) = h/2 * sum1;

    % ===================================================================
    % ====================== ELEMENTS M(i,i) ===========================+
    % ===================================================================

    M(i,i)   = h * rho(x(i)) ;  

    % ===================================================================
    % ====================== ELEMENTS M(i,i+1) ==========================
    % ===================================================================
    
    % Points of rference for [x_now,x_next] = int_next
    x_j_p_next = ( x_next + x_now )/2 + h/(2*sqrt(3));
    x_j_n_next = ( x_next + x_now )/2 - h/(2*sqrt(3));

    % x_j_next = [x_j_n_next; x_j_p_next];

    prod3 = rho(x_j_p_next) * phi_i_next(x_j_p_next) * phi_i_now(x_j_p_next);
    prod4 = rho(x_j_n_next) * phi_i_next(x_j_n_next) * phi_i_now(x_j_n_next);
    sum2 = prod3 + prod4;

    M(i,i+1) = h/2 * sum2;
   
end
    

%Making M symmetric
M(1,2) = M(2,1);
for i = 1 : num_pts-1
    M(i+1,i) = M(i,i+1);
end

end
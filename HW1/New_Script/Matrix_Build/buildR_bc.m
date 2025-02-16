function [R_bc] = buildR_bc(Data,Femregion)


num_pts = length(Femregion.dof);

h =  Femregion.h;
beta = Data.mu(0) / h;


% M of size (N_h+1)x(N_h+1)
R_bc=zeros( num_pts , 1);

R_bc(1) = -beta * Data.gD;
R_bc(end) = - 1i* Data.alpha * Data.gR;

end

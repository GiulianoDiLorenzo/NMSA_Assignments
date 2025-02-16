function [R_star] = buildR_star(Data,Femregion)


num_pts = length(Femregion.dof);

h =  Femregion.h;
beta = Data.mu(0) / h;


% R_star of size (N_h+1)x(N_h+1)
R_star=zeros( num_pts , num_pts);

R_star(1,2) =  beta;
R_star(end,end) = 1i * Data.alpha;
end

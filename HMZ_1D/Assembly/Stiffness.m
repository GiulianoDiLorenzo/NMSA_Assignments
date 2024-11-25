function [K_loc] = Stiffness(Grad,w_1D,nln,BJ)
%% [K_loc] = Stiffness(Grad,w_1D,nln,BJ)
%==========================================================================
% Build the local stiffness matrix for the term u'*v'
%==========================================================================
%    called in C_matrix1D.m
%
%    INPUT:
%          Grad        : (array real) evaluation of the derivative on
%                        quadrature nodes
%          w_1D        : (array real) quadrature weights
%          nln         : (integer) number of local unknowns
%          BJ          : (array real) Jacobian of the map 
%
%    OUTPUT:
%          K_loc       :  (array real) Local stiffness matrix


K_loc = zeros(nln,nln);

%% General implementation -- to be used with general finite element spaces
for i=1:nln
    for j=1:nln
        for k=1:length(w_1D)
            Binv = 1./BJ;    % inverse
            Jdet = BJ;       % determinant 
            K_loc(i,j) = K_loc(i,j) + (Jdet.*w_1D(k)) .* ( (Grad(k,:,i) * Binv) * (Grad(k,:,j) * Binv )');
        end
    end
end



                                              
                                              


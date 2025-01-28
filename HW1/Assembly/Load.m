function [f] = Load(force, dphiq, BJ, w_1D, pphys_1D, nln)
%% [f] = Load(force,dphiq,BJ,w_1D,pphys_1D,nln)
%==========================================================================
% Build the right hand side vector (fv)
%==========================================================================
%    called in C_matrix1D.m
%
%    INPUT:
%          force       : (string) expression of the forcing term
%          dphiq       : (array real) basis functions evaluated at q.p.
%          BJ          : (array real) Jacobian of the map 
%          w_1D        : (array real) quadrature weights
%          pphys_1D    : (array real) quadrature nodes in the physical
%                                     space
%          nln         : (integer) number of local unknowns
%
%    OUTPUT:
%          f           : (array real) Local right hand side


f = zeros(nln,1);
% evaluation of the right hand side on the physical nodes
x = pphys_1D;
F = force(x);

% Evaluation of the local r.h.s.
for s = 1:nln
    for k = 1:length(w_1D)
        Jdet = BJ;  % determinant 
        f(s) = f(s) + w_1D(k)*Jdet*F(k)*dphiq(1,k,s);
    end    
end

function [F] = buildF(Data)
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

F = Data.force(Data.x)';
end
% quadrature nodes and whieghs for faces integrals
function  [node_1D,w_1D] = Quadrature(nqn)
%% function [node_1D,w_1D] = Quadrature(Data)
%==========================================================================
% Compute quadrature nodes and weights triangular and quadrilateral elements
%==========================================================================
%    called in C_matrix1D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%
%    OUTPUT:
%          node2D      : (struct) .num  (int) number of basis functions
%                                 .nedge (int) number of edges
%                                 .fbases (string) basis functions
%                                 .Gbases  (string) df/dx
%

% Number of quadrature nodes for each element
switch nqn
    
    case{1}
        %===================================================
        % mid point rule
        %===================================================
        
        xnod = 0.5;
        w = 1;
        
    case{2}
        
        %===================================================
        % Trapezi rule
        %===================================================
        
        w(1) = 0.5;
        w(2) = 0.5;
        
        xnod(1) = 0;
        xnod(2) = 1;
        
        
    case{3}
        %===================================================
        % Cavalieri-Simpson rule
        %===================================================
        
        w(1) = 1/3;
        w(2) = 4/3;
        w(3) = 1/3;
        
        xnod(1) = 0;
        xnod(2) = 0.5;
        xnod(3) = 1;
        
end



node_1D=[xnod'];
w_1D=w;

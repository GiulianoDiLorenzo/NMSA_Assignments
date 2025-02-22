function [F] = buildF(Data)
% ========================================================================
%   OUTPUT : Force vector Array of size (num_pts,1)

%   INPUTS : 
%       - Data --> Structure of Galerkin formulation of the pb
% ========================================================================

F = Data.force(Data.x)';


end
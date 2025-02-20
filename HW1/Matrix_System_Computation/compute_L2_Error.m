function [L2_err] = compute_L2_Error(error, h)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
L2_err = 0;

for i = 1:length(error)-1
    sum_sq_err = error(i).^2 + error(i+1).^2;
    L2_err = L2_err + sum_sq_err;
end

L2_err = sqrt( h/2 * L2_err) ;

end


function PointwiseSol(femregion, uh, u_ex)
%% PointwiseSol(femregion, uh, u_ex)
%==========================================================================
% PLOT THE EXACT SOLUTION ON THE DOFS
%==========================================================================
%    called in PostProcessing.m
%
%    INPUT:
%          femregion   : (struct)  see CreateFemregion.m
%          uh          : (sparse(ndof,1) real) solution vector
%          u_ex        : (sparse(ndof,1) real) exact solution vector
%          Data        : (struct) see DataTest.m
%


M = max(uh);
m = min(uh);
if (abs(m - M) < 0.1)
    M = m + 1;
end

figure;
plot(femregion.coord(:,1),full(uh));
title('u_h(x)'); xlabel('x-axis'); ylabel('y-axis');

if(min(u_ex) ~= 0 && max(u_ex) ~= 0)
    hold on;
    plot(femregion.coord(:,1),u_ex,'r*');
    legend('u_h', 'u', 'Interpreter', 'Latex');
end


function [uh] = Snapshot(femregion, uh, string)
%% Snapshot(femregion, uh, string)
%==========================================================================
% PLOT THE EXACT SOLUTION ON THE DOFS
%==========================================================================
%    called in MainHMZ.m
%
%    INPUT:
%          femregion   : (struct)  see CreateFemregion.m
%          uh          : (sparse(ndof,1) real) solution vector
%          string      : (string) title of the figure
%


% x1 = femregion.domain(1,1);
% x2 = femregion.domain(1,2);
% M  =  2; % max(uh)
% m  = -2; % min(uh);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% PLOT OF SOLUTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(femregion.coord(:,1),full(uh));
title( string ); xlabel('x-axis'); ylabel('y-axis');
% axis([x1,x2,m,M]); 

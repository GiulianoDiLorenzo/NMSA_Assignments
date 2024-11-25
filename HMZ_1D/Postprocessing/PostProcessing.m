function [solutions]=PostProcessing(Data,femregion,uh)
%% [solutions]=PostProcessing(Data,femregion,uh)
%==========================================================================
% POST PROCESSING OF THE SOLUTION
%==========================================================================
%    called in Main.m
%
%    INPUT:
%          Data        : (struct)  see DataTest.m
%          femregion   : (struct)  see CreateFemregion.m
%          uh          : (sparse(ndof,1) real) solution vector
%
%    OUTPUT:
%          solutions   : (struct) containg solution vector uh and
%                        analytical solution u_ex
%

fprintf('Plot the solution ... \n');


%==========================================================================
% EVALUATION OF THE EXACT SOLUTION
%==========================================================================

dof = femregion.dof;
x = dof(:,1);
u_ex = Data.uex(x,Data.omega,Data.ro,Data.vel);


%==========================================================================
% PLOT SOLUTION
%==========================================================================

if(Data.visual_graph)
    PointwiseSol(femregion, uh, u_ex);
end

%==========================================================================
% SAVE SOLUTIONS
%==========================================================================

solutions = struct('u_ex',u_ex,'uh',uh);

fprintf('============================================================\n')

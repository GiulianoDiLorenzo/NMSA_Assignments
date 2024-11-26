function [errors] = ComputeErrors(Data, femregion, solutions)
%% [errors] = C_compute_errors(femregion, solutions)
%==========================================================================
% Compute L2, semi-H1, H1 and L-inf errors
%==========================================================================
%    called in Main.m
%
%    INPUT:
%          Data        : (struct)  see DataTest.m
%          femregion   : (struct)  see CreateFemregion.m
%          solutions   : (struct)  see PostProcessing.m
%
%    OUTPUT:
%          errors      : (struct)  


err = solutions.uh - solutions.u_ex;

[E_L2, E_SEMI_H1] = Calc_L2_H1_errors(femregion, solutions.uh, Data);

errors = struct('L2',   E_L2,...
                'H1',   sqrt(E_L2.^2 + E_SEMI_H1.^2),...
                'H1_S', E_SEMI_H1,...
                'inf',  norm (err,inf));
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


err = real(solutions.u_num) - solutions.u_ex;

[E_L2, E_H1] = Calc_L2_H1_errors(femregion, err, Data);

errors = struct('L2_err',   E_L2,...
                'H1_err',   E_H1,...
                'inf_err',  norm (err,inf));
function [basis] = ShapeBasis
%% function [basis] = ShapeBasis
%==========================================================================
% Compute the shape functions for intervals
%==========================================================================
%    called in Matrix1D.m
%
%    INPUT:
%
%    OUTPUT:
%          basis       : (struct) .num  (int) number of basis functions
%                                 .n_edge (int) number of end-points
%                                 .fbases (string) basis functions
%                                 .Gbases (string) df/dx



nln = 2;
basis = struct('num',nln,...
               'n_edge',2,...
               'fbases',{'1-csi',...
                         'csi',...
                        },...
               'Gbases',{'-1+0.*csi',...
                         ' 1+0.*csi',...
                        });
        
        

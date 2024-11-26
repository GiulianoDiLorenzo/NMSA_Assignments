function [femregion] = CreateFemregion(Data,Region) 
%% [femregion] = CreateFemregion(Data,Region)
%==========================================================================
% Creates finite element space
%==========================================================================
%    called in C_main1D.m
%
%    INPUT:
%          Data        : (struct)  see DataTest.m
%          Region      : (struct)  see CreateMesh.m
%
%    OUTPUT:
%          femregion    : (struct) finite element space

fprintf('Creating finite element space ... \n');


nln = 2;  
        
bound_pts = ones(length(Region.coord),1);
bound_pts = find(Region.coord(:,1)== Data.domain(1,1) | Region.coord(:,1) == Data.domain(1,2));


%==========================================================================
% COSTRUZIONE STRUTTURA FEMREGION
%==========================================================================
femregion=struct('fem',1,...
                'domain',Region.domain,...
                'h', Region.h,...
                'nln',nln,...
                'ndof',length(Region.coord),...
                'ne',Region.ne,...
                'dof',Region.coord,...
                'nqn_1D',2,...
                'coord',Region.coord,...
                'connectivity',Region.connectivity,...
                'boundary_points',bound_pts);
            
            
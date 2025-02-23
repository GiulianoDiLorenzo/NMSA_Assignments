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


nln = 2;   %localo degrees of freedom = 2 in 1D as we consider 2 points
        
bound_pts = ones(length(Region.coord),1); %Region.coord contains all the points coordinates
bound_pts = find(Region.coord(:,1)== Data.domain(1,1) | Region.coord(:,1) == Data.domain(1,2)); %Finding the points coreesponding to the boundary 


%==========================================================================
% COSTRUZIONE STRUTTURA FEMREGION
%==========================================================================
femregion=struct('fem',1,...
                'domain',Region.domain,...
                'h', Region.h,...
                'nln',nln,...
                'ndof',length(Region.coord),...
                'ne',Region.ne,...          %Number of elements in the mesh
                'dof',Region.coord,...      %Points coordinates in the mesh
                'nqn_1D',2,...
                'coord',Region.coord,...    %Points coordinates
                'connectivity',Region.connectivity,...
                'boundary_points',bound_pts);
            
            
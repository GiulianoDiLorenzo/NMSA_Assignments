function [Mesh] = CreateMesh(Data,nEl)
%% [Region] = C_create_mesh(Data, nEl)
%=========================================================================
% Creates regular mesh
%==========================================================================
%    called in Main.m
%
%    INPUT:
%          Data    : (struct)  see DataTest.m
%          nEl     : (int)    Number of mesh elements  
%
%    OUTPUT:
%          Region  : (struct) having fields: dimension
%                                            domain 
%                                            mesh size
%                                            number of vertices
%                                            number of elements
%                                            coordinates
%                                            boundary points
%                                            connectivity



x0 = Data.domain(1);
xL = Data.domain(2);


%================================================
% Geometrical info
 n_pts =nEl;                   % vertices of the mesh
 p = linspace(x0,xL,n_pts);         % points of the mesh
 t = [[1:n_pts-1]' [2:n_pts]']';    % intervals of the mesh, in terms of indexes
 MeshSize = (xL-x0)./(nEl-1);           % mesh size of the domain
%================================================

% Mesh data structure
Mesh = struct('dim',1,...
               'domain',Data.domain,...
               'h',MeshSize,...
               'n_pts',n_pts,...
               'coord',p',...               %coordinates of the mesh points
               'boundary_points',[x0,xL]);
           
           

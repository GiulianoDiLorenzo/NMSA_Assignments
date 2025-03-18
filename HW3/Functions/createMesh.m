function [Mesh] = createMesh(L, T, Nx, Nt)

Mesh = struct();

Mesh.L = L;         % Length of the road
Mesh.T = T;         % Duration of the simulation    

Mesh.domain = [[0,L]; [0,T]];


Mesh.Nx = Nx;       % Number of spatial points
Mesh.Nt = Nt;       % Number of time steps


Mesh.dx = L / Nx;                % Spatial step
Mesh.dt = T / Nt;                % Time step (CFL condition should be checked)

Mesh.x = linspace(0, L, Nx).';   % Spatial grid
Mesh.t = linspace(0,T,Nt+1).';     % time grid


end
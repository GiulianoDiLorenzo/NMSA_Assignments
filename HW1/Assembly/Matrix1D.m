function [A, f]=Matrix1D(Data,femregion)
%% [A,f] = Matrix1D(Data,Femregion)
%==========================================================================
% Assembly of the stiffness matrices A and rhs f
%==========================================================================
%    called in Main.m
%
%    INPUT:
%          Data        : (struct)  see DataTest.m
%          femregion   : (struct)  see CreateFemregion.m
%
%    OUTPUT:
%          A           : (sparse(ndof,ndof) real) stiffnes matrix
%          f           : (sparse(ndof,1) real) rhs vector


fprintf('Assembling matrices and right hand side ... \n');


% connectivity infos
ndof         = femregion.ndof;         % degrees of freedom
nln          = femregion.nln;          % local degrees of freedom
ne           = femregion.ne;           % number of elements
connectivity = femregion.connectivity; % connectivity matrix, list of intervals


% shape functions
[basis] = ShapeBasis;

% quadrature nodes and weights for integrals
[nodes_1D, w_1D] = Quadrature(2);

% evaluation of shape bases on quadrature nodes
[Phi,GradPhi] = EvalShapeBasis(basis,nodes_1D);


% Assembly begin ...
A = sparse(ndof,ndof);  % Global Stiffness matrix
M = sparse(ndof,ndof);  % Global mass matrix

f = sparse(ndof,1);     % Global Load vector

for ie = 1 : ne
     
    % Local to global map --> To be used in the assembly phase
    iglo = connectivity(1:nln,ie); %current interval of interest, in terms of indexes
  
    [BJ, nodes_1D_phys] = GetJacobian(femregion.coord(iglo,:), nodes_1D);
    % BJ        = Jacobian of the elemental map 
    % nodes_1D_phys  = vertex coordinates in the physical domain, i.e. in
    % [0,L], x_ie, and x_{ie+1}
   
    %=============================================================%
    % STIFFNESS MATRIX
    %=============================================================%
    
    % Local stiffness matrix 
    [A_loc] = Stiffness(GradPhi, w_1D, nln, BJ);

    % Assembly phase for stiffness matrix
    A(iglo,iglo) = A(iglo,iglo) + A_loc;    %[Data.mu(ie), Data.mu(ie+1); Data.mu(ie), Data.mu(ie+1)].*A_loc; 
    
    %=============================================================%
    % MASS MATRIX
    %=============================================================%
    
    % Local mass matrix 
    [M_loc] = Mass(Phi, w_1D, nln, BJ);         %Data.omega*Data.rho(ie)*BJ);

    % Assembly phase for mass matrix
    M(iglo,iglo) = M(iglo,iglo) + Data.omega*Data.rho(ie)*M_loc;   
    
    %==============================================
    % FORCING TERM --RHS
    %==============================================

    % Local load vector
    [load] = Load(Data.force, Phi, BJ, w_1D, nodes_1D_phys, nln);    

    % Assembly phase for the load vector
    f(iglo) = f(iglo) + load;

end


A = M-A;

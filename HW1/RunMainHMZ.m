clear; close all;

%% Path
addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing


%% Data for Test
Data = DataTest('Test1');

%% Options
Data.visual_graph = 1;
Data.calc_errors = 1;
Data.eig_analysis = 0;
Data.plot_eigvct = 0;


%% Main routine
%[err1,sol1,fem1,D1] = MainHMZ(Data,50);
[err2,sol2,fem2,D2] = MainHMZ(Data,100);
%[err3,sol3,fem3,D3] = MainHMZ(Data,200);
%[err4,sol4,fem4,D4] = MainHMZ(Data,400);


%% Plot Errors

if Data.calc_errors
    
    hVec   = [fem1.h, fem2.h, fem3.h, fem4.h];
    eVecL2 = [err1.L2, err2.L2, err3.L2, err4.L2];
    eVecH1 = [err1.H1, err2.H1, err3.H1, err4.H1];
    
    hs = subplot(2,1,1);
    loglog(hVec,hVec.^2,'-+b','Linewidth',2); hold on;
    loglog(hVec,eVecL2,'-or','Linewidth',2);
    legend(sprintf('h^%i',2),'||u-u_h||_{L^2}');
    ylabel('L^2-error');
    xlabel('h');
    hs.FontSize = 12;
    
    hs = subplot(2,1,2);
    loglog(hVec,hVec,'-+b','Linewidth',2); hold on;
    loglog(hVec,eVecH1,'-or','Linewidth',2);
    legend(sprintf('h^%i',1),'||u-u_h||_{H^1}');
    ylabel('H^1-error')
    xlabel('h');
    hs.FontSize = 12;
    
end
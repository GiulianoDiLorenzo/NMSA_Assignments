clear; close all;

%% Path
addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing


%% Data for Test

TestName = 'TestHW1_1';

Data = DataTest(TestName);
%Test2 solves the easy version of the problem with Dirichlet conditions put
%to 0 and f = 0

%% Options
Data.visual_graph = 1;
Data.calc_errors = 1;

%% Question 3_a
if strcmp(TestName,'TestHW1_3a')
    lambda = 1/(3*Data.f_0);
    nEl = ceil(10*Data.L/lambda);
    [err3a,sol3a,fem3a,D3a] = Main(Data,nEl);
end

%% Main routine
[err1,sol1,fem1,D1] = Main(Data,100);
[err2,sol2,fem2,D2] = Main(Data,200);
[err3,sol3,fem3,D3] = Main(Data,400);
[err4,sol4,fem4,D4] = Main(Data,800);

%% Plot Errors
hVec   = [fem1.h, fem2.h, fem3.h, fem4.h]; 
eVecL2 = [err1.L2, err2.L2, err3.L2, err4.L2];
eVecH1 = [err1.H1, err2.H1, err3.H1, err4.H1];

figure()
sgtitle('Graph of the error, and mesh size')
hs = subplot(2,1,1);
loglog(hVec,hVec.^2,'-b','Linewidth',1); hold on;
loglog(hVec,eVecL2,'-+r','Linewidth',2);
legend(sprintf('h^%i',2),'||u-u_h||_{L^2}');
ylabel('L^2-error');
xlabel('h');
hs.FontSize = 12;
grid on;

hs = subplot(2,1,2);
loglog(hVec,hVec,'-b','Linewidth',1); hold on;
loglog(hVec,eVecH1,'-+r','Linewidth',2);
legend(sprintf('h^%i',1),'||u-u_h||_{H^1}');
ylabel('H^1-error')
xlabel('h');
hs.FontSize = 12;
grid on;
 

clear; close all;
clc;
%% Path

addpath Errors
addpath FESpace
addpath Postprocessing
addpath Matrix_Build



%% Question 3_1
TestName = 'TestHW1_1';
omega = 1;



nEl = 10+1; % Number of points in the mesh [0,L]
omega = 1;

[Data] = NewDataTest(TestName, nEl, omega);
[err0, sol0, fem0, data0] = Pipeline(Data, nEl);


nEl = 100+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err1, sol1, fem1, data1] = Pipeline(Data, nEl);


nEl = 200+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err2, sol2, fem2, data2] = Pipeline(Data, nEl);


nEl = 400+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err3, sol3, fem3, data3] = Pipeline(Data, nEl);


nEl = 800+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err4, sol4, fem4, data4] = Pipeline(Data, nEl);


nEl = 1e3+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err5, sol5, fem5, data5] = Pipeline(Data, nEl);

nEl = 2*1e3+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err6, sol6, fem6, data6] = Pipeline(Data, nEl);

nEl = 4*1e3+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err7, sol7, fem7, data7] = Pipeline(Data, nEl);

nEl = 8*1e3+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err8, sol8, fem8, data8] = Pipeline(Data, nEl);
%%

L2_errors = [err0.L2_err , err1.L2_err , err2.L2_err , err3.L2_err , ...
             err4.L2_err , err5.L2_err , err6.L2_err , err7.L2_err , err8.L2_err , ];

H1_errors = [err0.H1_err , err1.H1_err , err2.H1_err , err3.H1_err , ...
             err4.H1_err , err5.H1_err , err6.H1_err , err7.H1_err ,err8.H1_err];

stepsizes = [fem0.h, fem1.h, fem2.h, fem3.h,fem4.h, fem5.h , fem6.h , fem7.h , fem8.h];

figure()

loglog(stepsizes, L2_errors, '-+r');
grid on;

title('Graph of the L2 error $||u_{ex}-u_{num}||$, as a function of the stepszize $h$', Interpreter='latex');
ylabel('$|| u_{ex} - u_{num}||_{L^2}$', Interpreter='latex')
xlabel('$h$ [m]', Interpreter='latex')

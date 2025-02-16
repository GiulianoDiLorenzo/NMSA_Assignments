clear; close all;
clc;
%% Path

addpath Errors
addpath FESpace
addpath Postprocessing
addpath Matrix_Build



%% Question 2a
TestName = 'TestHW1_2a';
omega = 1;



nEl = 10+1; % Number of points in the mesh [0,L]
omega = 1;

[Data] = NewDataTest(TestName, nEl, omega);
[err0, sol0, fem0, data0] = Pipeline(Data, nEl);


nEl = 100+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err1, sol1, fem1, data1] = Pipeline(Data, nEl);


nEl = 500+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err2, sol2, fem2, data2] = Pipeline(Data, nEl);


nEl = 1e3+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err3, sol3, fem3, data3] = Pipeline(Data, nEl);


nEl = 5e3+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err4, sol4, fem4, data4] = Pipeline(Data, nEl);


nEl = 1e4+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err5, sol5, fem5, data5] = Pipeline(Data, nEl);

%% PLots of the L2-norm error

L2_errors = [err0.L2_err , err1.L2_err , err2.L2_err , err3.L2_err , ...
             err4.L2_err , err5.L2_err];

H1_errors = [err0.H1_err , err1.H1_err , err2.H1_err , err3.H1_err , ...
             err4.H1_err , err5.H1_err];

stepsizes = [fem0.h, fem1.h, fem2.h, fem3.h,fem4.h, fem5.h];

figure()

loglog(stepsizes, L2_errors, '-+r');
grid on;

title('Graph of the L2 error $||u_{ex}-u_{num}||$, as a function of the stepszize $h$', Interpreter='latex');
ylabel('$|| u_{ex} - u_{num}||_{L^2}$', Interpreter='latex')
xlabel('$h$ [m]', Interpreter='latex')


%% Question 2b
TestName = 'TestHW1_2b';
omega = 1;



nEl = 10+1; % Number of points in the mesh [0,L]
omega = 1;

[Data] = NewDataTest(TestName, nEl, omega);
[err0, sol0, fem0, data0] = Pipeline(Data, nEl);


nEl = 100+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err1, sol1, fem1, data1] = Pipeline(Data, nEl);


nEl = 500+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err2, sol2, fem2, data2] = Pipeline(Data, nEl);


nEl = 1e3+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err3, sol3, fem3, data3] = Pipeline(Data, nEl);


nEl = 5e3+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err4, sol4, fem4, data4] = Pipeline(Data, nEl);


nEl = 1e4+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err5, sol5, fem5, data5] = Pipeline(Data, nEl);

%% PLots of the L2-norm error

L2_errors = [err0.L2_err , err1.L2_err , err2.L2_err , err3.L2_err , ...
             err4.L2_err , err5.L2_err];

H1_errors = [err0.H1_err , err1.H1_err , err2.H1_err , err3.H1_err , ...
             err4.H1_err , err5.H1_err];

stepsizes = [fem0.h, fem1.h, fem2.h, fem3.h,fem4.h, fem5.h];

figure()

loglog(stepsizes, L2_errors, '-+r');
grid on;

title('Graph of the L2 error $||u_{ex}-u_{num}||$, as a function of the stepszize $h$', Interpreter='latex');
ylabel('$|| u_{ex} - u_{num}||_{L^2}$', Interpreter='latex')
xlabel('$h$ [m]', Interpreter='latex')

good_h = fem4.h;


%% Question 3a
TestName = 'TestHW1_3a';
omega = 1;



nEl = 10+1; % Number of points in the mesh [0,L]
omega = 1;

[Data] = NewDataTest(TestName, nEl, omega);
[err0, sol0, fem0, data0] = Pipeline(Data, nEl);


nEl = 100+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err1, sol1, fem1, data1] = Pipeline(Data, nEl);


nEl = 500+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err2, sol2, fem2, data2] = Pipeline(Data, nEl);


nEl = 1e3+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err3, sol3, fem3, data3] = Pipeline(Data, nEl);


nEl = 5e3+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err4, sol4, fem4, data4] = Pipeline(Data, nEl);


nEl = 1e4+1;
[Data] = NewDataTest(TestName, nEl, omega);
[err5, sol5, fem5, data5] = Pipeline(Data, nEl);

%% PLots of the L2-norm error

L2_errors = [err0.L2_err , err1.L2_err , err2.L2_err , err3.L2_err , ...
             err4.L2_err , err5.L2_err];

H1_errors = [err0.H1_err , err1.H1_err , err2.H1_err , err3.H1_err , ...
             err4.H1_err , err5.H1_err];

stepsizes = [fem0.h, fem1.h, fem2.h, fem3.h,fem4.h, fem5.h];

figure()

loglog(stepsizes, L2_errors, '-+r');
grid on;

title('Graph of the L2 error $||u_{ex}-u_{num}||$, as a function of the stepszize $h$', Interpreter='latex');
ylabel('$|| u_{ex} - u_{num}||_{L^2}$', Interpreter='latex')
xlabel('$h$ [m]', Interpreter='latex')

good_h = fem4.h;

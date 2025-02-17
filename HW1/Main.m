clear; close all;
clc;
%% Path

addpath Pipeline
addpath Matrix_Build
addpath Results

% % Define the path to the 'results' folder
resultsFolder = fullfile('HW1', 'Results');
% 
% % Make sure the folder exists (if not, create it)
% if ~exist(resultsFolder, 'dir')
%     mkdir(resultsFolder);
% end


%% Question 2a
TestName = 'TestHW1_2a';
omega = 1;

nEl = 4+1;
[Data] = createData(TestName, nEl, omega);
[err1, sol1, mesh1, data1] = getResults(Data, nEl);


nEl = 50+1;
[Data] = createData(TestName, nEl, omega);
[err2, sol2, mesh2, data2] = getResults(Data, nEl);


nEl = 200+1;
[Data] = createData(TestName, nEl, omega);
[err3, sol3, mesh3, data3] = getResults(Data, nEl);


nEl = 400+1;
[Data] = createData(TestName, nEl, omega);
[err4, sol4, mesh4, data4] = getResults(Data, nEl);


nEl = 800+1;
[Data] = createData(TestName, nEl, omega);
[err5, sol5, mesh5, data5] = getResults(Data, nEl);


%% Question 2a - Results plotting

h_vals = [mesh1.h , mesh2.h , mesh3.h , mesh4.h , mesh5.h];

sol_ex = Data.uex;

L = {};
figure();
plot(data5.x, sol_ex(data5.x), "LineWidth",4, "Color",'r', LineStyle=':'); % Exact solution
hold on;
grid on
L{end+1} =  'exact solution';

% Approximate solutions
plot(data1.x , sol1 , '-x',Color="#7E2F8E");
plot(data2.x , sol2 , 'xm');
plot(data3.x , sol3 , 'xk');
plot(data4.x , sol4 , 'b', LineWidth=2.5);
plot(data5.x , sol5 , 'c', LineWidth=1.5);


for i = 1:length(h_vals)
    L{end+1} =  ['h =  ' num2str(h_vals(i))];
end
title('Comparison of exact and numerical solutions for different mesh sizes.')
legend(L);
xlabel('x [m]');
ylabel('amplitude');

print('Q2a_mesh_size_comp.png', '-r300'); % 300 DPI resolution

%% Question 2a - PLots of the L2-norm error

L2_errors = [err1 , err2 , err3 , err4 , err5];

h_vals = [mesh1.h, mesh2.h, mesh3.h,mesh4.h, mesh5.h];

h_vals_str = sprintf('%.3f, ', h_vals);  % Convert to string with 3 decimals
h_vals_str = h_vals_str(1:end-2);           % Remove last comma and space
titleText = sprintf('Graph of the L2 error ||u_{ex}-u_{num}||, as a function of the stepszize \n h = %s', h_vals_str);


figure();

semilogy(h_vals, L2_errors, '-+r', MarkerSize=8);
grid on;

title(titleText);
ylabel('$log \left( || u_{ex} - u_{num}||_{L^2} \right)$', Interpreter='latex')
xlabel('$h$ [m]', Interpreter='latex');


print('Q2a_L2_error', '-dpng', '-r300'); % 300 DPI resolution

slopes = zeros(length(h_vals) -1,1);
for i = 2 :length(h_vals)
    slopes(i-1) = ( log(L2_errors(1) - L2_errors(i)) ) / log((h_vals(1) - h_vals(i) ));
end
avg_slope = mean(slopes);
% slope =  ( log(L2_errors(1) - L2_errors(5)) ) / log((h_vals(1) - h_vals(5) ))


%% Question 2b - Computing routine
TestName = 'TestHW1_2b';
omega = 1;

nEl = 5+1;
[Data] = createData(TestName, nEl, omega);
[err1, sol1, mesh1, data1] = getResults(Data, nEl);


nEl = 100+1;
[Data] = createData(TestName, nEl, omega);
[err2, sol2, mesh2, data2] = getResults(Data, nEl);


nEl = 200+1;
[Data] = createData(TestName, nEl, omega);
[err3, sol3, mesh3, data3] = getResults(Data, nEl);


nEl = 300+1;
[Data] = createData(TestName, nEl, omega);
[err4, sol4, mesh4, data4] = getResults(Data, nEl);


nEl = 400+1;
[Data] = createData(TestName, nEl, omega);
[err5, sol5, mesh5, data5] = getResults(Data, nEl);


%% Question 2b - Results plotting

h_vals = [mesh1.h , mesh2.h , mesh3.h , mesh4.h , mesh5.h];

sol_ex = Data.uex;

L = {};
figure();
plot(data5.x, sol_ex(data5.x), "LineWidth",4, "Color",'r', LineStyle=':'); % Exact solution
hold on;
grid on
L{end+1} =  'exact solution';

% Approximate solutions
plot(data1.x , sol1 , '-x',Color="#7E2F8E");
plot(data2.x , sol2 , 'xm');
plot(data3.x , sol3 , 'xk');
plot(data4.x , sol4 , 'b', LineWidth=2.5);
plot(data5.x , sol5 , 'c', LineWidth=1.5);


for i = 1:length(h_vals)
    L{end+1} =  ['h =  ' num2str(h_vals(i))];
end
title('Comparison of exact and numerical solutions for different mesh sizes.')
legend(L);
xlabel('x [m]');
ylabel('amplitude');

print('Q2a_mesh_size_comp', '-dpng', '-r300'); % 300 DPI resolution



%% Question 2b - PLots of the L2-norm error

L2_errors = [err1 , err2 , err3 , ...
             err4 , err5];


stepsizes = [mesh1.h, mesh2.h, mesh3.h,mesh4.h, mesh5.h];

figure()

loglog(stepsizes, L2_errors, '-+r');
grid on;

title('Graph of the L2 error $||u_{ex}-u_{num}||$, as a function of the stepszize $h$', Interpreter='latex');
ylabel('$|| u_{ex} - u_{num}||_{L^2}$', Interpreter='latex')
xlabel('$h$ [m]', Interpreter='latex')

good_h = mesh4.h;



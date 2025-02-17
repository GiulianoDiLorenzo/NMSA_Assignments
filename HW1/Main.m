clear; 
close all;
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

% Define the list of nEl values
nElList = [10+1, 30+1, 50+1, 100+1, 150+1];

% Initialize a structure to store results
datas = struct();

% Loop over each nEl value
for i = 1:length(nElList)

    % Generate data
    [Data] = createData(TestName, nElList(i), omega);
    
    % Compute results
    [sol, mesh, data] = getResults(Data, nElList(i));
    % Compute error
    error = Data.uex(Data.x).' - sol;
    L2_err = compute_L2_Error(error, mesh.h);
    
    % Store results in the datas structure
    datas(i).n_pts = nElList(i);  % Store nEl value
    datas(i).err = L2_err;  % Store error
    datas(i).sol = sol;  % Store solution
    datas(i).mesh = mesh; % Store mesh structure
    datas(i).data = data; % Store data structure
end

%% Question 2a - Results plotting

x_axis = datas(end).data.x;
sol_ex = datas(end).data.uex(x_axis);

markOpts = ["-x", "-x", "-x" , "-b" , "-c"];
colorOpts = ["#7E2F8E" , "m" , "k" , "b" , "c"];
lineWidths = [1 , 1 , 1 , 2.5 , 1.5];

L = {};
figure();

plot( x_axis, sol_ex, "LineWidth",6, "Color",'r', LineStyle=':'); % Exact solution
hold on;
grid on
L{end+1} =  'exact solution';

for i = 1:length(datas)
    x_axis = datas(i).data.x;
    sol_num_i = datas(i).sol;
    % disp([' i =', num2str(i), ' size(sol_num) = ', num2str(size(sol_num_i)) ]);

    plot(x_axis , sol_num_i , markOpts(i), 'Color', colorOpts(i), 'LineWidth', lineWidths(i));

    L{end+1} =  ['h =  ' , num2str(datas(i).mesh.h)];
end

title('Comparison of exact and numerical solutions for different mesh sizes.')
legend(L);
xlabel('x [m]');
ylabel('amplitude');

print('Q2a_mesh_size_comp.png', '-r300'); % 300 DPI resolution

%% Question 2a - PLots of the L2-norm error
L2_errors = zeros(length(data),1);
h_vals = zeros(length(data),1);

for i = 1:length(datas)
    L2_errors(i) = datas(i).err;
    h_vals(i) = datas(i).mesh.h;
end

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

slopes = zeros(length(datas) -1,1);
for i = 2 :length(datas)
    slopes(i-1) = ( log(L2_errors(1) - L2_errors(i)) ) / log((h_vals(1) - h_vals(i) ));
end
avg_slope = mean(slopes);
% slope =  ( log(L2_errors(1) - L2_errors(5)) ) / log((h_vals(1) - h_vals(5) ))


%% Question 2b - Computing routine
TestName = 'TestHW1_2b';

omega = 1;
% Define the list of nEl values
nElList = [2+1, 6+1, 12+1, 50+1, 100+1];

% Initialize a structure to store results
datas = struct();

% Loop over each nEl value
for i = 1:length(nElList)

    % Generate data
    [Data] = createData(TestName, nElList(i), omega);
    
    % Compute results
    [sol, mesh, data] = getResults(Data, nElList(i));
    % Compute error
    error = Data.uex(Data.x).' - sol;
    L2_err = compute_L2_Error(error, mesh.h);
    
    % Store results in the datas structure
    datas(i).n_pts = nElList(i);  % Store nEl value
    datas(i).err = L2_err;  % Store error
    datas(i).sol = sol;  % Store solution
    datas(i).mesh = mesh; % Store mesh structure
    datas(i).data = data; % Store data structure
end

%% Question 2b - Results plotting

x_axis = datas(end).data.x;
sol_ex = datas(end).data.uex(x_axis);

markOpts = ["-x", "-x", "-x" , "-b" , "-c"];
colorOpts = ["#7E2F8E" , "m" , "k" , "b" , "c"];
lineWidths = [1 , 1 , 1 , 2.5 , 1.5];

L = {};
figure();

plot( x_axis, sol_ex, "LineWidth",6, "Color",'r', LineStyle=':'); % Exact solution
hold on;
grid on
L{end+1} =  'exact solution';

for i = 1:length(datas)
    x_axis = datas(i).data.x;
    sol_num_i = datas(i).sol;

    plot(x_axis , sol_num_i , markOpts(i), 'Color', colorOpts(i), 'LineWidth', lineWidths(i));

    L{end+1} =  ['h =  ' , num2str(datas(i).mesh.h)];
end

title('Comparison of exact and numerical solutions for different mesh sizes.')
legend(L);
xlabel('x [m]');
ylabel('amplitude');

print('Q2b_mesh_size_comp.png', '-r300'); % 300 DPI resolution


%% Question 2b - PLots of the L2-norm error

L2_errors = zeros(length(data),1);
h_vals = zeros(length(data),1);

for i = 1:length(datas)
    L2_errors(i) = datas(i).err;
    h_vals(i) = datas(i).mesh.h;
end

h_vals_str = sprintf('%.3f, ', h_vals);  % Convert to string with 3 decimals
h_vals_str = h_vals_str(1:end-2);           % Remove last comma and space
titleText = sprintf('Graph of the L2 error ||u_{ex}-u_{num}||, as a function of the stepszize \n h = %s', h_vals_str);

figure();

semilogy(h_vals, L2_errors, '-+r', MarkerSize=8);
grid on;

title(titleText);
ylabel('$log \left( || u_{ex} - u_{num}||_{L^2} \right)$', Interpreter='latex')
xlabel('$h$ [m]', Interpreter='latex');


print('Q2b_L2_error', '-dpng', '-r300'); % 300 DPI resolution

slopes = zeros(length(datas) -1,1);
for i = 2 :length(datas)
    slopes(i-1) = ( log(L2_errors(1) - L2_errors(i)) ) / log((h_vals(1) - h_vals(i) ));
end

avg_slope = mean(slopes);

%% Question 3a
L = 5;
f0 = 3;
N_pts = 1 + 30 * f0 * L;

omegaList = [0.1 , 1 , 2 , 5];

TestName = 'TestHW1_3a';

% Initialize a structure to store results
datas = struct();

for i = 1 : length(omegaList)
    % Generate data
    [Data] = createData(TestName, N_pts , omegaList(i));
    % Compute results
    [sol, mesh, data] = getResults(Data, N_pts);

    % Store results in the datas structure
    datas(i).n_pts = N_pts;  % Store N_pts value
    datas(i).sol = sol;      % Store solution
    datas(i).mesh = mesh;    % Store mesh structure
    datas(i).data = data;    % Store data structure
end


%% Question 3a - Plotting
L = {};
figure();

for i = 1:length(datas)
    x_axis = datas(i).data.x;
    sol_num_i = datas(i).sol;


    disp(size(sol_num_i));

    plot(x_axis , sol_num_i ); % , markOpts(i), 'Color', colorOpts(i), 'LineWidth', lineWidths(i));
    hold on;
    L{end+1} =  ['omega =  ' , num2str(datas(i).data.omega)];
end

title('Comparison of exact and numerical solutions for different omega.')
legend(L);
xlabel('x [m]');
ylabel('amplitude');

print('Q3a_mesh_omega_shift.png', '-r300'); % 300 DPI resolution

%% Question 3b
L = 5;
f0 = 3;
N_pts = 1 + 30 * f0 * L;

omegaList = [0.5 , 1, 1.5, 5 , 20];

TestName = 'TestHW1_3b';

% Initialize a structure to store results
datas = struct();

for i = 1 : length(omegaList)
    % Generate data
    [Data] = createData(TestName, N_pts , omegaList(i));
    % Compute results
    [sol, mesh, data] = getResults(Data, N_pts);

    % Store results in the datas structure
    datas(i).n_pts = N_pts;  % Store nEl value
    datas(i).sol = sol;  % Store solution
    datas(i).mesh = mesh; % Store mesh structure
    datas(i).data = data; % Store data structure
end


%% Quetion 3b - Results Comparison
L = {};
figure();

for i = 1:length(datas)
    x_axis = datas(i).data.x;
    sol_num_i = datas(i).sol;

    plot(x_axis , sol_num_i ); % , markOpts(i), 'Color', colorOpts(i), 'LineWidth', lineWidths(i));
    hold on;
    L{end+1} =  ['\omega =  ' , num2str(datas(i).data.omega)];
end

title('Comparison of exact and numerical solutions for different omega.')
legend(L);
xlabel('x [m]');
ylabel('amplitude');

%% Question 4 
T = 2;
dt = 1e-3;
N_inst = T / dt + 1;
t = linspace(0,T, N_inst );

L = 5;
f0 = 3;
omegaList = [0.5 , 1, 1.5, 5 , 10];

N_pts = 1 + 30 * f0 * L;

v = zeros(N_pts , N_inst, length(omegaList));

for i = 1 : length(omegaList)
    % Generate data
    [Data] = createData(TestName, N_pts , omegaList(i));
    % Compute results
    [sol, mesh, data] = getResults(Data, N_pts);

    % Store results in the datas structure
    datas(i).n_pts = N_pts;  % Store nEl value
    datas(i).sol = sol;  % Store solution
    datas(i).mesh = mesh; % Store mesh structure
    datas(i).data = data; % Store data structure

    v(:,:,i) = exp(1i*omegaList(i)*t) .* datas(i).sol;
end

%% Question 4 - Surface plots

% Create figure
figure;

for i = 1:length(omegaList)
    % Extract the i-th slice of v
    dataSlice = real(v(:,:,i));
    
    % Create surface plot
    subplot(ceil(length(omegaList)/2), 2,i);
    surf(dataSlice);
    
    % Formatting
    title(['Surface Plot for \omega = ', num2str(omegaList(i))]);
    xlabel('t-axis [s]');
    ylabel('x-axis [m]');
    zlabel('Values');
    colorbar; % Add color bar for reference
    shading interp; % Smooth shading
    view(3); % Set 3D view
    axis tight;
    
    % Pause or update figure to visualize each slice
    pause(0.5); % Adjust pause duration as needed
end


clear; 
close all;
clc;


%% Path

addpath Pipeline
addpath Matrix_Build
addpath Results


%% Question 2a 

TestName = 'TestHW1_2a';
omega = 1;

% Define the list of N_pts values
N_pts_list = [10+1, 30+1, 50+1, 100+1, 150+1];

% Question 2a - Numerical Run
for i = 1 : length(N_pts_list)
    [result] = runNumericalSolution(TestName, omega, N_pts_list(i), true);
    results(i) = result;
end

% Question 2a - Results plotting
titleText = ['Comparison of exact and numerical solutions for different' ...
             'mesh sizes'];
plotNumericalSolution(results, titleText, 'N_pts');
print('Q2a_mesh_size_comp.png', '-r300'); % 300 DPI resolution

% Question 2a L2 Error Computation and Plot
[slope, L2_errors, h_vals] = computeAndPlotL2Error(results);
print('Q2a_L2_error', '-dpng', '-r300'); % 300 DPI resolution

%% Question 2b - Computing routine
TestName = 'TestHW1_2b';
omega = 1;

% Define the list of N_pts values
N_pts_list = [10+1, 30+1, 50+1, 100+1, 150+1];

% Question 2b - Numerical Run
for i = 1 : length(N_pts_list)
    [result] = runNumericalSolution(TestName, omega, N_pts_list(i), true);
    results(i) = result;
end

% Question 2b - Results plotting
titleText = {['Comparison of exact and numerical solutions for different ...' ...
              'mesh sizes'], ...
             ['\mu_1 = ', num2str(results(1).data.mu1), ...
              ' \mu_2 = ', num2str(results(1).data.mu2), ...
              ' \rho_1 = ', num2str(results(1).data.rho1), ...
              ' \rho_2 = ', num2str(results(1).data.rho2)]};

plotNumericalSolution(results, titleText, 'N_pts');
print('Q2b_mesh_size_comp.png', '-r300'); % 300 DPI resolution

% Question 2b L2 Error Computation and Plot
[slope, L2_errors, h_vals] = computeAndPlotL2Error(results);
print('Q2b_L2_error', '-dpng', '-r300'); % 300 DPI resolution

%% Question 3a
TestName = 'TestHW1_3a';
L = 5;
f0 = 3;
N_pts = 1 + 30 * f0 * L;

% omegaList = [1 , 5];
omegaList = 1 : 0.5 : 19.5;


% Question 3a - Numerical Run
for i = 1 : length(omegaList)
    [result] = runNumericalSolution(TestName, omegaList(i), N_pts, false);
    resultsOmega(i) = result;
end

% Question 3a - Results plotting
titleText = 'u(x) for sweeping \omega \in ]0,20[';
plotNumericalSolution(resultsOmega, titleText, 'omega');
print('Q2b_mesh_size_comp.png', '-r300'); % 300 DPI resolution

peaks = zeros(length(omegaList) , 1);
phase = zeros(length(omegaList) , 1);


% Question 3a - Omega influence on magnitude and phase
for i = 1 : length(omegaList)
    peaks(i) = max( resultsOmega(i).sol);

    [pks1,locs1] = findpeaks(resultsOmega(i).sol, NPeaks=1);
    [pks2,locs2] = findpeaks( - resultsOmega(i).sol, NPeaks=1);
    scatter(resultsOmega(i).data.x(locs1) , pks1, 'v', 'filled' , 'HandleVisibility', 'off')
    scatter(resultsOmega(i).data.x(locs2) ,  - pks2, 'v', 'filled' , 'HandleVisibility', 'off')
    hold on;
   
    dT_pks = mean(resultsOmega(i).data.x(locs1) - resultsOmega(i).data.x(locs2));
    f = 1 / (2 * mean(resultsOmega(i).data.x(locs1) - resultsOmega(i).data.x(locs2))); % Frequency in Hz
    omega = 2 * pi * f; % Angular frequency

    % Compute phase shift
    comp =  min(locs1, locs2) ;
    phases = zeros(length(omegaList),1);
    if comp == locs1
        sgn = 0;
        loc = locs1;
    else
        sgn = 1;
        loc = locs2;
    end
    phases(i) =  sgn * pi/2 -  omega * resultsOmega(i).data.x(loc); % If first extremum is a peak
end
%%
figure();
subplot(2,1,1);
plot(omegaList, peaks);
xlabel('$\omega$', Interpreter='latex');
ylabel('$max ( u_{num} ( \omega ) )$', Interpreter='latex');
grid on

subplot(2,1,2);
plot(omegaList, phases);
grid on;

xlabel('$\omega$' , Interpreter='latex');
ylabel('$\phi ( u_{num} ( \omega ) )$' , Interpreter='latex');

print('Q3a_omega_influence.png', '-r300'); % 300 DPI resolution

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
t = linspace(0,T, N_inst);

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


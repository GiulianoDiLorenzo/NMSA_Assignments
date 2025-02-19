clear; 
close all;
clc;

%% Path
addpath Matrix_System_Computation
addpath Output_Functions
addpath Matrix_Build


%% Question 2a  - Setup
TestName = 'TestHW1_2a';
omega = 1;

N_pts_list = [10+1, 30+1, 50+1, 100+1, 150+1]; % Define the list of N_pts values

%% Question 2a - Numerical Run
for i = 1 : length(N_pts_list)
    [result2a] = runNumericalSolution(TestName, omega, N_pts_list(i), true);
    results2a(i) = result2a;
end

disp('================================================');
disp('===================== DONE =====================');
disp('================================================'); 

%%  Question 2a - Numerical Solution Results Plots

titleText = ['Comparison of exact and numerical solutions for different' ...
             'mesh sizes'];
saveName = 'Q2a_mesh_size_comp';

plotNumericalSolution(results2a, titleText, saveName, 'N_pts', true);


%% Question 2a - L2 Error Computation and Plot

saveName = 'Q2a_L2_error';
[lin_coeff, L2_errors, h_vals] = plotL2Error(results2a, saveName);


idx = (h_vals > min(h_vals) ) & (h_vals < max(h_vals)); % Choose an appropriate range
p = polyfit( log10(L2_errors(idx)) , log10(h_vals(idx)) , 1 );
slope = p(1);

%% Question 2b - Set up
TestName = 'TestHW1_2b';
omega = 1;

N_pts_list = [50+1, 100+1, 500+1, 1000+1, 1500+1, 3000 + 1]; % Define the list of N_pts values

%% Question 2b - Numerical Run
for i = 1 : length(N_pts_list)
    [result2b] = runNumericalSolution(TestName, omega, N_pts_list(i), true);
    results2b(i) = result2b;
end

disp('================================================');
disp('===================== DONE =====================');
disp('================================================'); 

%% Question 2b - Numerical Solution Results Plots

titleText = {['Comparison of exact and numerical solutions for different' ...
              'mesh sizes'], ...
             ['\mu_1 = ', num2str(results2b(1).data.mu1), ...
              ' \mu_2 = ', num2str(results2b(1).data.mu2), ...
              ' \rho_1 = ', num2str(results2b(1).data.rho1), ...
              ' \rho_2 = ', num2str(results2b(1).data.rho2)]};

saveName = 'Q2b_mesh_size_comp';

plotNumericalSolution(results2b, titleText, saveName , 'N_pts', false);

%% Question 2b -  L2 Error Computation and Plot
saveName = 'Q2b_L2_error';
[lin_coeff, L2_errors, h_vals] = plotL2Error(results2b, saveName);

idx = (h_vals > min(h_vals) ) & (h_vals < max(h_vals)); % Choose an appropriate range
p = polyfit( log10(L2_errors(idx)) , log10(h_vals(idx)) , 1 );
slope = p(1); 


%% Question 3a - Set up
TestName = 'TestHW1_3a';
L = 5;
f0 = 3;

N_pts = 1 + 30 * f0 * L; % Number of points to respect 10 points per wavelength condition
omegaList = 1 : 0.1 : 19.5;


%% Question 3a - Numerical Run
for i = 1 : length(omegaList)
    [result3a] = runNumericalSolution(TestName, omegaList(i), N_pts, false);
    results3a(i) = result3a;
end

disp('================================================');
disp('===================== DONE =====================');
disp('================================================'); 

%% Question 3a - Numerical Solution Plots

titleText = 'u(x) for sweeping \omega \in ]0,20[';

saveName = 'Q3a_u(x)_omega_comp';

plotNumericalSolution(results3a, titleText, saveName , 'omega', false);


%% Question 3a - Omega influence on magnitude and phase
[peaks, phases] = computePeakPhase(omegaList, results3a); 

figure();
sgtitle('Influence of \omega on magnitude and phase');

subplot(2,1,1);
plot(omegaList, peaks);
xlabel('$\omega$', Interpreter='latex');
ylabel('$max ( u_{num} ( \omega ) )$', Interpreter='latex');
title('Magnitude');
grid on

subplot(2,1,2);
plot(omegaList, phases);
grid on;

xlabel('$\omega$' , Interpreter='latex');
ylabel('$\phi ( u_{num} ( \omega ) )$' , Interpreter='latex');
title('Phase');

print('Q3a_omega_influence', '-dpng',  '-r300');

%% Question 3b - Set up
TestName = 'TestHW1_3b';

L = 5;
f0 = 3;

N_pts = 1 + 300 * f0 * L;           % Number of points to respect the 10 points per wavelength condition
omegaList = [0.5 , 1, 1.5, 5 , 20]; % Omega values we will sweep

%% Question 3b - Numerical Run

for i = 1 : length(omegaList)
    [result3b] = runNumericalSolution(TestName, omegaList(i), N_pts, false);
    results3b(i) = result3b;
end

disp('================================================');
disp('===================== DONE =====================');
disp('================================================'); 


%% Question 3b -  Omega influence on magnitude and phase ??

titleText = 'u(x) for sweeping \omega \in ]0,20[';
saveName = 'Q3b_u(x)_omega_comp';

plotNumericalSolution(results3b, titleText, saveName , 'omega', false);

%% Question 4 - Set Up
% TestName = 'TestHW1_4';
T = 2;
dt = 1e-3;
N_inst = T / dt + 1;
t = linspace(0,T, N_inst);

L = 5;
f0 = 3;

omegaList = [1 , 5, 10 , 15];
N_pts = 1 + 300 * f0 * L;

%% Question 4 - Generating v(x,t) = u(x) * exp(i * w * t)

v = zeros(N_pts , N_inst, length(omegaList));

for i = 1 : length(omegaList)
    [result4] = runNumericalSolution(TestName, omegaList(i), N_pts, false);
    results4(i) = result4;

    v(:,:,i) = exp(1i*omegaList(i)*t) .* results4(i).sol;
end

disp('================================================');
disp('===================== DONE =====================');
disp('================================================'); 

%% Question 4 - Surface plots

figure();
for i = 1:length(omegaList)
    
    dataSlice = real(v(:,:,i));      % Extract the i-th slice of v

    [xAx , tAx] = meshgrid(results4(i).data.x , t);
    
    % Create surface plot
    subplot(ceil(length(omegaList)/2), 2,i);
    surf(xAx , tAx , dataSlice');
    
    % Formatting
    title(['Surface Plot for \omega = ', num2str(omegaList(i))]);
    xlabel('x-axis [m]');
    ylabel('t-axis [s]');
    zlabel('Values');
    colorbar; % Add color bar for reference
    shading interp; % Smooth shading
    view(3); % Set 3D view
    axis tight;
end

print('Q4_v(x,t)', '-dpng' , '-r300');


%% Clearing Workspace
clear; 
close all;
clc;

%% Path

addpath Matrix_Build
addpath Matrix_System_Computation
addpath Output_Functions

%% Question 2a - Setup
TestName = 'TestHW1_2a';
omega = 1;

N_pts_list = [50+1, 100+1, 150 + 1 , 200 + 1, 500+1]; % Define the list of N_pts values

%% Question 2a - Numerical Run
mu_vals =  [1 , 1];
rho_vals = [1 , 1];

for i = 1 : length(N_pts_list)
    [result2a] = runNumericalSolution(TestName, omega, N_pts_list(i), true, mu_vals, rho_vals);
    results2a(i) = result2a;
end

fprintf('\n================================================');
fprintf('\n===================== DONE =====================');
fprintf('\n================================================\n'); 

%%  Question 2a - Numerical Solution Results Plots

titleText = ['Comparison of exact and numerical solutions for different ' ...
             'mesh sizes'];
saveName = 'Plots/Q2a_mesh_size_comp';

plotNumericalSolution(results2a, titleText, saveName, 'N_pts', true);

%% Question 2a - L2 Error Computation and Plot

saveName = 'Plots/Q2a_L2_error';
[lin_coeff, L2_errors, h_vals] = plotL2Error(results2a, saveName);


idx = (h_vals > min(h_vals) ) & (h_vals < max(h_vals)); % Choose an appropriate range
p = polyfit( log10(L2_errors(idx)) , log10(h_vals(idx)) , 3 );
slope = p(1);

%% Question 2b - Set up
TestName = 'TestHW1_2b';
omega = 1;

% N_pts_list = [50+1, 100+1, 200+1, 500+1, 1000+1, 1500+1, 3000 + 1]; % Define the list of N_pts values

N_pts = 500 + 1;
%% Question 2b - Numerical Run for mu sweep
mu_vals =  [ 1 , 1; 
             1 , 5 ; 
             1 , 10 ; 
             1 , 15 ;  
             1 , 20 ];
rho_vals = [1 , 1];

for i = 1 : length(mu_vals)
    [result2b_mu] = runNumericalSolution(TestName, omega, N_pts, true, mu_vals(i,:), rho_vals);
    results2b_mu(i) = result2b_mu;
end

fprintf('\n================================================');
fprintf('\n===================== DONE =====================');
fprintf('\n================================================\n'); 


%% Question 2b - Numerical Solution Results Plots for mu sweep

titleText = ['Studying influence of $\mu$ with $N_{pts}$ = ' , num2str(N_pts) , ' and mu_1 = ' , num2str(results2b_mu(1).data.mu1)];

saveName = 'Plots/Q2b_mu_comp';

plotNumericalSolution(results2b_mu, titleText, saveName , 'mu', true);

%% Question 2b - Numerical Run for rho sweep
rho_vals =  [ 1 , 1; 
              1 , 50 ; 
              1 , 100 ; 
              1 , 150 ;  
              1 , 200 ];
mu_vals = [1 , 1];

for i = 1 : length(rho_vals)
    [result2b_rho] = runNumericalSolution(TestName, omega, N_pts, true, mu_vals, rho_vals(i,:));
    results2b_rho(i) = result2b_rho;
end

fprintf('\n================================================');
fprintf('\n===================== DONE =====================');
fprintf('\n================================================\n'); 

%% Question 2b - Numerical Solution Results Plots for rho sweep

titleText = ['Studying influence of \rho, with N_{pts} = ' , num2str(N_pts) , ' and \rho_1 = ' , num2str(results2b_rho(1).data.rho1)];

saveName = 'Plots/Q2b_rho_comp';

plotNumericalSolution(results2b_rho, titleText, saveName , 'rho', true);

% %% Question 2b -  L2 Error Computation and Plot
% saveName = 'Q2b_L2_error';
% [lin_coeff, L2_errors, h_vals] = plotL2Error(results2b_mu, saveName);
% 
% idx = (h_vals > min(h_vals) ) & (h_vals < max(h_vals)); % Choose an appropriate range
% p = polyfit( log10(L2_errors(idx)) , log10(h_vals(idx)) , 1 );
% slope = p(1); 

%% Question 3a - Set up
TestName = 'TestHW1_3a';
L = 5;
f0 = 3;

N_pts = 1 + 30 * f0 * L; % Number of points to respect 10 points per wavelength condition

omegaList = 1 : 0.1 : 19.5; % PUT PLOTME = TRUE FOR computePeakPhase
% omegaList =  [1 , 2 , 3 , 4 , 5, 10 , 12 , 15 , 20];    % PUT PLOTME = FALSE FOR computePeakPhase

mu_vals =  [1 , 1];
rho_vals = [1 , 1];



%% Question 3a - Numerical Run
for i = 1 : length(omegaList)
    [result3a] = runNumericalSolution(TestName, omegaList(i), N_pts, false , 1 , 1);
    results3a(i) = result3a;
end

fprintf('\n================================================');
fprintf('\n===================== DONE =====================');
fprintf('\n================================================\n'); 

%% Question 3a - Numerical Solution Plots

titleText = ['$u_{num}(x)$ with $N_{pts} = $', num2str(N_pts), ' for sweeping $\omega \in ]0,20[$'];

saveName = 'Plots/Q3a_u(x)_omega_comp';

plotNumericalSolution(results3a, titleText, saveName , 'omega', false);


%% Question 3a - Computing peaks and phases for each omega with plots

saveName = 'Plots/Q3a_omega_influence_dB';
[peaks, locs, phases] = computePeakPhase(omegaList, results3a, saveName , true); 


%% Question 3b - Set up
TestName = 'TestHW1_3b';

L = 5;
f0 = 3;

N_pts = 1 + 300 * f0 * L;           % Number of points to respect the 10 points per wavelength condition

omegaList = 1 : 0.1 : 19.5; % PUT PLOTME = TRUE FOR computePeakPhase
% omegaList =  [1 , 2 , 3 , 4 , 5, 10 , 12 , 15 , 20];    % PUT PLOTME = FALSE FOR computePeakPhase

mu_vals =  [1 , 1];
rho_vals = [1 , 1];

%% Question 3b - Numerical Run

for i = 1 : length(omegaList)
    fprintf(['omega = ' , num2str(omegaList(i)) , '\n']);
    [result3b] = runNumericalSolution(TestName, omegaList(i), N_pts, false, mu_vals , rho_vals);
    results3b(i) = result3b;
end

fprintf('\n================================================');
fprintf('\n===================== DONE =====================');
fprintf('\n================================================\n'); 


%% Question 3b -  Omega influence on magnitude and phase ??

titleText = 'u(x) for sweeping \omega \in ]0,20[';
saveName = 'Plots/Q3b_u(x)_omega_comp';

plotNumericalSolution(results3b, titleText, saveName , 'omega', false);

%% Question 3b - Computing peaks and phases for each omega

saveName = 'Plots/Q3b_omega_influence_dB';
[peaks, locs, phases] = computePeakPhase(omegaList, results3b, saveName , true); 


%% Question 3b - Plotting gains and phases as functions of omega
figure();
sgtitle('Influence of \omega on magnitude and phase, piece-wise constant velocity');

subplot(2,1,1);
plot(omegaList, db(abs(peaks)),  LineWidth=2);
xlabel('$\omega$', Interpreter='latex');
ylabel('$max | u_{num} ( \omega ) | [dB]$', Interpreter='latex');
title('Magnitude');
xline(omegaList(locs))
grid on

subplot(2,1,2);
plot(omegaList, phases, LineWidth=2);
hold on;
yline(pi/2 , LineStyle="--");
grid on;

xlabel('$\omega$' , Interpreter='latex');
ylabel('$\angle u_{num} ( \omega )$ [rad]' , Interpreter='latex');
title('Phase');

print('Plots/Q3b_omega_influence_dB', '-dpng',  '-r300');

%% Question 4 - Set Up
TestName = 'TestHW1_2a';


L = 5;
f0 = 3;

omegaList = [1.2 , 2.6 , 4.7, 10 , 15];
N_pts = 1 + 300 * f0 * L;

mu_vals =  [1 , 1];
rho_vals = [1 , 1];

T = 2*pi/omegaList(1);
dt = 1e-3;
N_inst = round(T / dt) + 1;
t = linspace(0,T, N_inst);

%% Question 4 - Generating v(x,t) = u(x) * exp(i * w * t)

v = zeros(N_pts , N_inst, length(omegaList));


for i = 1 : length(omegaList)
    

    [result4] = runNumericalSolution(TestName, omegaList(i), N_pts, false, mu_vals, rho_vals);
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

print('Plots/Q4_v(x,t)', '-dpng' , '-r300');


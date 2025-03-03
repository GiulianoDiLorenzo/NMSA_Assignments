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
L = 1;

N_pts_list = [4+1, 10 + 1 , 100 + 1, 500+1, 1000+1]; % Define the list of N_pts values

%% Question 2a - Numerical Run
mu_vals =  [1 , 1];
rho_vals = [1 , 1];

for i = 1 : length(N_pts_list)
    [result2a] = runNumericalSolution(TestName, L, omega, N_pts_list(i), true, mu_vals, rho_vals);
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
[complexity_2a, L2_errors_2a, h_vals] = plotL2Error(results2a, saveName, 'N_pts');


%% Question 2b - Set up
TestName = 'TestHW1_2b';
omega = 1;
L = 1;
N_pts = 5000 + 1;   % Extra fine mesh

%% Question 2b - mu sweep
% =======================================
% ============ NUMERICAL RUN ============
% =======================================

mu_vals =  [ 1 , 1; 
             1 , 5 ; 
             1 , 10 ; 
             1 , 15 ;  
             1 , 20 ];
rho_vals = [1 , 1];

for i = 1 : length(mu_vals)
    [result2b_mu] = runNumericalSolution(TestName, L, omega, N_pts, true, mu_vals(i,:), rho_vals);
    results2b_mu(i) = result2b_mu;
end

fprintf('\n================================================');
fprintf('\n===================== DONE =====================');
fprintf('\n================================================\n'); 

%% Question 2b - mu sweep
% =======================================
% ====== NUMERICAL SOLUTION PLOTS =======
% =======================================

titleText = ['Studying influence of $\mu$ with $N_{pts} = $' , num2str(N_pts) , ' and $\mu_1 = $' , num2str( results2b_mu(1).data.mu1 )];

saveName = 'Plots/Q2b_mu_comp';

plotNumericalSolution(results2b_mu, titleText, saveName , 'mu', true);

%% Question 2b - mu sweep
% =======================================
% ==== L2 ERROR COMPUTATION AND PLOT ====
% =======================================

saveName = 'Plots/Q2b_L2_error_mu_sweep';
[complexity_2b_mu, L2_errors_2b_mu, mu_diff_vals] = plotL2Error(results2b_mu, saveName, 'mu');


%% Question 2b - rho sweep
% =======================================
% ============ NUMERICAL RUN ============
% =======================================

rho_vals =  [ 1 , 0.1; 
              1 , 1 ; 
              1 , 10 ; 
              1 , 50 ;  
              1 , 80 ];
mu_vals = [1 , 1];

for i = 1 : length(rho_vals)
    [result2b_rho] = runNumericalSolution(TestName, L , omega, N_pts, true, mu_vals, rho_vals(i,:));
    results2b_rho(i) = result2b_rho;
end

fprintf('\n================================================');
fprintf('\n===================== DONE =====================');
fprintf('\n================================================\n'); 

%% Question 2b - rho sweep
% =======================================
% ====== NUMERICAL SOLUTION PLOTS =======
% =======================================

titleText = ['Studying influence of $\rho$, with $N_{pts} = $' , num2str(N_pts) , ' and $\rho_1 = $' , num2str(results2b_rho(1).data.rho1)];

saveName = 'Plots/Q2b_rho_comp';
plotNumericalSolution(results2b_rho, titleText, saveName , 'rho', true);

%% Question 2b -  rho sweep
% =======================================
% ==== L2 ERROR COMPUTATION AND PLOT ====
% =======================================

saveName = 'Plots/Q2b_L2_error_rho_sweep';
[complexity_2b_rho, L2_errors_2b_rho, rho_diff_vals] = plotL2Error(results2b_rho, saveName, 'rho');


%% Question 3a - Set up
TestName = 'TestHW1_3a';
L = 5;
f0 = 3;

N_pts = 1 + 30 * f0 * L; % Number of points to respect 10 points per wavelength condition

omegaList = 1 : 0.1 : 19.5; % PUT PLOTME = TRUE FOR computePeakPhase
% omegaList =  [1 , 2 , 3 , 4 , 5, 10 , 12 , 15 , 19];    % PUT PLOTME = FALSE FOR computePeakPhase

mu_vals =  [1 , 1];
rho_vals = [1 , 1];



%% Question 3a - Numerical Run
for i = 1 : length(omegaList)
    [result3a] = runNumericalSolution(TestName, L , omegaList(i), N_pts, false , mu_vals , rho_vals);
    results3a(i) = result3a;
end

fprintf('\n================================================');
fprintf('\n===================== DONE =====================');
fprintf('\n================================================\n'); 

%% Question 3a - Numerical Solution Plots

titleText = ['Numerical solution $u_{num}(x)$ with $N_{pts} = $', num2str(N_pts), ' and $\omega \in ]0,20[$'];

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

omegaList = 1 : 0.1 : 19.9; % PUT PLOTME = TRUE FOR computePeakPhase
% omegaList =  [1 , 5, 10, 19];    % PUT PLOTME = FALSE FOR computePeakPhase

mu_vals =  [1 , 1];
rho_vals = [1 , 1];

%% Question 3b - Numerical Run

for i = 1 : length(omegaList)
    fprintf(['omega = ' , num2str(omegaList(i)) , '\n']);
    [result3b] = runNumericalSolution(TestName, L, omegaList(i), N_pts, false, mu_vals , rho_vals);
    results3b(i) = result3b;
end

fprintf('\n================================================');
fprintf('\n===================== DONE =====================');
fprintf('\n================================================\n'); 


%% Question 3b -  Omega influence on magnitude and phase ??

titleText = ['Numerical solution $u_{num}(x)$ with $N_{pts}$ = ', num2str(N_pts), ' and $\omega \in ]0,20[$'];

saveName = 'Plots/Q3b_u(x)_omega_comp_4_values';

plotNumericalSolution(results3b, titleText, saveName , 'omega', false);



%% Question 3b - Computing peaks and phases for each omega

saveName = ['Plots/Q3b_omega_influence_dB'];
[peaks, locs, phases] = computePeakPhase(omegaList, results3b, saveName , true); 


%% Question 4 - Set Up
TestName = 'TestHW1_3b';

L = 5;
f0 = 3;

omegaList = [1.2 , 2.6 , 4.7, 10 , 15];
N_pts = 1 + 300 * f0 * L;

mu_vals =  [1 , 1];
rho_vals = [1 , 1];

dt = 1e-3;

%% Question 4 - Computing and Plotting results

figure();
for i = 1 : length(omegaList)
    T = 2*pi/omegaList(i);
    N_inst = round(T / dt) + 1;
    t = linspace(0,T, N_inst);

    v = zeros(N_pts , N_inst, length(omegaList));

    [result_i] = runNumericalSolution(TestName, L, omegaList(i), N_pts, false, mu_vals, rho_vals);
    v(:,:,i) = exp(1i*omegaList(i)*t) .* result_i.sol;

    [x_axis , t_axis] = meshgrid(result_i.data.x , t);

    subplot(round(length(omegaList)/2), floor(length(omegaList)/2),i);
    surf(x_axis , t_axis , real(v(:,:,i))');
    
    title(['$\omega = $', num2str(omegaList(i)), ', T = ', num2str(T), ' s'], 'Interpreter','latex');
    xlabel('x-axis [m]', Interpreter='latex');
    ylabel('t-axis [s]', Interpreter='latex');
    zlabel('v(x,t)', Interpreter='latex');
    colorbar; 
    shading interp; 
    view(3); % Set 3D view
    axis tight;
end
sgtitle('$v(x,t)$ for different values of $\omega$', 'Interpreter', 'latex');



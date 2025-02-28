function plotNumericalSolution(results, titleText, saveName , sweepType, isExSol)
% ========================================================================
%   OUTPUT : Plots for the discrete solution
%
%   INPUTS : 
%       - results   --> Structure with all data : problem formulation and 
%                       solution
%       - titleText --> String for the title of the graph
%       - saveName  --> String for saving the graph
%       - sweeType  --> String to indicate which sweep we perform
%       - isExSol   --> Boolean to also plot the exact solution (==True)
% ========================================================================

L = {};
figure();

if ( isExSol )
    x_axis = results(end).data.x;
    sol_ex = results(end).data.uex(x_axis);
    plot( x_axis, sol_ex, "LineWidth",3.5, "Color",'r', LineStyle='-'); % Exact solution
    L{end+1} =  'exact solution';
end

if strcmp(sweepType , 'N_pts')
   markOpts = ["-o", "-o", "-x" , "-" , "--" ];
   colorOpts = [ "k" , "#7E2F8E" , "g" , "c" , "b" ];
   lineWidths = [1 , 1 , 1 , 1 , 1.5 ];
end

hold on; grid on

for i = 1:length(results)
    x_axis = results(i).data.x;
    sol_num_i = results(i).sol;
    % disp([' i =', num2str(i), ' size(sol_num) = ', num2str(size(sol_num_i)) ]);

    if strcmp(sweepType , 'N_pts')    
        L{end+1} =  ['h =  ' , num2str(results(i).mesh.h)];
        plot(x_axis , sol_num_i , markOpts(i), 'Color', colorOpts(i), 'LineWidth', lineWidths(i));

    elseif strcmp(sweepType , 'omega')
        L{end+1} = ['\omega = ' , num2str(results(i).data.omega)];
        plot(x_axis , sol_num_i);

    elseif strcmp(sweepType , 'mu')
        L{end+1} = ['\mu_2 = ' , num2str(results(i).data.mu2)];
        plot(x_axis , sol_num_i, 'LineWidth', 2);

    elseif strcmp(sweepType , 'rho')
        L{end+1} = ['\rho_2 = ' , num2str(results(i).data.rho2)];
        plot(x_axis , sol_num_i, 'LineWidth', 2);
    end
end


title(titleText, Interpreter="latex");
legend(L);
xlabel('x [m]');
ylabel('amplitude');

print(saveName,'-dpng', '-r300');
end
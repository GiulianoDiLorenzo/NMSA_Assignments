function [complexity, L2_errors, h_vals] = plotL2Error(results, saveName)
% ========================================================================
%   OUTPUT : 
%       - lin_coeff --> Coefficient of the slope
%       - L2_errors --> Array containing the L2 error of each run
%       - h_vals    --> Array containing the mesh size for each run
%
%   INPUTS : 
%       - results   --> Structure with all data : problem formulation and 
%                       solution
%       - saveName  --> String for saving the graph
% ========================================================================

    L2_errors = zeros(length(results),1);
    h_vals = zeros(length(results),1);

    for i = 1:length(results)
        L2_errors(i) = results(i).err;
        h_vals(i) = results(i).mesh.h;
    end
    
    h_vals_str = sprintf('%.3f, ', h_vals);  % Convert to string with 3 decimals
    h_vals_str = h_vals_str(1:end-2);           % Remove last comma and space
    titleText = sprintf(['graph of the error $||u_{ex}-u_{num}||_{L^2}$ ' ...
                        ' as a function of the stepszize \n h = %s'], h_vals_str);
    

    figure();
    % semilogy(h_vals, L2_errors, '-+r', MarkerSize=8);
    hold on ;
    plot(log10(h_vals), log10(L2_errors), '-ob', MarkerSize=8);

    
    grid on;
    title(titleText, Interpreter="latex");
    ylabel('$\log_{10} ( || u_{ex} - u_{num}||_{L^2} ) $', Interpreter='latex')
    xlabel('$\log_{10} (h)$ [m]', Interpreter='latex');
    
    slopes = zeros(length(results) -1,1);
    slopes = ( log10(L2_errors(1)) - log10(L2_errors(2:end)) ) ./ ( log10(h_vals(1)) - log10(h_vals(2:end)) );

    complexity = (mean(slopes));

    print(saveName,'-dpng', '-r300');

end
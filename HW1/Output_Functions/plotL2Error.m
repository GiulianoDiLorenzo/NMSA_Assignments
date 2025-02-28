function [complexity, L2_errors, param_vals] = plotL2Error(results, saveName, sweepType)
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
    slopes = zeros(length(results) -1,1);
    param_vals = zeros(length(results),1);
    

    for i = 1:length(results)
            L2_errors(i) = results(i).err;
    end
    
    figure();
    hold on ;

    if strcmp(sweepType , 'N_pts')
        h_vals = zeros(length(results),1);
        for i = 1:length(results)
            h_vals(i) = results(i).mesh.h;
            
        end
        h_vals_str = sprintf('%.3f, ', h_vals);     % Convert to string with 3 decimals
        h_vals_str = h_vals_str(1:end-2);           % Remove last comma and space
        titleText = sprintf(['Graph of the log of the error $||u_{ex}-u_{num}||_{L^2}$ ' ...
                             ' as a function of the stepszize \n with $h = %s $'], h_vals_str);
        plot(log10(h_vals), log10(L2_errors), '-ob', MarkerSize=8);

        ylabel('$\log_{10} ( || u_{ex} - u_{num}||_{L^2} ) $', Interpreter='latex')
        xlabel('$\log_{10} (h)$ [m]', Interpreter='latex');

        slopes = ( log10(L2_errors(1)) - log10(L2_errors(2:end)) ) ./ ( log10(h_vals(1)) - log10(h_vals(2:end)) );
        param_vals = h_vals;

    elseif strcmp(sweepType , 'mu')
        mu_2_vals = zeros(length(results),1);

        for i = 1:length(results)
            mu_2_vals(i) = results(i).data.mu2;
        end

        slopes = ( log10(L2_errors(end)) - log10(L2_errors(1:end-1)) ) ./ ...
                 ( log10(mu_2_vals(end)) - log10(mu_2_vals(1:end-1)) );

        mu_2_vals_str = sprintf('%.2f, ', mu_2_vals);  % Convert to string with 3 decimals
        mu_2_vals_str = mu_2_vals_str(1:end-2);        % Remove last comma and space
        % 
        % titleText = sprintf(['graph of the error $||u_{ex} - u_{num}||_{L^2} $' ...
        %                      ' as a function of the difference in $ \\mu $ values,'...
        %                      '\n $ \\mu_1 = %d $ and $ \\mu_2 = $ %s'], results(i).data.mu1, mu_2_vals_str);

        titleText = ['Graph of the log of the error $||u_{ex} - u_{num}||_{L^2} $' ...
                     ' as a function of the difference in $ \\mu $ values,' ...
                     '\n $ \\mu_1 = ' num2str(results(i).data.mu1) ...
                     ' $ and $ \\mu_2 = ' mu_2_vals_str ' $'];

        plot(log10(mu_2_vals), log10(L2_errors), '-ob', MarkerSize=8);
        ylabel('$\log_{10} ( || u_{ex} - u_{num}||_{L^2} )$', Interpreter='latex')
        xlabel('$\log_{10} ( \mu_2 ) $', Interpreter='latex');

        param_vals = mu_2_vals;

    elseif strcmp(sweepType , 'rho')
        rho_2_vals = zeros(length(results),1);
        for i = 1:length(results)
            rho_2_vals(i) = results(i).data.rho2;
        end

        rho_2_vals_str = sprintf('%.2f, ', rho_2_vals);     % Convert to string with 3 decimals
        rho_2_vals_str = rho_2_vals_str(1:end-2);               % Remove last comma and space

        titleText = ['Graph of the log of the error $||u_{ex} - u_{num}||_{L^2} $' ...
                     ' as a function of the difference in $ \\rho $ values,' ...
                     '\n $ \\rho_1 = ' num2str(results(i).data.rho1) ...
                     ' $ and $ \\rho_2 = ' rho_2_vals_str ' $'];

        plot(log10(rho_2_vals), log10(L2_errors), '-ob', MarkerSize=8);
        ylabel('$\log_{10} ( || u_{ex} - u_{num}||_{L^2} )$', Interpreter='latex')
        xlabel('$\log_{10} ( \rho_2) $', Interpreter='latex');
        
        slopes = ( log10(L2_errors(end)) - log10(L2_errors(1:end-1)) ) ./...
                 ( log10(rho_2_vals(end)) - log10(rho_2_vals(1:end-1)) );

        param_vals = rho_2_vals;

    else
        error('Unrecognized sweep type');
    end

    grid on;

    complexity = (mean(slopes));

    titleText = sprintf([titleText, '\n Computed slope $p = %.5f $'], complexity);
    title(titleText, Interpreter="latex");

    print(saveName,'-dpng', '-r300');

    

end


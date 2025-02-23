function [peaks, locs, phases] = computePeakPhase(omegaList, results, saveName, plotMe)
% ========================================================================
%   OUTPUT : 
%       - peaks      --> Maximum amplitude of the solution, array of size
%                        (num_omega , 1)
%       - phases     -->  Phase of the solution, array of size
%
%   INPUTS : 
%       - omegaList  --> List of omega values to try
%       - results    --> Structure with all data : problem formulation and 
%                        solution
% ========================================================================


    peaks = zeros(length(omegaList), 1);
    phases = zeros(length(omegaList), 1);
    locs = zeros(length(omegaList), 1);

    for i = 1 : length(omegaList)
        peaks(i) = max( (results(i).sol) );
        phases(i) = asin(results(i).sol(1) / peaks(i) );
    end

    [~, locs] = findpeaks(db(abs(peaks)));


if plotMe

    figure();
    sgtitle('Influence of \omega on magnitude and phase, constant velocity');
    
    subplot(2,1,1);
    
    plot(omegaList, db(abs(peaks)),  LineWidth=2);
    hold on;
    xline(omegaList(locs), 'r--', 'LineWidth', 2); 
    
    xlabel('$\omega$', Interpreter='latex');
    ylabel('$max | u_{num} ( \omega ) | [dB]$', Interpreter='latex');
    title('Magnitude');
    grid on
    
    subplot(2,1,2);
    
    plot(omegaList, phases, LineWidth=2);
    hold on;
    xline(omegaList(locs), 'r--', 'LineWidth', 2)
    
    grid on;
    
    xlabel('$\omega$' , Interpreter='latex');
    ylabel('$\angle u_{num} ( \omega )$ [rad]' , Interpreter='latex');
    title('Phase');
    
    print(saveName, '-dpng',  '-r300');
end


end

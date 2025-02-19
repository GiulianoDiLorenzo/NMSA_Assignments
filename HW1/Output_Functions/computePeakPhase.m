function [peaks, phases] = computePeakPhase(omegaList, results)

peaks = zeros(length(omegaList), 1);
phases = zeros(length(omegaList), 1);


for i = 1 : length(omegaList)
    peaks(i) = max( results(i).sol);

    [pks1,locs1] = findpeaks(results(i).sol, NPeaks=1);
    [pks2,locs2] = findpeaks( - results(i).sol, NPeaks=1);
    scatter(results(i).data.x(locs1) , pks1, 'v', 'filled' , 'HandleVisibility', 'off')
    scatter(results(i).data.x(locs2) ,  - pks2, 'v', 'filled' , 'HandleVisibility', 'off')
    hold on;
   
    dT_pks = mean(results(i).data.x(locs1) - results(i).data.x(locs2));
    f = 1 / (2 * mean(results(i).data.x(locs1) - results(i).data.x(locs2))); % Frequency in Hz
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
    phases(i) =  sgn * pi/2 -  omega * results(i).data.x(loc); % If first extremum is a peak
end


end
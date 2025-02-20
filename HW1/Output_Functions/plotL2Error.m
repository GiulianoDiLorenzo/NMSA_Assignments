function [lin_coeff, L2_errors, h_vals] = plotL2Error(datas, saveName)

    L2_errors = zeros(length(datas),1);
    h_vals = zeros(length(datas),1);

    for i = 1:length(datas)
        L2_errors(i) = datas(i).err;
        h_vals(i) = datas(i).mesh.h;
    end
    
    h_vals_str = sprintf('%.3f, ', h_vals);  % Convert to string with 3 decimals
    h_vals_str = h_vals_str(1:end-2);           % Remove last comma and space
    titleText = sprintf(['Log scale graph of the error $||u_{ex}-u_{num}||_{L^2}$' ...
                        ' as a function of the stepszize \n h = %s'], h_vals_str);
    

    figure();
    % semilogy(h_vals, L2_errors, '-+r', MarkerSize=8);
    hold on ;
    semilogy(h_vals, L2_errors, '-ob', MarkerSize=8);

    
    grid on;
    title(titleText, Interpreter="latex");
    ylabel('$|| u_{ex} - u_{num}||_{L^2}$', Interpreter='latex')
    xlabel('$h$ [m]', Interpreter='latex');
    
    slopes = zeros(length(datas) -1,1);
    slopes2 = zeros(length(datas) -1,1);
    slopesEx = ( L2_errors(1:end-1) - L2_errors(2:end))  ./ ( h_vals(1:end-1) - h_vals(2:end));

    for i = 2 :length(datas)
         slopes(i-1) = ( log10(L2_errors(1)) - log10(L2_errors(i)) ) / ( log10(h_vals(1)) - log10(h_vals(i)) ) ;
         slopes3(i-1) = ( log10(L2_errors(1)) - log10(L2_errors(i)) ) / ( h_vals(1) - h_vals(i)) ;
         slopes2(i-1) = ( L2_errors(1) - L2_errors(i) )  / ( h_vals(1) - h_vals(i) ) ;
     end

    lin_coeff = log10(mean(slopes));

    print(saveName,'-dpng', '-r300');

end
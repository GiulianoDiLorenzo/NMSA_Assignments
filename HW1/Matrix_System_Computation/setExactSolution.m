function [uex , uex_x, uex_xx] = setExactSolution(Data)
    % Define options
    options = {'2 sinus of frequency 2*pi , 8*pi' , '2 sinus of frequency 2*pi , 4*pi', 'simple sine' };
    
    % Create selection dialog
    [selection, ok] = listdlg('PromptString', 'Select an option:', ...
                              'SelectionMode', 'single', ...
                              'ListString', options);
   
    % Check if user made a selection
    if ok
        selectedOption = options{selection};
        fprintf('You selected: %s\n', selectedOption);
        
        % Set variables based on selection
        switch selection
            case 1
                L = Data.L
                uex = @(x)  (x >= 0  & x <= L/2) .*  sin(2*pi*x) ...
                          + (x >= L/2 & x <= L)   .*  sin(8*pi*x);

                uex_x = @(x) (x >= 0  & x <= L/2) .*  2*pi .* cos(2*pi*x) ...
                           + (x > L/2 & x <= L)   .*  8*pi .* cos(8*pi*x);

                uex_xx = @(x) (x >= 0  & x <= L/2) .* (-(2*pi)^2) .* sin(2*pi*x) ...
                            + (x > L/2 & x <= L)   .* (-(8*pi)^2)  .* sin(8*pi*x);
        
            case 2
                L = Data.L
                uex = @(x)  (x >= 0  & x <= L/2) .*  sin(2*pi*x) ...
                          + (x > L/2 & x <= L)   .*  sin(4*pi*x);

                uex_x = @(x) (x >= 0  & x <= L/2) .*  2*pi .* cos(2*pi*x) ...
                           + (x > L/2 & x <= L)   .*  4*pi .* cos(4*pi*x);

                uex_xx = @(x) (x >= 0  & x <= L/2) .* (-(2*pi)^2) .* sin(2*pi*x) ...
                            + (x > L/2 & x <= L)   .* (-(4*pi)^2)  .* sin(4*pi*x);
            case 3
                L = Data.L
                uex = @(x)  (x >= 0  & x <= L) .*  sin(2*pi*x);

                uex_x = @(x) (x >= 0  & x <= L) .*  2*pi .* cos(2*pi*x); 

                uex_xx = @(x) (x >= 0  & x <= L) .* (-(2*pi)^2) .* sin(2*pi*x) ;
        end
        
    else
        selectedOption = 'No selection';
        disp('No option selected. Exiting...');
    end
end
function [scenario] = setScenario(rho_max, x, Nx)

    scenario = struct();
    options = {'Traffic jam' , 'Green light', 'Traffic flow' };
        
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
                    scenario.name = "Traffic jam";
                    scenario.rho_L = 0.8 * rho_max;
                    scenario.rho_R = rho_max;
                    scenario.rho_0 = zeros(Nx, 1);
                    scenario.rho_0(x < 0.5) = scenario.rho_L;
                    scenario.rho_0(x >= 0.5) = scenario.rho_R;
                case 2
                    scenario.name = "Green light";
                    scenario.rho_L = 0.7 * rho_max;
                    scenario.rho_R = 0.2 * rho_max;
                    scenario.rho_0 = zeros(Nx, 1);
                    scenario.rho_0(x < 0.5) = scenario.rho_L;
                    scenario.rho_0(x >= 0.5) = scenario.rho_R;
                case 3
                    scenario.name = "Traffic flow";
                    scenario.rho_L = rho_max;
                    scenario.rho_R = 0.5 * rho_max;
                    scenario.rho_0 = zeros(Nx, 1);
                    scenario.rho_0(x < 0.5) = scenario.rho_L;
                    scenario.rho_0(x >= 0.5) = scenario.rho_R;
            end
            
        else
            selectedOption = 'No selection';
            disp('No option selected. Exiting...');
        end

end
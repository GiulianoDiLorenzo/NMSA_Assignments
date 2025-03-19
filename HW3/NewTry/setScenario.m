function [scenario] = setScenario(sce, rho_max, rho_c, x, Nx)

    scenario = struct();

    if strcmp(sce, 'Traffic jam')
        scenario.name = "Traffic jam";
        % scenario.rho_L = 0.2 * rho_max;
        % scenario.rho_R = rho_max;
        % scenario.rho_0 = zeros(Nx, 1);
        % scenario.rho_0(x < 0.5) = scenario.rho_L;
        % scenario.rho_0(x >= 0.5) = scenario.rho_R;
        % 
        scenario.rho_L = 0.8 * rho_max;
        scenario.rho_R = rho_max;
        scenario.rho_0 = zeros(Nx, 1);
        scenario.rho_0(x < 0.5) = scenario.rho_L;
        scenario.rho_0(x >= 0.5) = scenario.rho_R;
    elseif strcmp(sce, 'Green light')
        scenario.name = "Green light";
        % scenario.rho_L = rho_max;
        % scenario.rho_R = 0.7 * rho_max;
        % scenario.rho_0 = zeros(Nx, 1);
        % scenario.rho_0(x < 0.5) = scenario.rho_L;
        % scenario.rho_0(x >= 0.5) = scenario.rho_R;
        scenario.rho_L = 0.9 * rho_max;
        scenario.rho_R = 0.6 * rho_max;
        scenario.rho_0 = zeros(Nx, 1);
        scenario.rho_0(x < 0.5) = scenario.rho_L;
        scenario.rho_0(x >= 0.5) = scenario.rho_R;


    elseif strcmp(sce, 'Traffic flow')
        scenario.name = "Traffic flow";

        scenario.rho_L = rho_max;
        scenario.rho_R = 0.5 * rho_max;
        scenario.rho_0 = zeros(Nx, 1);
        scenario.rho_0(x < 0.5) = scenario.rho_L;
        scenario.rho_0(x >= 0.5) = scenario.rho_R;
        % scenario.rho_L = rho_max;
        % scenario.rho_R = 0.2 * rho_max;
        % scenario.rho_0 = zeros(Nx, 1);
        % scenario.rho_0(x < 0.5) = scenario.rho_L;
        % scenario.rho_0(x >= 0.5) = scenario.rho_R;
    end
               

    % options = {'Traffic jam' , 'Green light', 'Traffic flow' };
    % 
    %     % Create selection dialog
    %     [selection, ok] = listdlg('PromptString', 'Select an option:', ...
    %                               'SelectionMode', 'single', ...
    %                               'ListString', options);
    % 
    %     % Check if user made a selection
    %     if ok
    %         selectedOption = options{selection};
    %         fprintf('You selected: %s\n', selectedOption);
    % 
    %         % Set variables based on selection
        % 
        % 
        % else
        %     selectedOption = 'No selection';
        %     disp('No option selected. Exiting...');
        % end


        scenario.rho_c = rho_c;

        if(scenario.rho_L < rho_c)
            fprintf('Information travels right (-->) on the left part\n');
        elseif  (scenario.rho_L == rho_c)
            fprintf('Information stands still on the left part\n');
        else  
            fprintf('Information travels left (<--) on the left part\n');
        end

        if (scenario.rho_R < rho_c)
            fprintf('Information travels right (<--) on the right part\n')
        elseif  (scenario.rho_R == rho_c)
            fprintf('Information stands still on the right part\n');
        else  
            fprintf('Information travels left (<--) on the right part\n');
        end

end
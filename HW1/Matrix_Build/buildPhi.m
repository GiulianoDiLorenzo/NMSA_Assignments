function [Phi] = buildPhi(Mesh, plotMe)
% ========================================================================
%   OUTPUT : Structure with handler for the N_pts basis functions

%   INPUTS : 
%       - Mesh   --> Structure of the Mesh representation 
%       - plotMe --> Boolean for plotting (==True)
% ========================================================================

Phi = cell(Mesh.n_pts, 1);

h = Mesh.h;
x = Mesh.coord;


% First basis function, only in [x_1, x_2]
phi_1 = @(x) (x >= x(1) & x <= (x(1) + h)) .* ((x(1) + h - x) / h);

%First column of Phi matrix
Phi{1} = phi_1;

for j = 2: Mesh.n_pts - 1
    
    % Storing points of interest
    x_prev = x(j-1);
    x_now = x(j);
    x_next = x(j+1);
    
    % Handler creation, following the linear basis functions definitions
    phi_j = @(x) (x >= x_prev & x < x_now)  .* (x - x_prev) / h  + ...
                 (x >= x_now & x <= x_next) .* (x_next - x) / h;

    Phi{j} = phi_j;
    
     % % Some prints for checks
     % fprintf('x_prev = %f \n', x_prev);
     % fprintf('x_now = %f \n', x_now);
     % fprintf('x_next = %f \n', x_next);
     % fprintf('dx = %f', x_now - x_prev);
     % fprintf('\n\n');
end

% Last basis function, only in [x_{N_h}, x_{N_h+1}]
phi_end = @(x) (x >= (x(end)-h) & x <= x(end)) .*  (x - (x(end)-h) ) / h;
Phi{end} = phi_end;

if plotMe
    figure();
    
    L = {}; % Initialize an empty cell array
    for i = 1  : Mesh.n_pts
        phi_i = Phi{i};
    
        phi_i_val = phi_i(x);
        
        plot(x, phi_i_val);
        hold on
        grid on;
        tag_i = sprintf('$\\varphi_{%d}(x)$', i);   % Proper LaTeX format with
        L{end+1} = tag_i;
    end
    legend(L, 'Interpreter','latex');
    title(['Linear basis function graph, for N_{pts} = ', num2str(Mesh.n_pts)]);
    print('Plots/Basis_functions', '-dpng',  '-r300');
end

end
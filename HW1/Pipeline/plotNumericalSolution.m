function plotNumericalSolution(datas, titleText, sweepType)
L = {};
figure();

if strcmp(sweepType , 'N_pts')
    x_axis = datas(end).data.x;
    sol_ex = datas(end).data.uex(x_axis);
    plot( x_axis, sol_ex, "LineWidth",4, "Color",'r', LineStyle=':'); % Exact solution
    L{end+1} =  'exact solution';
end

if strcmp(sweepType , 'N_pts')
   % markOpts = ["-x", "-x", "-x" , "-b" , "-c"];
   % colorOpts = ["#7E2F8E" , "m" , "k" , "b" , "c"];
   % lineWidths = [1 , 1 , 1 , 2.5 , 1.5];
end

hold on; grid on

for i = 1:length(datas)
    x_axis = datas(i).data.x;
    sol_num_i = datas(i).sol;
    % disp([' i =', num2str(i), ' size(sol_num) = ', num2str(size(sol_num_i)) ]);

    
    if strcmp(sweepType , 'N_pts')
        plot(x_axis , sol_num_i ); % markOpts(i), 'Color', colorOpts(i), 'LineWidth', lineWidths(i));
        L{end+1} =  ['h =  ' , num2str(datas(i).mesh.h)];
    elseif strcmp(sweepType , 'omega')
        plot(x_axis , sol_num_i);
        L{end+1} = ['\omega = ' , num2str(datas(i).data.omega)];
    end
end


title(titleText);
legend(L);
xlabel('x [m]');
ylabel('amplitude');
end
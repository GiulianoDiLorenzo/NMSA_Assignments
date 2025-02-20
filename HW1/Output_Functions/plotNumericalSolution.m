function plotNumericalSolution(datas, titleText, saveName , sweepType, isExSol)
L = {};
figure();

if ( isExSol )
    x_axis = datas(end).data.x;
    sol_ex = datas(end).data.uex(x_axis);
    plot( x_axis, sol_ex, "LineWidth",3, "Color",'r', LineStyle='-.'); % Exact solution
    L{end+1} =  'exact solution';
end

if strcmp(sweepType , 'N_pts')
   markOpts = ["-o", "-o", "-x" , "-" , "--" ];
   colorOpts = [ "k" , "#7E2F8E" , "m" , "c" , "b" ];
   lineWidths = [1 , 1 , 1 , 1 , 1.5 ];
end

hold on; grid on

for i = 1:length(datas)
    x_axis = datas(i).data.x;
    sol_num_i = datas(i).sol;
    % disp([' i =', num2str(i), ' size(sol_num) = ', num2str(size(sol_num_i)) ]);

    if strcmp(sweepType , 'N_pts')    
        L{end+1} =  ['h =  ' , num2str(datas(i).mesh.h)];
        plot(x_axis , sol_num_i , markOpts(i), 'Color', colorOpts(i), 'LineWidth', lineWidths(i));

    elseif strcmp(sweepType , 'omega')
        L{end+1} = ['\omega = ' , num2str(datas(i).data.omega)];
        plot(x_axis , sol_num_i);

    elseif strcmp(sweepType , 'mu')
        L{end+1} = ['\mu_2 = ' , num2str(datas(i).data.mu2)];
        plot(x_axis , sol_num_i, 'LineWidth', 2, 'LineStyle','-.');

    elseif strcmp(sweepType , 'rho')
        L{end+1} = ['\rho_2 = ' , num2str(datas(i).data.rho2)];
        plot(x_axis , sol_num_i);
    end
end


title(titleText, Interpreter="latex");
legend(L);
xlabel('x [m]');
ylabel('amplitude');

print(saveName,'-dpng', '-r300');
end
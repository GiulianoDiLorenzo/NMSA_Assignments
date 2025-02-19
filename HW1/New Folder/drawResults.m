function drawResults(x , sol_ex, solutions, x_vals, h_vals)
L = {};
figure();
plot(x, sol_ex, "LineWidth",2, "Color","r");
hold on;

for i = 1:length(solutions)
    plot(x_vals(i) , solutions(i), '-');
    L{end+1} =  ['h =  ' num2str(h_vals(i))];
end

legend(L);


end
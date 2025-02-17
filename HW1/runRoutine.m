function [err_mat, sol_mat, fem_mat, data_mat] = runRoutine(TestName, omega, numElemList)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Preallocate struct arrays
err_mat(length(numElemList)) = struct();  
sol_mat(length(numElemList)) = struct();
fem_mat(length(numElemList)) = struct();
data_mat(length(numElemList)) = struct();

for i = 1 : length(numElemList) 
    [Data] = NewDataTest(TestName, numElemList(i), omega);
    [err, sol, fem, data] = Pipeline(Data, numElemList(i));
    err_mat{i} = err;
    sol_mat{i} = sol;
    fem_mat{i} = fem;
    data_mat{i} = data;
end

end
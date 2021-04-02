load solutions_smooth

[~,abso] = config_coefficient_smooth(u{1}.Mesh.Nodes(1,:),u{1}.Mesh.Nodes(2,:));
numSol = length(u);
lenSol = length(u{1}.NodalSolution);
H = zeros(lenSol,numSol);
parfor i = 1: numSol
    H(:,i) = abso'.*u{i}.NodalSolution;
end

save data_smooth H
disp('Data generated!')
pause(.1)

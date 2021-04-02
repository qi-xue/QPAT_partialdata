load solutions_discrete

[~,abso] = config_coefficient_discrete(u{1}.Mesh.Nodes(1,:),u{1}.Mesh.Nodes(2,:));
numSol = length(u);
lenSol = length(u{1}.NodalSolution);
H = zeros(lenSol,numSol);
for i = 1: numSol
    H(:,i) = abso'.*u{i}.NodalSolution;
end

save data_discrete H

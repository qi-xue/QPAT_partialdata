% Demonstrate the solution of -div(diff*∇(u))+abso*u = f
%   with suitable BC,
%   using PDE toolbox.
%% Step 1: Create PDEModel with N equations.
model = createpde(1);

%% Step 2: Create geometry.
% create a rectangle
circ = [1;0;0;1]; % basic shape -- a circle;
gd   = [circ]; % combine the shapes into one matrix;
ns   = char('circle'); % give each shape a name;
ns   = ns';
sf   = 'circle'; % the union method of shapes;
[dl,bt] = decsg(gd,sf,ns); % combine shapes using set formula;
[dl,bt] = csgdel(dl,bt); % remove face boundaries;
geometryFromEdges(model,dl);
% % % pdegplot(model,'EdgeLabels','on')
% % % xlim([-1.1,1.1])
% % % axis equal

%% Step 4: Create the PDE coefficients.
% m ∂^2u/∂t^2+d ∂u/∂t−∇·(c∇u)+au=f
specifyCoefficients(model,'m',0,'d',0,'c',@coef_c,'a',@coef_a,'f',0);

%% Step 5: Generate mesh.
generateMesh(model,'Hmax',.01);

center = 4/9*pi:pi/18:5/9*pi;
u = cell(length(center),1);
parfor i = 1: length(center)
    %% Step 3: Create boundary conditions.
    % Dirichlet hu = r
    applyBoundaryCondition(model,'dirichlet','edge',1:4,'h',1,'r',@(region,state)coe_r(region,state,center(i)));
    %% Step 6: Solve the PDE.
    u{i} = solvepde(model);
    disp("The "+num2str(i)+"th boundary condition is finished.")
    pause(.05)
end
save solutions_discrete u center

for i = 1: length(center)
    %% Step 7: Plot the solution
    uR = real(u{i}.NodalSolution);
    pdeplot(model,'XYData',uR,'ZData',uR,'Mesh','on')
    colormap jet
    shading interp
    xlabel('x')
    ylabel('y')
    view([0 0 1])
    axis equal
    pause(1)
end

function a = coef_a(region,state)
    [~,a] = config_coefficient_discrete(region.x,region.y);
end

function c = coef_c(region,state)
    [c,~] = config_coefficient_discrete(region.x,region.y);
end

function r = coe_r(region,state,center)
    r = config_boundary(region.x,region.y,center);
end


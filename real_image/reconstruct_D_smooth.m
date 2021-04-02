% Solve the Cauchy problem (-\Delta+b)\sqrt D = H/\sqrt(sigma)

%% Load data and interpolate
load coefficient_D_smooth.mat
indSin = find(H(:,ref)<pre/60);
coeD(indSin) = 20/.2;
intD = scatteredInterpolant(x,y,coeD,'natural','linear');
source = H(:,ref)./sqrt(sigma);
source(indSin)= 20/sqrt(.2);
intS = scatteredInterpolant(x,y,source,'natural','linear');

%% Construct and solve PDE model for \sqrt(D)
model = createpde(1);
circ = [1;0;0;1]; % basic shape -- a circle;
gd   = [circ]; % combine the shapes into one matrix;
ns   = char('circle'); % give each shape a name;
ns   = ns';
sf   = 'circle'; % the union method of shapes;
[dl,bt] = decsg(gd,sf,ns); % combine shapes using set formula;
[dl,bt] = csgdel(dl,bt); % remove face boundaries;
geometryFromEdges(model,dl);
generateMesh(model,'Hmax',.01);

specifyCoefficients(model,'m',0,'d',0,'c',1,...
    'a',@(region,state)coeFunD(region,state,intD),...
    'f',@(region,state)coeFunS(region,state,intS));

applyBoundaryCondition(model,'dirichlet','edge',1:4,'h',1,'r',@coeFunB);

u = solvepde(model);

%% Plot reconstructed D
xx = model.Mesh.Nodes(1,:)';
yy = model.Mesh.Nodes(2,:)';
indNor = find(yy>0.1);
ind = setdiff(1:length(xx),indNor);

% D
D  = real(u.NodalSolution).^2;
D(ind) = .2;
fig = figure('Position',[50,50,400,400]);
pdeplot(model,'XYData',D,'ZData',D,'Mesh','on')
shading interp
view([0 0 1])
colormap jet
caxis([.1,.3])
axis equal tight
ylim([.1,1])
xlim([-1,1])
hc = colorbar;
hc.FontSize = 15;
ax = gca;
set(ax,'FontSize',15)
print(fig,'-depsc','diffusion_real_reconstruct')

fig = figure('Position',[50,50,400,400]);
DTru = config_coefficient_smooth(xx,yy);
pdeplot(model,'XYData',DTru,'ZData',DTru,'Mesh','on')
shading interp
colormap jet
caxis([.1,.3])
view([0 0 1])
axis equal tight
ylim([.1,1])
xlim([-1,1])
hc = colorbar;
hc.FontSize = 15;
ax = gca;
set(ax,'FontSize',15)
print(fig,'-depsc','diffusion_real_true')

fig = figure('Position',[50,50,400,400]);
pdeplot(model,'XYData',DTru,'ZData',DTru,'Mesh','on')
shading interp
colormap jet
caxis([.1,.3])
view([0 0 1])
axis equal tight
hc = colorbar;
hc.FontSize = 15;
ax = gca;
set(ax,'FontSize',15)
print(fig,'-depsc','diffusion_real')

% mu
HNew     = scatteredInterpolant(x,y,H(:,ref),'natural','linear');
sigmaNew = scatteredInterpolant(x,y,sigma   ,'natural','linear');
mu = HNew(xx,yy)./sqrt(sigmaNew(xx,yy)./D);
mu(ind) = 20;
fig = figure('Position',[50,50,400,400]);
pdeplot(model,'XYData',mu,'ZData',mu,'Mesh','on')
shading interp
view([0 0 1])
colormap jet
caxis([10,30])
axis equal tight
ylim([.1,1])
xlim([-1,1])
hc = colorbar;
hc.FontSize = 15;
ax = gca;
set(ax,'FontSize',15)
print(fig,'-depsc','absorption_real_reconstruct')

fig = figure('Position',[50,50,400,400]);
[~,muTru] = config_coefficient_smooth(xx,yy);
pdeplot(model,'XYData',muTru,'ZData',muTru,'Mesh','on')
shading interp
colormap jet
caxis([10,30])
view([0 0 1])
axis equal tight
ylim([.1,1])
xlim([-1,1])
hc = colorbar;
hc.FontSize = 15;
ax = gca;
set(ax,'FontSize',15)
print(fig,'-depsc','absorption_real_true')

fig = figure('Position',[50,50,400,400]);
pdeplot(model,'XYData',muTru,'ZData',muTru,'Mesh','on')
shading interp
colormap jet
caxis([10,30])
view([0 0 1])
axis equal tight
hc = colorbar;
hc.FontSize = 15;
ax = gca;
set(ax,'FontSize',15)
print(fig,'-depsc','absorption_real')

close all

%% Compare with background truth
disp("Relative error for diffusion coefficient:")
disp(norm(D(indNor)-DTru(indNor))/norm(DTru(indNor))) % 4.74%
disp("Relative error for absorption coefficient:")
disp(norm(mu(indNor)-muTru(indNor))/norm(muTru(indNor))) % 2.47%

function f = coeFunD(region,state,intD)
    f = intD(region.x,region.y);
end
function f = coeFunS(region,state,intS)
    f = intS(region.x,region.y);
end
function f = coeFunB(region,state)
    f = sqrt(config_coefficient_smooth(region.x,region.y));
end

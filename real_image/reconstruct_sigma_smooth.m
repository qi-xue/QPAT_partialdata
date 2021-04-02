load data_smooth
load solutions_smooth

% For complex medium, we don't need smooth at all!!


ref = round(numel(center)/2);
D1  = H(:,ref-1)./H(:,ref);
D2  = H(:,ref+1)./H(:,ref);
pre = .007;
tri = u{1}.Mesh.Elements;
tri = tri';
x   = u{1}.Mesh.Nodes(1,:);
x   = x';
y   = u{1}.Mesh.Nodes(2,:);
y   = y';

%% compute \nabla\ln\sigma
[dxD1,dyD1] = trigradient(tri,x,y,D1);
[dxD2,dyD2] = trigradient(tri,x,y,D2);
dxD1        = trismooth(x,y,dxD1,pre/2);
dyD1        = trismooth(x,y,dyD1,pre/2);
dxD2        = trismooth(x,y,dxD2,pre/2);
dyD2        = trismooth(x,y,dyD2,pre/2);
[dxdxD1,~]  = trigradient(tri,x,y,dxD1);
[~,dydyD1]  = trigradient(tri,x,y,dyD1);
[dxdxD2,~]  = trigradient(tri,x,y,dxD2);
[~,dydyD2]  = trigradient(tri,x,y,dyD2);
lapD1       = dxdxD1+dydyD1;
lapD2       = dxdxD2+dydyD2;
lapD1       = trismooth(x,y,lapD1,pre/2);
lapD2       = trismooth(x,y,lapD2,pre/2);
lenDat      = size(H,1);
graLnSig = zeros(lenDat,2);
for i = 1: lenDat
    temp = -[dxD1(i),dyD1(i);dxD2(i),dyD2(i)]\[lapD1(i);lapD2(i)];
    graLnSig(i,:) = temp';
end

%% compute \ln\sigma & \sigma
% Start the integration from different points and average.
staAng  = linspace(center(ref)-pi/9,center(ref)+pi/9,10);
numSta  = numel(staAng);
lnSig   = zeros(lenDat,numSta);
sigma   = zeros(lenDat,numSta);
staX    = cos(staAng);
staY    = sin(staAng);
staDif  = config_coefficient_smooth(staX,staY);
staSol  = config_boundary(staX,staY,center(ref));
staSig  = log(staDif.*staSol.*staSol);
% Construct truth boundary value of \sigma
indBou = find(x.^2+y.^2>=.998^2);
difBou = config_coefficient_smooth(x(indBou),y(indBou));
solBou = config_boundary(x(indBou),y(indBou),center(ref));
sigma(indBou,:) = difBou.*solBou.*solBou*ones(1,numSta);
lnSig(indBou,:) = log(sigma(indBou,:));
funGraX = scatteredInterpolant(x,y,graLnSig(:,1),'natural','linear');
funGraY = scatteredInterpolant(x,y,graLnSig(:,2),'natural','linear');
for indSta = 1: numSta
    x0     = staX(indSta);
    y0     = staY(indSta);
    lnSig0 = staSig(indSta);
    ind    = setdiff(1:lenDat,indBou);
    for i = ind
        xt     = x(i);
        yt     = y(i);
        dir    = [xt-x0,yt-y0];
        lenDir = norm(dir);
        t      = 0:pre/2:lenDir;
        dir    = dir/lenDir;
        xInt = ones(numel(t),2)*[x0,0;0,y0]+t'*dir;
        xMid = (xInt(2:end,1)+xInt(1:end-1,1))/2;
        yMid = (xInt(2:end,2)+xInt(1:end-1,2))/2;
        M = funGraX(xMid,yMid);
        N = funGraY(xMid,yMid);
        lnSig(i,indSta) = lnSig0+M'*(xInt(2:end,1)-xInt(1:end-1,1))+N'*(xInt(2:end,2)-xInt(1:end-1,2));
        sigma(i,indSta) = exp(lnSig(i,indSta));
    end
    disp("Starting point "+num2str(indSta)+" is finished.")
    pause(.1)
end
sigma = sum(sigma,2)/numSta;
% Verify with background truth
D      = config_coefficient_smooth(x,y);
sigTru = D.*u{ref}.NodalSolution.*u{ref}.NodalSolution;
disp("The relative error for \sigma")
indSin = find(H(:,ref)<pre/30);
ind = setdiff(1:size(H,1),indSin);
disp(norm(sigma(ind)-sigTru(ind))/norm(sigTru))
% Plot
datPlo = sigma;
datPlo(indSin) = 0;
figure
trisurf(tri,x,y,datPlo,'edgecolor','w','facecolor','interp')
shading interp
colorbar
view([0 0 1])
datPlo = sigTru;
datPlo(indSin) = 0;
figure
trisurf(tri,x,y,datPlo,'edgecolor','w','facecolor','interp')
shading interp
colorbar
view([0 0 1])



%% compute \Delta\sqrt(\sigma) / \sqrt(\sigma))
% = |\nabla\ln\sigma|^2/4+\Delta\ln\sigma/2
temp1 = graLnSig(:,1);
temp2 = graLnSig(:,2);
temp1 = trismooth(x,y,graLnSig(:,1),pre/2);
temp2 = trismooth(x,y,graLnSig(:,2),pre/2);
[dxTemp1,~] = trigradient(tri,x,y,temp1);
[~,dyTemp2] = trigradient(tri,x,y,temp2);
coeD = (dxTemp1+dyTemp2)/2+(temp1.^2+temp2.^2)/4;
% verify with background truth
[dxSig,dySig] = trigradient(tri,x,y,sqrt(sigTru));
dxSig         = trismooth(x,y,dxSig,pre/2);
dySig         = trismooth(x,y,dySig,pre/2);
[dxdxSig,~]   = trigradient(tri,x,y,dxSig);
[~,dydySig]   = trigradient(tri,x,y,dySig);
coeDTru       = (dxdxSig+dydySig)./sqrt(sigTru);
disp("The relative error for coefficient of D")
indSin = find(H(:,ref)<pre/30);
ind    = setdiff(1:lenDat,indSin);
disp(norm((coeDTru(ind)-coeD(ind)))/norm(coeDTru(ind)))
datPlo  = coeD;
datPlo(indSin) = 100;
figure
trisurf(tri,x,y,datPlo,'edgecolor','w','facecolor','interp')
shading interp
colorbar
view([0,0,1])
datPlo  = coeDTru;
datPlo(indSin) = 100;
figure
trisurf(tri,x,y,datPlo,'edgecolor','w','facecolor','interp')
shading interp
colorbar
view([0,0,1])


save coefficient_D_smooth x y tri ref H pre sigma sigTru coeD coeDTru D1 dxD1 lapD1 graLnSig


%%% This file is only used to debug and verify the numerical computation!

%% verify D1 dxD1 \Delta D1
% % % indSin = find(H(:,ref)<pre);
% % % datPlo = D1;
% % % datPlo(indSin) = 0;
% % % figure
% % % trisurf(tri,x,y,datPlo,'edgecolor','w','facecolor','interp')
% % % shading interp
% % % colorbar
% % % view([0 0 1])
% % % datPlo = dxD1;
% % % datPlo(indSin) = 0;
% % % figure
% % % trisurf(tri,x,y,datPlo,'edgecolor','w','facecolor','interp')
% % % shading interp
% % % colorbar
% % % view([0 0 1])
% % % datPlo = lapD1;
% % % datPlo(indSin) = 0;
% % % figure
% % % trisurf(tri,x,y,datPlo,'edgecolor','w','facecolor','interp')
% % % shading interp
% % % colorbar
% % % view([0 0 1])

%% verify \nabla \ln \sigma
close all
indSin = find(H(:,ref)<pre*2);
indBou = find(x.^2+y.^2>=.999^2);
indBou = setdiff(indBou,indSin);
indRin = find(x.^2+y.^2>=.99^2);
indRin = setdiff(indRin,indSin);
ind    = setdiff(indRin,indBou);
xBou   = x(ind);
yBou   = y(ind);
figure
plot(xBou,graLnSig(ind,1),'bo')
hold on
datPlo = trigradient(tri,x,y,sigTru)./sigTru;
plot(xBou,datPlo(ind),'r*')
plot(x(indBou),graLnSig(indBou,1),'kx')
plot(x(indBou),datPlo(indBou),'g.')
legend('Approx. Int.','True Int.','Approx. Bound.','True Bound.')

datPlo = graLnSig(:,1);
datPlo(indSin) = 0;
figure
trisurf(tri,x,y,datPlo,'edgecolor','w','facecolor','interp')
shading interp
view([0 0 1])
grid off
colorbar
axis tight

datPlo = trigradient(tri,x,y,sigTru)./sigTru;
datPlo(indSin) = 0;
figure
trisurf(tri,x,y,datPlo,'edgecolor','w','facecolor','interp')
shading interp
view([0 0 1])
grid off
colorbar
axis tight


%% verify sigma
% % % close all
% % % indSin = find(H(:,ref)<pre*2);
% % % indBou = find(x.^2+y.^2>=.999^2);
% % % indBou = setdiff(indBou,indSin);
% % % indRin = find(x.^2+y.^2>=.99^2);
% % % indRin = setdiff(indRin,indSin);
% % % ind    = setdiff(indRin,indBou);
% % % figure
% % % plot(x(ind),sigma(ind),'bo')
% % % hold on
% % % plot(x(ind),sigTru(ind),'r*')
% % % legend('Approx. Near Bound.','True Near Bound.')
% % % 
% % % datPlo = sigma;
% % % datPlo(indSin) = 0;
% % % figure
% % % trisurf(tri,x,y,datPlo,'edgecolor','w','facecolor','interp')
% % % shading interp
% % % colorbar
% % % view([0 0 1])
% % % 
% % % datPlo = sigTru;
% % % datPlo(indSin) = 0;
% % % figure
% % % trisurf(tri,x,y,datPlo,'edgecolor','w','facecolor','interp')
% % % shading interp
% % % colorbar
% % % view([0 0 1])
% % % 
% % % ind = setdiff(1:size(H,1),indSin);
% % % disp(norm(sigma(ind)-sigTru(ind))/norm(sigTru))
% % % 
% % % datPlo = trismooth(x,y,sqrt(sigma),.01);
% % % datPlo = trilaplace(tri,x,y,datPlo);
% % % indBou = find(x.^2+y.^2>=.97^2);
% % % [difBou,absBou] = config_coefficient_smooth(x(indBou),y(indBou));
% % % datPlo(indBou)  = absBou./difBou.*sqrt(sigma(indBou));
% % % indBou = find(y>=.9);
% % % [difBou,absBou] = config_coefficient_smooth(x(indBou),y(indBou));
% % % datPlo(indBou)  = absBou./difBou.*sqrt(sigma(indBou));
% % % datPlo(indSin) = 0;
% % % figure
% % % trisurf(tri,x,y,datPlo,'edgecolor','w','facecolor','interp')
% % % shading interp
% % % colorbar
% % % view([0 0 1])
% % % 
% % % datPlo = trismooth(x,y,sqrt(sigTru),.01);
% % % datPlo = trilaplace(tri,x,y,datPlo);
% % % datPlo(indSin) = 0;
% % % figure
% % % trisurf(tri,x,y,datPlo,'edgecolor','w','facecolor','interp')
% % % shading interp
% % % colorbar
% % % view([0 0 1])

%% verify coefficient of D
% % % close all
% % % indSin = find(H(:,ref)<pre);
% % % datPlo = coeD;
% % % datPlo(indSin) = 100;
% % % figure
% % % trisurf(tri,x,y,datPlo,'edgecolor','w','facecolor','interp')
% % % shading interp
% % % colorbar
% % % view([0,0,1])
% % % 
% % % datPlo  = coeDTru;
% % % datPlo(indSin) = 100;
% % % figure
% % % trisurf(tri,x,y,datPlo,'edgecolor','w','facecolor','interp')
% % % shading interp
% % % colorbar
% % % view([0,0,1])
% % % 
% % % ind = setdiff(1:size(H,1),indSin);
% % % disp(norm((coeDTru(ind)-coeD(ind)))/norm(coeDTru(ind)))
% % % 


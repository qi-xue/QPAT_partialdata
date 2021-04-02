x     = -1:.005:1;
lenX  = numel(x);
y     = -1:.005:1;
lenY  = numel(y);
[X,Y] = meshgrid(x,y);
[diff,abso] = config_coefficient_discrete(X(:),Y(:));
diff(X.^2+Y.^2>1) = nan;
abso(X.^2+Y.^2>1) = nan;

surf(X,Y,reshape(diff,lenY,lenX))
shading interp
view([0 0 1])
axis equal tight
caxis([.08,.35])
colorbar
title("Diffusion coefficent")


figure
surf(X,Y,reshape(abso,lenY,lenX))
shading interp
view([0 0 1])
axis equal tight
caxis([10,35])
colorbar
title("Absorption coefficent")

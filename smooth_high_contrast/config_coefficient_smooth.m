function [diff,abso] = config_coefficient_smooth(x,y)

pre  = .003;
xx   = -1.3:pre:1.3;
numx = length(xx);
yy   = -1.3:pre:1.3;
numy = length(yy);
[X,Y]= meshgrid(xx,yy);
temp = 0.2+.13*incircle (X(:),Y(:),[.5,.3],.15)...
          -.10*inpolygon(X(:),Y(:),[.05,.3,.3,.05],[.55,.55,.75,.75])...
          -.12*inpolygon(X(:),Y(:),[-.65,-.45,-.45,-.65],[-.5,-.5,.45,.45])...
          +.15*inpolygon(X(:),Y(:),[-.3,.3,0,-.3],[.05,.05,.35,.35])...
          +.08*incircle (X(:),Y(:),[-.2,-.6],.15)...
          -.08*inpolygon(X(:),Y(:),[.2,.35,.85,.7],[-.4,-.6,-.2,0]);
sigm = .05;
t    = -2.5*sigm:pre:2.5*sigm;
[XX,YY]= meshgrid(t,t);
gaus = exp(-(XX.^2+YY.^2)/2/sigm^2)/2/pi/sigm^2;
gaus = gaus/sum(sum(gaus));
temp = conv2(reshape(temp,numy,numx),gaus,'same');
diff = interp2(X,Y,temp,x,y);

temp = 20+15*incircle (X(:),Y(:),[.5,.3],.15)...
         +13*inpolygon(X(:),Y(:),[-.4,-.25,-.1],[.65,.4,.7])...
         -10*inpolygon(X(:),Y(:),[-.65,-.45,-.45,-.65],[-.5,-.5,.45,.45])...
         +14*inpolygon(X(:),Y(:),[-.3,.3,0,-.3],[.05,.05,.35,.35])...
         +15*incircle (X(:),Y(:),[-.2,-.6],.15)...
         -08*inpolygon(X(:),Y(:),[.2,.35,.85,.7],[-.4,-.6,-.2,0]);
temp = conv2(reshape(temp,numy,numx),gaus,'same');
abso = interp2(X,Y,temp,x,y);
end

function y = incircle(x,y,center,radius)

y = (x-center(1)).^2+(y-center(2)).^2 <= radius*radius;

end


function [diff,abso] = config_coefficient_smooth(x,y)

pre  = .003;
xx   = -1.3:pre:1.3;
numx = length(xx);
yy   = -1.3:pre:1.3;
numy = length(yy);
[X,Y]= meshgrid(xx,yy);
temp = 2+0.30*incircle (X(:),Y(:),[.5,.3],.15)...
        -0.15*inpolygon(X(:),Y(:),[.05,.3,.3,.05],[.55,.55,.75,.75])...
        -0.20*inpolygon(X(:),Y(:),[-.65,-.45,-.45,-.65],[-.5,-.5,.45,.45])...
        +0.20*inpolygon(X(:),Y(:),[-.3,.3,0,-.3],[.05,.05,.35,.35])...
        +0.10*incircle (X(:),Y(:),[-.2,-.6],.15)...
        -0.10*inpolygon(X(:),Y(:),[.2,.35,.85,.7],[-.4,-.6,-.2,0]);
sigm = .05;
t    = -2.5*sigm:pre:2.5*sigm;
[XX,YY]= meshgrid(t,t);
gaus = exp(-(XX.^2+YY.^2)/2/sigm^2)/2/pi/sigm^2;
gaus = gaus/sum(sum(gaus));
temp = conv2(reshape(temp,numy,numx),gaus,'same');
diff = interp2(X,Y,temp,x,y);

temp = 20+05*incircle (X(:),Y(:),[.5,.3],.15)...
         +03*inpolygon(X(:),Y(:),[-.4,-.25,-.1],[.65,.4,.7])...
         -03*inpolygon(X(:),Y(:),[-.65,-.45,-.45,-.65],[-.5,-.5,.45,.45])...
         +04*inpolygon(X(:),Y(:),[-.3,.3,0,-.3],[.05,.05,.35,.35])...
         +02*incircle (X(:),Y(:),[-.2,-.6],.15)...
         -02*inpolygon(X(:),Y(:),[.2,.35,.85,.7],[-.4,-.6,-.2,0]);
temp = conv2(reshape(temp,numy,numx),gaus,'same');
abso = interp2(X,Y,temp,x,y);
end

function y = incircle(x,y,center,radius)

y = (x-center(1)).^2+(y-center(2)).^2 <= radius*radius;

end


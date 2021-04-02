function zNew = trismooth(x,y,z,pre)
% First interpolate on a regular grid.
% Then smooth the interpolated function.
% At last interpolate the value on the original nodes.

%% Interpolate
xMin = min(x);
xMax = max(x);
xSpa = xMax-xMin;
yMin = min(y);
yMax = max(y);
ySpa = yMax-yMin;
xInt = xMin-xSpa/20:pre:xMax+xSpa/20;
yInt = yMin-ySpa/20:pre:yMax+ySpa/20;
[X,Y]= meshgrid(xInt,yInt);
fInt = scatteredInterpolant(x,y,z,'natural','linear');
zInt = fInt(X(:),Y(:));
zInt = reshape(zInt,size(X));

%% Smooth
sigm = pre*2;
t    = -3*sigm:pre:3*sigm;
[XX,YY]= meshgrid(t,t);
gaus = exp(-(XX.^2+YY.^2)/2/sigm^2)/2/pi/sigm^2;
gaus = gaus/sum(sum(gaus));
mean = ones(size(gaus));
mean = mean/sum(sum(mean));
zInt = conv2(zInt,gaus,'valid');
X    = conv2(X,mean,'valid');
Y    = conv2(Y,mean,'valid');

%% Interpolate
zNew = interp2(X,Y,zInt,x,y,'spline');

end
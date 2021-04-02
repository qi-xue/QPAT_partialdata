function [diff,abso] = config_coefficient_smooth(x,y)

imagRead = imread('photoacoustic_color.jpeg');
imagRead = rgb2gray(imagRead);
imagRead = imagRead(130:250,180:300); % crop
imagRead = flipud(imagRead);
[numRow,numCol] = size(imagRead);
imagRead(imagRead>180) = 180;

xx   = linspace(-1.1,1.1,numCol);
yy   = linspace(-1.1,1.1,numRow);
[X,Y]= meshgrid(xx,yy);
temp = .2+.2*(double(imagRead)/180-.5);
diff = interp2(X,Y,temp,x,y,'linear');

temp = 20+20*(double(imagRead)/180-.5);
abso = interp2(X,Y,temp,x,y,'linear');

end

function bou = config_boundary(x,y,center)

norX  = x./sqrt(x.^2+y.^2);
theta = acos(norX).*(sign(y)>=0)+...
        (2*pi-acos(norX)).*(sign(y)<0);
minThe= min(abs(theta-center),2*pi-abs(theta-center));
sigma = .3;
bou   = exp(-minThe.^2/2/sigma^2)/sqrt(2*pi)/sigma;

end
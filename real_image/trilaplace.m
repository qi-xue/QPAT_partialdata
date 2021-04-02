function lap = trilaplace(tri,x,y,z)
% compute \Delta z based on trigradient.

[dxZ,dyZ] = trigradient(tri,x,y,z);
[dxdxZ,~] = trigradient(tri,x,y,dxZ);
[~,dydyZ] = trigradient(tri,x,y,dyZ);
lap       = dxdxZ+dydyZ;

end
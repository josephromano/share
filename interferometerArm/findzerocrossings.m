function [xc, ic] = findzerocrossings(x, y, slope)
% 
% find zero crossings and indices for y(x) with 'increasing' or 'decreasing' slope
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(y);
ind = 1:(n-1);

switch slope
  case 'increasing'
    ic = find((y(ind)<=0) & (y(ind+1)>0));

  case 'decreasing'
    ic = find((y(ind)>=0) & (y(ind+1)<0));

  otherwise
    disp('unknown slope');
end

% perform nearest-neighborhood interpolation to get x-value of
%  zero crossings
nzeros = length(ic);
xc = zeros(nzeros,1);
for j=1:nzeros
  i1 = ic(j);
  i2 = ic(j)+1;
  x1 = x(i1);
  x2 = x(i2);
  y1 = y(i1);
  y2 = y(i2);
  xc(j) = x1 - y1*(x2-x1)/(y2-y1);
end

return


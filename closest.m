function [i,j,dist]=closest(X,Y,xi,yi)
% $Id: closest.m 358 2008-04-07 14:15:03Z zhang $
%       [I,J]=CLOSEST(X,Y,XI,YI)
%	finds I,J indicies of coordinate arrays X and Y that give the 
%	point(s) closest to the input points(s) XI,YI.
%	The point(s) of interest XI,YI may be specified as a pair
%	points, a pair of vectors, or a matrix XY with two columns
%	i.e. [I,J]=CLOSEST(X,Y,XY).  This last option allows the
%	direct output of GINPUT to be used as the input XY,
%	e.g. [I,J]=CLOSEST(X,Y,GINPUT)
%
%	John Wilkin, 4 November 93 
% mods 9 Oct 2007 to allow roms grid structure as first input (to save 
% typing lon,lat all the time) 

if isstruct(X)
  % assume this is a roms grd structure
  g = X;
  if nargin == 2
    xi = Y;
    yi = xi(:,2);
    xi = xi(:,1);
  elseif nargin == 3
    yi = xi;
    xi = Y;
  end
  X = g.lon_rho;
  Y = g.lat_rho;
else  
  if nargin == 3
    yi = xi(:,2);
    xi = xi(:,1);
  end
end

for k=1:length(xi)
  dist = abs( (xi(k)-X) + sqrt(-1)*(yi(k)-Y));
  [ii,jj] = find(dist==min(dist(:)));
  i(k) = ii;
  j(k) = jj;
  dist = min(dist(:));
end

function [han,data] = roms_quivergrd(u,v,grd,d,uscale,varargin)
% Makes a quiver plot of u,v data on the ROMS C-grid
% [han,data] = roms_quivergrd(u,v,grd,d,uscale,varargin{:})
%                                                     ^^^^^
%                                    important to use this syntax if
%                                    embedding in other functions
%
% The vectors are averaged and plotted at the RHO points grid
%
% Inputs:
%
% u,v    velocity components on their natural postions of the 
%        staggered C-grid
% grd    ROMS grd structure (see function ROMS_GET_GRID)
%
% WARNING: THE ANGLE RECORDED IN SOME ROMS GRID FILES MAY BE IN DEGREES
% THIS FUNCTION ASSUMES RADIANS
% 
% Optional inputs:
%
% d              causes the vectors to be decimated (if a scalar, both 
%                  dimensions are decimated by this factor, if a vector [nx ny]
%                  then each direction is decimated accordingly
% uscale         scale factor for vectors
%                  quiver is called with quiver(x,y,uscale*u,uscale*v,0,...)
%                    see help quiver for an explanation
% varargin       any list of valid quiver linetype/style options
%
%
% Outputs:
%   
% han = handle to graphics objects from quiver
% data = structure of ue,vn,x,y data passed to the quiver call (so that 
%        users can retain the values interpolated from the staggered C-grid
%        to rho-grid locations, and rotated to east-north coordinates)
% 
% John Wilkin
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_quivergrd.m 556 2020-01-19 01:08:38Z robertson $

% get plot state
nextplotstatewas = get(gca,'nextplot');

% hold whatever is already plotted
set(gca,'nextplot','add')

if length(size(u))>2
  % to handle singleton dimensions inherited from the time or vertical slicing
  u = squeeze(u);
  v = squeeze(v);
end

% rho points coordinates
x = grd.lon_rho;
y = grd.lat_rho;
angle = grd.angle;

if max(abs(angle(:)))>2*pi
    warning('there are angles greater than 2*pi')
end

% average vector components to rho points
u = av2(u')';
v = av2(v);

% pad to correct dimension with NaNs at edges
u = u(:,[1 1:end end]);
u(:,[1 end]) = NaN;
v = v([1 1:end end],:);
v([1 end],:) = NaN;

% decimate the vectors
if nargin < 4
  d = 1;
end
if isempty(d)
  d = 1;
end
if length(d)~=1
  dx = d(1);
  dy = d(2);
else
  dx = d;
  dy = d;
end
if dx>1 || dy>1
  sx = 1:dx:size(x,1);
  sy = 1:dy:size(x,2);
  x = x(sx,sy);
  y = y(sx,sy);
  u = u(sx,sy);
  v = v(sx,sy);
  angle = angle(sx,sy);
end

% rotate coordinates if required
uveitheta = (u+sqrt(-1)*v).*exp(sqrt(-1)*angle);
u = real(uveitheta);
v = imag(uveitheta);

% clip the zeros and NaNs (this eliminates zero-length vectors (dots) from
% the quiver plot
if 1
  % keep = find(isnan(u)~=1 & u~=0 & isnan(v)~=1 & v~=0 ...
  %  & isnan(x)~=1 & isnan(y)~=1);
  clip = find(isnan(x)==1 | isnan(y)==1 | isnan(u)==1 | isnan(v)==1 ...
   | (u==0 & v==0) ); 
  x(clip) = [];
  y(clip) = [];
  u(clip) = [];
  v(clip) = [];
end

% flag the vectors outside the axes to avoid the dots at the vector
% origin that occur with some X displays
if 1
  ax = axis;
  clip = find(x<ax(1)|x>ax(2)|y<ax(3)|y>ax(4));
  x(clip) = NaN;
end

h = quiver(x,y,uscale*squeeze(u),uscale*squeeze(v),0,varargin{:});

% restore nextplotstate to what it was
set(gca,'nextplot',nextplotstatewas);

if nargout > 0
  han = h;
end
if nargout > 1
  data.ue = u;
  data.vn = v;
  data.x = x;
  data.y = y;
end

function a = av2(a)
%AV2	grid average function.  
%       If A is a vector [a(1) a(2) ... a(n)], then AV2(A) returns a 
%	vector of averaged values:
%	[ ... 0.5(a(i+1)+a(i)) ... ]  
%
%       If A is a matrix, the averages are calculated down each column:
%	AV2(A) = 0.5*(A(2:m,:) + A(1:m-1,:))
%
%	TMPX = AV2(A)   will be the averaged A in the column direction
%	TMPY = AV2(A')' will be the averaged A in the row direction
%
%	John Wilkin 21/12/93

[m,n] = size(a);
if m == 1
	a = 0.5 * (a(2:n) + a(1:n-1));
else
	a = 0.5 * (a(2:m,:) + a(1:m-1,:));
end

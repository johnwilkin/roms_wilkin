function [han,data] = roms_curquivergrd(u,v,grd,lon0,lat0,td,nd,varargin)
% Makes a curvy track plot of u,v data on the ROMS C-grid
% [han,data] = roms_curquivergrd(u,v,grd,lon0,lat0,td,nd,...
%                  [LineStyle/'noplot'],Marker,varargin{:})
%                                               ^^^^^
%                                        important to use this syntax if
%                                           embedding in other functions
% Inputs:
%        u,v  velocity components on their natural postions of the
%                 staggered C-grid
%        grd  ROMS grd structure (see function ROMS_GET_GRID)
% lon0, lat0  vectors of starting coordinates for the tracks
%         td  duration of the tracks in days
%
% Optional inputs:
%         nd  number of points in the track (default 20)
%
%  If the next (8th) input is 'noplot' the tracks are calculated 
%         but not plotted, otherwise it is interpretted as ...
%  LineStyle  string describing line format (to plot function)
%     Marker  string describing symbol to plot at origin - use [] to
%                 disable default 'k.'
%   varargin  linetype/style options to pass to plot function
%
% Outputs:
%
% han = handle to graphics objects from plot command
% data = lon,lat coordinates of the tracks plotted
%
% John Wilkin - September 2019
% Created to support curvy vectors in roms_zview/roms_sview
% $Id: roms_curquivergrd.m 563 2020-03-30 15:22:38Z robertson $

% there is an option to compute tracks but not plot them
noplot = false;

% parse the plot line options
if nargin > 9
  if strncmp(varargin{1},'no',2)
    noplot = true;
  end
  OtherOpts = {varargin{3:end}};
  Marker = varargin{2};
  LineStyle = varargin{1};
elseif nargin == 9
  if strncmp(varargin{1},'no',2)
    noplot = true;
  end
  OtherOpts = [];
  Marker = varargin{2};
  LineStyle = varargin{1};
elseif nargin == 8
  if strncmp(varargin{1},'no',2)
    noplot = true;
  end
  OtherOpts = [];
  Marker = 'k.';
  LineStyle = varargin{1};
elseif nargin < 8
  OtherOpts = [];
  Marker = 'k.';
  LineStyle = 'k-';
end
if nargin == 6
  nd = [];
end
if isempty(nd)
  nd = 20;
end

% to handle singleton dimensions inherited from the time or vertical slicing
if length(size(u))>2
  u = squeeze(u);
  v = squeeze(v);
end

% convert starting coordinates lon0,lat0 to i,j space
[I1,J1] = ndgrid(1:size(grd.h',1),1:size(grd.h',2));
Fi = scatteredInterpolant(squash(grd.lon_rho'),squash(grd.lat_rho'),I1(:));
Fj = scatteredInterpolant(squash(grd.lon_rho'),squash(grd.lat_rho'),J1(:));
x0 = Fi(lon0,lat0);
y0 = Fj(lon0,lat0);

% average vector components to rho points
u = av2(u')';
v = av2(v);

% pad to correct dimension with NaNs at edges
u = u(:,[1 1:end end]);
u(:,[1 end]) = NaN;
v = v([1 1:end end],:);
v([1 end],:) = NaN;

% normalize by grid metrics so that velocity is fractional i,j per second
ui = u.*grd.pm;
vi = v.*grd.pn;

% work in grid index i,j coordinates in meshgrid convention that the
% velocity is given in - John Wilkin's way of doing things
[I,J] = meshgrid(1:size(grd.pm,2),1:size(grd.pm,1));

% rk4 wants velocity at two time levels, so just hold the velocity steady
% by duplicating u,v in time
U(:,:,1) = ui;
U(:,:,2) = ui;
V(:,:,1) = vi;
V(:,:,2) = vi;
T(1) = 0;
T(2) = td*86400;

% times along the track
ts = linspace(T(1),T(2),nd);

% starting points
xp(1,:) = x0;
yp(1,:) = y0;
for i=2:length(ts)
  [xp(i,:),yp(i,:)] = rk4(I,J,T,U,V,xp(i-1,:),yp(i-1,:),...
    ts(i-1),ts(i)); % ,size(X,2),size(X,1));
end

%% convert xp,yp (in i,j space) to lon/lat
[I2,J2] = ndgrid(1:size(grd.h',1),1:size(grd.h',2));
Fx = griddedInterpolant(I2,J2,grd.lon_rho');
Fy = griddedInterpolant(I2,J2,grd.lat_rho');
lonp = Fx(xp,yp);
latp = Fy(xp,yp);

if ~noplot
  
  % get plot state
  nextplotstatewas = get(gca,'nextplot');
  
  % hold whatever is already plotted
  set(gca,'nextplot','add')
  
  %% plot the tracks
  if isempty(OtherOpts) && isempty(Marker)
    han = plot(lonp,latp,LineStyle);
  elseif isempty(OtherOpts)
    han = plot(lonp,latp,LineStyle,...
      lonp(1,:),latp(1,:),Marker);
  else
    han = plot(lonp,latp,LineStyle,...
      lonp(1,:),latp(1,:),Marker,OtherOpts{:});
  end
  
  % restore nextplotstate to what it was
  set(gca,'nextplot',nextplotstatewas);
  
else
  
  han = [];
  
end

if nargout > 1
  data.lon = lonp;
  data.lat = latp;
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

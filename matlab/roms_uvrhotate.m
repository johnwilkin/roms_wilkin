function [uout,vout,xout,yout] = roms_uvrhotate(u,v,grd,d)
% Place ROMS vector component data given on the u,v c-grid positions 
% onto common rho points (simple 2-point average) and rotate to 
% east-north coordinates (for plotting with quiver). Optionally 
% decimate the data in space.
%
% [ueast,vnorth,x,y] = roms_uvrhotate(u,v,grd,d)
%
% Inputs:
%
% u,v    velocity (or stress) components on their natural postions on 
%        the staggered C-grid
%        Function works for [u,v], [ubar,vbar], [sustr,svstr], ...
% grd    grid structure (see function ROMS_GET_GRID) 
%        OR filename/url from which grid coordinates can be read
%
% Optional inputs:
%
% d      causes the vectors to be decimated. If a scalar, both
%        dimensions are decimated by this factor, if a vector [nx ny]
%        then each direction is decimated accordingly
%
% Outputs
%
% u,v, averaged to rho points and rotated to geographic coordinates
%
% If extra outputs are requested, x,y (lon/lat) are returned decimated
% the same way as the vector components. 
%
% John Wilkin - March 2010
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_uvrhotate.m 550 2020-01-18 23:25:36Z robertson $

if ~isstruct(grd)
  % assume filename was given for grd
  grd = roms_get_grid(grd,grd);
end
if isfield(grd,'angle')
  angle = grd.angle;
else
  angle = zeroes(size(grd.h));
  warning('ROMS_UVRHOTATE:noangle',...
    'ANGLE not found in grd - assuming zero rotation')
end

% decimate?
if nargin < 4
  d = 1;
end
if isempty(d)
  d = 1;
end
if d ~= 1
  d1 = d(1);
  d2 = d1;
  if length(d) > 1
    d2 = d(2);
  end
end

if nargout > 2
  % rho points coordinates
  x = grd.lon_rho;
  y = grd.lat_rho;
end

if ndims(u)>2
  % handle singleton dimensions inherited from time or vertical slicing
  u = squeeze(u);
  v = squeeze(v);
end

% Now if ndims(u) > 2 then we have data at multiple times
% or depths. A convenient way to handle this is loop over the added
% dimension applying the averaging to rho points and rotation in turn.

switch ndims(u)
  
  case 4
    [nt,nz,nlatlon] = size(u);
    uout = zeros([nt nz size(grd.lon_rho)]);
    vout = uout;
    for t=1:nt
      for k=1:nz
        [uout(t,k,:,:),vout(t,k,:,:)] = ...
          uv2rho(squeeze(u(t,k,:,:)),squeeze(v(t,k,:,:)),angle);
      end
    end
    % decimate
    if d~=1
      uout = uout(:,:,1:d2:end,1:d1:end);
      vout = vout(:,:,1:d2:end,1:d1:end);
    end
    
  case 3
    [nz,nlatlon] = size(u);
    uout = zeros([nz size(grd.lon_rho)]);
    vout = uout;
    for k=1:nz % this still works if the first index is actually times
      [uout(k,:,:),vout(k,:,:)] = ...
        uv2rho(squeeze(u(k,:,:)),squeeze(v(k,:,:)),angle);
    end
    % decimate
    if d~=1
      uout = uout(:,1:d2:end,1:d1:end);
      vout = vout(:,1:d2:end,1:d1:end);
    end
    
  case 2
    % input is 2-d data
    [uout,vout] = uv2rho(u,v,angle);
    % decimate
    if d~=1
      uout = uout(1:d2:end,1:d1:end);
      vout = vout(1:d2:end,1:d1:end);
    end
    
  otherwise
    error('wrong number of dimensions in velocity')
    
end

if nargout > 2
  % decimate the coordinates
  if d ~= 1
    xout = x(1:d2:end,1:d1:end);
    yout = y(1:d2:end,1:d1:end);
  else
    xout = x;
    yout = y;
  end
end


function [u,v] = uv2rho(u,v,angle)

% average vector components to rho points
u = av2(u')';
v = av2(v);

% pad to correct dimension with NaNs at edges
u = u(:,[1 1:end end]);
u(:,[1 end]) = NaN;
v = v([1 1:end end],:);
v([1 end],:) = NaN;

% rotate coordinates
uveitheta = (u+sqrt(-1)*v).*exp(sqrt(-1)*angle);
u = real(uveitheta);
v = imag(uveitheta);


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

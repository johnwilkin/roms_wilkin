function [ueast,vnorth] = roms_surface_geostrophic_velocity(zeta,grd)
% [ueast,vnorth,x,y] = roms_surface_geostrophic_velocity(zeta,grd)
%
% Compute surface geostrophic velocity from sea level (zeta) and convert to
% geographic east-north coordinates.
%
% Inputs:
%
% zeta   sea level on the ROMS rho points grid 
% grd    grid structure (see function ROMS_GET_GRID) with grid metrics pm,
%        pn, angle (for curvilinear or rotated grids) and Coriolis f (or 
%        lat_rho from which to calculate it)
%        OR ROMS grid file or output filename/url or from which grid 
%        coordinates can be read
%
% Outputs
%
% ueast,vnorth 
%        components of surface geostrophic velocity on the rho points grid
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
%
% Obtain an up-to-date version of this code from 
% https://github.com/johnwilkin/roms_wilkin

if ischar(grd)
  % assume filename/url was given for grd
  grd = roms_get_grid(grd,grd);
end
if isfield(grd,'angle')
  angle = grd.angle;
else
  angle = zeroes(size(grd.h));
  warning('ROMS_UVRHOTATE:noangle',...
    'ANGLE not found in grd - assuming zero rotation')
end
if ~isfield(grd,'f')
  try 
    grd.f = sw_f(grd.lat_rho);
  catch
    error('Input grd does not contain f or lon_rho to set Coriolis')
  end
end

% grid metrics at respective u and v cell faces
pm_u = 0.5*(grd.pm(:,1:end-1)+grd.pm(:,2:end));
pn_v = 0.5*(grd.pn(1:end-1,:)+grd.pn(2:end,:));

% grad zeta at u,v cell faces
dZetaDx_u = (zeta(:,2:end)-zeta(:,1:end-1)).*pm_u;
dZetaDy_v = (zeta(2:end,:)-zeta(1:end-1,:)).*pn_v;

% grad zeta at rho points, rotated to east-north 
% (recycling code from roms_uvrhotate and roms_quivergrd)
[dZdx_r,dZdy_r] = uv2rho(dZetaDx_u,dZetaDy_v,angle);

% times g/f and choose the correct vector component
g = 9.81;
ueast = -g./grd.f.*dZdy_r;
vnorth = g./grd.f.*dZdx_r;

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







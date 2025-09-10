function g = roms_get_ijgrid(g)
% grd = roms_get_ijgrid(grd)
%
% Convert lon/lat in a GRD structure to notional i,j space for plotting
% with roms_zview etc. in i,j coordinates. Set angle to zero so that 
% vectors plot aligned with i,j, axes.
%
% Designed to help visualize curvilinear grids created with Delft grid
% generation tool that does not sensibly fill the masked land points and
% messes up roms_zview plots that use pcolorjw. Useful for NOAA NOS 
% CBOFS and DBOFS grids and Phil Duzinski's Delaware River grid.
%
% John Wilkin - December 2016
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_get_ijgrid.m 596 2020-12-29 16:46:14Z wilkin $

[nj,ni] = size(g.lon_rho);

[I,J] = meshgrid(0:(ni-1),0:(nj-1));
g.lat_rho = J;
g.lon_rho = I;

[I,J] = meshgrid(0:(ni-1),0.5:1:(nj-1.5));
g.lat_v = J;
g.lon_v = I;

[I,J] = meshgrid(0.5:1:(ni-1.5),0:(nj-1));
g.lat_u = J;
g.lon_u = I;

[I,J] = meshgrid(0.5:1:(ni-1.5),0.5:1:(nj-1.5));
g.lat_psi = J;
g.lon_psi = I;

% angle
g.angle_orig = g.angle;
g.angle = zeros(size(g.angle));
g.nolatlon = 1; 
g.merc = false;


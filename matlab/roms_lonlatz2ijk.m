function [Fi,Fj,Fk] = roms_lonlatz2ijk(grd)
% [Fi,Fj] = roms_lonlat2ij(grd)
%
% Create an Interpolant class that will convert lon,lat to fractional
% ROMS i,j coordinates on the FORTRAN rho points grid
%
% Input can be a ROMS grid file (or any file with lon_rho,lat_rho) or
% a roms_wilkin style grd structure from roms_get_grid
%
% Once the Interpolant is created it can be used repeatledy with ...
% Xgrid = Fi(obs_lon,obs_lat)
% Ygrid = Fj(obs_lon,obs_lat)
%
% John Wilkin - Nov 2018
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_lonlat2ij.m 548 2020-01-18 23:23:36Z robertson $


work in progress



if isstruct(grd)
  lon = grd.lon_rho';
  lat = grd.lat_rho';
else
  lon = ncread(grd,'lon_rho');
  lat = ncread(grd,'lat_rho');
end

% set up the i,j,k, to lon,lat,z transformation
[I3,J3,K3] = ndgrid(1:size(data,1),1:size(data,2),1:size(data,3));
[I2,J2] = ndgrid(1:size(data,1),1:size(data,2));

Fx = griddedInterpolant(I2,J2,permute(g.(lo),[2 1]));
Fy = griddedInterpolant(I2,J2,permute(g.(la),[2 1]));
Fz = griddedInterpolant(I3,J3,K3,permute(g.(zz),[3 2 1]));



% set up a griddedInterpolant class
[I2,J2] = ndgrid(1:size(lon,1),1:size(lon,2));
I2 = I2-1; % shift to 0-based coordinates
J2 = J2-1;
Fi = scatteredInterpolant(lon(:),lat(:),I2(:));
Fj = scatteredInterpolant(lon(:),lat(:),J2(:));

% Now Fi(lonp,latp) will be the fractional i-index, starting from 0,
% for the position of lonp,latp

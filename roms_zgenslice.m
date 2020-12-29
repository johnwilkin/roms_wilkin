function [Var4D,Coord4D] = ...
  roms_zgenslice(file,varname,zTrk,lonTrk,latTrk,timTrk,varargin)
% [Var4D,Coord4D] = ...
%   roms_zgenslice(file,varname,zTrk,lonTrk,latTrk,timTrk,varargin)
%                               ****
%
% Get ROMS output at a set of arbitrary z, lon, lat, time coordinates,
% e.g. 4-D glider trajectory, or a scattered set of in situ observation
% points such as associated with moving platforms like animals or
% profiling vehicles.
%
% Input arguments are as for roms_genslice but with an additional input,
% zTrk, specifying the depth.
%     if scalar then use this depth at all points
%     if vector (same length as lon,lat) then use each respective value
%     if > 1 treat as a vertical s-coord k index
%     if < 0 treat as a z depth (meters) below sea surface
%
% John Wilkin - May 2019
% $Id: roms_zgenslice.m 580 2020-09-08 17:17:33Z wilkin $

if length(zTrk)==1
  zTrk = zTrk*ones(size(lonTrk));
end

% first slice along lon, lat and time
[slice,geo] = roms_genslice(file,varname,lonTrk,latTrk,timTrk,varargin{:});

% If giving a vector of z values to go along with lon,lat then it doesn't
% make sense to have interpolated the track, so trap this condition

if isfield(geo,'was_interpolated')
  disp('lon,lat track was interpolated so can''t match up to z input')
  disp('Do not give a Ntrk number of points for this case')
  disp('Interpolation to z not executed')
  return
end

% now interpolate in z
Var4D = zeros(size(lonTrk));
zoffset = zeros(size(lonTrk));
for j=1:length(lonTrk)
  if zTrk(j)>=1 
    % assume vertical k index
    Var4D(j) = slice(zTrk(j),j);
  else
    % z<0 is a water column depth below sea surface
    % use surface and bottom zw because so as to bracket the entire water column
    z = geo.z(:,j);
    z(1) = geo.zw(1,j); % zw(1,:) is actually z=-h
    z(end) = geo.zw(end,j);
    if ~any(isnan(z))
      if zTrk(j) < z(1) % point is below the seafloor
        zoffset(j) = z(1)-zTrk(j);
        z(1) = -10000;
      end
      Var4D(j) = interp1(z,slice(:,j),zTrk(j));
    else
      Var4D(j) = NaN;
      zoffset(j) = NaN;
    end
  end
end

% copy the relevant coordinate information
for field = {'dis','time','dislen'}
  Coord4D.(char(field)) = geo.(char(field))(1,:);
end
for field = {'h','zeta','tindex','lonTrk','latTrk','en','ep','angle'}
  Coord4D.(char(field)) = geo.(char(field));
end
Coord4D.zTrk = zTrk;
Coord4D.zoffset = zoffset;




function [R,button] = roms_get_river_source_locations(g,lon,lat)
% [R,button] = roms_get_river_source_locations(g,[lon lat])
%
% Output structure R contains the
%    river_direction (R.rdir)
%    river_Xposition (R.xpos)
%    river_Eposition (R.ypos) 
% for a ROMS river file point source located nearest input to lon,lat. 
% 
% The sign of the river transport (R.rsgn) is negative if the direction 
% from land to sea is in the negative ROMS coordinate direction. 
%
% Rivers entering the domain at R.rsgn < 0 locations must be given a
% negative volume transport (this is because ROMS allows for volume
% sinks as well as volume sources - the sign indicates the direction of
% flow across the cell face). 
%
% If R.rsgn = 0 the selected xpos,ypos i,j index is not a coastline point
% according to the masks in the grid structure. 
%
% Inputs:
%    g is a grid structure from roms_get_grid
%    lon,lat is the desired source location
%        if not given ginput will be invoked to select the point manually
%        if lon is a 2-elemnt vector then lat is not required
%
% Example usage:
%
% Often you need a reasonably high level of zoom to use this effectively,
% so you can combine interactive pan and point selection in turn 
% like this ...
% 
%   clear R Xpos Epos Rdir Rsgn
%   for k=1:10
%     pan on
%     disp(' ')
%     disp('Pan the plot ... press <enter> when ready to select point')
%     pause()
%     [Rselect,button] = roms_get_river_source_locations(g);
%     if button==3
%       break
%     end
%     R.xpos(k) = Rselect.xpos;
%     R.epos(k) = Rselect.epos;
%     R.rdir(k) = Rselect.rdir;
%     R.rsgn(k) = Rselect.rsgn;
%     pan off
%   end 
%
% John Wilkin - August 27, 2016 - jwilkin@rutgers.edu
%
% Copyright (c) 2020 - John L. Wilkin
% $Id: roms_get_river_source_locations.m 595 2020-12-28 22:03:45Z wilkin $

switch nargin
  case 1
    % click to get position
    disp(' select source position interactively ...')
    disp('            ... right mouse click to exit')
    [lon,lat,button] = ginput(1);
    if isempty(lon)
      error('no point selected')
    end
  case 2
    % input was a 2-element vector
    lon = lon(1);
    lat = lon(2);
  case 3
    % lon and lat were given individually
  otherwise
    error('check input argument list')
end

% find closest u point
glon = g.lon_u;
glat = g.lat_u;
[j,i] = closest(glon,glat,lon,lat);
udist = sw_dist([lat glat(j,i)],[lon glon(j,i)],'km');

% find closest v point
glon = g.lon_v;
glat = g.lat_v;
[j,i] = closest(glon,glat,lon,lat);
vdist = sw_dist([lat glat(j,i)],[lon glon(j,i)],'km');

% find if source closest to u or v face
if vdist < udist
  disp('v-face point')
  % set ROMS river_direction, river_Xposition, river_Eposition
  R.rdir = 1;
  R.xpos = i-1;
  R.epos = j;
  % check that point is on the coast
  R.rsgn = diff(g.mask_rho([j j+1],i));
  if R.rsgn == 0
    warning('selected point was not on a coastline')
    return
  end
else
  disp('u-face point')
  glon = g.lon_u;
  glat = g.lat_u;
  [j,i] = closest(glon,glat,lon,lat);
  % set ROMS river_direction, river_Xposition, river_Eposition
  R.rdir = 0;
  R.xpos = i;
  R.epos = j-1;
  % check that point is on the coast
  R.rsgn = diff(g.mask_rho(j,[i i+1]));
  if R.rsgn == 0
    warning('selected point was not on a coastline')
    return
  end
end


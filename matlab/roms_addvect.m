function [thedata,thegrid,han] = roms_addvect(file,var,time,depth,grd,vec_d,uscale,varargin)
% Adds vectors from a ROMS file to the current axes
% [thedata,thegrid,han] = roms_addvect(file,var,time,depth,grd,...
%                         vec_d,uscale,varargin)
% 
% file = roms his/avg/rst etc nc file
% var = 'u' or 'v'
%     = 'ubar' or 'vbar'
%     = 'sustr', 'svstr', 'stress', 'windstress'
%     = 'Uwind', 'Vwind', 'wind', 'winds'
%     = 'u_eastward' or 'v_northward' 
% time = time index into nc file
% depth = z depth of horizontal slice (m)
% grd can be 
%       grd structure (from roms_get_grid)
%       grd_file name
%       [] (will attempt to get grid from roms file)
% vec_d = density (decimation factor) of velocity vectors to plot over 
%       if 0 no vectors are plotted
% varargin are quiver arguments passed on to roms_quiver
%
% This needs a little work to generalize the distinction between
% plotting data from s levels or z depths
%
% John Wilkin jwilkin@rutgers.edu
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_addvect.m 601 2020-12-29 19:20:10Z wilkin $

if nargin < 5
  grd = [];
end

if isinf(depth)
  depth = grd.N;
end

% get plot state
nextplotstatewas = get(gca,'nextplot');

% hold whatever is already plotted
set(gca,'nextplot','add')

ROMScoordinates = false;

if vec_d
  switch lower(var)
    case { 'sustr','svstr','stress','windstress'}
      [u,x,y,t,grd] = roms_2dslice(file,'sustr',time,grd);
      v = roms_2dslice(file,'svstr',time,grd);
    case { 'bustr','bvstr','botstress','bstress','drag'}
      [u,x,y,t,grd] = roms_2dslice(file,'bustr',time,grd);
      v = roms_2dslice(file,'bvstr',time,grd);
    case { 'uwind','vwind','wind','winds'}
      warning('the logic for winds might be wrong w.r.t. rotation - need to check')
      [u,x,y,t,grd] = roms_2dslice(file,'Uwind',time,grd);
      v = roms_2dslice(file,'Vwind',time,grd);
    case { 'ubar','vbar'}
      [u,x,y,t,grd] = roms_2dslice(file,'ubar',time,grd);
      v = roms_2dslice(file,'vbar',time,grd);
    case { 'u','v'}
      if depth > 0
        % assume an s-level
        u = nc_varget(file,'u',[time-1 depth-1 0 0],[1 1 -1 -1]);
        v = nc_varget(file,'v',[time-1 depth-1 0 0],[1 1 -1 -1]);
      else
        [u,x,y,t,grd] = roms_zslice(file,'u',time,depth,grd);
        v = roms_zslice(file,'v',time,depth,grd);
      end
    case {'u_eastward','v_northward'}
      % case that only east-north velocities are saved (e.g.
      % MoanaProject)
      if depth > 0
        % assume an s-level
        disp('try plotting at a z-level if this fails')
        u = nc_varget(file,'u_eastward',[time-1 depth-1 0 0],[1 1 -1 -1]);
        v = nc_varget(file,'v_northward',[time-1 depth-1 0 0],[1 1 -1 -1]);
        x = grd.lon_rho;
        y = grd.lat_rho;
      else
        [u,x,y,t,grd] = roms_zslice(file,'u_eastward',time,depth,grd);
        v = roms_zslice(file,'v_northward',time,depth,grd);
      end
      ROMScoordinates = false;
    otherwise
      % assume the vector data will be for a pair of 3D variables named
      % usomething and vsomething (e.g. utemp, vtemp from the
      % quadratic averages)
      disp([ 'Plotting vector data for u' var(2:end)])
      if depth > 0
        % assume an s-level
        u = nc_varget(file,[ 'u' var(2:end)],[time-1 depth-1 0 0],[1 1 -1 -1]);
        v = nc_varget(file,[ 'v' var(2:end)],[time-1 depth-1 0 0],[1 1 -1 -1]);
      else
        % assume a z-depth
        [u,x,y,t,grd] = roms_zslice(file,[ 'u' var(2:end)],time,depth,grd);
        v = roms_zslice(file,[ 'v' var(2:end)],time,depth,grd);
      end
  end
  u = squeeze(u);
  v = squeeze(v);
  if ROMScoordinates
     hanq = roms_quivergrd(u,v,grd,vec_d,uscale,varargin{:});
  else
     hanq = quiver(x,y,uscale*u,uscale*v,0,varargin{:});
  end
end

% restore nextplotstate to what it was
set(gca,'nextplot',nextplotstatewas);

if nargout > 0
  thedata.u = u;
  thedata.v = v;
  thedata.t = t;
end
if nargout > 1
  thegrid = grd;
end
if nargout > 2
  han = hanq;
end


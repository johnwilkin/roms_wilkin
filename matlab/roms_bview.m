function [Data,han] = roms_bview(file,varname,time,bndy,grd,xcoord)
% [data,han] = roms_bview(file,var,time,bndy,grd,xcoord)
%
% file   = roms his/avg/rst etc nc file
% var    = variable to plot
% time   = time index into nc file
% bndy   = 'north','south','east','west'
% grd can be
%       grd structure (from roms_get_grid)
%       grd_file name
%       [] (will attempt to get grid from roms file)
% xcoord = 'lon','lat', or 'dist' (default) to specify plot abscissa
%
% John Wilkin jwilkin@rutgers.edu
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_bview.m 592 2020-12-28 21:56:28Z wilkin $

if nargin < 5
  grd = [];
end
if nargin < 6
  xcoord = 'dist';
end

if iscell(varname) % in case called in a cell list loop of varnames
  varname = char(varname);
end
if iscell(bndy)
  bndy = char(bndy);
end

if isempty(grd)
  try
    grd = roms_get_grid(file,file);
  catch
    error([ 'Unable to generate grd structure from ' file])
  end
end

% trap case that varname was given with boundary included
nsew_check = strfind(varname,'_');
if ~isempty (nsew_check)
  bndy = varname((nsew_check+1):end);
else
  varname = [varname '_' bndy];
end

% time coordinate variable
tvarname = roms_get_time_varname(file,varname);

% time information
if isinf(time)
  time = 'latest';
end
if ischar(time)
  time = roms_get_time_index(file,tvarname,time);
end

% get appropriate z data for u, v or rho points
switch varname(1)
  case 'u'
    pos = '_u';
  case 'v'
    pos = '_v';
  otherwise
    pos = '_rho';
end
z   = grd.(['z' pos(1:2)]);
lon = grd.(['lon' pos]);
lat = grd.(['lat' pos]);
m   = grd.(['mask' pos]);

% choose z data for correct boundary
switch bndy(1)
  case 'w'
    z = z(:,:,1);
    lon = lon(:,1);
    lat = lat(:,1);
    m = m(:,1);
  case 'e'
    z = z(:,:,end);
    lon = lon(:,end);
    lat = lat(:,end);
    m = m(:,end);
  case 's'
    z = z(:,1,:);
    lon = lon(1,:);
    lat = lat(1,:);
    m = m(1,:);
  case 'n'
    z = z(:,end,:);
    lon = lon(end,:);
    lat = lat(end,:);
    m = m(end,:);
end
z = squeeze(z);
lon = lon(:);
lat = lat(:);
m = m(:);
m(m==0) = NaN;

% compute approximate distance on sphere between lon/lat coordinate pairs
rearth = 6370.800; % km
dy = rearth*pi/180*diff(lat);
dx = rearth*pi/180*diff(lon).*cos(pi/180*0.5*(lat(2:end)+lat(1:end-1)));
dist = cumsum([0; sqrt(dx(:).^2+dy(:).^2)]);

dist = repmat(dist',[size(z,1) 1]);
lon = repmat(lon',[size(z,1) 1]);
lat = repmat(lat',[size(z,1) 1]);
m = repmat(m',[size(z,1) 1]);

data = nc_varget(file,varname,[time-1 0 0],[1 -1 -1]);
data = squeeze(data);

% get the time/date for plot label
t = roms_get_time(file,time);
if isdatetime(t)
  tdate = ['on day ' datestr(t,0)];
else
  tdate = ['on day ' num2str(t,'%8.2f')];
end

% pcolor plot of the variable
titlestr = ...
  {['file: ' strrep_(file) ],tdate,upper(strrep_(varname))};

switch xcoord
  case 'lon'
    hant = pcolorjw(lon,z,m.*data);
    xlabel('longitude')
    title(titlestr)
  case 'lat'
    hant = pcolorjw(lat,z,m.*data);
    xlabel('latitude')
    title(titlestr)
  case 'dist'
    hant = pcolorjw(dist,z,m.*data);
    xlabel('distance (m)')
    title(titlestr)
  otherwise
    % hack to skip the plot and only return the data to use with the
    % multiple boundary plot roms_bviewm.m
end

if nargout > 0
  Data.data = data;
  Data.lon = lon;
  Data.lat = lat;
  Data.dist = dist;
  Data.mask = m;
  Data.z = z;
  Data.t = t;
  Data.tstr = titlestr;
end

if nargout > 1
  han = hant;
end

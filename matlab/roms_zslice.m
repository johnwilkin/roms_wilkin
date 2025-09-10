function [data,x,y,t,grd] = roms_zslice(file,var,time,depth,grd)
% [data,x,y] = roms_zslice(file,var,time,depth,grd)
% Get a constant-z slice out of a ROMS history, averages or restart file
%
% Inputs
%    file = his or avg nc file
%    var = variable name
%    time = time index in nc file, a datestr format string, or 'last'
%    depth = depth in metres of the required slice
%    grd (optional) is the structure of grid coordinates from roms_get_grid
%
% Outputs
%
%    data = the 2d slice at requested depth
%    x,y = horizontal coordinates
%    t = time in days for the data
%
% John Wilkin
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_zslice.m 564 2020-03-30 15:23:14Z robertson $

% accommodate idiots who give positive z when it should be negative
depth = -abs(depth);

try 
  ncinfo(file,var);
catch
  error([ 'No ' var ' in ' file])
end

% allow for time in datestr format
if ischar(time)
  fdnums = nj_time(file,var);
  if ~any(strcmp(time,{'end','last','latest'}))
    dnum = datenum(time);
    if dnum >= fdnums(1) && dnum <= fdnums(end)
      % date is in the file
      [~,time] = min(abs(dnum-fdnums));
      time = time(1); % in case request falls exactly between output times
    else
      % date string was not in the file
      disp(['Date ' time ' is not between the dates in the file:'])
      disp([datestr(fdnums(1),0) ' to ' datestr(fdnums(end),0)])
      return
    end
  else
    % date string was logically 'latest'
    time = length(fdnums);
  end
  t = fdnums(time);
else
  % get the time  
  Vinfo = ncinfo(file,var);
  time_variable = Vinfo.Dimensions(end).Name;
  % t = ncread(file,time_variable,time,1);
  t = roms_get_time(file,time);
end

% check the grid 
if nargin<5 || (nargin==5 && isempty(grd))
  % no grd input given so try to get grd_file name from the history file
  grd_file = file;
  grd = roms_get_grid(grd_file,file);
else
  if ischar(grd)
    grd = roms_get_grid(grd,file);
  else
    % input was a grd structure but check that it includes the z values
    if ~isfield(grd,'z_r')
      error('grd does not contain z values');
    end
  end
end

% get the 3D chunk of data to be zsliced
data = ncread(file,var,[1 1 1 time],[Inf Inf Inf 1]);
data = double(permute(data,ndims(data):-1:1));
data = squeeze(data);

% slice at requested depth
[data,x,y] = roms_zslice_var(data,1,depth,grd);

% Apply mask to catch shallow water values where the z interpolation does
pos = roms_cgridpos(size(data),grd);
mask = grd.(['mask_' pos]);
data(mask==0) = NaN;


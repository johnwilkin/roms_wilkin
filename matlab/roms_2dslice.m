function [data,x,y,t,grd] = roms_2dslice(file,var,time,grd)
% $Id: roms_2dslice.m 479 2017-08-01 18:12:18Z wilkin $
% Get a horiztonal slice of a 2-D variable out of a ROMS history, averages or restart file
% [data,x,y,t] = roms_2dslice(file,var,time,grd)
%
% Inputs
%    file = his or avg nc file
%    var = variable name
%    time = time index in nc file
%    grd (optional) is the structure of grid coordinates from roms_get_grid 
%
% Outputs
%    
%    data = the 2d slice at requested depth 
%    x,y = horizontal coordinates
%    t = time in matlab datenum convention
%
% John Wilkin
%
% converted to use snctools 2006-10-23
% output t should now be a matlab datenum

if ~nc_isvar(file,var)
  error([ 'Variable ' var ' is not present in file ' file])
end

  if ischar(time)
    fdnums = roms_get_date(file,-1);
    if strcmp(time,'latest')
      time = length(fdnums);
    else
      dnum = datenum(time);
      if dnum >= fdnums(1) && dnum <= fdnums(end)
        [~,time] = min(abs(dnum-fdnums));
        time = time(1);
      else
        warning(' ')
        disp(['Requested date ' time ' is not between the dates in '])
        disp([file ' which are ' datestr(fdnums(1),0) ' to ' ])
        disp(datestr(fdnums(end),0))
        thedata = -1;
        return
      end
    end
  end
  
% get the time
info = nc_vinfo(file,var);
time_variable = info.Dimensions(end).Name;

% time_variable = nc_attget(file,var,'time');
% if isempty(time_variable)
%   time_variable = 'scrum_time'; % doubt this is really of any use any more 
% end

% support for old uses of roms_2dslice to load variables that are not
% progostic model output (i.e. have no time dimension)
isprognostic = 0;
dimension = getfield(nc_getvarinfo(file,var),'Dimension');
if any( strcmp(dimension,time_variable) )
  isprognostic = 1;
end

if isprognostic
  if nc_varsize(file,time_variable)<time
    disp(['Requested time index ' int2str(time) ' not available'])
    disp(['There are ' int2str(nc_varsize(file,time_variable)) ...
      ' time records in ' file])
      error(' ')
  end
  data = nc_varget(file,var,[time-1 0 0],[1 -1 -1]);
  t = roms_get_date(file,time); % gets output in matlab datenum convention
else
  data = nc_varget(file,var);
  t = [];
end
data = squeeze(data);

% check the grid information
if nargin<4 || (nargin==4 && isempty(grd))
  % no grd input given so get grd from the roms file
  grd = roms_get_grid(file,file);
else
  if ischar(grd)
    % assume input is the name of the grid file
    grd = roms_get_grid(grd,file);
  else
    % do nothing - grd structure has been input
    % 
    % might want to trap the case that grid info is invalid ...
  end
end

switch roms_cgridpos(size(data),grd)
  case 'u'
    x = grd.lon_u;
    y = grd.lat_u;
    mask = grd.mask_u;
  case 'v'
    x = grd.lon_v;
    y = grd.lat_v;
    mask = grd.mask_v; 
  case 'psi'
    x = grd.psi_v;
    y = grd.psi_v;
    mask = grd.mask_psi; 
  otherwise    
    % for zeta, Hsbl etc
    x = grd.lon_rho;
    y = grd.lat_rho;
    mask = grd.mask_rho; 
end

% Apply mask
dry = find(mask==0);
mask(dry) = NaN;
data = data.*mask;

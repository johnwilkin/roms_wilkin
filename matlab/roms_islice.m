function [data,z,lon,lat,t] = roms_islice(file,var,time,iindex,grd)
% [data,z,lon,lat,t] = roms_islice(file,var,time,iindex,grd)
% Get a constant-i slice out of a ROMS history, averages or restart file
%
% Inputs
%    file = his or avg nc file
%    var = variable name
%    time = time index in nc file
%    iindex = i index of the required slice
%    grd (optional) is the structure of grid coordinates from roms_get_grid
%
% Outputs
%
%    data = the 2d slice at requested depth
%    z (2d) matrix of depths
%    lon,lat = horizontal coordinates along the slice
%    t = time in days for the data
%
% John Wilkin
%
% July 2016 (JLW) fixed bug in the padding/extrapolation to sea surface
%           added padding/extrapolation to sea floor
%           Replaced snctools calls with native netcdf
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_islice.m 530 2019-05-15 19:35:53Z wilkin $


% get the data
% data = nc_varget(file,var,[time-1 0 0 iindex-1],[1 -1 -1 1]);

% special handling for derived variables
switch lower(var)
  case {'omegaca','omegaar','ph','pphotal'}
    temp = ncread(file,'temp',[iindex 1 1 time],[1 Inf Inf 1]);
    salt = ncread(file,'salt',[iindex 1 1 time],[1 Inf Inf 1]);
    alk = ncread(file,'alkalinity',[iindex 1 1 time],[1 Inf Inf 1]);
    tic = ncread(file,'TIC',[iindex 1 1 time],[1 Inf Inf 1]);
    press = -permute(grd.z_r,[3 2 1]);
    press = press(iindex,:,:);
    data = roms_co2sys_var(var,temp,salt,alk,tic,press);
    var = 'salt'; % Vinfo won't return coordinates for derived variables
  otherwise
    data = ncread(file,var,[iindex 1 1 time],[1 Inf Inf 1]);
end
data = double(permute(data,ndims(data):-1:1));
data = squeeze(data);

% Forecast Model Run Collection (FMRC) changes time coordinate to "time" 
% but leaves the coordinates attribute pointed to ocean_time
% Vinfo = ncinfo(file,var);
% time_variable = Vinfo.Dimensions(end).Name;
% t = nc_varget(file,time_variable,time-1,1);
% t = ncread(file,time_variable,time,1);
time_variable = roms_get_time_varname(file,var);
t = roms_get_time(file,time_variable,time);

% determine where on the C-grid these values lie
try
  varcoords = ncreadatt(file,var,'coordinates');
  if contains(varcoords,'_u')
    pos = 'u';
  elseif contains(varcoords,'_v')
    pos = 'v';
  elseif contains(varcoords,'_rho')
    pos = 'rho';
  else
    error('Unable to parse the coordinates attribute to know where the data fall on C-grid')
  end
catch
  switch var
    case 'u'
      pos = 'u';
    case 'v'
      pos = 'v';
    otherwise
      pos = 'rho';
  end
  warning('Trouble parsing coordinates attribute to know where the data fall on C-grid')
end

% check the grid information
if nargin<5 || (nargin==5 && isempty(grd))
  % no grd input given so try to get grd_file name from the file
  grd = roms_get_grid(file,file);
else
  if ischar(grd)
    grd = roms_get_grid(grd,file);
  else
    % input was a grd structure but check that it includes the z values
    if ~isfield(grd,'z_r')
      try
        grd = roms_get_grid(grd,file,0,1);
      catch
        error('grd does not contain z values');
      end
    end
  end
end

% get section depth coordinates
z_r = grd.z_r;
z_w = grd.z_w;
isw = false;

switch pos
  
  case 'u'
    % average z_r to Arakawa-C u points
    % this might be redundant if z u,v values are already in structure
    z = 0.5*(z_r(:,:,1:(end-1))+z_r(:,:,2:end));
    zw = 0.5*(z_w(:,:,1:(end-1))+z_w(:,:,2:end));
    x = grd.lon_u;
    y = grd.lat_u;
    mask = grd.mask_u;
    
  case 'v'
    % average z_r to Arakawa-C v points
    z = 0.5*(z_r(:,1:(end-1),:)+z_r(:,2:end,:));
    zw = 0.5*(z_w(:,1:(end-1),:)+z_w(:,2:end,:));
    x = grd.lon_v;
    y = grd.lat_v;
    mask = grd.mask_v;
    
  otherwise
    % for temp, salt, rho, w
    z = z_r;
    zw = z_w;
    x = grd.lon_rho;
    y = grd.lat_rho;
    mask = grd.mask_rho;
    if size(data,1) ~= size(z,1)
      % trap the var=='omega' case
      % but omega can be N or N+1 depending on whether a rst or his file
      z = grd.z_w;
      isw = true;
    end
    
end

% extract the i slice of the z coordinates
z = squeeze(z(:,:,iindex));

if ~isw
  % pad z to sea surface and seafloor for plotting
  z_surf = squeeze(zw(end,:,iindex));
  z_surf = z_surf(:)';
  z_bot = squeeze(zw(1,:,iindex));
  z_bot = z_bot(:)';
  z = [z_bot; z; z_surf];
  data = [data(1,:); squeeze(data); data(end,:)];
end

lon = repmat(x(:,iindex)',[size(z,1) 1]);
lat = repmat(y(:,iindex)',[size(z,1) 1]);

% land/sea mask
mask(mask==0) = NaN;
mask = repmat(mask(:,iindex)',[size(z,1) 1]);

% remove singleton dimensions
z = squeeze(z);
data = mask.*squeeze(data);

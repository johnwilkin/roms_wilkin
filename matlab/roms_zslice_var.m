function [data,x,y] = roms_zslice_var(data,time,depth,grd)
% Get a constant-z slice out of a 4-D ROMS variable 
% [data,x,y] = roms_zslice_var(data,time,depth,grd)
%
% Bugs/features:
% 'time' isn't used
% input must be squeezed to dimensions  Z ETA_var XI_var
%
% John Wilkin
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_zslice_var.m 477 2017-08-01 18:11:05Z wilkin $

sizexy = fliplr(size(data));
sizexy = fliplr(sizexy([1 2]));

depth = -abs(depth);

% interpolate to requested depth 

% make 2d 'ribbon' out of the data
[N L M] = size(data);
data = data(:,1:L*M);

% make 2d 'ribbon' out of the depth coordinates
z_r = grd.z_r;

if ~any(sizexy-size(grd.lon_u))
  var = 'u';
elseif ~any(sizexy-size(grd.lon_v))
  var = 'v';
else
  var = 'a rho-shaped variable';
end

switch var
  case 'u'    
    % average z_r to Arakawa-C u points
    z = 0.5*(z_r(:,:,1:(end-1))+z_r(:,:,2:end));
    x = grd.lon_u;
    y = grd.lat_u;
    mask = grd.mask_u;
    
  case 'v'
    % average z_r to Arakawa-C v points
    z = 0.5*(z_r(:,1:(end-1),:)+z_r(:,2:end,:));
    x = grd.lon_v;
    y = grd.lat_v;
    mask = grd.mask_v; 
    
  otherwise    
    % for temp, salt, rho, w
    z = z_r;
    x = grd.lon_rho;
    y = grd.lat_rho;
    mask = grd.mask_rho; 
    if size(data,1) ~= size(z,1)
      % trap the var='omega' case
      % but omega can be N or N+1 depending on whether a rst or his file
      z = grd.z_w;
    end
    
end

z = reshape(z,[size(z,1) size(z,2)*size(z,3)]);

% code lifted from omviz/scrum_zslice:

% pad out bottom and surface z values with -Inf and 0 respectively
z = [-Inf*ones([1 L*M]); z; zeros([1 L*M])];

% pad out bottom and surface data values
data = [NaN*ones([1 L*M]); data; data(N,:)];

z = flipud(z);
data = flipud(data);

% Find the indices of data values that have just greater depth than
% depth

zg_ind = find(diff(z<depth)~=0);
zg_ind = zg_ind + [0:1:length(zg_ind)-1]';
data_greater_z = data(zg_ind);
depth_greater_z = z(zg_ind);
        
% Find the indices of the data values that have just lesser depth
% than depth
zl_ind = find(diff(z>depth)~=0);
zl_ind = zl_ind + [1:1:length(zg_ind)]';
data_lesser_z = data(zl_ind);
depth_lesser_z = z(zl_ind);
        
% Interpolate between the data values.
alpha = (depth-depth_greater_z)./(depth_lesser_z-depth_greater_z);
data_at_depth = (data_lesser_z.*alpha)+(data_greater_z.*(1-alpha));
data = reshape(data_at_depth,[L M]);

% Apply mask to catch shallow water values where the z interpolation does
% not create NaNs in the data
% dry = find(mask==0);
% mask(dry) = NaN;
% data = data.*mask;
data(mask==0) = NaN;




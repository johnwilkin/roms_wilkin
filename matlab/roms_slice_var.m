function [data,z,lon,lat] = roms_slice_var(data,grd,index,ij,time)  
% $Id: roms_slice_var.m 358 2008-04-07 14:15:03Z zhang $
% Get an i or j direction slice of a roms variable 
% 
% Inputs:
%    data = 3D or 4D variable in roms coordinates
%    grd = grid structure including z_r
%    index = i or j index on the roms horizontal grid
%    ij = 'i' or 'j'
%    time = optional time index into 'data' 
%    
% Outputs:
%    data, z, lon, lat suitable for pcolor(lon,z,data) 
%
% Usage:
%    [data,z,lon,lat] = roms_slice_var(data,grd,index,ij,time)  

if length(size(data))==4
  data = squeeze(data(time,:,:,:));
end

% get section depth coordinates
z_r = grd.z_r;
z_w = grd.z_w;
Np = grd.N+1;

switch roms_cgridpos(size(data),grd)
  case 'u'    
    % average z_r to Arakawa-C u points
    zM = size(z_r,2);
    zMm = zM-1;
    zL = size(z_r,3);
    zLm = zL-1;
    zh = 0.5*(z_w(1,:,1:zLm)+z_w(1,:,2:zL));
    z = 0.5*(z_r(:,:,1:zLm)+z_r(:,:,2:zL));
    z0 = 0.5*(z_w(Np,:,1:zLm)+z_w(Np,:,2:zL));
    z = [zh; z; z0];
    x = grd.lon_u;
    y = grd.lat_u;
    mask = grd.mask_u;
    
  case 'v'
    % average z_r to Arakawa-C v points
    zM = size(z_r,2);
    zMm = zM-1;
    zL = size(z_r,3);
    zLm = zL-1;
    zh = 0.5*(z_w(1,1:zMm,:)+z_w(1,2:zM,:));
    z = 0.5*(z_r(:,1:zMm,:)+z_r(:,2:zM,:));
    z0 = 0.5*(z_w(Np,1:zMm,:)+z_w(Np,2:zM,:));
    z = [zh; z; z0];
    x = grd.lon_v;
    y = grd.lat_v;
    mask = grd.mask_v; 
    
  otherwise    
    % for temp, salt, rho, w
    zh = grd.z_w(1,:,:);
    z = z_r;
    z0 = grd.z_w(Np,:,:);
    z = [zh; z; z0];
    x = grd.lon_rho;
    y = grd.lat_rho;
    mask = grd.mask_rho; 
    % might have to do something special to handle omega
    
end

% extract the slices of the coordinates
switch ij
  case 'i'
    % slice along constant i
    z = squeeze(z(:,:,index));
    lon = repmat(x(:,index)',[size(z,1) 1]);
    lat = repmat(y(:,index)',[size(z,1) 1]);
    data = squeeze(data(:,:,index));
    % pad surface and bottom
    kindex = [1 1:size(data,1) size(data,1)];
    data = data(kindex,:);
    % land/sea mask
    dry = find(mask==0);
    mask(dry) = NaN;
    mask = repmat(mask(:,index)',[size(z,1) 1]);
    data = data.*mask;
  case 'j'
    % slice along constant j
    z = squeeze(z(:,index,:));
    lon = repmat(x(index,:),[size(z,1) 1]);
    lat = repmat(y(index,:),[size(z,1) 1]);
    data = squeeze(data(:,index,:));
    % pad surface and bottom
    kindex = [1 1:size(data,1) size(data,1)];
    data = data(kindex,:);
    % land/sea mask
    dry = find(mask==0);
    mask(dry) = NaN;
    mask = repmat(mask(index,:),[size(z,1) 1]);
    data = data.*mask;
  otherwise
    error([ 'Specify either ''i'' or ''j'' for input ''ij'' '])
end

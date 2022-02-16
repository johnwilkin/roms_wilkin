function mld = roms_mld(file,time,grd,opt)
% mld = roms_mld(file,time,grd,opt)
% Plot Mixed Layer Depth (MLD) from a ROMS  history or averages file
%     or THREDDS aggregation
%
% Depth where density is 0.125 kg/m3 greater than rho(z=-10 m) 
%
% John Wilkin jwilkin@rutgers.edu
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu

if isinf(time)
  time = 'latest';
end
if ischar(time) || isdatetime(time)
  time = roms_get_time_index(file,time);
end

% MLD criterion
% Add support for other criteria later
% if nargin < 4
%  opt = 'rho10m_plus_p125';
% end
opt = 'rho10m_plus_p125';

% read temp and salt
START = [time-1 0  0  0];
COUNT = [1     -1 -1 -1];
temp = nc_varget(file,'temp',START,COUNT);
salt = nc_varget(file,'salt',START,COUNT);
salt(salt<0) = 0;

% potential density from salt, temp, pressure
rho = sw_dens(salt(:,:),temp(:,:),-grd.z_r(:,:));
rho = reshape(rho,size(temp));

switch opt
  
  case {'rho10m_plus_p125','aristizabal'}
    % depth where rho is 0.0125 kg/m^3 greater than rho at -10 m
    
    % density at z = -10
    rho10m = roms_zslice_var(rho,1,-10,grd);
    
    % density minus 10 m value
    data = rho-reshape(rho10m,[1 size(rho10m)]);
    data2d = data(:,:);
    
    % find where this 0.125 kg/m3
    Rvalue = 0.125;
    
    % only compute on wet cells
    water = 1:length(data2d(1,:));
    water(grd.mask_rho==0) = [];
    
    % prealllocate for speed
    mld = nan(size(data2d(1,:)));
    
    for i = water
      % find shallowest point with data less than "value"
      k = find(data2d(:,i)>Rvalue,1,'last');
      if isempty(k)
        % all points were greater than value
        % want the whole water column
        mld(i) = -grd.h(i);
      else
        if k < length(data2d(:,1))-1
          mld(i) = interp1(data2d([k k+1],i),grd.z_r([k k+1],i),Rvalue);
        else
          % all points less than value - the surface is not in the profile
          mld(i) = 0;
        end
      end
    end
    mld = -reshape(mld,size(grd.lon_rho));
   
end

    
    
    
    
    
    
    

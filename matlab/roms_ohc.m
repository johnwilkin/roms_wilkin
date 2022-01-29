function ohc = roms_ohc(file,time,grd)
% ohc = roms_ohc(file,time,grd,opt)
% Calculate Ocean Heat Content from a ROMS history or averages file
%     or THREDDS aggregation
%
% Ocean Heat Content is rho0*Cp * integral_T=26^z=0 (T-26) dz in kJ/cm^2
% 
% Whitaker, W. D., 1967: Quantitative determination of heat transfer from
%     sea to air during passage of Hurricane Betsy. Ph.D. thesis, 
%     Texas A&M University.
% Leipper, D. F., and D. Volgenau, 1972: Hurricane heat potential of the 
%     Gulf of Mexico. Journal of Physical Oceanography, 2 (3), 218–224.
% Aristizabal, M., T. Miles, S. Glenn, Evaluation of Upper Ocean Metrics
%     Relevant to Air-Sea Heat Fluxes in the HWRF-POM and HWRF-HYCOM 
%     Forecasting Models during Hurricane Dorian (2019), unpublished
%     manuscript
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

% read temp and salt
START = [time-1 0  0  0];
COUNT = [1     -1 -1 -1];
temp = nc_varget(file,'temp',START,COUNT);

data = temp;
data2d = data(:,:);

% Find depth of the 26 C isotherm
% If all water is T>26, use the bottom depth
% If all water is T<26, depth of isosurface is zero
Tvalue = 26;

% only compute for wet cells
water = 1:length(data2d(1,:));
water(grd.mask_rho==0) = [];

% prealllocate for speed
zohc = nan(size(data2d(1,:)));

for i = water 
  % find shallowest point with data less than "value"
  k = find(data2d(:,i)<Tvalue,1,'last');
  if isempty(k)
    % all points were greater than value
    % want the whole water column
    zohc(i) = -grd.h(i);  
  else
    if k < length(data2d(:,1))-1
      zohc(i) = interp1(data2d([k k+1],i),grd.z_r([k k+1],i),Tvalue);
    else
      % all points less than value - the surface is not in the profile
      zohc(i) = 0;
    end
  end
end
zohc = reshape(zohc,size(grd.lon_rho));

% OHC calculation
rho0 = 1025;   % assume value, but also in output netcdf
Cp = 3985;     % Joules/kg/degC from mod_scalars.F

% vertical integral of T-26 from isosurface z = zohc to sea surface z = 0 
Tm26_int = roms_zint(data-Tvalue,grd,zohc,zeros(size(zohc)));

% convert J/m^2 to kJ/cm^2
fac = 1e-7; % (1/1000)/(100^2)
ohc = fac*rho0*Cp*Tm26_int;


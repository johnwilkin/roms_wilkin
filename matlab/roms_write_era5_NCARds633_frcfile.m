% Create a ROMS meteorology forcing file from ERA5 reanalysis extracted
% from NCAR dataset ds633.0 using roms_get_era5_ncar_ds633.m, i.e. ...
%
% First run e.g. E = roms_get_era5_ncar_ds633(year,month,bounding_box...)
% and then run this script to create the forcimng netcdf file.
%
% See the help on roms_get_era5_ncar_ds633 regarding obtaining login
% credentials to access the ERA5 archive at NCAR Research Data Archive.
%
% -------------------------------------------------------------------------
%
% NOTE: This is an example script. The user needs to configure the ROMS
% time coordinate basedate (Time0) and the output filename (ncname) and
% title and any other pertinent metadata.
%
% This script uses routines from the myroms.org Matlab tools
% (roms_metadata, nc_create, nc_write, etc.).
%
% John Wilkin - December 2020
%
% Copyright (c) - 2021 John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_write_era5_NCARds633_frcfile.m 600 2020-12-29 19:19:51Z wilkin $

%% ------------------------------------------------------------------------

% USER SETS PARAMETERS IN THIS BLOCK

% shift time to the ROMS basedate you want to use
Time0 = datenum(2011,1,1);
time = E.time.data - Time0;
Tname = 'time'; % need to reset some parameters from roms_metadata

% yyyy and mm are inherited from the call to roms_get_era5_ncar_ds633
YYYY = upper(int2str(E.yyyy));
MM = upper(sprintf('%02d',E.mm));

% Set the output file name prefix
ncname = strcat('frc_watl_ERA5_bulkflux_',YYYY,MM,'.nc');
titlestr = 'ERA-5 meteorology forcing (from NCAR ds633.0) for WATL grid ';
citation = E.citation;

spherical = true;

% Coordinates and all the data necessary to write this file are in
% structure E loaded by roms_get_era5_ncar_ds633 for a requested year,
% month, and lon/lat bounding box.
lon = E.lon.data;
lat = E.lat.data;

%% ------------------------------------------------------------------------

% Code below adapted from d_ecmwf2roms.m by H. Arango and J. Wilkin in 
% the Matlab codes distributed at myroms.org. 

% Create surface forcing NetCDF files. Build creation parameters in
% structure, S. Combine all variables in monthly files.

nctype    = 'nc_float';  
Unlimited = true;  % time dimension is unlimited

ReplaceValue = NaN;
PreserveType = true;

mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET'));

fprintf('\n')
disp(['Creating ROMS NetCDF forcing file: ' ncname]);
fprintf('\n')

clear S
S.Filename = ncname;
Outfile = fullfile(pwd,S.Filename);

Im = length(lon);
Jm = length(lat);

S.Attributes(1).Name = 'type';
S.Attributes(1).Value = 'FORCING file';

S.Attributes(2).Name = 'title';
S.Attributes(2).Value = titlestr;

S.Attributes(3).Name = 'history';
S.Attributes(3).Value = ['Forcing file created with ',...
  which(mfilename) ' on ' datestr(now)];

S.Attributes(4).Name = 'citation';
S.Attributes(4).Value = citation;

% Any attributes added here (following the same syntax) become global 
% attributes of the output netcdf file
S.Attributes(5).Name = 'geospatial_lat_max';
S.Attributes(5).Value = max(lat);
S.Attributes(6).Name = 'geospatial_lat_min';
S.Attributes(6).Value = min(lat);
S.Attributes(7).Name = 'geospatial_lat_resolution';
S.Attributes(7).Value = diff(lat(1:2));
S.Attributes(8).Name = 'geospatial_lon_max';
S.Attributes(8).Value = max(lon);
S.Attributes(9).Name = 'geospatial_lon_min';
S.Attributes(9).Value = min(lon);
S.Attributes(10).Name = 'geospatial_lon_resolution';
S.Attributes(10).Value = diff(lon(1:2));
S.Attributes(11).Name = 'time_coverage_start';
S.Attributes(11).Value = datestr(E.time.data(1),31);
S.Attributes(12).Name = 'time_coverage_end';
S.Attributes(12).Value = datestr(E.time.data(end),31);
S.Attributes(13).Name = 'roms_application';
S.Attributes(13).Value = 'Geospatial extent covers domain of WATL (West Atlantic) grid';

% nc file dimensions
S.Dimensions(1).Name = 'lon';
S.Dimensions(1).Length = Im;
S.Dimensions(1).Unlimited = false;

S.Dimensions(2).Name = 'lat';
S.Dimensions(2).Length = Jm;
S.Dimensions(2).Unlimited = false;

S.Dimensions(3).Name = Tname;
S.Dimensions(3).Length = nc_constant('nc_unlimited');
S.Dimensions(3).Unlimited = true;

% add standard ROMS metadata using roms_metadata.m (myroms.org Arango)
S.Variables(1) = roms_metadata('spherical');
S.Variables(2) = roms_metadata('lon');
S.Variables(3) = roms_metadata('lat');
S.Variables(4) = roms_metadata('ocean_time', [], [], Unlimited);
S.Variables(4).Name = Tname;
S.Variables(4).Dimensions.Name = Tname;

% Edit the time variable "units" attribute for the correct reference
% time and add calendar attribute
ivar = findstrinstruct(S.Variables,'Name','time');
iatt  = strcmp({S.Variables(ivar).Attributes.Name}, 'units');
natts = length(S.Variables(ivar).Attributes);
S.Variables(ivar).Attributes(iatt).Value = ['days since ' datestr(Time0,31)];
S.Variables(ivar).Attributes(natts+1).Name  = 'calendar';
S.Variables(ivar).Attributes(natts+1).Value = 'gregorian';

% Check ROMS metadata structure.  Fill unassigned fields.

S = check_metadata(S);

% Create forcing NetCDF files. Write lon/lat coordinates

ncid = nc_create(Outfile, mode, S);
lon = repmat(lon,[1 Jm]);
lat = repmat(lat',[Im 1]);
nc_write(Outfile,'spherical',int32(spherical));
nc_write(Outfile,'lon',lon);
nc_write(Outfile,'lat',lat);

% ---------------------------
% Process a list of ROMS forcing variables

vlist = roms_varlist('bulkflux');
first = true;

for v = vlist
  
  Vroms = char(v);
  
  if first
    % add time to the file if that hasn't been done yet
    nc_write(Outfile,'time',time);
    first = false;
  end
  
  % Add the ROMS variables to the forcing file
  S.Variables(5) = roms_metadata(Vroms, spherical, nctype, Unlimited);
  S.Variables(5).Dimensions(3).Name = 'time';
  S.Variables(5).Attributes(3).Value = 'time';
  iattc = findstrinstruct(S.Variables(5).Attributes,'Name','coordinates');
  S.Variables(5).Attributes(iattc).Value = 'lon lat time';
  
  % Prepare to update some long names to be more descriptive of ERA5
  ilongname = findstrinstruct(S.Variables(5).Attributes,'Name','long_name');
  
  % write the data ----------------------------------------
  switch Vroms
    case 'Tair'
      field = E.t2.data - 273.15; % Kelvin to Celsius
      S.Variables(5).Attributes(ilongname).Value = ...
        'surface air temperature at 2 m';
    case 'Qair'
      % Compute relative humidity using Clausius-Clapeyron equation
      tsur  = E.t2.data;          % 2m temperature
      tdew  = E.d2.data;          % 2m dewpoint
      tsur  = tsur - 273.15;
      tdew  = tdew - 273.15;
      VP   = 6.11 .* 10.0 .^ (7.5 .* tdew ./ (237.7 + tdew));
      VPsat = 6.11 .* 10.0 .^ (7.5 .* tsur ./ (237.7 + tsur));
      field = 100.0 .* (VP ./ VPsat);
      S.Variables(5).Attributes(ilongname).Value = ...
        'surface air relative humidity at 2 m';
    case 'Pair'
      field = E.msl.data * 0.01;  % Pa to millibar
    case 'Uwind'
      field = E.u10.data;
      S.Variables(5).Attributes(ilongname).Value = ...
        'surface u-wind component (east) at 10 m';
    case 'Vwind'
      field = E.v10.data;
      S.Variables(5).Attributes(ilongname).Value = ...
        'surface v-wind component (north) at 10 m';
    case 'swrad'
      field = E.msnswrf.data;
      S.Variables(5).Attributes(ilongname).Value = ...
        'net solar shortwave radiation flux';
    case 'lwrad_down'
      field = E.msdwlwrf.data;
    case 'lwrad'
      field = E.msnlwrf.data;
    case 'rain'
      field = E.mtpr.data;
    otherwise
      disp(['Skipping ' Vroms])
  end
  
  % Add the variables metadata and write the data
  nc_append(Outfile,S)
  nc_write(Outfile,Vroms,field);
end

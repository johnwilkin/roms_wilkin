% Script roms_write_gfs_NCARds083_frcfile.m
%
% Create a ROMS meteorology forcing file from GFS reanalysis extracted
% from NCAR dataset ds083.3 (0.25 deg lon/lat 3-hourly data)
%
% First run e.g. E = roms_get_gfs_NCARds083_bulkflux
% and then run this script to create the forcing netcdf file.
%
% See the help on roms_get_gfs_NCARds083_bulkflux regarding obtaining login
% credentials to access the archive at NCAR Research Data Archive.
%
% -------------------------------------------------------------------------
%
% NOTE: This is an example script. The user needs to configure the ROMS
% time coordinate basedate (Time0) and the output filename (ncname) and
% title and any other pertinent metadata.
%
% This script uses routines from the myroms.org Matlab tools
% (roms_metadata, nc_create, nc_write, nc_constant etc.).
%
% John Wilkin - May 2021
%
% Copyright (c) - 2021 John L. Wilkin - jwilkin@rutgers.edu
%
% Obtain an up-to-date version of this code from 
% https://github.com/johnwilkin/roms_wilkin
%
% See also roms_get_gfs_NCARds083_bulkflux

%% ------------------------------------------------------------------------

% USER SETS PARAMETERS IN THIS BLOCK

% shift time to the ROMS basedate you want to use
Time0 = datetime(2011,1,1);
time = E.time - Time0; % this is a datetime duration
time = days(time);     % convert duration to days since ...
Tname = 'time'; % need to reset some parameters from roms_metadata

% % yyyy and mm are deduced from the time coordinate
YYYYMMDD = datestr(E.time(1),'yyyymm');

% Output file name prefix
ROMS_APP = 'watl'; % ROMS_APP = 'MYAPPLICATION'
% Output path. If not set data are written to the working directory.
Outdir = '/Volumes/home/om/roms/watl';
ncname = strcat('frc_',ROMS_APP,'_GFS_bulkflux_',YYYYMMDD,'.nc');
if ~exist('Outdir','var')
  Outdir = pwd;
end

titlestr = "GFS 0.25 3-hourly meteorology forcing (from NCAR ds083.3) for " +ROMS_APP;
citation = E.citation;

spherical = true;

% Coordinates and all the data necessary to write this file are in
% structure E loaded by roms_get_gfs_NCARds083_bulkflux for a requested 
% datetime interval and lon/lat grid
lon = E.lon(:);
lat = E.lat(:);

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
Outfile = fullfile(Outdir,S.Filename);

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
S.Attributes(11).Value = datestr(E.time(1),31);
S.Attributes(12).Name = 'time_coverage_end';
S.Attributes(12).Value = datestr(E.time(end),31);
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
  iattT = findstrinstruct(S.Variables(5).Attributes,'Name','time');
  S.Variables(5).Attributes(iattT).Value = 'time';
  iattC = findstrinstruct(S.Variables(5).Attributes,'Name','coordinates');
  S.Variables(5).Attributes(iattC).Value = 'lon lat time';
  
  % Prepare to update some long names to be more descriptive of ERA5
  ilongname = findstrinstruct(S.Variables(5).Attributes,'Name','long_name');
  
  % write the data ----------------------------------------
  switch Vroms
    case 'Tair'
      field = E.Tair;
      S.Variables(5).Attributes(ilongname).Value = ...
        'surface air temperature at 2 m';
    case 'Qair'
      % Compute relative humidity using Clausius-Clapeyron equation
      field = E.Qair;
      S.Variables(5).Attributes(ilongname).Value = ...
        'surface air relative humidity at 2 m';
    case 'Pair'
      field = E.Pair;
    case 'Uwind'
      field = E.Uwind;
      S.Variables(5).Attributes(ilongname).Value = ...
        'surface u-wind component (east) at 10 m';
    case 'Vwind'
      field = E.Vwind;
      S.Variables(5).Attributes(ilongname).Value = ...
        'surface v-wind component (north) at 10 m';
    case 'swrad'
      field = E.swrad;
      S.Variables(5).Attributes(ilongname).Value = ...
        'net solar shortwave radiation flux';
    case 'lwrad_down'
      field = E.lwrad_down;
    case 'lwrad'
      field = E.lwrad;
    case 'rain'
      field = E.rain;
    otherwise
      disp(['Skipping ' Vroms])
  end
  
  % Add the variables metadata and write the data
  nc_append(Outfile,S)
  nc_write(Outfile,Vroms,field);
end

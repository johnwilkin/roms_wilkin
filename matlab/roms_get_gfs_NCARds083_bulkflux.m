function G = roms_get_gfs_NCARds083_bulkflux(romsvname,outcoords,userpass)
% E = roms_get_gfs_NCARds083_bulkflux(romsvname,outcoords,userpass)
%
% Read GFS data via OPeNDAP from NCAR ds083.3 (0.25 deg lon/lat 3-hourly
% data) and interpolate to requested lon/lat grid and time range.
%
% Inputs:
%
%   ROMSVNAME - ROMS forcing variable name(s) to process as a cell array
%     selection {'Pair,'Tair'}, a single variable string 'Pair'
%     or the string 'all' to do them all
%   OUTCOORDS - a structure with target coordinates named lon, lat
%     and t (2-element vector start/end datetime)
%   USERPASS (string) - RDA authentication in the format:
%     'username:password' (notice the colon between username and password)
%     This string augments the OPeNDAP data URL thus:
%     url = 'https://username:password@rda.ucar.edu/thredds/dodsC/...
%                    ^^^^^^^^^^^^^^^^^
%     If your username or password includes text that would be interpretted
%     by the http protocol it must be URL encoded. In particular, since it
%     is common practise to use an email address as RDA username, the @
%     must be encoded as %40, e.g. my username and password would be
%     userpass = 'jwilkin%40rutgers.edu:mypassword' (not my real password!)
%     If you need help on this conversion, go to:
%     https://www.w3schools.com/tags/ref_urlencode.ASP
%     and enter your username or password string to find out how to URL
%     encode them to build the username:password string. DON'T URL encode
%     the colon - that will throw an error
%
% John Wilkin - May 2021
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
%
% Obtain an up-to-date version of this code from
% https://github.com/johnwilkin/roms_wilkin
%
% See also roms_write_gfs025_NCARds083_frcfile

% ERA5 analysis is username/password restricted
if nargin < 3
  try
    userpass = userpass_from_netrc;
    disp('Using username:password for rda.ucar.edu from $HOME/.netrc')
  catch
    disp(['Unsuccessful parsing username and password credentials' ...
      'for rda.ucar.edu from .netrc']);
    error('Try giving userpass string as input to this function')
  end
end
dataurl = fullfile('https://USERPASS@rda.ucar.edu',...
  'thredds/dodsC/aggregations/g/ds083.3/1/Best');
dataurl = strrep(dataurl,'USERPASS',userpass);

%% determine which WRF variables to process for ROMS

% vars.atmos_name   = {'Temperature_height_above_ground',...
%   'Relative_humidity_height_above_ground',...
%   'Pressure_surface',...
%   {'Downward_Short-Wave_Radiation_Flux_surface',...
%   'Upward_Short-Wave_Radiation_Flux_surface'},...
%   {'Downward_Long-Wave_Radp_Flux_surface',...
%   'Upward_Long-Wave_Radp_Flux_surface'},...
%   'Downward_Long-Wave_Radp_Flux_surface',...
%   'u-component_of_wind_height_above_ground',...
%   'v-component_of_wind_height_above_ground',...
%   'Precipitation_rate_surface'};

% In ds084.4/Best
% https://rda.ucar.edu/thredds/dodsC/aggregations/g/ds084.4/1/Best.html
%
% Y 'Temperature_height_above_ground',...@height_above_ground1=2 @reftime
% N 'Relative_humidity_height_above_ground',...
% N 'Pressure_reduced_to_MSL_msl',...
% Y {'Downward_Short-Wave_Radiation_Flux_surface',
% Y 'Upward_Short-Wave_Radiation_Flux_surface'},...
% Y {'Downward_Long-Wave_Radp_Flux_surface',...
% Y 'Upward_Long-Wave_Radp_Flux_surface'},...
% Y 'Downward_Long-Wave_Radp_Flux_surface',...
% Y 'u-component_of_wind_height_above_ground',... @height_above_ground=10
% Y 'v-component_of_wind_height_above_ground',... @height_above_ground=10
% N 'Precipitation_rate_surface'};
%
% But has these ...
%   'Pressure_surface'
%          do we need to also use Land_cover_0__sea_1__land_surface
%   'Specific_humidity_height_above_ground' @height_above_ground1=2
%   'Precipitation_rate_surface_Mixed_intervals_Average @time6

% Landing page for dataset ds083.3 https://rda.ucar.edu/datasets/ds083.3
% NCEP GDAS/FNL 0.25 Degree Global Tropospheric Analyses and Forecast Grids
% In ds083.3/Best
% https://rda.ucar.edu/thredds/dodsC/aggregations/g/ds083.3/1/Best.html
%
% Y 'Temperature_height_above_ground',...@height_above_ground1=2 @reftime
% Y 'Relative_humidity_height_above_ground',...
% Y 'Pressure_reduced_to_MSL_msl',...
% N {'Downward_Short-Wave_Radiation_Flux_surface',
% N 'Upward_Short-Wave_Radiation_Flux_surface'},...
% N {'Downward_Long-Wave_Radp_Flux_surface',...
% Y 'Upward_Long-Wave_Radp_Flux_surface'},...
% N 'Downward_Long-Wave_Radp_Flux_surface',...
% Y 'u-component_of_wind_height_above_ground',... @height_above_ground=10
% Y 'v-component_of_wind_height_above_ground',... @height_above_ground=10
% N 'Precipitation_rate_surface'};
%
% But has these ...
% Y 'Downward_Short-Wave_Radiation_Flux_surface_Mixed_intervals_Average
% Y 'Upward_Short-Wave_Radiation_Flux_surface_Mixed_intervals_Average
% Y 'Downward_Long-Wave_Radp_Flux_surface_Mixed_intervals_Average
% Y 'Upward_Long-Wave_Radp_Flux_surface_Mixed_intervals_Average
% Y 'Precipitation_rate_surface_Mixed_intervals_Average'

vars.atmos_name   = {...
  'Temperature_height_above_ground',...
  'Relative_humidity_height_above_ground',...
  'Pressure_reduced_to_MSL_msl',...
  {'Downward_Short-Wave_Radiation_Flux_surface_Mixed_intervals_Average',...
  'Upward_Short-Wave_Radiation_Flux_surface_Mixed_intervals_Average'},...
  {'Downward_Long-Wave_Radp_Flux_surface_Mixed_intervals_Average',...
  'Upward_Long-Wave_Radp_Flux_surface_Mixed_intervals_Average'},...
  'Downward_Long-Wave_Radp_Flux_surface_Mixed_intervals_Average',...
  'u-component_of_wind_height_above_ground',...
  'v-component_of_wind_height_above_ground',...
  'Precipitation_rate_surface_Mixed_intervals_Average'};

% match-up table of variable names
vars.roms_name = {'Tair','Qair','Pair','swrad','lwrad','lwrad_down',...
  'Uwind','Vwind','rain'};

%% loop over variable names expects a cell array
if ~iscell(romsvname)
  if strcmpi(romsvname,'all')
    romsvname = vars.roms_name;
  else
    romsvname = {romsvname};
  end
end

%% target grid to which we will interpolate the weather model
lon = outcoords.lon;
lat = outcoords.lat;
if isvector(lon)
  [lon,lat] = ndgrid(lon,lat);
end

%% open weather model from which we read the meteorology data
disp([' Reading data from ' dataurl])
U = ncinfo(dataurl);
G.dataurl = dataurl;
G.I = U;

%% weather model coordinates

% Define bounding box [west east south north] that encompasses the target
% grid (slightly inflated by lp)
lp = 0.75;
bbox = [min(lon(:)) max(lon(:)) min(lat(:)) max(lat(:))] + lp*[-1 1 -1 1];

% catch request with negative west longitudes
if all(bbox(1:2)<0)
  x = ncread(dataurl,'lon')-360;
elseif all(bbox(1:2)>0)
  x = ncread(dataurl,'lon');
else
  error('Requested longitude range straddles the prime meridian')
end
y = ncread(dataurl,'lat');

% Find minimal data set to pass to griddedInterpolant
I = find(x>=bbox(1) & x<=bbox(2));
J = find(y>=bbox(3) & y<=bbox(4));
[lon_atmos,lat_atmos] = ndgrid(double(x(I)),double(y(J)));

% data is arranged from north to south so flip latitudes to get ascending
% inputs for griddedInterpolant (no need to flip lon)
lat_atmos = flip(lat_atmos,2);


%% Process the list of requested ROMS met forcing variables

disp(' Processing input ROMS forcing variable list')

timecheck = false;

for k = 1:length(romsvname)
  rname = romsvname{k};
  disp(['  Doing ' rname])
  
  % get corresponding WRF name(s)
  kk = strcmp(rname,vars.roms_name);
  aname = vars.atmos_name{kk};
  
  % find entry in the available variables
  if iscell(aname)
    anamechk = aname{1};
  else
    anamechk = aname;
  end
  m = findstrinstruct(U.Variables,'Name',anamechk);
  
  % Parse coordinates attribute to get time coordinate name. Since this is
  % an FMRC the second entry should be the relevant coordinate variable
  c = findstrinstruct(U.Variables(m).Attributes,'Name','coordinates');
  coordinates = U.Variables(m).Attributes(c).Value;
  p = split(coordinates); % separate the words in the coordinates string
  tvarname = char(p{2});
  time = roms_get_time(dataurl,tvarname);
  tindex = find(isbetween(time,outcoords.t(1),outcoords.t(end))==1);
  time = roms_get_time(dataurl,tvarname,tindex);
  
  if timecheck
    if ~isempty(setdiff(time,time_last))
      error('time coordinate does not match previous variable')
    end
  else
    timecheck = true;
    G.time = time;
  end
  time_last = time;
  
  switch length(U.Variables(m).Size)
    case 3
      start = [I(1) J(1) tindex(1)];
      count = [length(I) length(J) length(tindex)];
    case 4 % singleton altitude dimension
      start = [I(1) J(1) 1 tindex(1)];
      count = [length(I) length(J) 1 length(tindex)];
    otherwise
      error('check number of dimensions of data')
  end
  
  if iscell(aname)
    % two variables are required for up/down radiation to form net
    disp(['   Reading ' aname{1}])
    down = ncread(dataurl,aname{1},start,count);
    disp(['   Reading ' aname{2}])
    up = ncread(dataurl,aname{2},start,count);
    data_atmos = down-up;
  else
    disp(['   Reading ' aname])
    data_atmos = ncread(dataurl,aname,start,count);
  end
  
  % special handling for units conversion
  switch rname
    case 'Tair'
      if strcmp(U.Variables(m).Attributes(2).Value,'K')
        disp('    changing units K to C')
        data_atmos = data_atmos - 273.15; % Convert K to deg C
      else
        error('Was expecting temperature in K')
      end
    case 'Pair'
      if strcmp(U.Variables(m).Attributes(2).Value,'Pa')
        disp('    changing units Pa to millibar')
        data_atmos = data_atmos*0.01; % Convert Pa to millibars
      else
        error('Was expecting pressure in Pa')
      end
    otherwise
      disp(['    got units ' U.Variables(m).Attributes(2).Value])
  end
  data_atmos = squeeze(double(data_atmos));
  data_atmos = flip(data_atmos,2); % to match new ascending order for latitude
  
  first = true; % first time index for this variable
  for tn = 1:length(tindex)
    if first
      Fn = griddedInterpolant(lon_atmos,lat_atmos,data_atmos(:,:,tn),...
        'linear','none');
      first = false;
    else
      Fn.Values = data_atmos(:,:,tn);
    end
    data_out(:,:,tn) = Fn(lon,lat);
  end
  
  G.(rname) = data_out;
  
end

G.lon = lon(:,1);
G.lat = lat(1,:);
G.description = 'https://rda.ucar.edu/datasets/ds083.3';
G.citation = [ ...
  'National Centers for Environmental Prediction/National Weather ',...
  'Service/NOAA/U.S. Department of Commerce. 2000, updated daily. ',...
  'NCEP FNL Operational Model Global Tropospheric Analyses, ',...
  'continuing from July 1999. Research Data Archive at the National ',...
  'Center for Atmospheric Research, Computational and Information ',...
  'Systems Laboratory. https://doi.org/10.5065/D6M043C6. Accessed ',...
  datestr(now)];

function index = findstrinstruct(S,field,string)
% find INDEX into a structure S for which S.FIELD matches STRING
% John Wilkin - Nov 2018
index = find(arrayfun(@(n) strcmp(S(n).(field),string), 1:numel(S)));



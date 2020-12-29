function [data,geo,info] = podaac_get_viirsL3U(sat,dnum,bbox)
% [data,geo,info] = podaac_get_viirsL3U(satellite,dnum,BoundingBox)
%
% Read VIIRS L3C (Level-3 Uncollated) SST data from PO-DAAC THREDDS
%    L3 = fixed lon/lat grid (not L2 on a changing swath).
%
% Inputs:
%
%    SATELLITE - string 'N20' or 'NPP'
%
%      There are two PO.DAAC VIIRS L3 datasets from two satellites and you
%      must choose which one to download.
%        SUOMI-NPP
%          Short Name: VIIRS_NPP-OSPO-L3U-v2.61
%          doi: 10.5067/GHVRS-3UO61
%          https://podaac.jpl.nasa.gov/dataset/VIIRS_NPP-OSPO-L3U-v2.61
%        NOAA-20
%          Short Name: VIIRS_N20-OSPO-L3U-v2.61
%          doi: 10.5067/GHV20-3UO61
%          https://podaac.jpl.nasa.gov/dataset/VIIRS_N20-OSPO-L3U-v2.61
%
%    DNUM:
%      A 2-element vector of datenums, for which all data in the interval
%         sst_time <= DNUM < sst_time is returned
%      A single value DNUM is not permitted because the data is stored in
%      numerous granules at short intervals of time.
%
% Optional inputs:
%
%    BoundingBox: 4-element vector in the format of Matlab axis
%     [minlon maxlon minlat maxlat] to extract a subregion of data.
%     Default is to return the global dataset. If there are no data
%     inside BoundingBox the result DATA is empty
%
% Outputs:
%
%   DATA returned is ** sea_surface_temperature MINUS sses_bias ***
%     Quoting the documentation:
%        sses_bias is derived against Piecewise Regression SST produced by
%        local regressions with buoys. Subtracting sses_bias from
%        sea_surface_temperature produces more accurate estimate of SST at
%        the depth of buoys. Further information at Petrenko et al., JTECH,
%        2016; doi:10.1175/JTECH-D-15-0166.1
%     PO.DAAC THREDDS only returns quality_level = 5 data as valid which is
%        recommended for operational applications
%
%   GEO is a structure with the data coordinates
%
%   INFO gives aspects of the documentation of this dataset.
%
% John Wilkin - July 2020
% $Id: podaac_get_viirsL3U.m 582 2020-10-12 15:12:20Z wilkin $
%
% EXAMPLE USAGE:
%
% bbox = [-80.46 -59.75 32.28 46.57];
% dnum = datenum(2020,7,1)+[0 1];
% [data1,geo1,info1] = podaac_get_viirsL3U('N20',dnum,bbox);
% [data2,geo2,info2] = podaac_get_viirsL3U('NPP',dnum,bbox);

if nargin < 1
  
  disp('First input must select which VIIRS satellite')
  disp('  Options are:')
  disp('   ''N20'' or ''NOAA-20'' for NOAA-20/JPSS-1 or ...')
  disp('   ''NPP'' or ''SUOMI'' for Suomi-NPP/NPOESS')
  disp('  PO.DAAC documentation is here: ')
  podaac_search = [ ...
    'https://podaac.jpl.nasa.gov/datasetlist?ids=',...
    'ProcessingLevel:SatelliteSpatialResolution&values=',...
    '*3*:%5B0.0%20TO%205.0%5D&search=VIIRS&view=list'];
  web(podaac_search)
  return
end
  
switch lower(sat(1:3))
  
  case {'npp','suo'}
    
    dataurl = fullfile('https://thredds.jpl.nasa.gov/thredds/dodsC',...
      'OceanTemperature/VIIRS_NPP-OSPO-L3U-v2.61.nc');
    infourl = 'https://podaac.jpl.nasa.gov/dataset/VIIRS_NPP-OSPO-L3U-v2.61';
    ShortName = 'VIIRS_NPP-OSPO-L3U-v2.61';
    doi = '10.5067/GHVRS-3UO61';
    
  case {'n20','noa'}
    
    dataurl = fullfile('https://thredds.jpl.nasa.gov/thredds/dodsC',...
      'OceanTemperature/VIIRS_N20-OSPO-L3U-v2.61.nc');
    infourl = 'https://podaac.jpl.nasa.gov/dataset/VIIRS_N20-OSPO-L3U-v2.61';
    ShortName = 'VIIRS_N20-OSPO-L3U-v2.611';
    doi = '10.5067/GHV20-3UO61';
    
  otherwise
    
    disp('First input must select which VIIRS satellite')
    disp('  Options are:')
    disp('   ''N20'' or ''NOAA-20'' for NOAA-20/JPSS-1 or ...')
    disp('   ''NPP'' or ''SUOMI'' for Suomi-NPP/NPOESS')
    return
    
end

if nargin < 2
  % View the OPeNDAP Data Access Form in browser window
  weburl = [dataurl '.html'];
  web(weburl)
  % View the dataset information page in a browser window
  web(infourl)
  return
end

if numel(dnum)==1
  disp('  Input DNUM must be a 2-element vector range of dates')
  return
end

lon = ncread(dataurl,'lon');
lat = ncread(dataurl,'lat');
time = double(ncread(dataurl,'time'))/86400+datenum(1981,1,1);

% find the index for time coordainte in the request date range
tindex = find(time>=dnum(1) & time<dnum(2));
if isempty(tindex)
  data = NaN;
  geo.time = dnum;
  disp('  Requested DNUM is not in time range of available data')
  return
end

if nargin < 3
  iindex = 1:length(lon);
  jindex = 1:length(lat);
else
  iindex = find(lon>=bbox(1) & lon<=bbox(2));
  jindex = find(lat>=bbox(3) & lat<=bbox(4));
end

% report
disp([' Reading ' int2str(length(tindex)) ' time records to find data'])
disp([' in time interval ' datestr(dnum(1)) ' to ' datestr(dnum(2))])
if nargin ==3
  disp([' in bounding box ' mat2str(bbox)])
end 

kount = 1;
han = waitbar(0,'Progress reading time records ...');
han.Children.FontSize = 16;

for k=1:length(tindex)
  
  waitbar(k/length(tindex),han,...
    ['Reading time record ' int2str(k) ' of ' int2str(length(tindex))])
  
  % read sea_surface_temperature in Kelvin
  sea_surface_temperature = ncread(dataurl,'sea_surface_temperature',...
    [iindex(1) jindex(1) tindex(k)],[length(iindex) length(jindex) 1]);
  
  if any(isfinite(sea_surface_temperature(:)))
    
    geo.time(kount) = time(tindex(k));
    disp(['  Got data for ' datestr(geo.time(kount))])
    
    % apply sses_bias to adjust skin temperature to bulk temperature at
    % notional buoy sensor depth
    sses_bias = ncread(dataurl,'sses_bias',...
      [iindex(1) jindex(1) tindex(k)],[length(iindex) length(jindex) 1]);
    
    % adjust for sses_bias and convert to Celsius
    data(:,:,kount) = sea_surface_temperature -272.15 - sses_bias;
    
    kount = kount+1;
    
  end
end

[geo.lon,geo.lat] = meshgrid(lon(iindex),lat(jindex));

citation = ['NOAA Office of Satellite and Product Operations (OSPO). ',...
  '2019. GHRSST Level 3U OSPO dataset v2.61 from VIIRS on S-NPP ',...
  'Satellite (GDS v2). Ver. 2.61. PO.DAAC, CA, USA. Dataset accessed ',...
  datestr(now,'YYYY-MM-DD') ' at https://doi.org/' doi];

info.These_Data_Are_Adjusted_To_Buoy_Depth_Using = ...
  'sea_surface_temperature minus sses_bias';
info.infoURL = infourl;
info.ShortName = ShortName;
info.doi = doi;
info.citation = citation;
gatlist = {'acknowledgement','creator_name','creator_email','history',...
  'institution','platform','sensor','id'};
for gat = gatlist
  attname = char(gat);
  info.(char(attname)) = ncreadatt(dataurl,'/',attname);
end

close(han)




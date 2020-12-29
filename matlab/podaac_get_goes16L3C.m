function [data,geo,info] = podaac_get_goes16L3C(dnum,bbox)
% [data,geo,info] = podaac_get_goes16L3C(dnum,BoundingBox)
%
% Read GOES-16 L3C (Level-3 Collated) SST data from PO-DAAC THREDDS
%    Level-3 means on a fixed lon/lat grid (not a changing swath).
%    Collated means gather together passes from a single sensor - here
%    passes meaning multiple full disk scans from geostationary orbit.
%    This is PO.DAAC dataset with Short Name ABI_G16-STAR-L3C-v2.70 and
%    doi: 10.5067/GHG16-3UO27. Full documentation at:
%    https://podaac.jpl.nasa.gov/dataset/ABI_G16-STAR-L3C-v2.70
%
% Inputs:
%    DNUM:
%      A datenum for which to extract data. The data returned is for the
%         closest observed time, with the time offset reported in GEO to
%         allow a check on whether the offset is out of desired scope.
%      A 2-element vector of datenums, for which all data in the interval
%         sst_time <= DNUM < sst_time is returned
%      DNUM = Inf will the latest set of data.
%      The DNUM input is required because we don't want the entire THREDDS
%      catalog of 20,000 scanes of data.
%
% Optional inputs:
%    BoundingBox: 4-element vector in the format of Matlab axis
%     [minlon maxlon minlat maxlat] to extract a subregion of data.
%     Default is to return the global dataset. If there are no data
%     inside BoundingBox the result DATA is empty
%
% Outputs:
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
%   INFO gives aspects of the documentation of this dataset.
%
% John Wilkin - June 2020
% $Id: podaac_get_goes16L3C.m 586 2020-10-12 20:21:53Z wilkin $
%
% EXAMPLE USAGE: 
%
% bbox = [-80.46 -59.75 32.28 46.57];
% dnum = datenum(2020,7,1)+[0 1];
% [data,geo,info] = podaac_get_goes16L3C(dnum,bbox);

dataurl = fullfile('https://thredds.jpl.nasa.gov/thredds/dodsC',...
  'OceanTemperature/ABI_G16-STAR-L3C-v2.70.nc');

infourl = 'https://podaac.jpl.nasa.gov/dataset/ABI_G16-STAR-L3C-v2.70';
ShortName = 'ABI_G16-STAR-L3C-v2.70';
doi = '10.5067/GHG16-3UO27';

if nargin < 1
  % View the OpenDAP Data Access Form in browser window
  weburl = [dataurl '.html'];
  web(weburl)
  % View the dataset information page in a browser window
  web(infourl)
  return
end

lon = ncread(dataurl,'lon');
lat = ncread(dataurl,'lat');
time = double(ncread(dataurl,'time'))/86400+datenum(1981,1,1);

if numel(dnum) > 1
  tindex = find(time>=dnum(1) & time<dnum(2));
  if isempty(tindex)
    data = NaN;
    geo.time = dnum;
    geo.message = 'Requested DNUM is not in time range of data';
  else
    geo.time = time(tindex);
  end
else
  [~,tindex] = min(abs(time-dnum));
  geo.time = time(tindex);
  geo.time_offset = dnum-time(tindex);
  if (geo.time_offset)>3
    disp('Warning: Data are more than 3 days older than requested date')
  end
end

if nargin < 2
  iindex = 1:length(lon);
  jindex = 1:length(lat);
else
  iindex = find(lon>=bbox(1) & lon<=bbox(2));
  jindex = find(lat>=bbox(3) & lat<=bbox(4));
end

% sea_surface_temperature in Kelvin
sea_surface_temperature = ncread(dataurl,'sea_surface_temperature',...
  [iindex(1) jindex(1) tindex(1)],[length(iindex) length(jindex) ...
  length(tindex)]);

% sses_bias to adjust skin temperature to bulk temperature at notional buoy
% sensor depth
sses_bias = ncread(dataurl,'sses_bias',...
  [iindex(1) jindex(1) tindex(1)],[length(iindex) length(jindex) ...
  length(tindex)]);

data = sea_surface_temperature -272.15 - sses_bias;

[geo.lon,geo.lat] = meshgrid(lon(iindex),lat(jindex));

citation = ['NOAA/NESDIS USA, 5200 Auth Rd, Camp Springs, MD, 20746. ',...
  '2019. GHRSST NOAA/STAR GOES-16 ABI L3C America Region SST v2.70 ',...
  'dataset in GDS2. Ver. 2.70. PO.DAAC, CA, USA. Dataset accessed ',...
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





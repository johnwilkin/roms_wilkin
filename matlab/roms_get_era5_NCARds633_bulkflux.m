function E = roms_get_era5_NCARds633_bulkflux(yyyy,mm,bbox,userpass)
% E = roms_get_era5_NCARds633_bulkflux(yyyy,mm,bbox,userpass)
%
% Read ECMWF ERA5 meteorological reanalysis from the NCAR Research Data 
% Archive (RDA) dataset ds633.0 https://rda.ucar.edu/datasets/ds633.0
%
% Use this function to extract data for a chosen month, then use 
% roms_write_era5_NCARds633_frcfile to create and write the ROMS format 
% surface forcing netcdf file. 
%
% This function reads everything ROMS needs to build a surface forcing 
% file for the BULK_FLUX option, i.e. marine boundary layer pressure, 2-m 
% air temperature, 2-m dew point (for conversion to relative humidity), 
% 10-m winds, net shortwave radiation, and both net longwave and downward 
% longwave radiation from which ROMS selects depending on whether the ROMS
% user has #define LONGWAVE_OUT or not.
%
% ERA5 data are freely available but you must be a registered user at
% rda.ucar.edu to obtain access. See the "Register Now" link at the top 
% left of https://rda.ucar.edu to obtain a login). 
%
% Login credentials must be passed to this function as the input USERPASS.
% If they not given as an input, the function attempts to parse the 
% information from your .netrc file. 
%
% This NCAR ds633 version of ERA5 is global 0.25 by 0.25 deg lon/lat at 
% hourly intervals. It is a rolling archive that is updated monthly and 
% has data available from 1979 up until now minus approximetely 4 months.  
%
% Inputs:
%
%   yyyy (scalar) - the year number
%
%   mm (scalar) - the month number
%
%   bbox - A 4-element vector in the format of Matlab axis defining the
%      lon/lat bounding box region to subset with the OPeNDAP query
%      e.g. bbox = [-110 -30 0 55] for West Atlantic Model. 
%      NOTE: The longitude coordinate in ERA5 breaks at the prime meridian. 
%      This function detects whether in input longitudes in BBOX are 
%      negative (west of prime meridian) or positive (east of prime 
%      and adjusts accordingly. But the query cannot stradle the prime 
%      meridian. The user will have to make two files and merge them. If I
%      ever have a project in the east Atlantic I might code this merger
%      automatically, until then you are on your own.
%
%   userpass (string) - RDA authentication in the format:
%     'username:password' (notice the colon between username and password)
%     This string augments the OPeNDAP data URL thus:
%     url = 'https://username:password@rda.ucar.edu/thredds/dodsC/...
%                    ^^^^^^^^^^^^^^^^^
%     and therefore if your username or password includes text that would 
%     be interpretted by the http protocol it must be URL encoded. In
%     particular, since it is common practise to use an email address as
%     RDA username, the @ must be encoded as %40, e.g. my username and 
%     password would be entered as 
%     userpass = 'jwilkin%40rutgers.edu:mypassword' (not my real password!)
%
% Outputs:
%
%   E (structure) - contains all the information required by script
%   roms_write_era5_NCARds633_frcfile.m to generate a ROMS format forcing
%   netcdf file.
%
%   Note: Data for entire year/month requested here are returned at hourly
%   intervals. MIDNIGHT ON THE LAST DAY OF THE MONTH IS NOT INCLUDED. That
%   will be in the next file in the year/month sequence passed to ROMS
%   input.
%
% NOTE on Matlab versions:
%   This code runs successfully with Matlab Release 2020b (9.9.0.1467703).
%   Some Matlab versions throw an error in ncread when reading flux 
%   variables (longwave radiation etc.). Circumstantial evidence links this
%   to ncread built with NetCDF 4.6.1. You can check the version of NetCDF 
%   that your Matlab is using with: >> libvers = netcdf.inqLibVers
%    
% John Wilkin - December 2020
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_get_era5_NCARds633_bulkflux.m 589 2020-12-28 21:50:00Z wilkin $
%
% Obtain an up-to-date version of this code from 
% https://github.com/johnwilkin/roms_wilkin
%
% See also roms_write_era5_NCARds633_frcfile

% ERA5 analysis is username/password restricted
if nargin < 4
  try
  userpass = userpass_from_netrc;
  disp('Using username:password for rda.ucar.edu from $HOME/.netrc')
  catch 
    disp(['Unsuccessful parsing username and password credentials' ...
      'for rda.ucar.edu from .netrc']); 
    error('Try giving userpass string as input to this function')
  end
end

urlbase = 'rda.ucar.edu/thredds/dodsC/files/g/ds633.0';
server = strcat('https://',userpass,'@',urlbase,'/');

% ERA5 data in this archive use time since 01-01-1900
epoch = datenum(1900,1,1);
YYYYMM = [sprintf('%d',yyyy) sprintf('%02d',mm);];
YYYY = sprintf('%d',yyyy);

% All the ERA5 variables needed to ROMS BULK_FLUX option
ecmwf_bulkflux_vars = {'msl','t2','d2','u10','v10',...
  'msdwlwrf','msnlwrf','msnswrf','mtpr'};

% ------------------------------------------------------------------------
%
% Structure E describes variable names and other information
% for building the data access query

clear E

E.u10.long = '10-m wind (u)';
E.u10.name = 'VAR_10U';
E.u10.code = '165';
E.u10.units = 'm s-1';
E.u10.set = 'e5.oper.an.sfc'; % 1 file
E.u10.v = '10u';

E.v10.long = '10-m wind (v)';
E.v10.name = 'VAR_10V';
E.v10.code = '166';
E.v10.units = 'm s-1';
E.v10.set = 'e5.oper.an.sfc'; % 1 file
E.v10.v = '10v';

E.d2.long = '2 meter dewpoint temperature';
E.d2.name = 'VAR_2D';
E.d2.code = '168';
E.d2.units = 'K';
E.d2.set = 'e5.oper.an.sfc'; % 1 file
E.d2.v = '2d';

E.t2.long = '2 meter temperature';
E.t2.name = 'VAR_2T';
E.t2.code = '167';
E.t2.units = 'K';
E.t2.set = 'e5.oper.an.sfc'; % 1 file
E.t2.v = '2t';

E.msl.long = 'mean sea-level pressure';
E.msl.name = 'MSL';
E.msl.code = '151';
E.msl.units = 'Pa';
E.msl.set = 'e5.oper.an.sfc'; % 1 file
E.msl.v = 'msl';

E.msshf.long = 'mean surface sensible heat flux';
E.msshf.name = 'MSSHF';
E.msshf.code = '033';
E.msshf.units = 'W m-2';
E.msshf.set = 'e5.oper.fc.sfc.meanflux'; % 3 files
E.msshf.v = 'msshf';

E.mslhf.long = 'mean surface latent heat flux';
E.mslhf.name = 'MSLHF';
E.mslhf.code = '034';
E.mslhf.units = 'W m-2';
E.mslhf.set = 'e5.oper.fc.sfc.meanflux';  % 3 files
E.mslhf.v = 'mslhf';

E.msdwlwrf.long = 'mean surface downward long-wave radiation flux';
E.msdwlwrf.name = 'MSDWLWRF';
E.msdwlwrf.code = '036';
E.msdwlwrf.units = 'W m-2';
E.msdwlwrf.set = 'e5.oper.fc.sfc.meanflux';  % 3 files
E.msdwlwrf.v = 'msdwlwrf';

E.msnlwrf.long = 'mean surface net long-wave radiation flux';
E.msnlwrf.name = 'MSNLWRF';
E.msnlwrf.code = '038';
E.msnlwrf.units = 'W m-2';
E.msnlwrf.set = 'e5.oper.fc.sfc.meanflux';  % 3 files
E.msnlwrf.v = 'msnlwrf';

E.msnswrf.long = 'mean surface net short-wave radiation flux';
E.msnswrf.name = 'MSNSWRF';
E.msnswrf.code = '037';
E.msnswrf.units = 'W m-2';
E.msnswrf.set = 'e5.oper.fc.sfc.meanflux';  % 3 files
E.msnswrf.v = 'msnswrf';

E.mtpr.long = 'mean total precipitation rate';
E.mtpr.name = 'MTPR';
E.mtpr.code = '055';
E.mtpr.units = 'kg m-2 s-1';
E.mtpr.set = 'e5.oper.fc.sfc.meanflux';  % 3 files
E.mtpr.v = 'mtpr';

%% ------------------------------------------------------------------------

for vname = ecmwf_bulkflux_vars
  
  v = char(vname);
  
  % Build data urls -------------------------------------------------------
  
  disp(' ')
  disp(['Processing variable ' E.(v).long])
  
  clear files
  
  % ERA5 THREDSS data organization:
  %
  % Marine Boundary Layer variables are from analysis in monthly files. 
  %    One file contains the entire month of 1-hourly data.
  %
  % Flux variables are from forecast files organized in half-month files
  %    but with 6 hours of the 1st of the month in the previous month's
  %    ddile. A full month of data therefore requires access to three 
  %    files. Flux variables are given as mean rates during each forecast
  %    hour, as opposed to time accumulations used at other ECMWF sources.
  
  switch E.(v).set
    
    case 'e5.oper.fc.sfc.meanflux'
      
      % Flux forecasts are in two files for the 1st to 16th, and 16th to 
      % end of month, but the first 6 hours are actually in the preceding 
      % month file. To read a complete calendar month therefore requires 
      % reading from 3 files
      
      varf = strcat(E.(v).set,'.235_',E.(v).code,'_',E.(v).v,'.ll025sc.');
      
      % build string for the time span for each of the 3 files
      switch mm
        
        case 1
          
          % start December previous year
          drange = [sprintf('%d',yyyy-1) '121606_' YYYYMM '0106'];
          files(1).url = strcat(server,E.(v).set,'/',...
            sprintf('%d',yyyy-1),'12','/',varf,drange,'.nc');
 
          drange = [YYYYMM '0106_' YYYYMM '1606'];
          files(2).url = strcat(server,E.(v).set,'/',YYYYMM,'/',...
            varf,drange,'.nc');
          
          drange = [YYYYMM '1606_' YYYY sprintf('%02d',mm+1) '0106'];
          files(3).url = strcat(server,E.(v).set,'/',YYYYMM,'/',...
            varf,drange,'.nc');
          
        case 12
          
          drange = [YYYY sprintf('%02d',mm-1) '1606_' YYYYMM '0106'];
          files(1).url = strcat(server,E.(v).set,'/',YYYY,...
            sprintf('%02d',mm-1),'/',varf,drange,'.nc');
          
          drange = [YYYYMM '0106_' YYYYMM '1606'];
          files(2).url = strcat(server,E.(v).set,'/',YYYYMM,'/',...
            varf,drange,'.nc');
          
          % through January next year
          drange = [YYYYMM '1606_' sprintf('%d',yyyy+1) '010106']; 
          files(3).url = strcat(server,E.(v).set,'/',YYYYMM,'/',...
            varf,drange,'.nc');
          
        otherwise
          
          % previous month
          drange = [YYYY sprintf('%02d',mm-1) '1606_' YYYYMM '0106'];
          files(1).url = strcat(server,E.(v).set,'/',YYYY,...
            sprintf('%02d',mm-1),'/',varf,drange,'.nc');
          
          drange = [YYYYMM '0106_' YYYYMM '1606'];
          files(2).url = strcat(server,E.(v).set,'/',YYYYMM,'/',...
            varf,drange,'.nc');
          
          drange = [YYYYMM '1606_' YYYY sprintf('%02d',mm+1) '0106'];
          files(3).url = strcat(server,E.(v).set,'/',YYYYMM,'/',...
            varf,drange,'.nc');
          
      end
      
    case 'e5.oper.an.sfc'
      
      varf = strcat(E.(v).set,'.128_',E.(v).code,'_',E.(v).v,'.ll025sc.');
      
      % build string defining the month time span
      ld = datenum(yyyy,mm+1,1)-datenum(yyyy,mm,1); % last day of month     
      drange = [YYYYMM '0100_' YYYYMM sprintf('%02d',ld) '23'];
      files(1).url = strcat(server,E.(v).set,'/',YYYYMM,'/',...
        varf,drange,'.nc');
      
  end
  
  % -----------------------------------------------------------------------
  % read spatial coordinates subset
  
  url = files(1).url;
  
  lon = ncread(url,'longitude')-360;
  Is = find(lon>=bbox(1),1,'first');
  Ie = find(lon<=bbox(2),1,'last');
  if isempty(Is) || isempty(Ie)
    disp('No longitudes -360 to 0 found in bounding box limits')
    warning('Adding 360 to longitude coordinate and trying again')
    lon = lon+360;
    Is = find(lon>=bbox(1),1,'first');
    Ie = find(lon<=bbox(2),1,'last');
    if isempty(Is) || isempty(Ie)
      error('Failed to find longitude points in bounding box limits')
    end
  end
  Ilen = Ie-Is;
  
  lat = ncread(url,'latitude');
  Js = find(lat<=bbox(4),1,'first');
  Je = find(lat>=bbox(3),1,'last');
  Jlen = Je-Js;
  
  lon = ncread(url,'longitude',Is,Ilen)-360;
  lat = ncread(url,'latitude',Js,Jlen);
  lat = flip(lat);
 
  % read time coordinates and data subset
  
  switch E.(v).set
    
    case 'e5.oper.fc.sfc.meanflux'
      
      for k=1:3 % 3 files to read
        
        url = files(k).url;
        
        switch k
          
          case 1
            
            I = ncinfo(url);
            dim = findstrinstruct(I.Dimensions,...
              'Name','forecast_initial_time');
            Last = I.Dimensions(dim).Length;
            
            % Get last 7 hours 
            itime = double(ncread(url,'forecast_initial_time',Last,1))/24 ...
              + epoch;
            ni = length(itime);
            fhour = double(ncread(url,'forecast_hour',6,7))/24;
            nf = length(fhour);
            time = repmat(itime',[nf 1])+repmat(fhour,[1 ni]);
            disp(' Reading file 1 of 3 ...  ')
            tic
            data = ncread(url,E.(v).name,[Is Js 6 Last],[Ilen Jlen 7 1]);
            data = flip(data,2);
            fprintf('\b')
            toc
            
            TIME = time(:);
            DATA = data(:,:,:);
            
          case 2
            
            % Get the full 16 days
            % I have tested getting this in smaller chunks, but the fastest
            % method is a single query loading all initial times and forecast
            % hours at once.
            itime = double(ncread(url,'forecast_initial_time'))/24 + epoch;
            ni = length(itime);
            fhour = double(ncread(url,'forecast_hour'))/24;
            nf = length(fhour);
            time = repmat(itime',[nf 1])+repmat(fhour,[1 ni]);
            disp('  Reading file 2 of 3 ...  ')
            tic
            data = ncread(url,E.(v).name,[Is Js 1 1],[Ilen Jlen Inf Inf]);
            data = flip(data,2);
            fprintf('\b')
            toc
            
            TIME = cat(1,TIME,time(:));
            DATA = cat(3,DATA,data(:,:,:));
            
          case 3
            
            % Get the full 16 days
            % Then trim the extra 6 hours of the next month
            itime = double(ncread(url,'forecast_initial_time'))/24 + epoch;
            ni = length(itime);
            fhour = double(ncread(url,'forecast_hour'))/24;
            nf = length(fhour);
            time = repmat(itime',[nf 1])+repmat(fhour,[1 ni]);
            disp('   Reading file 3 of 3 ...  ')
            tic
            data = ncread(url,E.(v).name,[Is Js 1 1],[Ilen Jlen Inf Inf]);
            data = flip(data,2);
            fprintf('\b')
            toc
            
            % Trim off the extra 7 hours that are from the next month so
            % data ends at 23:00 and there is no duplication of timne in a
            % multifile sequence
            time(end+(-6:0)) = [];
            data(:,:,end+(-6:0)) = [];
            
            TIME = cat(1,TIME,time(:));
            DATA = cat(3,DATA,data(:,:,:));
            
        end % switch
      end % files
      
      E.(v).time = TIME;
      E.(v).data = DATA;
      
    case 'e5.oper.an.sfc'
      
      % 1 file to read
      url = files(1).url;
      I = ncinfo(url);
      time = epoch + double(ncread(url,'time'))/24;
      TIME = time(:);
      disp(' Reading 1 file ...')
      Nt = I.Dimensions(1).Length;
      clear data
      tic
      ndays = 16;
      fprintf(['  in ' int2str(length(1:(ndays*24):Nt)) ' ' ...
        int2str(ndays) '-day chunks: '])
      count = 1;
      for ch=1:(ndays*24):Nt
        fprintf('%d',count)
        Tlen = min(ndays*24,Nt-ch+1);
        data(:,:,ch-1+(1:Tlen)) = ncread(url,E.(v).name,[Is Js ch],...
          [Ilen Jlen Tlen]);
        count = count+1;
        fprintf('\b')
      end
      data = flip(data,2);
      fprintf('  ')
      toc
      
      DATA = data;
      E.(v).time = TIME;
      E.(v).data = DATA;
      
  end
end

E.time.data = E.(v).time;
E.lon.data = lon;
E.lat.data = lat;
E.yyyy = yyyy;
E.mm = mm;
E.description = 'https://rda.ucar.edu/datasets/ds633.0';

E.citation = [ ...
  'European Centre for Medium-Range Weather Forecasts, 2019,' ,...
  'updated monthly. ERA5 Reanalysis (0.25 Degree Latitude-Longitude ',...
  'Grid). Research Data Archive at the National Center for ',...
  'Atmospheric Research, Computational and Information Systems ',...
  'Laboratory. https://doi.org/10.5065/BH6N-5N20. Accessed ' datestr(now)];

function index = findstrinstruct(S,field,string)
% index = findstrinstruct(S,field,string)
%
% find INDEX into a structure S for which S.FIELD matches STRING
%
% John Wilkin - Nov 2018
% from roms_wilkin toolbox
index = find(arrayfun(@(n) strcmp(S(n).(field),string), 1:numel(S)));

function userpass = userpass_from_netrc
% attempt to parse username:password for rad.ucar./edu from .netrc
[~,HOME] = system('echo $HOME');
fid = fopen(fullfile(strip(HOME),'.netrc'));
while 1  
  str = fgetl(fid);
  if ~ischar(str)
    break
  else
    if contains(str,'machine') && contains(str,'rda.ucar.edu')
      % the next 2 lines of .netrc are login and password
      str = strrep(fgetl(fid),'login','');
      user = strrep(str,'@','%40');
      pass = strrep(fgetl(fid),'password','');
      userpass = strcat(strip(user),':',strip(pass));
    end
  end
end

  




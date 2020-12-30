function url = maracoos(varargin)
% url = maracoos(options)
% Get the THREDDS Data Server dataset URL for output from the MARACOOS.org
% DOPPIO data assimilative ocean analysis and forecast system
%
% Not all of these options are yet available, or openly accessible
%
% Options
%    A sequence of string arguments to specify the catalog
%    If no options are given, default is 'rt' and 'his' (see below)
%
%    'rt'   Forecast Model Run Collection (FMRC) best time series
%    'ra' or 'reanalysis' V1R3 historical reanalysis
%
%    'his'  1-hourly shapshots (ROMS 'history' files)
%    'avg'  1-day averages
%
%    'mo' monthly means of reanalysis
%    'an' annual means of reanalysis
%    'me' or 'em' monthly ensemble means of reanalysis
%
% Example usage
%    dataurl = maracoos('rt','avg'); % for daily averages from NRT system
%
% John Wilkin - March 2019
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu

root = 'http://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio';
instance = '2017_da/his/History_Best';
url = fullfile(root,instance);
avg = false;

for k=1:nargin
  opt = varargin{k};
  switch lower(opt(1:2))
    case 'rt'
      instance = '2017_da/his/History_Best';
      url = fullfile(root,instance);
      disp('MARACOOS real-time')
    case {'re','ra'}
      instance = 'restricted_access/DopAnV1R3-ini2007_da/his';
      url = fullfile(root,instance);
      disp('V1R3 reanalysis')
    case 'mo'
      instance = 'restricted_access/DopAnV1R3-ini2007_da/mon_avg';
      url = fullfile(root,instance);
      disp('V1R3 reanalysis month by month averages')
    case {'me','em'}
      instance = 'restricted_access/DopAnV1R3-ini2007_da/mon_ens_means';
      url = fullfile(root,instance);
      disp('V1R3 reanalysis monthly ensemble average')
    case {'an','ye'}
      instance = 'restricted_access/DopAnV1R3-ini2007_da/year_avg';
      url = fullfile(root,instance);
      disp('V1R3 reanalysis annual averages')
    case 'av'
      avg = true;
    case 'hi'
      avg = false;
    otherwise
      disp(['Huh? Did not understand option ' opt])
  end
end
if avg
  disp('averages')
else
  disp('history')
end
if avg
  url = strrep(url,'his','avg');
  url = strrep(url,'History','Averages');
end



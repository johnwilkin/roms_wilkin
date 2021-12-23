function [dnum,tindex] = roms_get_time(file,varargin)
% [dtime,tindex] = roms_get_time(file,[tvarname],tindex,[tlen],[tstride])
%
% Read time from netcdf file and convert to a DATETIME obect
% Parses the UNITS string to get base date and conversion factor to days.
%
% Inputs:
%    FILE      filename or data URL
%
% Optional inputs:
%    TVARNAME  name of the time variable
%              If not given, scan the file for variables named
%                 'time' then 'ocean_time' then 'bry_time' then 'frc_time'
%    TINDEX    If not given, get all times
%              If TINDEX is a scalar, get this value
%              If TINDEX is a vector, get all these values ... but the
%                 STRIDE in TINDEX must be constant (start:N:end)
%    TLEN      If TLEN is present, then get TLEN values starting at TINDEX
%                 but if TLEN < 0 get TLEN values prior to TINDEX
%    TSTRIDE   If TSTRIDE is present, get TLEN values skipping TSTRIDE
%              
% Examples:
%
%    dtime = roms_get_time(file);                   % get all dates
%    dtime = roms_get_time(file,Inf);               % get all dates
%    dtime = roms_get_time(file,'ocean_time',10,3); % get dates 10, 11, 12
%    dtime = roms_get_time(file,10,3,2);            % get dates 10, 12, 14
%    dtime = roms_get_time(file,'last',-7);         % get last 7 dates
%
%    Auxiliary output TINDEX is the index of the first output DTIME
%
% THIS FUNCTION replaces roms_get_date
%
% John Wilkin - February 2021
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
%
% See also ROMS_GET_TIME_INDEX

if isempty(varargin) % get all times in file
  tvarname = [];
  tindex = Inf;
else
  if ~isempty(strfind(varargin{1},'time')) % time name was given
    tvarname = varargin{1};
    varargin = varargin(2:end); % varargin now tindex, [tlen], [tstride]
  else % time name not given - must figure it out
    tvarname = [];
  end
end

if isempty(tvarname)
  timenamelist = {'time','ocean_time','bry_time','frc_time','sea_time'};
  % scan file for time coordinate variable
  I = ncinfo(file);
  for tn = timenamelist
    tnstr = char(tn);
    if findstrinstruct(I.Variables,'Name',tnstr)
      tvarname = tnstr;
      break
    end
  end
end
if isempty(tvarname)
  error("Unable to deduce time coordinate variable name in" +newline+file)
end

switch length(varargin)
  case 0
    tindex = Inf;
  case 1
    tindex = varargin{1};
    tlen = 1;
  case 2
    tindex = varargin{1};
    tlen = varargin{2};
  case 3
    tindex = varargin{1};
    tlen = varargin{2};
    tstride = varargin{3};
end

if ischar(tindex)
  switch tindex
    case {'last','latest','end'}
      I = ncinfo(file);
      tindex = I.Variables(findstrinstruct(I.Variables,'Name',tvarname)).Size;
    otherwise
      error("Unable to comprehend tindex string: " + tindex)
  end
end
if isinf(tindex)
  tindex = 1;
  tlen = Inf;
end
if numel(tindex) > 1 % a vector of time indices was input
  % but it must have a constant stride for this to work
  tlen = length(tindex);
  dt = unique(diff(tindex));
  if numel(dt)~=1
    error('TINDEX values must step by a constant STRIDE')
  end
  if dt~=1
    tstride = dt;
  end
end

% tlen < 0 is interpretted as get the values PRIOR to tindex
if tlen < 0
  tindex(1) = tindex(1)+tlen+1;
  tlen = abs(tlen);
end

if exist('tstride','var')
  time = ncread(file,tvarname,tindex(1),tlen,tstride);
else
  time = ncread(file,tvarname,tindex(1),tlen);
end
time = double(time);

units = ncreadatt(file,tvarname,'units');
switch lower(units(1:4))
  case 'mill'
    fac = 86400e3;
  case 'seco'
    fac = 86400;
  case 'minu'
    fac = 1440;
  case 'hour'
    fac = 24;
  case 'days'
    fac = 1;
end

% convert to datenum
try
  tsince = units(6+strfind(units,'since'):end);
  % some date formats aren't automatically parsed by datenum
  tsince = strrep(tsince,'UTC','');
  tsince = strrep(tsince,'Z','');
  tsince = strrep(tsince,'T',' ');
  if strncmp(tsince,'0001-01-01',10)||strncmp(tsince,'0000-00-00',10)
    % assume not a real calendar date and just convert model time to days
    dnum = time/fac;
  else
    dnum = time/fac + datenum(tsince);
    dnum = datetime(dnum,'ConvertFrom','datenum');
  end
catch
  warning("Unable to parse base date from units string:" + newline + units)
  dnum = time/fac;
end

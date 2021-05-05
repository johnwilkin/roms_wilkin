function tindex = roms_get_time_index(file,varargin)
% tindex = roms_get_time_index(file,[tvarname],date)
%
% Read time from netcdf file and find the index closest to requested DATE
% Parses the UNITS string to get base date and conversion factor to days.
%
% Inputs:
%    FILE      filename or data URL
%
% Optional inputs:
%    TVARNAME  name of the time variable
%              If not given, scan the file for variables named
%                 'time' then 'ocean_time' then 'bry_time' then 'frc_time'
%    DATE      a DATE STRING that will be parsed to a DATENUM, or
%              'last', 'latest' or 'end', or
%              a numeric value, in which case the index to the nearest
%                 time value is returned (with no regard for dates)
%
% John Wilkin - May 2021
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
%
% See also ROMS_GET_TIME

switch length(varargin)
  case 1
    tvarname = [];
    dstr = varargin{1};
  case 2
    tvarname = varargin{1};
    dstr = varargin{2};
end

if isempty(tvarname)
  timenamelist = {'time','ocean_time','bry_time','frc_time'};
  % scan file for time coordinate variable
  I = ncinfo(file);
  for tn = timenamelist
    tnstr = char(tn);
    if findstrinstruct(I.Variables,'Name',tn)
      tvarname = tnstr;
      break
    end
  end
end
if isempty(tvarname)
  error("Unable to deduce time coordinate variable name in" +newline+file)
end

% read the time coordinate variable
time = double(ncread(file,tvarname));

if ~ischar(dstr)
  % input is a time coordinate VALUE
  if dstr < min(time) || dstr > max(time)
    error("Requested value " + num2str(dstr,'%8.2f') + ...
      " is not in the range " + num2str(time(1),'%8.2f') + " to " ...
      + num2str(time(end),'%8.2f'))
  else
    [~,tindex] = min(abs(dstr-time));
    tindex = tindex(1);
  end
else
  switch dstr
    case {'last','latest','end'}
      I = ncinfo(file);
      tindex = I.Variables(findstrinstruct(I.Variables,'Name',tvarname)).Size;
    otherwise
      try
        dnumin = datenum(dstr);
      catch
        error("Failed to convert requested date '" + dstr + "' into a DATENUM")
      end
      % convert times in file to datenum
      units = ncreadatt(file,tvarname,'units');
      switch units(1:4)
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
        end
      catch
        warning("Unable to parse base date from units string:" + newline + units)
        dnum = time/fac;
      end
      if dnumin < min(dnum) || dnumin > max(dnum)
        error("Requested date '" + dstr + "' is not in the range " ...
          + datestr(min(dnum)) + " to " + datestr(max(dnum)))
      else
        [~,tindex] = min(abs(dnum-dnumin));
        tindex = tindex(1);
      end
  end
end

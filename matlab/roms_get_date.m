function [Dnum,Dstr] = roms_get_date(file,tindex,dateformat,timevarname)
% $Id: roms_get_date.m 583 2020-10-12 20:19:38Z wilkin $
% [dnum,dstr] = roms_get_date(file,tindex,dateformat,timevarname)
%
% Read ocean_time value from FILE for time index TINDEX and try to parse
% this into a Matlab datenum value. TINDEX starts at 1.
%
% If TINDEX == -1 get all times in the FILE
%
% If a 3rd input is given, interpret this as a format option to datestr and 
% convert the datenum value into a string with this format. If the
% dateformat < 0 the string returned will simply be "day ...."
%
% John Wilkin
% Uses snctools and function parsetnc
%
% WOULD BE GOOD TO DO AWAY WITH PARSETNC <<<<<<<<<< !!!!!!! <<<<<<<<<<<<<
% and allow time = Inf to get the latest value in the file
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu

if nargin < 2
  tindex = -1;
end

if nargin < 3
  dateformat = 0;
  basedate = NaN;
end

if nargin < 4
  if nc_isvar(file,'ocean_time')
    % for most roms output
    timevarname = 'ocean_time';
    if nc_isvar(file,'time')
      % if both ocean_time and time then probably a fmrc
      timevarname = 'time';
    end
  elseif nc_isvar(file,'time')
    % probably fmrc
    timevarname = 'time';
  else
    error('Can''t determine what time variable to use for date')
  end
end

if tindex == -1
  % get all times
  ocean_time  = nc_varget(file,timevarname,0,-1);
else
  %ocean_time  = nc_varget(file,timevarname,0,-1);
  %ocean_time  = ocean_time(tindex);
  ocean_time  = nc_varget(file,timevarname,tindex-1,1);
end
tunits = nc_attget(file,timevarname,'units');

switch tunits(1:3)
  case 'day'
    fac = 1;
  case 'hou'
    fac = 1/24;
  case 'sec'
    fac = 1/86400;
  otherwise
    warning('Cannot interpret tunits string:')
    disp(tunits)
    disp('Making no time units assumption - returning ocean_time from file')
    fac = 0;
end

% convert to days
if fac == 0
  dnum = ocean_time;
  if nargout > 1
    dstr = [num2str(ocean_time) ' ' tunits];
  end
else
  % try to parse the time units string
  try
    basedate = datenum(parsetnc(tunits));
  catch
    try
      tunits = tunits((strfind(tunits,'since')+6):end);
      basedate = datenum(tunits);
    catch
      try
        % this for FMRC
        tunits = strrep(tunits,'T',' ');
        basedate = datenum(tunits,'yyyy-mm-dd HH:MM:SS');
      catch
        % can't determine the basedate
        basedate = NaN;
      end
    end
  end
  if isnan(basedate)
    dnum = fac*ocean_time;
  else
    dnum = basedate + fac*ocean_time;
  end
end

if nargout > 1
  if dateformat < 0 || isnan(basedate)
    % not a date format
    for i=1:length(dnum)
      dstr(i,:) = ['day ' num2str(dnum(i),'%8.2f')];
    end
  else
    for i=1:length(dnum)
      dstr(i,:) =  datestr(dnum(i),dateformat);
    end
  end
end   

if nargout == 0
  disp(datestr(dnum))
end
if nargout > 0
  Dnum = dnum;
end
if nargout > 1
  Dstr = dstr;
end

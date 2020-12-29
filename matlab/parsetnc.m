function [gregorian_base, rescale_serial_rel, serial_base_jd, serial_base] = parsetnc(time_units_str)
% $Id: parsetnc.m 423 2014-01-13 17:13:38Z wilkin $
% PARSETNC: parses the COARDS string that specifies time units
%
%function [gregorian_base, rescale_serial_rel, serial_base_jd, ...
%          serial_base] = parsetnc(time_units_str)
%
%        INPUT:
% time_units_str: The string stored as the units attribute of the time
%    variable in a COARDS standard netcdf file. It is like:
%    'seconds since 1992-10-8 15:15:42.5 -6:00'
%
%        OUTPUT:
% gregorian_base: a 6-vector giving the year, month, day, hour, minute,
%    second of the base time as specified in time_units_str. This is in UT.
% rescale_serial_rel: number used to convert the time vector in the netcdf
%    file to days. For example if time_units_str is
%    'seconds since 1992-10-8 15:15:42.5 -6:00' then rescale_serial_rel is
%    1/(24*60*60).
% serial_base_jd: the Julian day number of the base time, in UT, as
%    determined by get_julian_day. Thus gregorian_base =
%    get_calendar_date(serial_base_jd).
% serial_base: the serial time of the base time, in UT, as determined by
%    matlab's datenum function. Thus gregorian_base = datevec(serial_base).
%    serial_base will be problematic for times before October 15 1582, when
%    the Gregorian calendar was adopted, since datenum is proleptic gregorian
%    and so not what you might expect.
%
%        Notes:
% In a COARDS standard netcdf file time is specified relative to a standard
% time by the time variable having a units attribute of the form:
%      'seconds since 1992-10-8 15:15:42.5 -6:00'
% This indicates seconds since October 8th, 1992 at 3 hours, 15
% minutes and 42.5 seconds in the afternoon in the time zone
% which is six hours to the west of Coordinated Universal Time
% (i.e. Mountain Daylight Time). The time zone specification can
% also be written without a colon using one or two-digits
% (indicating hours) or three or four digits (indicating hours
% and minutes).  Instead of 'seconds' the string may contain 'minutes',
% 'hours', 'days' and 'weeks' and all of these may be singular or plural
% and have capital first letters.  I also allow the letters 'UTC' or
% 'UT' at the end of the string, but these are ignored.
%
% All calculations are done using the functions get_calendar_date and
% get_julian_day which know about both the Julian and Gregorian calendars and
% so work back to julian day 0 in the year -4712.

% This function calls: get_calendar_date.m, get_julian_day.m
% This function is called by: get_serial_time.m, timenc.m

% Copyright J. V. Mansbridge, CSIRO, Thursday September 16 11:12:35 EST 1999

% Note 1) I had trouble using strtok with ':' and '-' and so I converted
%         the first 2 instances of each of them to spaces.  This leaves
%         the correction to Universal time (if there is one) untouched.

fm = find(time_units_str == '-');
fc = find(time_units_str == ':');
str = time_units_str;

if length(fm) >= 2
  str(fm(1:2)) = ' ';
else
  error('The year-month-day string is wrongly specified in the base date')
end

if length(fc) == 1
  str(fc(1)) = ' ';
elseif length(fc) >= 2
  str(fc(1:2)) = ' ';
end
str = char(str);

%Convert the time vector to serial days since the base time.

[str, remainder] = strtok(str);
if ( strcmp('week', str) | strcmp('weeks', str) | ...
      strcmp('Week', str) | strcmp('week', str) )
  rescale_serial_rel = 7;
elseif ( strcmp('day', str) | strcmp('days', str) | ...
      strcmp('Day', str) | strcmp('day', str) )
  rescale_serial_rel = 1;
elseif ( strcmp('hour', str) | strcmp('hours', str) | ...
      strcmp('Hour', str) | strcmp('Hours', str) )
  rescale_serial_rel = 1/24;
elseif ( strcmp('minute', str) | strcmp('minutes', str) | ...
      strcmp('Minute', str) | strcmp('Minutes', str) )
  rescale_serial_rel = 1/(24*60);
elseif ( strcmp('second', str) | strcmp('seconds', str) | ...
      strcmp('Second', str) | strcmp('Seconds', str) )
  rescale_serial_rel = 1/(24*60*60);
else
  error(['bad time baseline string = ' time_units_str ])
end

%Find the serial base time (initially without paying attention to any
%reference to Universal time).

[str, remainder] = strtok(remainder);
[str, remainder] = strtok(remainder);
if isempty(str)
  error('The base time has no string for the year')
end
year_base = str2num(str);
[str, remainder] = strtok(remainder);
if isempty(str)
  error('The base time has no string for the month')
end
month_base = str2num(str);
[str, remainder] = strtok(remainder);
if isempty(str)
  error('The base time has no string for the day')
end
day_base = str2num(str);
[str, remainder] = strtok(remainder);
if isempty(str)
  % disp('The base time has no string for the hour, assume hour = min = sec = 0')
  hour_base = 0;
  minute_base = 0;
  second_base = 0;
else
  hour_base = str2num(str);
  [str, remainder] = strtok(remainder);
  if isempty(str)
    % disp('The base time has no string for the minute, assume min = sec = 0')
    minute_base = 0;
    second_base = 0;
  else
    minute_base = str2num(str);
    [str, remainder] = strtok(remainder);
    if isempty(str)
      % disp('The base time has no string for the second, assume sec = 0')
      second_base = 0;
    else
      second_base = str2num(str);
    end
  end
end
gregorian_base = [year_base month_base day_base hour_base ...
        minute_base second_base];
serial_base = datenum(year_base, month_base, day_base, hour_base, ...
    minute_base, second_base);
serial_base_jd = get_julian_day(gregorian_base);

% Strip off some strings whose meaning we understand and also any trailing
% blanks from remainder.

string_list = {'UTC', 'utc', 'UT', 'ut'};
for ii = 1:length(string_list)
  xx = findstr(remainder, string_list{ii});
  if ~isempty(xx)
    remainder = remainder(1:(xx-1));
    break
  end
end
remainder = deblank(remainder);

% If the remainder of the string is not empty (or filled with blanks)
% then we assume that there is information about the conversion to
% Universal time.  This is parsed and serial_base and gregorian_base are
% then modified appropriately.

if ~isempty(remainder)

  % Find the number of hours and minutes that the time is offset from
  % Coordinated Universal Time.
  
  fc = find(remainder == ':');
  if length(fc) == 0
    intxx = str2num(remainder);
    if ( (-99 < intxx) & ( intxx < 99) )
      hour_extra = intxx;
      min_extra = 0;
    elseif ( (-9999 < intxx) & ( intxx < 9999) )
      hour_extra = fix(0.01*intxx);
      min_extra = intxx - 100*hour_extra;
    else
      error(['1:Universal time offset is faulty in ' time_units_str])
    end
  elseif length(fc) == 1
    [str, remainder] = strtok(remainder, ':');
    hour_extra = str2num(str);
    [str, remainder] = strtok(remainder, ':');
    min_extra = sign(hour_extra)*str2num(str);
  else
    error(['2:Universal time offset is faulty in ' time_units_str])
  end
 
  % Error checks
  
  if ((hour_extra < -12) | (hour_extra > 12))
    error(['3:Universal time offset is faulty in ' time_units_str])
  end
  if ((min_extra < -59) | (min_extra > 59))
    error(['4:Universal time offset is faulty in ' time_units_str])
  end

  % Convert the Universal time correction to days.
  
  time_extra = (hour_extra + min_extra/60)/24;

  % Correct serial_base from the local time, as specified in the early
  % part of the string, to Universal time.  Thus in the example
  % 'seconds since 1992-10-8 15:15:42.5 -6:00' we will have
  % time_extra = -6/24 days and this value must be subtracted from
  % serial_base.

  serial_base = serial_base - time_extra;
  serial_base_jd = serial_base_jd - time_extra;
  gregorian_base = get_calendar_date(serial_base_jd);

end


function [jd] = get_julian_day(y, m, d, h)
% GET_JULIAN_DAY   Converts calendar dates to corresponding Julian day numbers.
%
%     Usage: [jd] = get_julian_day(y, m, d, h)
%        or
%            [jd] = get_julian_day([y m d hour min sec])
%     ************************************************************
%
%        jd... decimal Julian day number
%        d.... day (1-31) component of Gregorian date
%        m.... month (1-12) component
%        y.... year (e.g., 1979) component
%        h.... decimal hours (can be a fraction, assumed 0 if absent)
%
%     ************************************************************
%
%    NOTES
%    1) Formally, Julian days start and end at noon. In this convention,
% Julian day 2440000 begins at 1200 hours, May 23, 1968.
%    2) gtime may be a matrix containing columns of years, months,
% etc. Alternatively, y, m, d, h may be vectors and must all be of the same
% type, i.e. a row or column vector. In any case jd will be a column vector
% of julian days.
%    3) The algorithm is taken from Astronomical Algorithms by Jean Meeus and
% gives the correct dates on the Gregorian and Julian calenders. Thus it is
% accurate back to the beginning of the year -4712.
%    4) Because the Christian calendar does not have a year zero then what
% historians call 10 BC is actually the year -9.
%    5) The standard matlab functions datenum, datevec and datestr assume
% that we are using a Gregorian calendar and so cannot be used before 15
% October 1582. The same restriction applies to the locally written matlab
% routines julian and gregorian.

% Copyright J. V. Mansbridge, CSIRO, Tuesday June 27 15:15:51 EST 2000

if nargin == 3,
  h=0.;
elseif nargin == 1,
  h = y(:,4) + (y(:,5) + y(:,6)/60)/60; % hour,min,sec as fraction of hour
  d=y(:, 3);
  m=y(:, 2);
  y=y(:, 1);
elseif nargin == 4
  
else
  error('get_julian_day must have either 1, 3 or 4 input arguments')
end

d = d + h/24; % Make d non-integer as used by Meeus

greg_lim = 15 + 31*(10 + 12*1582);
ff = find(m <= 2);
if ~isempty(ff)
  m(ff) = m(ff) + 12;
  y(ff) = y(ff) - 1;
end

ff_greg = find((d + 31*(m + 12*y)) >= greg_lim);
b = zeros(size(y));
if ~isempty(ff_greg)
  a = floor(y(ff_greg)/100);
  b(ff_greg) = 2 - a + floor(a/4);
end
jd = floor(365.25*(y + 4716)) + floor(30.6001*(m + 1)) + d + b - 1524.5;
jd = jd(:);




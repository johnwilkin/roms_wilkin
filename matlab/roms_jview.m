function [Data,han] = roms_jview(file,var,time,jindex,grd,cplot,clev)
% [data,han] = roms_jview(file,var,time,jindex,grd,[cplot],[clev])
%
% file   = roms his/avg/rst etc nc file
%
% var    = variable to plot
%
% time  = time index into nc FILE
%      or date/time string (in DATESTR format 0) 'dd-mmm-yyyy HH:MM:SS'
%      in which case the function finds the closest time index to that time
%      If time is Inf, 'latest' or 'last' the last record is plotted.
%
% jindex = jindex for slice
%        if jindex < 0  the x-axis coordinate will be lat instead of lon
%
% grd can be 
%       grd structure (from roms_get_grid)
%       grd_file name
%       [] (will attempt to get grid from roms file)
%
% [cplot] optional
%       if true the result will be a contour plot, not pcolor
%       [clev] optional specifies the contours (see N or V in help contour)
%
% Hidden options
%       If there is a global variable logdata and it is true, the data
%       will be log10 transformed before plotting
%       If there is a global variable log_chl and it is true, a specific
%       log scaling suited to typical chlorophyll values is applied
%        
% John Wilkin jwilkin@rutgers.edu
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_jview.m 592 2020-12-28 21:56:28Z wilkin $
%
% See also ROMS_ZVIEW, ROMS_SVIEW, ROMS_IVIEW

if nargin == 0
  help roms_iview
  return
end

if isinf(time)
  time = 'latest';
end
if ischar(time)
  time = roms_get_time_index(file,time);
end

if iscell(var) % case function is called in a cell list loop of varnames
  var = char(var);
end

if nargin < 5
  grd = [];
end

xcoord = 'lon';
if jindex<0
  xcoord = 'lat';
  jindex = abs(jindex);
end

[data,z,lon,lat,t] = roms_jslice(file,var,time,jindex,grd);
data = double(data);

global logdata log_chl
if logdata
  data = log10(data);
else
  if log_chl
    data = max(0.01,data);
    data = (log10(data)+1.4)/0.012;
    ct = [0.01 .03 .1 .3 1 3 10 30 66.8834];
    logct = (log10(ct)+1.4)/0.012;
    cax = range(logct);
  end
end

% option to make a contour plot instead
if nargin > 5
  if nargin == 6
    clev = 30;
  end
else
  cplot = false;
end
if cplot
  % contour plot of the variable
  switch xcoord
    case 'lon'
      hant = contour(lon,z,data,clev);
      labstr = [' - MeanLat ' num2str(mean(lat(:)),4)];
    case 'lat'
      hant = contour(lat,z,data,clev);
      labstr = [' - MeanLon ' num2str(mean(lon(:)),4)];
  end
else
  % pcolor plot of the variable
  switch xcoord
    case 'lon'
      hant = pcolorjw(lon,z,data);
      labstr = [' - MeanLat ' num2str(mean(lat(:)),4)];
    case 'lat'
      hant = pcolorjw(lat,z,data);
      labstr = [' - MeanLon ' num2str(mean(lon(:)),4)];
  end
end

if log_chl
  caxis(cax);
  hancb = colorbar;
  set(hancb,'ytick',logct(1:end-1),'yticklabel',ct)
  set(get(hancb,'xlabel'),'string','mg m^{-3}')
end

% get the time/date for plot label
[t,tindex] = roms_get_time(file,time);
if isdatetime(t)
  tdate = ['on day ' datestr(t,0)];
else
  tdate = ['on day ' num2str(t,'%8.2f')];
end

titlestr = ...
    {['file: ' strrep_(file) ],...
    [(strrep_(var)) ' ' tdate labstr]};
title(titlestr)

if nargout > 0
  Data.var = data;
  Data.lon = lon;
  Data.lat = lat;
  Data.z = z;
  Data.t = t;
  Data.tindex = tindex;
end

if nargout > 1
  han = hant;
end

function [thedata,thegrid,han] = roms_sview(file,var,time,k,grd,vec_d,uscale,varargin)
% [theData,theGrid,theHan] = roms_sview(file,var,time,k,grd,vec_d,uscale,varargin)
% Plot a constant-s slice out of a ROMS history, averages or restart file
%
% Inputs:
%
% file  = roms his/avg/rst/dia etc netcdf file
%
% var = (1) name of the ROMS output variable to plot
%       (2) if isstruct(var) then
%            var.name is the variable name
%            var.cax  is the color axis range
%     See notes below about special variables not in the file
%
% time  = time index into nc FILE
%      or date/time string (in DATESTR format 0) 'dd-mmm-yyyy HH:MM:SS'
%      in which case the function finds the closest time index to that time
%      If time is Inf, 'latest' or 'last' the last record is plotted.
%
% k     = index of vertical (s-coordinate) level of horizontal slice 
%       if k==0 any vector plot will be for ubar,vbar
%       if k==Inf plot the surfacemost layer
%
% grd can be 
%       grd structure (from roms_get_grid)
%       grd_file name
%
% vec_d = density (decimation factor) of velocity vectors to plot over 
%       if 0 no vectors are plotted
%
% uscale = vector length scale
%       if uscale < 0 then pseudo particle tracks are plotted instead of 
%       quiver velocity vectors, and abs(uscale) is intepretted as the 
%       duration in days of the track length. See roms_curquivergrd.
%       This can be very slow on a large grid. If you have zoomed in the
%       view it will be faster to add curved vectors separately with
%       function roms_addvect.
%
% varargin are additional arguments passed on to roms_quivergrd or 
%       roms_curquivergrd to format the plot
%
% If wet/dry masks are present in the output file, then they will be 
% applied to the plot if this preference is set:
%      setpref('ROMS_WILKIN','USE_WETDRY_MASK',true);
%
% Requesting plots of special variables not actually in the file
%       varname = ...
%            'Chlorophyll' with a captial C plots chlorophyll log
%              transformed before pcolor
%            'ubarmag' or 'vbarmag' plots velocity magnitude computed
%              from vector (ubar,vbar)
%            'umag' or 'vmag' plots velocity magnitude computed
%              from vector (u,v)
%            'stress' or 'bstress' plot the magnitude of the surface or
%               bottom stress, respectively
%            'rvor' plot relative velocity
%            'ow' plot Okubo-Weiss parameter
%            'wind' plot the magnitude of vector (Uwind,Vwind) as from a
%               forcing file; both components must be in the same file 
%               and on the ROMS rho-points grid (no regrid option) 
%            'omegaca','omegaar','ph','phtotal' use CO2SYS to compute
%            constituents of the ocean carbon state - alkalinity, TIC must
%            be present in the ROMS file from the Fennel/BGC model
%
% Outputs:
% 
% thedata = structure of pcolored data and velocities
% thegrid = roms grid structure
% han = structure of handles for pcolor, quiver and title objects
%
% John Wilkin jwilkin@rutgers.edu
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_sview.m 592 2020-12-28 21:56:28Z wilkin $
%
% See also ROMS_ZVIEW, ROMS_IVIEW, ROMS_JVIEW

if nargin == 0
  error('Usage: roms_sview(file,var,time,k,grd,vec_d,uscale,varargin)');
end

if isinf(time)
  time = 'latest';
end
if ischar(time) || isdatetime(time)
  time = roms_get_time_index(file,time);
end

% % check if input TIME is in datestr format, and if so find the
% % time index in FILE that is the closest
% if isinf(time)
%   time = 'latest';
% end
% if ischar(time)
%   fdnums = roms_get_date(file,-1);
%   if strncmp(time,'latest',2)
%     time = length(fdnums);
%   else
%     dnum = datenum(time);
%     if dnum >= fdnums(1) && dnum <= fdnums(end)
%       [~,time] = min(abs(dnum-fdnums));
%       time = time(1);
%     else
%       warning(' ')
%       disp(['Requested date ' time ' is not between the dates in '])
%       disp([file ' which are ' datestr(fdnums(1),0) ' to ' ])
%       disp(datestr(fdnums(end),0))
%       thedata = -1;
%       return
%     end
%   end
% end

if isinf(k)
  k = grd.N;
end

% sneaky little trick to allow me to send a caxis range through 
% the input - should really be done with attribute/value pairs
if isstruct(var)
  cax = var.cax;  
  var = var.name;
else 
  cax = 'auto';
end
if iscell(var) % case function is called in a cell list loop of varnames
  var = char(var);
end

if strcmp(var,'temp')
  varlabel = 'temperature';
else
  varlabel = strrep_(var);
end

% sneaky trick to force log transformation of chlorophyll
% ... give input varname with a capital C
log_chl = 0;
if strcmp(var,'Chlorophyll')
  log_chl = 1;
  var = 'chlorophyll';
end

%% get the data

% figure out whether a 2-D or 3-D variable by testing the dimensions of the
% variable
switch var
  case {'ubarmag','vbarmag'}
    vartest = 'ubar';
  case 'stress'
    vartest = 'sustr';
  case 'bstress'
    vartest = 'bustr';
  case 'wind'
    vartest = 'Uwind';
  case {'umag','rvor','ow'}
    vartest = 'u';
  case {'omegaca','omegaar','ph','phtotal'}
    vartest = 'salt';
  otherwise
    vartest = var;
end
Vi = ncinfo(file,vartest);
dimNames = {Vi.Dimensions.Name};
if any(strcmp(dimNames,'s_rho')) || any(strcmp(dimNames,'s_w')) ...
    || any(strcmp(dimNames,'sc_r'))
  START = [time-1 1 k-1  0  0]; % 2nd dim is for perfect restart files
  COUNT = [1      1 1   -1 -1];
  depstr = [ ' at level ' int2str(k) ' '];
else
  START = [time-1 1  0  0];
  COUNT = [1      1 -1 -1];
  depstr = ' ';
end
% perfect restart files have an added dimension named 'two'
if ~any(strcmp(dimNames,'two')) && ~any(strcmp(dimNames,'three'))
  % not a perfect restart so take these out
  START(2) = [];
  COUNT(2) = [];
end

long_name = ' ';
units = ' ';
switch var
  % derived variables
  case { 'ubarmag','vbarmag'}
    datau = squeeze(nc_varget(file,'ubar',START,COUNT));
    datau = datau(:,[1 1:end end]);
    datau = av2(datau')';
    datav = squeeze(nc_varget(file,'vbar',START,COUNT));
    datav = datav([1 1:end end],:);
    datav = av2(datav);
    data = abs(datau+sqrt(-1)*datav);
    depstr =  ' depth average ';
    % var = 'ubar'; % for time handling
  case 'stress'
    warning('option not debugged yet')
    datau = squeeze(nc_varget(file,'sustr',START,COUNT));
    datau = datau(:,[1 1:end end]);
    datau = av2(datau')';
    datav = squeeze(nc_varget(file,'svstr',START,COUNT));
    datav = datav([1 1:end end],:);
    datav = av2(datav);
    data = abs(datau+sqrt(-1)*datav);
    depstr =  ' at surface ';
    % var = 'sustr'; % for time handling
  case 'bstress'
    warning('option not debugged yet')
    datau = squeeze(nc_varget(file,'bustr',START,COUNT));
    datau = datau(:,[1 1:end end]);
    datau = av2(datau')';
    datav = squeeze(nc_varget(file,'bvstr',START,COUNT));
    datav = datav([1 1:end end],:);
    datav = av2(datav);
    data = abs(datau+sqrt(-1)*datav);
    depstr =  ' at bottom ';
    % var = 'bustr'; % for time handling
  case 'wind'
    datau = squeeze(nc_varget(file,'Uwind',START,COUNT));
    datav = squeeze(nc_varget(file,'Vwind',START,COUNT));
    data = abs(datau+sqrt(-1)*datav);
    depstr =  ' 10 m above surface ';
    % var = 'Uwind'; % for time handling
  case 'umag'
    datau = squeeze(nc_varget(file,'u',START,COUNT));   
    datau(isnan(datau)==1) = 0; 
    datau = datau(:,[1 1:end end]);
    datau = av2(datau')';   
    datav = squeeze(nc_varget(file,'v',START,COUNT));
    datav(isnan(datav)==1) = 0;
    datav = datav([1 1:end end],:);
    datav = av2(datav);
    data = abs(datau+sqrt(-1)*datav);
    depstr = [ ' at level ' int2str(k) ' '];
    % var = 'temp'; % for time handling
  case 'rvor'
    datau = squeeze(nc_varget(file,'u',START,COUNT)); 
    % datau(isnan(datau)) = 0; 
    datav = squeeze(nc_varget(file,'v',START,COUNT));
    % datav(isnan(datav)) = 0;
    data = roms_vorticity(datau,datav,grd,'relative');
  case 'ow'
    datau = squeeze(nc_varget(file,'u',START,COUNT)); 
    % datau(isnan(datau)) = 0; 
    datav = squeeze(nc_varget(file,'v',START,COUNT));
    % datav(isnan(datav)) = 0;
    data = roms_vorticity(datau,datav,grd,'okubo-weiss');
  case {'omegaca','omegaar','ph','phtotal'} 
    temp = nc_varget(file,'temp',START,COUNT);
    salt = nc_varget(file,'salt',START,COUNT);
    alk = nc_varget(file,'alkalinity',START,COUNT);
    tic = nc_varget(file,'TIC',START,COUNT);
    press = -grd.z_r(k,:,:);
    data = roms_co2sys_var(var,temp,salt,alk,tic,press); 
    data = squeeze(data);
  otherwise
    data = squeeze(nc_varget(file,var,START,COUNT));
    try
      units = ncreadatt(file,var,'units');
      long_name = ncreadatt(file,var,'long_name');
    catch
    end
end

%% get the appropriate land/sea or wet/dry mask
usewetdry = false;
haswetdry = nc_isvar(file,'wetdry_mask_rho');
if haswetdry % wet dry mask exists in file
  try % override with preference
    usewetdry = getpref('ROMS_WILKIN','USE_WETDRY_MASK');
  catch % no preference - assume not
  end
else
  % averages don't have wetdry mask
  avg_wetdry = false;
  try 
    avg_wetdry = getpref('ROMS_WILKIN','AVG_WETDRY_MASK');
  catch
  end
  if avg_wetdry
    data(data==0) = NaN;
  end
end
pos = roms_cgridpos(data,grd);
ma = ['mask_' pos];
lo = ['lon_' pos];
la = ['lat_' pos];
if usewetdry
  mask = squeeze(nc_varget(file,['wetdry_' ma],[time-1 0 0],[1 -1 -1]));
else
  mask = grd.(ma);
end
x = grd.(lo);
y = grd.(la);
mask(mask==0) = NaN;

%%
if log_chl
  data = max(0.01,data,'includenan');
  data = (log10(data)+1.4)/0.012;
  ct = [0.01 .03 .1 .3 1 3 10 30 66.8834];
  logct = (log10(ct)+1.4)/0.012;
  cax = range(logct);
end

%% special handling for some grids to blank out regions
if isfield(grd,'special')
  if iscell(grd.special)
    % potentially several special options
    vlist = grd.special;
  else
    % single option but copy to cell for handling below
    vlist{1} = grd.special;
  end
  for k=1:length(vlist)
    opt = char(vlist{k});
    switch char(opt)
      case 'jormask'
        %           for opt = vlist
        %     opt = char
        %     switch char(opt)
        %       case 'jormask'
        % apply Jay O'Reilly's mask to trim the plotted nena data
        xpoly = [-82 -79.9422 -55.3695 -55.3695 -82];
        ypoly = [24.6475 24.6475 44.0970 46 46];
        ind = inside(x,y,xpoly,ypoly);
        mask(ind==0) = NaN;
      case 'nestedge'
        xpoly = vlist{k+1}(:,1);
        ypoly = vlist{k+1}(:,2);
        ind = inside(x,y,xpoly,ypoly);
        mask(ind==1) = NaN;
        break
      case 'logdata'
        % this would be a better place to log transform data before
        % plotting
        data = max(0.01,data);
        data = log10(data);
      case 'say'
        unix(['say ' var]);
    end
  end
end

%% make the plot
hanpc = pcolorjw(x,y,data.*mask);
caxis(cax)
if isfield(grd,'nolatlon')
  if ~grd.nolatlon
    xtickformat('degrees');
    ytickformat('degrees');
  end
else
  xtickformat('degrees');
  ytickformat('degrees');
end
hancb = colorbar;

%%
if log_chl
  set(hancb,'ytick',logct(1:end-1),'yticklabel',ct)
  set(get(hancb,'xlabel'),'string','mg m^{-3}')
end

%% add vectors
if nargin > 5
  if vec_d
    
    % nc = netcdf(file);
    % add vectors
    % ! sorry, this doesn't allow for {u,v}bar vectors on a 3d variable
    if k > 0
      if numel(START)<4 % field was 2D but k>0 wants 3D velocity
        START = START([1 1 2 3]);
        START(2) = k-1;
        COUNT = COUNT([1 1 2 3]);
        COUNT(2) = 1;
      end
      u = nc_varget(file,'u',START,COUNT);
      v = nc_varget(file,'v',START,COUNT);
      depstr = [depstr ' Vectors at level ' int2str(k) ' '];
    else
      u = nc_varget(file,'ubar',START,COUNT);
      v = nc_varget(file,'vbar',START,COUNT);
      % a forcing file won't have u,v ...
      if isempty(u)
        u = nc_varget(file,'sustr',START,COUNT);
        v = nc_varget(file,'svstr',START,COUNT);
        depstr = [depstr ' Wind stress vectors '];
      else
        depstr = [depstr ' Depth average velocity vectors '];
      end
    end
    if nargin < 7
      uscale = 1;
    end
    u = squeeze(u);
    v = squeeze(v);
    u(isnan(u)==1) = 0;
    v(isnan(v)==1) = 0;
    if uscale > 0
      % quiver plot
      [hanquiver,dataq] = roms_quivergrd(u,v,grd,vec_d,uscale,varargin{:});
    else
      % curvy track plot
      lon0 = x(1:vec_d:end,1:vec_d:end);
      lat0 = y(1:vec_d:end,1:vec_d:end);
      dmask = data(1:vec_d:end,1:vec_d:end);
      lon0(isnan(dmask)) = [];
      lat0(isnan(dmask)) = [];
      [hanquiver,curdata] = ...
        roms_curquivergrd(u,v,grd,lon0(:),lat0(:),-uscale,10,varargin{:});
    end
  end
end

% change plotaspectratio to be approximately Mercator
% if you don't like this, add variable merc = false to the grd structure
% to disable this
if isfield(grd,'merc')
  if grd.merc 
    set(gca,'DataAspectRatio',[1 cos(mean(ylim)*pi/180) 1]);
  end
else
  % default is to assume lon/lat axes and scale them to be mercator-like
  set(gca,'DataAspectRatio',[1 cos(mean(ylim)*pi/180) 1]);
end

%% my trick to plot a coast if it knows how to do this from the grd_file
% name
try
  if strfind('leeuwin',grd.grd_file)
    gebco_eez(0,'k')
  elseif strfind('eauc',grd.grd_file)
    plotnzb
  elseif strfind('nena',grd.grd_file)
    plotnenacoast(3,'k') 
  elseif strfind('sw06',grd.grd_file)
    plotnenacoast(3,'k')
  end
catch
end

% get the time/date for plot label
t = roms_get_time(file,time);
if isdatetime(t)
  tdate = ['on day ' datestr(t,0)];
else
  tdate = ['on day ' num2str(t,'%8.2f')];
end

% try
%   [t,tdate] = roms_get_date(file,time,0);
% catch
%   t = NaN;
%   tdate = 'unknown';
% end

% label
titlestr{1} = ['file: ' strrep_(file)];
titlestr{2} = [varlabel ' ' tdate ' ' depstr];
hantitle = title(titlestr);

% pass data to outputs
if nargout > 0
  thedata.x = x;
  thedata.y = y;
  thedata.data = data;
  thedata.t = t;
  thedata.tstr = tdate;
  thedata.tindex = time;
  thedata.varname = var;
  thedata.units = units;
  thedata.long_name = long_name;
  if nargin > 5
    if vec_d
      thedata.u = u;
      thedata.v = v;
      if uscale < 0
        thedata.xcurv = curdata.lon;
        thedata.ycurv = curdata.lat;
      end
    end
  end
end
if nargout > 1
  thegrid = grd;
end
if nargout > 2
  han.title = hantitle;
  han.pcolor = hanpc;
  han.colorbar = hancb;
  if exist('hanquiver','var')
    han.quiver = hanquiver;
  end
end

function str = caps(str)
str = lower(str);
str(1) = upper(str(1));

function s = strrep_(s)
s = strrep(s,'\','\\');
s = strrep(s,'_','\_');
s = strrep(s,'^','\^');

function a = av2(a)
%AV2	grid average function.  
%       If A is a vector [a(1) a(2) ... a(n)], then AV2(A) returns a 
%	vector of averaged values:
%	[ ... 0.5(a(i+1)+a(i)) ... ]  
%
%       If A is a matrix, the averages are calculated down each column:
%	AV2(A) = 0.5*(A(2:m,:) + A(1:m-1,:))
%
%	TMPX = AV2(A)   will be the averaged A in the column direction
%	TMPY = AV2(A')' will be the averaged A in the row direction
%
%	John Wilkin 21/12/93
[m,n] = size(a);
if m == 1
	a = 0.5 * (a(2:n) + a(1:n-1));
else
	a = 0.5 * (a(2:m,:) + a(1:m-1,:));
end

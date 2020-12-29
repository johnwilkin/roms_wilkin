function [thedata,thegrid,han] = roms_zview(file,var,time,depth,grd,vec_d,uscale,varargin)
% [Data,Grid,Han] = roms_zview(file,var,time,depth,grd,vec_d,uscale,varargin)
% Plot a constant-z slice out of a ROMS history, averages or restart file
%
% Inputs:
%
% file  = roms his/avg/rst/dia etc. netcdf file or OPeNDAP aggregation
%
% var = (1) name of the ROMS output variable to plot
%       (2) if isstruct(var) then
%            var.name is the variable name
%            var.cax  is the color axis range
%       (3) special variables not actually in the file
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
%
% time  = time index into nc FILE
%      or date/time string (in DATESTR format 0) 'dd-mmm-yyyy HH:MM:SS'
%      in which case the function finds the closest time index to that time
%      If time is Inf, 'latest' or 'last' the last record is plotted.
%
% depth = z depth of horizontal slice (m)
%       if depth==0 any vector plot will be for ubar,vbar
%
% grd can be
%       grd structure, output from roms_get_grid << RECOMMENDED
%       grd_file name
%       [] will attempt to get grid from the first input (file)
%
%     RECOMMENDED best practise is to load grd structure from the actual
%     output file you are using. This guarantees you get the correct
%     coordinates: grd = roms_get_grid(file,file)
%       See the help on roms_get_grid to understand about using zeta for 
%     to enable calculating time-varying vertical coordinates if this is 
%     important to you in plotting
%
% vec_d = density (decimation factor) of velocity vectors to plot over
%       if 0 no vectors are plotted
%
% uscale = vector length scale
%       if uscale < 0 then pseudo particle tracks are plotted instead of 
%       quiver velocity vectors, and abs(uscale) is intepretted as the 
%       duration in days of the track length. See roms_curquivergrd.m
%
% varargin are additional arguments passed on to roms_quivergrd or 
%       roms_curquivergrd to format the plot
%
% If this preference is set to true, and wet/dry masks are present in the
% output, then they will be applied to the plot:
%      setpref('ROMS_WILKIN','USE_WETDRY_MASK',true);
%
% Outputs:
%
% Data = structure of pcolored data and velocities
% Grid = roms grid structure
% Han = structure of handles for pcolor, quiver and title objects
%
% John Wilkin jwilkin@rutgers.edu
%
% Copyright (c) 2020 John L. Wilkin 
% $Id: roms_zview.m 592 2020-12-28 21:56:28Z wilkin $
%
% See also ROMS_SVIEW, ROMS_IVIEW, ROMS_JVIEW

if nargin == 0
  error('Usage: roms_zview(file,var,time,depth,grd,vec_d,uscale,varargin)');
end

% check if input TIME is in datestr format, and if so find the
% time index in FILE that is the closest
if isinf(time)
  time = 'latest';
end
if ischar(time)
  fdnums = roms_get_date(file,-1);
  if strncmp(time,'latest',2)
    time = length(fdnums);
  else
    dnum = datenum(time);
    if dnum >= fdnums(1) && dnum <= fdnums(end)
      [~,time] = min(abs(dnum-fdnums));
      time = time(1);
    else
      warning(' ')
      disp(['Requested date ' time ' is not between the dates in '])
      disp([file ' which are ' datestr(fdnums(1),0) ' to ' ])
      disp(datestr(fdnums(end),0))
      thedata = -1;
      return
    end
  end
end

if nargin < 5
  grd = [];
end

% sneaky trick to send a preset caxis range through the input - should 
% really be done with attribute/value pairs.
if isstruct(var)
  cax = var.cax;
  var = var.name;
else
  cax = 'auto';
end
if iscell(var) % case function is called in a cell list loop of varnames
  var = char(var);
end

varlabel = strrep_(var);
if strcmp(varlabel,'temp')
  varlabel = 'temperature';
end

% sneaky trick to force log transformation of chlorophyll
% ... give input varname with a capital C
log_chl = 0;
if strcmp(var,'Chlorophyll')
  log_chl = 1;
  var = 'chlorophyll';
end

%% pcolor plot of the variable
switch var
  case { 'ubar','vbar','zeta','Hsbl','h','f','pm','pn',...
      'swrad','SST','dqdsst','shflux','swflux','SSS',...
      'sustr','svstr','Uwind','Vwind','Tair','Pair',...
      'sensible','latent'}
    data = roms_2dslice(file,var,time,grd);
    depstr = [];
  case 'umag'
    datau = roms_zslice(file,'u',time,depth,grd);
    datav = roms_zslice(file,'v',time,depth,grd);
    % average to rho points
    datau(isnan(datau)==1) = 0;
    datav(isnan(datav)==1) = 0;    
    datau = datau(:,[1 1:end end]);
    datau = av2(datau')';
    datav = datav([1 1:end end],:);
    datav = av2(datav);
    data = abs(datau+sqrt(-1)*datav);
    var = 'temp'; % for mask and time handling
    x = grd.lon_rho;
    y = grd.lat_rho;
    datau = roms_zslice(file,var,time,depth,grd);
    depstr = [' - Vectors at depth ' num2str(abs(depth)) ' m '];
  case 'rvor'
    datau = roms_zslice(file,'u',time,depth,grd);
    datav = roms_zslice(file,'v',time,depth,grd);
    data = roms_vorticity(datau,datav,grd,'relative');
    depstr = [' - Depth ' num2str(abs(depth)) ' m '];
  case 'ow'
    datau = roms_zslice(file,'u',time,depth,grd);
    datav = roms_zslice(file,'v',time,depth,grd);
    data = roms_vorticity(datau,datav,grd,'okubo-weiss');
    depstr = [' - Depth ' num2str(abs(depth)) ' m '];
  case {'omegaca','omegaar','ph','pphotal'}
    temp = roms_zslice(file,'temp',time,depth,grd);
    salt = roms_zslice(file,'salt',time,depth,grd);
    alk = roms_zslice(file,'alkalinity',time,depth,grd);
    tic = roms_zslice(file,'TIC',time,depth,grd);
    press = -depth*ones(size(salt));
    data = roms_co2sys_var(var,temp,salt,alk,tic,press); 
    depstr = [' - Depth ' num2str(abs(depth)) ' m '];
  otherwise
    data = roms_zslice(file,var,time,depth,grd);
    depstr = [' - Depth ' num2str(abs(depth)) ' m '];
end

long_name = ' ';
units = ' ';
try
  units = ncreadatt(file,var,'units');
  long_name = ncreadatt(file,var,'long_name');
catch
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
  if grd.nolatlon ~= 1
    xtickformat('degrees');
    ytickformat('degrees');
  end
else
  xtickformat('degrees');
  ytickformat('degrees');
end
hancb = colorbar;

if log_chl
  set(hancb,'ytick',logct(1:end-1),'yticklabel',ct)
  set(get(hancb,'xlabel'),'string','mg m^{-3}')
end

if nargin > 5
  if vec_d
    % add vectors
    % ! sorry, this doesn't allow for {u,v}bar vectors on a 3d variable
    if depth
      u = roms_zslice(file,'u',time,depth,grd);
      v = roms_zslice(file,'v',time,depth,grd);
      if isempty(depstr)
        depstr = [' - Vectors at depth ' num2str(abs(depth)) ' m '];
      end
    else
      try
        % a forcing file won't have u,v ...
        u = roms_2dslice(file,'ubar',time,grd);
        v = roms_2dslice(file,'vbar',time,grd);
      catch
        % ... failing that look for wind stress
        % (should make this more general to allow for plotting
        % wind at 10m, but that would be on rho-points no u,v-points
        % so some extra checking is required -- maybe try to use
        % roms_addvect or roms_sview instead ... something for later)
        u = roms_2dslice(file,'sustr',time,grd);
        v = roms_2dslice(file,'svstr',time,grd);
        depstr = ' - Wind stress vectors ';
      end
      if isempty(depstr)
        depstr = ' - Vectors for depth-average velocity ';
      end
    end
    if nargin < 7
      uscale = 1;
    end
    u(isnan(u)==1) = 0;
    v(isnan(v)==1) = 0;
    if uscale > 0
      % quiver plot
      [hanquiver,dataq] = roms_quivergrd(u,v,grd,vec_d,uscale,varargin{:});
    else
      % curvy track plot
      lon0 = x(2:vec_d:end,2:vec_d:end);
      lat0 = y(2:vec_d:end,2:vec_d:end);
      dmask = data(2:vec_d:end,2:vec_d:end);
      hmask = grd.h(2:vec_d:end,2:vec_d:end);
      dmask(hmask>1000) = NaN;
      lon0(isnan(dmask)) = [];
      lat0(isnan(dmask)) = [];
      [hanquiver,curdata] = roms_curquivergrd(u,v,grd,lon0(:),lat0(:),...
        -uscale,20,varargin{:});
    end
  end
end

% change plotaspectratio to be approximately Mercator
% if you don't like this, add variable nolatlon (=1) to the grd structure
% to disable this
if isfield(grd,'nolatlon')
  if grd.nolatlon ~= 1
    set(gca,'DataAspectRatio',[1 cos(mean(ylim)*pi/180) 1]);
  end
else
  set(gca,'DataAspectRatio',[1 cos(mean(ylim)*pi/180) 1]);
end

% my trick to plot a coast if it knows how tofrom the grd_file name
try
  if strfind(grd.grd_file,'leeuwin')
    gebco_eez(0,'k')
  elseif strfind(grd.grd_file,'eauc')
    plotnzb
  elseif strfind(grd.grd_file,'nena')
    plotnenacoast(3,'k')
  elseif strfind(grd.grd_file,'sw06')
    plotnenacoast(3,'k')
  end
catch
end

% get the time/date
[t,tdate] = roms_get_date(file,time,0);

% label
titlestr{1} = ['file: ' strrep_(file)];
titlestr{2} = [varlabel ' ' tdate ' ' depstr];
%itlestr{2} = tdate; titlestr{3} = [varlabel ' ' depstr];

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
  if exist('dataq','var')==1
    thedata.ue = dataq.ue;
    thedata.vn = dataq.vn;
    thedata.xq = dataq.x;
    thedata.yq = dataq.y;
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

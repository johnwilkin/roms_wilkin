function [Tvar,Coord] = ...
  roms_genslice(file,varname,lonTrk,latTrk,timTrk,varargin)
% [slice,geo] = roms_genslice(file,varname,lonTrk,latTrk,timTrk,options)
%
% Get a vertical slice from ROMS output at coordinates specified by inputs
% of lon, lat and time, for example a glider track, an isobath, or the
% path that follows an estuary thalweg. Data is returned on the ROMS
% native vertical s-coordinate.
%
% To map to an arbitrary set of 4-D coordinates at user defined
% z-coordinates a subsequent vertical interpolation is required (see
% roms_zgenslice).
%
% Inputs:
%   file      ROMS output netcdf or THREDDS URL (aggregation across time)
%   varname   name of the variable we want along the track
%   lon, lat  coordinates along the track
%   timTrk    time values along the track
%             - default: Matlab datenum or with option 'date' (see below)
%             - string: 'yyyy-mm-dd HH:MM' will be converted to a datenum
%             - index: (integer) into file with option 'index' (see below)
%             - the dimension of timeTrk must match lon and lat, OR ... if
%               this is a scalar (index) or string then that single time is
%               used for all lon/lat
%
% Optional inputs:
%            integer Ntrk - interpolate to this many points to refine
%               resolution of output along trajectory. DEFAULT is return
%               values strictly at the requested lon,lat coordinates
%            string 'linear' for interpolation method [DEFAULT]
%            string 'nearest' for alternate interpolation method
%            string 'natural' for alternate interpolation method
%            string 'date' to interpret timTrk as datenum [DEFAULT]
%            string 'index' interpret timTrk as index into time
%            string 'verbose' to monitor progress
%            string 'zeta' to read zeta from file and compute time-varying
%               z coordinate. Without this option DEFAULT is to assume
%               zeta=0 in computing z (for speed)
%
% Output:
%   slice    data interpolated to the track in space and time
%   geo      structure of data coordinate information
%              lon  2D lon along the track suited to pcolor plotting
%              lat  2D lat along the track suited to pcolor plotting
%              time in datenum convention
%              tindex gives the record in the file for the time
%                   immediately before the interpolated timTrk value
%              z    2D depth along the track in vertical stretched coords
%                   suited to pcolor plotting
%              zw   corresponding depths of the layer interfaces
%                   (if the requested variable is a 'w' points variable
%                   then z and zw are the same)
%              dz   thickness of centered layers (so vertical sum of dz
%                   is (h+zeta)
%              dis  2D distance along the track (in kilometers) suited to
%                   pcolor plotting
%              dislen 2D track segment lengths (km) for alongtrack integral
%              en,ep unit vectors normal and parallel to the track for
%                   use by roms_velslice when computing across-track and
%                   along-track velocity
%
% Example usage of time options:
%
%  Plot time index 20 in the file and interpolate to 200 points along track
%   roms_genslice(file,varname,lonTrk,latTrk,20,'index',200)
%
%  Specify time as a string
%   roms_genslice(file,varname,lonTrk,latTrk,'2001-01-01')
%
%  Plot by interpolating in time to a given datenum
%   roms_genslice(file,varname,lonTrk,latTrk,datenum(1960,12,18))
%
%  Don't use zeta along the track to compute z. Interp to 500 points
%   roms_genslice(file,varname,lonTrk,latTrk,'1960-12-18','zero_zeta',500)
%
% History: Based on gslice originally by John Evans and John Wilkin
%   modified by Weifeng Gordon Zhang, and then adapted by Alexander
%   Crosby 2011-04 and entered into NCTOOLBOX as nc_genslice using
%   netcdf-java and CDM to add interoperability with other ocean models
%
% 2013: netcdf-java abandoned because nc_genslice would throw java heap
%       space errors in handling the time coordinate with grid_interop
%       which defaults to 4-D z for ROMS which is huge;
%       and geoij did not correctly extract the minimal i,j subset for
%       a grid rotated w.r.t. east-north coordinates
% 2015: Updates to allow timTrk be an index into the time dimension rather
%       than strictly a datenum, and to optionally recompute z as a
%       function of time varying zeta.
% 2016: Removed all nctoolbox functions for portability (they weren't
%           working right anyway)
%       Use Matlab native netcdf reads throughout instead of snctools
%       Corrected calculation of z coordinate for u,v slices
%       Added layer thicknesses, h and zeta to coordinates output
%       Corrected wrong logic for normal/tangential vectors components that
%           were used by roms_velslice
%       Companion roms_velslice.m updated
% 2018: Removed call to roms_get_date (which uses snctools) and processed
%           time units internally ... may not be robust for some ROMS
%           instances
%       Added waitbar to show progress, because this function is slow for a
%           long time series of data
% 2019: Handling of zeta in z coords calculation was confusing - clarified
%
% Future ... Matlab profiler shows that a significant amount of time is
% spent in the repeated ncread calls but not sure if anything much can be
% done about that.
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu 
% $Id: roms_genslice.m 597 2020-12-29 16:48:48Z wilkin $

%% Check inputs

if ~isequal(size(lonTrk(:)),size(latTrk(:)))
  error('The dimensions of lonTrk and latTrk do not match')
end
lonTrk = double(lonTrk);
latTrk = double(latTrk);
if ischar(timTrk)
  timTrk = datenum(timTrk);
end
if ~isscalar(timTrk)
  if ~isequal(size(lonTrk(:)),size(timTrk(:)))
    error('length of timTrk does not match dimensions of lonTrk,latTrk')
  end
end

% Default options
Ntrk = [];
verbose = false;
method = 'linear';
timMethod = 'date';
zero_zeta = true;

if nargin > 5
  for k = 1:length(varargin)
    opt = varargin{k};
    if ischar(opt)
      switch opt(1)
        case 'v'
          verbose = true;
        case {'n','l'} % nearest, natural or linear
          method = opt;
        case 'd'
          timMethod = 'date';
        case 'i'
          timMethod = 'index';
        case 'z'
          zero_zeta = false;
        otherwise
          error(['failed to comprehend option ' opt])
      end
      % elseif isstruct(opt)
      %  g = opt;
    else
      Ntrk = opt;
    end
  end
end

% some files (e.g. diagnostics) might not have zeta, in which case force
% the zero_zeta option but issue warning
if ~zero_zeta
  try
    ncinfo(file,'zeta');
  catch
    warning('no ''zeta'' variable in file - settinh zeta = 0 and proceeding')
    zero_zeta = true;
  end
end

%% Process input track coordinates ----------------------------------------

if isempty(Ntrk) && length(lonTrk)==2
  % only have section end points and no Ntrk to create viable slice
  Ntrk = 300;
  disp(['Interpolating to ' int2str(Ntrk) ' values between end points'])
  Coord.was_interpolated = true;
end

if ~isempty(Ntrk) && Ntrk > length(lonTrk)
  % Refine track resolution by adding more points
  dist = cumsum([0; sw_dist(latTrk(:),lonTrk(:),'km')]);
  lonTrk = interp1(dist,lonTrk(:),linspace(0,dist(end),Ntrk));
  latTrk = interp1(dist,latTrk(:),linspace(0,dist(end),Ntrk));
  if length(timTrk)~=1
    timTrk = interp1(dist,timTrk(:),linspace(0,dist(end),Ntrk));
  end
  Coord.was_interpolated = true;
end

% Alongtrack distance coordinate for output
lonTrk = lonTrk(:);
latTrk = latTrk(:);
timTrk = timTrk(:);
dist = cumsum([0; sw_dist(latTrk,lonTrk,'km')]);

% track segment lengths in case we want an alongtrack integral
% >>>>>> review this - might want a half dtmp at first/last points <<<<<<<<
dtmp = diff(dist);
dtmp = dtmp([1 1:end end]);
distlen = 0.5*(dtmp(1:(end-1))+dtmp(2:end));
dsm = distlen*1000; % alongtrack elemental distances in meters

if isscalar(timTrk)
  % slice at single time
  timTrk = timTrk*ones(size(lonTrk));
end

% Get time from file
if verbose
  disp('Reading full vector of all times from ')
  disp([' ' file])
end

%% Time coordinate
% The block below may not be robust for all possible descriptions of the
% time coordinate and its units. But this does appear to work with FMRC and
% recent regular ROMS files

try
  % if this is FMRC the coordinate is named time not ocean_time
  tunits = ncreadatt(file,'time','units');
  timename = 'time'; % probably a FMRC aggregation
catch
  try
    tunits = ncreadatt(file,'ocean_time','units');
    timename = 'ocean_time'; % default ROMS file
  catch
    error('Cannot determine the name of the time coordinate')
  end
end
switch tunits(1:3)
  case 'sec'
    fac = 86400;
  case 'min'
    fac = 3600;
  case 'hou'
    fac = 24;
  case 'day'
    fac = 1;
end

% basedate
tsince = tunits(6+strfind(tunits,'since'):end);
tsince = strrep(tsince,'UTC','');

% convert tim_mod to datenum
tim_mod = ncread(file,timename)/fac + datenum(tsince);

%% Time processing
% vector Ttrk is in fractional time index units regardless of how
% input timTrk was given (round(Ttrk) is used to index the data file)
switch timMethod
  case 'date'
    % interp the time index to the requested track times
    Ttrk = interp1(tim_mod,1:length(tim_mod),timTrk);
  case 'index'
    % time index was given, interp the model time to this index
    Ttrk = timTrk;
    if length(tim_mod)>1
      timTrk = interp1(1:length(tim_mod),tim_mod,Ttrk);
    else
      timTrk = tim_mod*ones(size(Ttrk));
    end
end

%% Begin reading ROMS data information ------------------------------------

if verbose
  disp('Reading full domain lon/lat grid coordinates from ')
  disp(['  ' file])
end

% determine what staggered grid location the requested variable is on
Vi = ncinfo(file,varname);
dimNames = {Vi.Dimensions.Name};
if any(strcmp(dimNames,'eta_rho'))
  igrid = 1;
  pos = '_rho';
elseif any(strcmp(dimNames,'eta_u'))
  igrid = 3;
  pos = '_u';
elseif any(strcmp(dimNames,'eta_v'))
  igrid = 4;
  pos = '_v';
else
  error(['could not determine c-grid position of ' varname])
end
w_points = false;
if any(strcmp(dimNames,'s_w'))
  w_points = true;
end

%%  Get coordinate information

ncid  = netcdf.open(file,'nc_nowrite'); % open

% lon, lat and mask
lonname = ['lon' pos];
latname = ['lat' pos];
mskname = ['mask' pos];
varid = netcdf.inqVarID(ncid,lonname);
lonfull = netcdf.getVar(ncid,varid)';
varid = netcdf.inqVarID(ncid,latname);
latfull = netcdf.getVar(ncid,varid)';

% s-coordinate parameters
for vlist = roms_varlist('s-param')
  vname = char(vlist);
  varid = netcdf.inqVarID(ncid,vname);
  S.(vname)= netcdf.getVar(ncid,varid);
end
N = length(ncread(file,'s_rho'));
S.N = N;

netcdf.close(ncid); % close

%% Determine the minimal model/data index range (Iax, Jax) required for the
%  bounding box encompassing the entire track

[nJ,nI] = size(lonfull);
[I,J] = meshgrid(1:nI,1:nJ);

% Fractional i,j indicies of the trajectory
F    = scatteredInterpolant(lonfull(:),latfull(:),I(:),method);
Itrk = F(lonTrk,latTrk);
Iax  = floor(min(Itrk)):ceil(max(Itrk));
F    = scatteredInterpolant(lonfull(:),latfull(:),J(:),method);
Jtrk = F(lonTrk,latTrk);
Jax  = floor(min(Jtrk)):ceil(max(Jtrk));

% Expand Iax,Jax by 1 in each dimension (if possible) else in some
% instances the round-off in fJ,fI gives trouble with space weights
% Might need to further inflate the Iax,Jax selection when slicing
% velocity because of the step that averages h to velocity points.
Iax = max(Iax(1)-1,1):min(Iax(end)+1,nI);
Jax = max(Jax(1)-1,1):min(Jax(end)+1,nJ);

if verbose
  disp('Spatial subset to be extracted from ')
  disp([' ' file])
  disp('  is for index limits ')
  disp(['  J = ' int2str(Jax(1)) ',' int2str(Jax(end))])
  disp(['  I = ' int2str(Iax(1)) ',' int2str(Iax(end))])
end

%%  Now get fractional track I,J grid on the SUBSET Iax,Jax grid
%  Emphasize: this I,J index is to the *subsetted* grid - not the full grid
%  and it is in i,j convention for the staggered grid for varname

lon_mod = ncread(file,lonname,[Iax(1) Jax(1)],[length(Iax) length(Jax)])';
lat_mod = ncread(file,latname,[Iax(1) Jax(1)],[length(Iax) length(Jax)])';
[I,J] = meshgrid(1:size(lon_mod,2),1:size(lon_mod,1));

if verbose
  disp('scatteredInterpolant is calculating i,j coords of the track ...')
end

F = scatteredInterpolant(lon_mod(:),lat_mod(:),J(:),method);
Jtrk = F(lonTrk,latTrk);
F = scatteredInterpolant(lon_mod(:),lat_mod(:),I(:),method);
Itrk = F(lonTrk,latTrk);

% Find track locations that are outside the model grid or available times
% Flag with valid=NaN. We will skip over these positions but leave NaNs
% in the output so that the space dimension of inputs and outputs match
I = ncinfo(file);
if ~isempty(findstrinstruct(I.Variables,'Name',mskname))
% if nc_isvar(file,mskname)
  msk_mod = ncread(file,mskname,[Iax(1) Jax(1)],[length(Iax) length(Jax)])';
  msk_mod(msk_mod~=1) = NaN;
  F = scatteredInterpolant(lon_mod(:),lat_mod(:),msk_mod(:),method);
  Mtrk = F(lonTrk,latTrk);
else
  Mtrk = ones(size(lonTrk));
end

valid = 1:length(lonTrk);
valid(isnan(Itrk)|isnan(Jtrk)|isnan(Ttrk)|isnan(Mtrk)) = NaN;
if ~any(isfinite(valid))
  error('No track points are encompassed by the data. Check lon/lat/time')
end

if verbose
  disp(' ... scatteredInterpolant done')
end

%% Grid and track metrics are required for across/along track velocity
%  Not sure how imprecise this method is for setting the grid metrics
%  on the track

pm = ncread(file,'pm')';
pn = ncread(file,'pn')';
lonr = ncread(file,'lon_rho')';
latr = ncread(file,'lat_rho')';
try
  angle = ncread(file,'angle')';
  gotangle = true;
catch
  gotangle = false;
end

% Can't use F from above because that might be for u or v
F = scatteredInterpolant(lonr(:),latr(:),pm(:),method);
pm = F(lonTrk,latTrk);
F = scatteredInterpolant(lonr(:),latr(:),pn(:),method);
pn = F(lonTrk,latTrk);
if gotangle
  F = scatteredInterpolant(lonr(:),latr(:),angle(:),method);
  angle = F(lonTrk,latTrk);
end

% 1/m dI/ds and 1/n dJ/ds
% These define the unit vector parallel and perpendicular to the track in
% i,j coordinates for computing tangential and normal velocity
dtmp = diff(Itrk);
dtmp = dtmp([1 1:end end]);
dIom = 0.5*(dtmp(1:(end-1))+dtmp(2:end))./(pm.*dsm);

dtmp = diff(Jtrk);
dtmp = dtmp([1 1:end end]);
dJon = 0.5*(dtmp(1:(end-1))+dtmp(2:end))./(pn.*dsm);

% force row vectors
dIom = dIom(:)';
dJon = dJon(:)';

% In principle, the vectors have unit magnitude, but if the track has major
% kinks in it the discrete calculation can lead to magnitudes that differ
% slightly from 1. For safety, renormalize the vectors

emag = abs(complex(dIom,dJon));
ep = complex( dIom,dJon)./emag; % parallel to track
en = complex(-dJon,dIom)./emag; % normal to track

%% h for z calculation ----------------------------------------------------

h = ncread(file,'h',[Iax(1) Jax(1)],[length(Iax) length(Jax)]);
% shift to same c-grid location as data
h = shift_h_to_uv(igrid,h);

% zeta for z calculation
if zero_zeta
  % Assume zeta=0 everywhere and calculate all z levels accordingly just
  % once here. In the more general case that the z sought is the true
  % varying coordainte, then zeta has to be read for every chunk of the
  % trajectory before z is calculated
  zeta = zeros(size(h));
  z = set_depth(S.Vtransform,S.Vstretching,S.theta_s,S.theta_b,...
    S.hc,S.N,1,h,zeta,0);
  z = permute(z,[3 2 1]);
  zw = set_depth(S.Vtransform,S.Vstretching,S.theta_s,S.theta_b,...
    S.hc,S.N,5,h,zeta,0);
  zw = permute(zw,[3 2 1]);
  z2 = z;
  zw2 = zw;
end

%% Now we are ready to interpolate the data to the track ------------------

% preallocate for speed
if w_points % vertical dimension is different for 'w'
  Tnan = nan([N+1 length(lonTrk)]);
  onez = ones([N+1 1]);
  nanz = nan([N+1 1]);
else
  Tnan = nan([N length(lonTrk)]);
  onez = ones([N 1]);
  nanz = nan([N 1]);
end
Tvar    = Tnan;
Tlon    = Tnan;
Tlat    = Tnan;
Tdis    = Tnan;
Tdislen = Tnan;
% z coordinates always have same vertical dimension
Tz   = nan([N length(lonTrk)]);
Tzw  = nan([N+1 length(lonTrk)]);
dzed = nan([N length(lonTrk)]);
% 1-D vectors of depth and zeta
Th = nan([1 length(lonTrk)]);
Ts = Th;

% Each chunk is the set of points of the track that fall within the bounds
% of a distinct pair of times
chunks = unique(floor(Ttrk));
chunks(isnan(chunks)) = [];

if verbose
  disp(['Processing this track requires ' int2str(length(chunks)) ...
    ' time intervals out of '])
  disp([' ' file])
  disp(' within time index limits ')
  disp([' T = ' int2str(min(Ttrk)) ':' int2str(max(Ttrk))])
end

%% Loop over time chunks --------------------------------------------------

time_t2 = NaN; % for first time through when we test recycling of t2 to t1

% Counter for alongtrack profiles
prof = 1;
nprofs = size(Tvar,2);

hanwaitbar = waitbar(0,['Processing ' int2str(nprofs) ' profiles ...']);
set(get(findobj(hanwaitbar,'type','axes'),'title'),'FontSize',18)

for nchunk = 1:length(chunks)
  
  %   waitbar(nchunk/length(chunks),hw1,...
  %     ['Processing ' int2str(nchunk) ' of ' int2str(length(chunks)) ...
  %     ' time chunks ...']);
  
  % Index of time at beginning of interval
  tindex0 = chunks(nchunk);
  
  % find the chunk of track coordinates in this time interval
  Section = find(floor(Ttrk)==tindex0);
  lonS = lonTrk(Section);
  latS = latTrk(Section);
  timS = timTrk(Section);
  distS = dist(Section);
  distlenS = distlen(Section);
  JtrkS = Jtrk(Section);
  ItrkS = Itrk(Section);
  
  %% Load time this tindex and the next -----------------------------------
  
  time_t2_last = time_t2;
  
  time_t1 = tim_mod(tindex0);
  if (tindex0+1) > length(tim_mod)
    time_t2 = tim_mod(tindex0);
  else
    time_t2 = tim_mod(tindex0+1);
  end
  
  if verbose
    disp([' Doing interval ' int2str(nchunk) ' with ' ...
      int2str(length(Section)) ' coordinates in time range:'])
    disp(['  ' datestr(time_t1) ' to ' datestr(time_t2)])
  end
  
  % Determine if the data for time_t2 on the previous chunk can be recycled
  % to time_t1 on this chunk. This is done to avoid an unnecessary ncread
  % operation
  if time_t2_last == time_t1
    time_shift = true;
  else
    time_shift = false;
  end
  if verbose && time_shift
    disp('  shifting time_t2 data from last chunk to time_t1')
  end
  
  %% Get the data subset - two time levels at a time ----------------------
  
  if time_shift     % move time_2 data to time_1
    data_t1 = data_t2;
  else              % must read time_1 data anew
    tmp = ncread(file,varname,[Iax(1) Jax(1) 1 tindex0],...
      [length(Iax) length(Jax) Inf 1]);
    data_t1 = permute(tmp,[3 2 1]);
  end
  if (tindex0+1) > length(tim_mod)
    data_t2 = data_t1;
  else
    tmp = ncread(file,varname,[Iax(1) Jax(1) 1 tindex0+1],...
      [length(Iax) length(Jax) Inf 1]);
    data_t2 = permute(tmp,[3 2 1]);
  end
  
  %% zeta subset (for z calculation) - two time levels --------------------
  if ~zero_zeta
    
    % zeta at time 1
    if time_shift     % move time_2 data to time_1
      z = z2;
      zw = zw2;
    else              % must read time_1 data anew
      zeta = ncread(file,'zeta',[Iax(1) Jax(1) tindex0],...
        [length(Iax) length(Jax) 1]);
      zeta = shift_h_to_uv(igrid,zeta); % zeta on same grid as data
      z = set_depth(S.Vtransform,S.Vstretching,S.theta_s,S.theta_b,...
        S.hc,S.N,1,h,zeta,0);
      z = permute(z,[3 2 1]);
      zw = set_depth(S.Vtransform,S.Vstretching,S.theta_s,S.theta_b,...
        S.hc,S.N,5,h,zeta,0);
      zw = permute(zw,[3 2 1]);
    end
    
    if (tindex0+1) > length(tim_mod)
      z2 = z;
      zw2 = zw;
    else
      % zeta at time 2
      zeta = ncread(file,'zeta',[Iax(1) Jax(1) tindex0+1],...
        [length(Iax) length(Jax) 1]);
      zeta = shift_h_to_uv(igrid,zeta);
      z2 = set_depth(S.Vtransform,S.Vstretching,S.theta_s,S.theta_b,...
        S.hc,S.N,1,h,zeta,0);
      z2 = permute(z2,[3 2 1]);
      zw2 = set_depth(S.Vtransform,S.Vstretching,S.theta_s,S.theta_b,...
        S.hc,S.N,5,h,zeta,0);
      zw2 = permute(zw2,[3 2 1]);
    end
  end
  
  %% loop over lon/lat ----------------------------------------------------
  
  for j = 1:length(Section)
    
    % each Section is the set of points of the track that fall
    % within the bounds of a distinct pair of times
    if verbose && rem(prof,10)==0
      disp(['   Profile ' int2str(prof) ' of ' int2str(length(lonTrk))])
    end
    if rem(prof-1,20)==0
      waitbar(nchunk/length(chunks),hanwaitbar,...
        ['Processing profile ' int2str(prof) ' of ' int2str(nprofs)]);
    end
    
    if isfinite(valid(prof))
      
      % Interpolate to the track positions at both times using linear
      % weights in I,J fractional coordinates. I guess this could be done
      % with another call to scatteredInterpolant
      
      % space weights
      fJ = floor(JtrkS(j));
      fI = floor(ItrkS(j));
      wJ = JtrkS(j)-fJ;
      wI = ItrkS(j)-fI;
      Dwgt(1,1) = (1-wJ)*(1-wI);
      Dwgt(1,2) = (1-wJ)*wI;
      Dwgt(2,1) = wJ*(1-wI);
      Dwgt(2,2) = wJ*wI;
      
      % time weights
      Twgt_t2 = (timS(j)-time_t1)/(time_t2-time_t1);
      if isnan(Twgt_t2)
        Twgt_t2 = 1;
      end
      Twgt_t1 = 1-Twgt_t2;
      
      % try/catch in case lon/lat point causes fJ,fI to be outside
      % the range of valid data points. For a track point in the cell
      % adjacent to the land mask some values will be invalid. This might
      % be accommodated by modifying weights to only use valid points.
      % Will revisit this some other time (J Wilkin).
      
      %% interpolate data to track ----------------------------------------
      
      try
        % spatial interpolation of data to track(time 1)
        data = data_t1;
        Dslice_t1 = Dwgt(1,1)*data(:,fJ,fI) + Dwgt(1,2)*data(:,fJ,fI+1)+...
          Dwgt(2,1)*data(:,fJ+1,fI) + Dwgt(2,2)*data(:,fJ+1,fI+1);
        
        % spatial interpolation of data to track (time 2)
        data = data_t2;
        Dslice_t2 = Dwgt(1,1)*data(:,fJ,fI) + Dwgt(1,2)*data(:,fJ,fI+1)+...
          Dwgt(2,1)*data(:,fJ+1,fI) + Dwgt(2,2)*data(:,fJ+1,fI+1);
        
        % time interpolation
        Dslice = Twgt_t1*Dslice_t1 + Twgt_t2*Dslice_t2;
        
      catch
        warning('lon/lat out of bounds of some data')
        Dslice = nan(size(data(:,1,1)));
      end
      
      % Enter this profile into the accumulated slice
      Tvar(:,prof) = Dslice;
      Tlat(:,prof) = latS(j)*onez;
      Tlon(:,prof) = lonS(j)*onez;
      Tdis(:,prof) = distS(j)*onez;
      Tdislen(:,prof) = distlenS(j)*onez;
      
      %% interpolate z and zw to track ------------------------------------
      
      try
        % spatial interpolation of z to track (time 1)
        data = z;
        Zslice_t1 = Dwgt(1,1)*data(:,fJ,fI)+Dwgt(1,2)*data(:,fJ,fI+1)+...
          Dwgt(2,1)*data(:,fJ+1,fI)+Dwgt(2,2)*data(:,fJ+1,fI+1);
        
        % spatial interpolation of zw to track (time 1)
        data = zw;
        Zwslice_t1 = Dwgt(1,1)*data(:,fJ,fI)+Dwgt(1,2)*data(:,fJ,fI+1)+...
          Dwgt(2,1)*data(:,fJ+1,fI)+Dwgt(2,2)*data(:,fJ+1,fI+1);
        
        % spatial interpolation of z to track (time 2)
        data = z2;
        Zslice_t2 = Dwgt(1,1)*data(:,fJ,fI)+Dwgt(1,2)*data(:,fJ,fI+1)+...
          Dwgt(2,1)*data(:,fJ+1,fI)+Dwgt(2,2)*data(:,fJ+1,fI+1);
        
        % spatial interpolation of zw to track (time 2)
        data = zw2;
        Zwslice_t2 = Dwgt(1,1)*data(:,fJ,fI)+Dwgt(1,2)*data(:,fJ,fI+1)+...
          Dwgt(2,1)*data(:,fJ+1,fI)+Dwgt(2,2)*data(:,fJ+1,fI+1);
        
        % time interpolation
        Zslice = Twgt_t1*Zslice_t1 + Twgt_t2*Zslice_t2;
        Zwslice = Twgt_t1*Zwslice_t1 + Twgt_t2*Zwslice_t2;
        
      catch
        warning('lon/lat out of bounds of some data')
        Zslice = nan(size(z(:,1)));
        Zwslice = nan(size(zw(:,1)));
      end
      
      % Enter this profile into the accumulated slice
      Tz(:,prof)   = Zslice;
      Tzw(:,prof)  = Zwslice;
      dzed(:,prof) = diff(Zwslice,1);  % layer thicknesses
      Th(prof)     = -Zwslice(1);      % interpolated h
      Ts(prof)     = Zwslice(end);     % interpolated zeta
      
    else % position not valid
      
      % track coordinates are outside the data so enter NaNs for this point
      Tz(:,prof)   = nan([N 1]);
      Tzw(:,prof)  = nan([N+1 1]);
      dzed(:,prof) = nan([N 1]);
      Tvar(:,prof) = nanz;
      Tlat(:,prof) = latS(j)*nanz;
      Tlon(:,prof) = lonS(j)*nanz;
      Tdis(:,prof) = distS(j)*nanz;
      Tdislen(:,prof) = distlenS(j)*nanz;
      
    end
    
    % update profile counter
    prof = prof + 1;
  end
end

close(hanwaitbar)

%% Output
if w_points
  Coord.z  = Tzw;
else
  Coord.z  = Tz;
  Coord.dz = dzed;
end
Coord.zw = Tzw;
Coord.lon = Tlon;
Coord.lat  = Tlat;
Coord.dis   = Tdis;
Coord.time   = repmat(timTrk(:)',[N 1]);
Coord.dislen  = Tdislen;
Coord.h   = Th;
Coord.zeta = Ts;
Coord.en    = en;
Coord.ep     = ep;
if gotangle
  Coord.angle    = angle;
end
Coord.lonTrk = lonTrk(:)';
Coord.latTrk  = latTrk(:)';
Coord.timTrk   = timTrk(:)';
Coord.tindex    = floor(Ttrk)';

function hs = shift_h_to_uv(igrid,h)
% Horizontal average to place extracted subset of h (at Iax,Jax)
% on same c-grid location as data
% This might have issues for tracks that are close to the perimeter
switch igrid
  case 1
    % don't need to shift for rho points
    hs = h;
  case 2
    error('have not debugged this code for psi points variables')
  case 3 % u
    hs = 0.5*(h(1:(end-1),:)+h(2:end,:));
  case 4 % v
    hs = 0.5*(h(:,1:(end-1))+h(:,2:end));
  otherwise
end


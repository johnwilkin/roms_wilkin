function [thedata,thegrid,thehan] = roms_nsview(file,var,ntime,k,grd,nvec_d,uscale,varargin)
% [thedata,thegrid,thehan] = roms_nsview(file,var,time,k,grd,vec_d,uscale,varargin)
%
% Overlay plots from a set of ROMS outputs.
% Useful for superimposing nested grid solutions
%
% Same input syntax as roms_sview except:
% FILE is a cell array of files names
% GRD is a cell array of grid structures, unless isempty(grd) in which
%    case grid structures will be retrieved from the respective files
% If TIME is a vector, each index is used in the respective grids. Input
%    a single datetime object or date string to have the function
%    determine which record to read to best match the times
% If VEC_D is a cell array or vector, the vector decimation values are
%    used in the respective grids.
% If VEC_D is a scalar that decimation ratio is used in all grids
%
% John Wilkin - October 2016
% $Id: roms_nsview.m 533 2019-09-13 13:41:54Z wilkin $

if iscell(file)
  fnames = file;
  nf = size(fnames,2);
end

get_grd = true;
if iscell(grd)
  gnames = grd;
  if size(gnames,2)~=nf
    error('FILE and GRD lists don''t match in number')
  end
  get_grd = false;
end

% preallocate outputs from roms_sview
thedata = cell([1 nf]);
thegrid = cell([1 nf]);
thehan  = cell([1 nf]);

% get plot state
nextplotstatewas = get(gca,'nextplot');

for i=1:nf
  
  file = fnames{i};
  if get_grd
    grd = roms_get_grid(file,file);
  else
    grd = gnames{i};
  end
  
  add_vectors = false;
  if nargin > 5
    add_vectors = true;
    if iscell(nvec_d)
      vec_d = nvec_d{i};
    elseif isscalar(nvec_d)
      vec_d = nvec_d;
    else
      vec_d = nvec_d(i);
    end
  end
  
  if iscell(ntime)
    time = ntime{i};
  elseif isdatetime(ntime)
    time = datestr(ntime);
  elseif ischar(ntime)
    time = ntime;
  elseif isscalar(ntime)
    time = ntime;
  else % vector of time indices one for each grid
    time = ntime(i);
  end
  
  if i>1
    % draw a white polygon to obscure the parent grid - otherwise where the
    % child is land the parent solution shows through, which looks terrible
    % where there is a lot of detail in the coastline and especially rivers
    X = [squash(grd.lon_psi(1:end,1)); ...
      squash(grd.lon_psi(end,1:end)); ...
      flipud(squash(grd.lon_psi(1:end,end))); ...
      flipud(squash(grd.lon_psi(1,1:end)))];
    Y = [squash(grd.lat_psi(1:end,1)); ...
      squash(grd.lat_psi(end,1:end)); ...
      flipud(squash(grd.lat_psi(1:end,end))); ...
      flipud(squash(grd.lat_psi(1,1:end)))];
    han = fill(X,Y,'w');
    han.EdgeColor = 'w';
  end
  
  try
    if getpref('ROMS_WILKIN','VERBOSE')
      disp(['plotting from ' file])
    end
  catch
  end
  if add_vectors
    [thedata{i},thegrid{i},thehan{i}] = ...
      roms_sview(file,var,time,k,grd,vec_d,uscale,varargin{:});
  else
    [thedata{i},thegrid{i},thehan{i}] = ...
      roms_sview(file,var,time,k,grd);
  end
  
  % hold so we can plot the next file over
  set(gca,'nextplot','add')
  
end

% restore nextplotstate to what it was
set(gca,'nextplot',nextplotstatewas);

function y = squash(x)
% SQUASH(X) = X(:) saves needing an intermediate variable when you want
% to squash an array subset, say X(I,J,K), into a vector to do something
% else to it (like make it the argument of MIN, MAX, ALL etc).
%
% John Wilkin
y = x(:);

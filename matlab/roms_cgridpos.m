function pos = roms_cgridpos(s,grd)
% function pos = roms_cgridpos(s,grd)
%
% From the size of a variable determine whether it is defined on the
% u, v, rho or psi location on the ROMS Arakawa-C grid
%
% John Wilkin
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_cgridpos.m 440 2016-02-18 20:34:17Z wilkin $

if ndims(s)~=2
  s = size(s);
else
  ss = size(s);
  if ss(1)~= 1
    s = size(s);
  end
end

% su = size(grd.lon_u);
% sv = size(grd.lon_v);
sr = size(grd.lon_rho);
if isfield(grd,'lon_psi')
  sp = size(grd.lon_psi);
end
if sr(1)==sr(2)
  warning([ 'The grid is square so the result of ' ...
    which(mfilename) ' may be unreliable']);
end

% only the right-most 2 dimensions matter
s = s([end-1 end]);
if all(~(s-sr))
  pos = 'rho';
elseif all(~(s-size(grd.lon_u)))
  pos = 'u';
elseif all(~(s-size(grd.lon_v)))
  pos = 'v';
elseif all(~(s-sp))
  pos = 'psi';
else
  error( 'Could not determine C-grid location')
end

function hanout = roms_plot_bathy(grd,cmap,clev,var)
% han = roms_plot_bathy(grd,cmap,clev,var);
% 
% defaults:
%    clev = [100 250 500 1000:1000:4000];
%    cmap = zebra(2,64);
%    var = 'h' (if var='r' r-value is plotted)
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_plot_bathy.m 387 2009-10-14 13:30:50Z wilkin $

if ischar(grd)
  grd_file = grd;
  grd = roms_get_grid(grd_file);
else
  grd_file = ' ';
end

if nargin < 2
  colormap(zebra(2,64))
else
  if isempty('cmap')
    colormap(zebra(2,64))
  else
    colormap(cmap)
  end
end

if nargin < 3
  clev = [100 250 500 1000:1000:4000];
else
  if isempty('clev')
    clev = [100 250 500 1000:1000:4000];
  end
end

if nargin < 4
  var = 'h';
end

switch var
  case 'h'
    hanp = pcolorjw(grd.lon_rho,grd.lat_rho,grd.h.*grd.mask_rho_nan);
  case 'r'
    hanp = pcolorjw(grd.lon_rho,grd.lat_rho,rvalue(grd.h).*grd.mask_rho_nan);
  case 'cfl'
    cfl = min(1./(grd.pm),1./grd.pn)./sqrt(9.81*grd.h);
    hanp = pcolorjw(grd.lon_rho,grd.lat_rho,cfl.*grd.mask_rho_nan);
    [worst,where] = sort(cfl(:));
    disp(num2str(worst(1)))
    ncheck = 30;
    worst = worst(1:ncheck);
    where = where(1:ncheck);
    hold on
    plot(grd.lon_rho(where),grd.lat_rho(where),'m*')
    [grd.lon_rho(where) grd.lat_rho(where) grd.h(where) cfl(where)]
    hold off
end
han.pcolor = hanp;

if numel(clev)==2
  caxis(clev)
else
  caxis([min(clev) max(clev)])
  if length(clev) > 2
    hold on
    [~,hanc] = contour(grd.lon_rho,grd.lat_rho,grd.h.*grd.mask_rho_nan,clev);
    set(hanc,'edgecolor',0.5*[1 1 1])
    hold off
    han.contours = hanc;
  end
end

set(gca,'tickdir','out')
titlestr{1} = ['Model bathymetry ' strrep(grd_file,'_','\_')];
switch var
  case 'r'
    titlestr{2} = 'r-value';
  case 'cfl'
    titlestr{2} = '\Deltat max < min(\Deltax,\Deltay)/sqrt(gh)';
end
title(titlestr,'fontsize',12)
han.colorbar = colorbar('h');

amerc
han.colorbar.Position = [0.45 0.16 0.35 0.02];

if nargout > 0
  hanout = han;
end


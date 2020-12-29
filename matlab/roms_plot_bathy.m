function han = roms_plot_bathy(grd,cmap,clev,var);
% $Id: roms_plot_bathy.m 387 2009-10-14 13:30:50Z wilkin $
% han = roms_plot_bathy(grd,cmap,clev,var);
% 
% defaults:
%    clev = [100 250 500 1000:1000:4000];
%    cmap = zebra(2,64);
%    var = 'h' (if var='r' r-value is plotted)

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
    pcolorjw(grd.lon_rho,grd.lat_rho,grd.h.*grd.mask_rho_nan)
  case 'r'
    pcolorjw(grd.lon_rho,grd.lat_rho,rvalue(grd.h).*grd.mask_rho_nan)
  case 'cfl'
    cfl = min(1./(grd.pm),1./grd.pn)./sqrt(9.81*grd.h);
    pcolorjw(grd.lon_rho,grd.lat_rho,cfl.*grd.mask_rho_nan)
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
if isstr(clev)
  caxis(clev)
else
  caxis([min(clev) max(clev)])
  if length(clev) > 2
    hold on
    [cs,hanc] = contour(grd.lon_rho,grd.lat_rho,grd.h.*grd.mask_rho_nan,clev);
    set(hanc,'edgecolor','k')
    hold off
  end
end
amerc
grid on
set(gca,'tickdir','out')
titlestr{1} = ['Model bathymetry ' strrep(grd_file,'_','\_')];
switch var
  case 'r'
    titlestr{2} = 'r-value';
  case 'cfl'
    titlestr{2} = '\Deltat max < min(\Deltax,\Deltay)/sqrt(gh)';
end
title(titlestr,'fontsize',12)
colorbar('h')

if nargout > 0
  han = hanc;
end


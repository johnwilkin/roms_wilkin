function [Snew,Pobj,Hobj] = roms_isosurface(file,varname,time,value,g,azel)
% [S,P,H] = roms_isosurface(file,varname,time,value,grd,azel)
%
% Plot an isosurface of a ROMS 3-d variable
% Presently configured to plot over a rendition of the model bathymetry
%
% Inputs:
%
% file    = roms his/avg/rst/dia etc netcdf file or THREDDS data URL
% varname = name of the 3-D ROMS output variable to plot
% time    = time index into FILE
% value   = value of the isosurface to render

% grd     = grd structure (from roms_get_grid)
%
% THIS CODE STILL UNDER DEVELOPMENT. PLAN IS TO GIVE IT CALLING
% FUNCTIONALITY SIMILAR TO ROMS_ZVIEW AND ROMS_SVIEW ETC. 
% See roms_isosurface_demo.m for examples on how to use this function.
%
% John Wilkin - Jan 2021
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu

% For now only support input 'time' as the tindex into the file. 
tindex = time;

% read the data
data = ncread(file,varname,[1 1 1 tindex],[Inf Inf Inf 1]);
data = squeeze(data);

% figure out where these data lie on the ROMS staggered C-grid
pos = roms_cgridpos(permute(data,[3 2 1]),g);
ma = ['mask_' pos];
lo = ['lon_' pos];
la = ['lat_' pos];
zz = ['z_' pos(1)];
    
% I haven't figured out an easy way to clip the isosurface after it is
% rendered, but you can hack things here by setting to NaN parts of DATA
% you don't want to visualize
% % exclude the shelf-break front etc. beyond the MAB 100-m isobath
% h = permute(g.h,[2 1]);
% m = h;
% m(m>100) = NaN;
% m(m<=100) = 1;
% M = repmat(m,[1 1 size(data,3)]);
% data = data.*M; 

% Isosurface works in plaid i,j,k space, so the surface will be rendered 
% in i,j,k coordinates that will then be mapped to lon,lat,z for display.
% This is by far preferable to regridding the data to a uniform lon/lat
% grid.

% set up the i,j,k, to lon,lat,z transformation
[I3,J3,K3] = ndgrid(1:size(data,1),1:size(data,2),1:size(data,3));
[I2,J2] = ndgrid(1:size(data,1),1:size(data,2));

Fx = griddedInterpolant(I2,J2,permute(g.(lo),[2 1]));
Fy = griddedInterpolant(I2,J2,permute(g.(la),[2 1]));
Fz = griddedInterpolant(I3,J3,K3,permute(g.(zz),[3 2 1]));

% get the isosurface in i,j,k space
S = isosurface(I3,J3,K3,data,value);
Snew = S;

% convert fractional i,j to lon/lat and fractional k to z by interpolation
Snew.vertices(:,1) = Fx(S.vertices(:,1),S.vertices(:,2));
Snew.vertices(:,2) = Fy(S.vertices(:,1),S.vertices(:,2));
Snew.vertices(:,3) = Fz(S.vertices(:,1),S.vertices(:,2),S.vertices(:,3));

%% plot -------------------------------------------------------------------

% bathymetry
hviz = -g.h; % bathymetry needs to be negative to plot correctly
hviz(g.mask_rho==0) = 10; % make land positive so it will get clipped
Hobj = surf(g.lon_rho,g.lat_rho,hviz,-hviz);
shading flat
caxis([-50 300])
colormap(demcmap(caxis,128,flipud(cmap_land),cmap_haxby))

% isosurface
hold on
Pobj = patch(Snew);
hold off

% set appearance of the isosurface
Pobj.FaceColor = 'b';
Pobj.EdgeColor = 'none';
Pobj.FaceAlpha = 0.5;
hanlightb = camlight(0,50);
hanlightb.Position = [-74 -50 500];
lighting phong

% 3-D view point
if nargin > 5
  view(azel)
end

% set axis limits according to the span of the isosurface
V = Snew.vertices;
xlim([min(V(:,1)) max(V(:,1))])
ylim([min(V(:,2)) max(V(:,2))])
zlim([min(V(:,3)) max(V(:,3))])






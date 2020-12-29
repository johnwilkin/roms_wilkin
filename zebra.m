function map = zebra(a,n,m);
% $Id: zebra.m 358 2008-04-07 14:15:03Z zhang $
% zebra palette colormap with NBANDS broad bands and NENTRIES rows in
% the color map - just try it, e.g. colormap(zebra) 
%
% The default is 4 broad bands
%
% MAP = ZEBRA(NBANDS,NENTRIES)
%
% see Hooker, S. B. et al, Detecting Dipole Ring Separatrices with Zebra 
% Palettes, IEEE Transactions on Geosciences and Remote Sensing, vol. 33,
% 1306-1312, 1995
%
% John Wilkin

if nargin < 3
  m = 0.5; % saturation and value go from m to 1
  % don't use m = 0
end
if nargin < 2
   n  = size(get(gcf,'colormap'),1);  % number of entries in the colormap
end
if nargin < 1
   a = 4; % there are this many large bands in the palette
end

x = 0:(n-1);
hue = exp(-3*x/n);
sat = m+(1-m)*(0.5*(1+sawtooth(2*pi*x/(n/a))));
val = m+(1-m)*0.5*(1+cos(2*pi*x/(n/a/2)));

map = [hue(:) sat(:) val(:)];
map = hsv2rgb(map);

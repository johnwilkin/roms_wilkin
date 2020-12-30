function han = roms_addvect_scale(pos,uv,uscale,label,varargin)
% $Id: roms_addvect_scale.m 480 2017-08-01 18:12:47Z wilkin $
%  han = roms_addvect_scale(pos,uv,uscale,label,varargin)
%
%  pos = [x y] position for scale vector (can use, e.g. ginput(1)
%  uv = [u v] scale vector
%  uscale = same scale as used in the roms_quiver command
%  label = string to label vector, default is 'm/s'
%  varargin = arguments to quiver, e.g. vector color
%
% John Wilkin
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id$
  
% add a scale vector to the plot

% get plot state
nextplotstatewas = get(gca,'nextplot');

% hold whatever is already plotted
set(gca,'nextplot','add')

% keyboard

% hq = quiver('v6',pos(1),pos(2),uscale*uv(1),uscale*uv(2),0,varargin{:});
hq = quiver(pos(1),pos(2),uscale*uv(1),uscale*uv(2),0,varargin{:});

if nargin < 4
  label = 'm/s';
end
if isempty(label)
  label = 'm/s';
end

% interactive text placement
% ht = gtext([num2str(uscale) ' m/s']);

umag = abs(uv(1)+sqrt(-1)*uv(2));
ht = text(pos(1),pos(2),[num2str(umag) ' ' label ' ']);
set(ht,'horizontalalignment','right')

if nargout > 0
  han = [hq];
  if exist('ht')==1
    han = [han; ht];
  end
end

% restore nextplotstate to what it was
set(gca,'nextplot',nextplotstatewas);


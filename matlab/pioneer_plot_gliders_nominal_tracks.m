function [Hanglider,Hanlegend] = pioneer_plot_gliders_nominal_tracks(leg)
% [hanglider,hanlegend] = pioneer_plot_gliders_nominal_tracks(leg)
% Pioneer glider sampling plan
%
% Input
%     leg = option to plot a legend
%
% Outputs
%     hanglider = handle to glider tracks
%     hanlegend = handle to legend object   
%
% adapted from diamondgliders.m from Pioneer Focus group COL April 2012 
% John Wilkin 13 Aug 2015
% $Id: pioneer_plot_gliders_nominal_tracks.m 585 2020-10-12 20:20:53Z wilkin $

% get plot state
nextplotstatewas = get(gca,'nextplot');

% hold whatever is already plotted
set(gca,'nextplot','add')

% restore nextplotstate to what it was
% set(gca,'nextplot',nextplotstatewas);

%% EB - pivotal
xe = -70;
xg1 = [xe xe -70.19 -70.19 xe];
yg1 = [40.67 39.83 39.83 40.09 40.67];
hanglider(1) = plot(xg1,yg1,'g-');
d1 = sum(sw_dist(yg1,xg1,'km'));

%% SS-P - pivotal
ys = 39.3333;
xw = -71.1667;
xg = linspace(xw,xe,5);
xg2 = xg([1 1 3 5 5 3 1]);
yg2 = [39.76 ys 39.82 ys 39.83 ys 39.76];
hanglider(2) = plot(xg2,yg2,'b-');
d2 = sum(sw_dist(yg2,xg2,'km'));

%% SS-D - default interlaced diamonds
%  aka "argyle"
xg = linspace(xw,xe,5);
xg3 = xg([1 2 4 5 4 2 1]);
ym = 39.5;
yg3 = [ym 39.82 ys ym 39.83 ys ym];
hanglider(3) = plot(xg3,yg3,'c-');
d3 = sum(sw_dist(yg3,xg3,'km'));

%% FZ - pivotal (and default second glider)
xg4 = [-70.3750  -71.1667  -71.1667  -70.3750  -70.3750];
yg4 = [40.0833   40.0833   39.8333   39.8333   40.0833];
hanglider(4) = plot(xg4,yg4,'r-');
d4 = sum(sw_dist(yg4,xg4,'km'));

%% SS-G - default SS to Gulf Stream
auvns = [
  -70.8833   39.9400
  -70.8833   40.4800
  -70.7071   40.4800
  -70.7071   39.9400
  -70.8833   39.9400 ];
xg5 = auvns([1 3 3 1 1]);
yg5 = [39.9 39.9 39.1 39.1 39.9];
hanglider(5) = plot(xg5,yg5,'w-');
d5 = sum(sw_dist(yg5,xg5,'km'));
set(hanglider(5),'color',0.2*[1 1 1])

%% legend
opt = 0;
if nargin > 1
  opt = leg;
end
switch opt
  case 3
    hanlegend = legend(hanglider,...
      [int2str(d1) ' km'],...
      [int2str(d2) ' km'],...
      [int2str(d3) ' km'],...
      [int2str(d4) ' km'],...
      [int2str(d5) ' km']);
    set(hanlegend,'position',[.35 .7 .03 .1]);
    legend boxoff
    set(hanglider,'linewidth',2)
  case 2
    % as time
    k2d = 1/25; % convert km to days at 25 km/day
    hanlegend = legend(hanglider,...
      [int2str(d1*k2d) ' days at 25 km/day'],...
      [int2str(d2*k2d) ' days'],...
      [int2str(d3*k2d) ' days'],...
      [int2str(d4*k2d) ' days'],...
      [int2str(d5*k2d) ' days']);
    set(hanlegend,'position',[.35 .7 .03 .1]);
    legend boxoff
    set(hanglider,'linewidth',2)
  case 1
    hanlegend = legend(hanglider,'EB','SS-1','SS-2','FZ','GS');
    set(hanlegend,'position',[.35 .5 .03 .1]);
    hanlegend.Position = [.25 .75 .1 .1];
    legend boxoff
    set(hanglider,'linewidth',1)    
end

% restore nextplotstate to what it was
set(gca,'nextplot',nextplotstatewas);

if nargout > 0
  Hanglider = hanglider;
end
if nargout > 1
  Hanlegend = hanlegend;
end



function [han,mooring,hant] = pioneer_plot_mooring_sites(varargin)
% [han,m,text] = pioneer_plot_mooring_sites(varargin)
%  where inputs arguments are any valid options to PLOT function
%
% Outputs
%     han  = handle to plot
%       m  = structure with coordinates, labels and long names
%    text  = handle to labels
%
% John Wilkin 13 August 2015
%
% Usage e.g.:
% plot_pioneer_mooring_sites('ko','markersize',15,'markerfacecolor','w')
% plot_pioneer_mooring_sites('r+','markersize',10)
%
% If the first option is one of these strings then text labels are added
% plot_pioneer_mooring_sites('label',...) for long name
% plot_pioneer_mooring_sites('LABEL',...) for designator e.g CP01CNPM
% plot_pioneer_mooring_sites('number',...) number  
% plot_pioneer_mooring_sites('Number',...) number and long name
% plot_pioneer_mooring_sites('NUMBER',...) number and designator
%
% Mooring sites
% location data corrected on Aug 7, 2018, from longitude/latitude data at
% at http://www.myroms.org:8080/erddap/info/index.html superceding nominal
% information at http://oceanobservatories.org/array/coastal-pioneer
% $Id: pioneer_plot_mooring_sites.m 585 2020-10-12 20:20:53Z wilkin $

% get plot state
nextplotstatewas = get(gca,'nextplot');

% hold whatever is already plotted
set(gca,'nextplot','add')

% at the end restore nextplotstate to what it was
% set(gca,'nextplot',nextplotstatewas);

%%

mlist = 1:7;

% change this to plot only a subset of moorings - here 2014 for PNI
% mlist = [1 3 5 6 7];

mlabel = false;
plotopts{1} = '^k';
if nargin > 0
  if strcmpi('label',varargin{1}) || strcmpi('number',varargin{1})
    mlabel = true;
    opt = varargin{1};
    if nargin > 1
      for k=2:nargin
        plotopts{k-1} = varargin{k};
      end
    end
  else
    plotopts = varargin;
  end
end

% The list of entries here is in the order at myroms.org ERDDAP
% http://www.myroms.org:8080/erddap/info/index.html
% The number (m) corresponds to the diagram
% https://oceanobservatories.org/wp-content/uploads/2018/03/PioneerArray_2017_labels.png
% which is the same numbering convention in the Levin et al. observation
% impacts papers (I and II)

m = 2;
mooring(m).lon = -70.77952;
mooring(m).lat =  40.13922;
mooring(m).dep =  133;
mooring(m).label = '2';
mooring(m).longname = 'Central Profiler Mooring';
mooring(m).designator = 'CP01CNPM';

m = 5;
mooring(m).lon = -70.88;
mooring(m).lat =  40.226;
mooring(m).dep =  127;
mooring(m).label = '5';
mooring(m).longname = 'Central Inshore Profiler';
mooring(m).designator = 'CP02PMCI';

m = 6;
%ooring(m).lon = -70.8818;
%ooring(m).lat =  40.1962;
mooring(m).lon = -70.8797;
mooring(m).lat =  40.0961;
mooring(m).dep =  148;
mooring(m).label = '6';
mooring(m).longname = 'Central Offshore Profiler';
mooring(m).designator = 'CP02PMCO';

m = 1;
mooring(m).lon = -70.775;
mooring(m).lat =  40.364;
mooring(m).dep =  95;
mooring(m).label = '1';
mooring(m).longname = 'Upstream Inshore Profiler';
mooring(m).designator = 'CP02PMUI';

m = 3;
mooring(m).lon = -70.785;
mooring(m).lat =  39.94;
mooring(m).dep =  452;
mooring(m).label = '3';
mooring(m).longname = 'Upstream Offshore Profiler';
mooring(m).designator = 'CP02PMUO';

m = 4;
mooring(m).lon = -70.8817;
mooring(m).lat =  40.36473;
mooring(m).dep =  92;
mooring(m).label = '4';
mooring(m).longname = 'Inshore Surface-piercing Profiler';
mooring(m).designator = 'CP03ISPM';

m = 7;
mooring(m).lon = -70.88;
mooring(m).lat =  39.94;
mooring(m).dep =  453;
mooring(m).label = '7';
mooring(m).longname = 'Offshore Profiler';
mooring(m).designator = 'CP04OSPM';

for m=mlist
  M.lon(m) = mooring(m).lon;
  M.lat(m) = mooring(m).lat;
end

han = plot(M.lon,M.lat,plotopts{:});
if mlabel
  for m=mlist
    % pad the string at left so it's offset from the symbol
    switch opt(1:2)
      case 'la'
        str = mooring(m).longname;
      case 'LA'
        str = mooring(m).designator;
      case 'nu'
        str = int2str(m);
      case 'Nu'
        str = ['(' int2str(m) ') ' mooring(m).longname];
      case 'NU'
        str = ['(' int2str(m) ') ' mooring(m).designator];
    end
    rot = 0;
    switch m
      case {1,2,3}
        hant(m) = text(mooring(m).lon,mooring(m).lat,['   ' str]);
        hant(m).HorizontalAlignment = 'left';
      case {4,5,6,7}
        hant(m) = text(mooring(m).lon,mooring(m).lat,[str '   ']);
        hant(m).HorizontalAlignment = 'right';
    end
    hant(m).Rotation = rot;
  end
end

% restore nextplotstate to what it was
set(gca,'nextplot',nextplotstatewas);


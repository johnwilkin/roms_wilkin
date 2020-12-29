function [han,mooring,hant] = pioneer_mplot_mooring_sites(varargin)
% [han,m,text] = plot_pioneer_mooring_sites(varargin)
%  where inputs arguments are any valid options to PLOT function
%
% Outputs
%     han  = handle to plot
%       m  = structure with coordinates, labels and long names
%    text  = handle to labels
%
% John Wilkin 13 August 2015
% converted for m_map Julia Levin
%
% Usage e.g.:
% plot_pioneer_mooring_sites('ko','markersize',15,'markerfacecolor','w')
% plot_pioneer_mooring_sites('r+','markersize',10)
%
% If the first option is the string 'label' then text labels are added
% plot_pioneer_mooring_sites('label',...) for long name
% plot_pioneer_mooring_sites('LABEL',...) for designator e.g CP01CNPM
%
% Mooring sites
% location data corrected on Aug 7, 2018, from longitude/latitude data at
% at http://www.myroms.org:8080/erddap/info/index.html superceding nominal
% information at http://oceanobservatories.org/array/coastal-pioneer
% $Id: pioneer_mplot_mooring_sites.m 573 2020-05-08 19:33:31Z wilkin $

% get plot state
nextplotstatewas = get(gca,'nextplot');

% hold whatever is already plotted
set(gca,'nextplot','add')

% at the end restore nextplotstate to what it was
% set(gca,'nextplot',nextplotstatewas);

%%

mlist = 1:7;

% change this to plot only a subset of moorings - here 2014 for PNI
%mlist = [1 3 5 6 7];

mlabel = false;
plotopts{1} = '^k';
if nargin > 0
  if strcmpi('label',varargin{1})
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

% list is in the order of entries at myroms.org ERDDAP
% number (m) corresponds to the diagram
% https://oceanobservatories.org/wp-content/uploads/2018/03/PioneerArray_2017_labels.png

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

han = m_plot(M.lon,M.lat,plotopts{:});
if mlabel
  for m=mlist
    % pad the string at left so it's offset from the symbol
    if strcmp(opt(1),'L')
      str = ['   ' mooring(m).designator];
    else
      str = ['   ' mooring(m).longname];
    end
    hant(m) = text(mooring(m).lon,mooring(m).lat,str);
    hant(m).HorizontalAlignment = 'left';
    hant(m).Rotation = 40;
  end
end

% restore nextplotstate to what it was
set(gca,'nextplot',nextplotstatewas);


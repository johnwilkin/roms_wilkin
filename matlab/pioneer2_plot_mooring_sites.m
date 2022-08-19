function [han,hant,mooring] = pioneer2_plot_mooring_sites(varargin)
% [hansym,hantext,mooringdata] = pioneer_plot_mooring_sites(varargin)
% where inputs arguments are a NUMBER to set the MARKERSIZE, and/or
% a string to indicate label type. any valid options to PLOT function
%
% If the string option is one of these then text labels are added
% pioneer2_plot_mooring_sites('label',...) for long name
% pioneer2_plot_mooring_sites('LABEL',...) for designator e.g CP01CNPM
% pioneer2_plot_mooring_sites('number',...) number
% pioneer2_plot_mooring_sites('Number',...) number and long name
% pioneer2_plot_mooring_sites('NUMBER',...) number and designator
%
% Outputs
%     han  = handle to plot
%       m  = structure with coordinates, labels and long names
%    text  = handle to labels
%
% John Wilkin 13 August 2015
%
% Usage e.g.:
% pioneer2_plot_mooring_sites(['k^','markersize',15,'markerfacecolor','w')
% pioneer2_plot_mooring_sites('r+','markersize',10)
%

% Mooring sites
%
% Copyright (c) 2022 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: pioneer_plot_mooring_sites.m 585 2020-10-12 20:20:53Z wilkin $

% get plot state
nextplotstatewas = get(gca,'nextplot');

% hold whatever is already plotted
set(gca,'nextplot','add')

% at the end restore nextplotstate to what it was
% set(gca,'nextplot',nextplotstatewas);


% change this to plot only a subset of moorings - here 2014 for PNI
mlist = 1:10;
% mlist = [1 3 5 6 7];

mlabel = false;
msize = 15; % default markersize
for k=1:nargin
  if strcmpi('label',varargin{k}) || strcmpi('number',varargin{k})
    mlabel = true;
    opt = varargin{k};
  else
    msize = varargin{k};
  end
end

% The list of entries here is in the order at myroms.org ERDDAP
% http://www.myroms.org:8080/erddap/info/index.html
% The number (m) corresponds to the diagram
% https://oceanobservatories.org/wp-content/uploads/2018/03/PioneerArray_2017_labels.png
% which is the same numbering convention in the Levin et al. observation
% impacts papers (I and II)

L = string({'SMW';'SME';'NPM';'CPM';'SPM';'NOPM';...
  'SOPM';'NSM';'CSM';'SSM'});

P = [...
  35	56.587	75	20.94  25
  35	56.587	75	7.5    35
  36	9.757	  74	50.351 94
  35	56.587	74	53.461 89
  35	42.86  	74	52.255 82
  36	3.172 	74	44.994 600
  35	49.724	74	43.533 600
  36	9.757 	74	49.701 100
  35	56.587	74	52.821 100
  35	42.86	  74	51.599 100 ];

for m = 1:10
  mooring(m).lon = -P(m,3)-P(m,4)/60;
  mooring(m).lat = P(m,1)+P(m,2)/60;
  mooring(m).dep = P(m,5);
  mooring(m).label = int2str(m);
  mooring(m).longname = L(m);
  mooring(m).designator = L(m);
end

for m=mlist
  M.lon(m) = mooring(m).lon;
  M.lat(m) = mooring(m).lat;
end

sym = char('rs','rs','b^','b^','b^','b^','b^','go','go','go');
for m=mlist
  han(m) = plot(M.lon(m),M.lat(m),sym(m,:));
  han(m).MarkerSize = msize;
  han(m).MarkerFaceColor = sym(m,1);
  switch m
    case {3,4,5,6,7}
      han(m).MarkerEdgeColor = 'w';
    otherwise
      han(m).MarkerEdgeColor = 'k';
  end
end

if ~mlabel
  hant = [];
end
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
      case {1,2}
        rot = 45;
        hant(m) = text(mooring(m).lon,mooring(m).lat,['    ' str]);
        hant(m).HorizontalAlignment = 'left';
      case {3,4,5,6,7}
        rot = 0;
        hant(m) = text(mooring(m).lon,mooring(m).lat,['    ' str]);
        hant(m).HorizontalAlignment = 'left';
      case {8,9,10}
        rot = 45;
        hant(m) = text(mooring(m).lon,mooring(m).lat,[str '    ']);
        hant(m).HorizontalAlignment = 'right';
    end
    hant(m).Rotation = rot;
  end
  set(hant,'FontSize',msize)
end

% restore nextplotstate to what it was
set(gca,'nextplot',nextplotstatewas);


function [Data,hax] = roms_bviews(file,varname,time,blist,grd)
% [DATA,HAX] = roms_bviews(file,var,time,blist,grd)
% $Id: roms_bviews.m 522 2018-12-12 19:32:06Z wilkin $
%
% Like roms_bview but stitches together several boundaries given in the
% boundary list BLIST in the order given using a set of axes adjacent to
% each other. (This could be modified to concatenate the data into a single
% plot view). 
%
% Inputs:
%   file   = roms his/avg/rst etc nc file
%   var    = variable to plot (without _west etc.)
%   time   = time index into nc file
%   blist e.g. {'west','south','east','north'}
%   grd can be
%       grd structure (from roms_get_grid)
%       grd_file name
%       [] (will attempt to get grid from roms file)
%
% Note that the xcoord argument ordinarily passed to roms_bview is ignored
% because the x-axis defaults to distance for each panel in the plot
%
% Output:
%   DATA is a structure of structures with the plotting data returned by
%     each call to roms_biview for the respective boundaries in the list
%   HAX is a vector of handles to the plot axes
%
% John Wilkin - Sept 2018

nax = length(blist);
axm = floor(nax/2)+1;

clf
hax = nfigaxes([1 nax],[0 0],[0.1 0.9],[0.2 0.8]);
ymin = Inf;
ymax = -Inf;

for k=1:nax
  
  bndy = char(blist{k});
  data = roms_bview(file,varname,time,bndy,grd,'stitch');
  
  switch bndy
    case {'west','north'}
      fn = fieldnames(data);
      for v = [1 2 3 5 6]
        data.(fn{v}) = fliplr(data.(fn{v}));
      end
      data.dist = fliplr(data.dist(:,end)-data.dist);
  end
  
  axes(hax(k))
  hant = pcolorjw(data.dist,data.z,data.mask.*data.var);
  xlim(range(data.dist));
  
  if k==axm
    hax(k).Title.String = data.tstr;
    hax(k).Title.FontWeight = 'normal';
    xlabel('distance (km) in each sector')
  end
  
  tmp = ylim;
  ymin = min(ymin,tmp(1));
  ymax = max(ymax,tmp(2));
  
  if k>1
    hold on
    plot([0 0],[ymin ymax],'w')
    hold off
  end
  
  Data(k).data = data;
  
end

hax(2).YTick = [];
hax(3).YTick = [];
linkprop(hax,'YLim');
ylim([ymin ymax])
linkprop(hax,'CLim');






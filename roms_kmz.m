function roms_kmz(filename,varargin)
% Converts the current figure (gcf) to a google earth kmz file
% Designed to work with e.g. roms_zview type plots but in principle 
% should work with just about any plot in lon/lat coordinates
%
% Usage:
%        roms_kmz(outputfilename)
%
%
% John Wilkin 7 March 2010
%
% need export_fig from
% http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig
% to create a png that preserves the transparency of the background
%
% need matlab googleearth toolbox from
% http://code.google.com/p/googleearthtoolbox/

if exist('export_fig','file')~=2
  warning('ROMS_KMZ:NeedExportFig',...
  'export_fig is not in the default matlabpath')
  disp('Trying to add export_fig from /home/om/matlab')
  if exist('/home/om/matlab/export_fig','dir')==7
    addpath /home/om/matlab/export_fig -end
    disp('Success')
  else
    error(['Need export_fig from http://www.mathworks.com/' ...
      'matlabcentral/fileexchange/23629-exportfig'])
  end
end
if exist('ge_kml','file')~=2
  warning('ROMS_KMZ:NeedGoogleearth',...
  'googleearth toolbox is not in the default matlabpath')
  disp('Trying to add googleearth from /home/om/matlab')
  if exist('/home/om/matlab/googleearth','dir')==7
    addpath /home/om/matlab/googleearth -end
    disp('Success')
  else
    error(['Need googleearth toolbox for matlab from ',...
      'http://code.google.com/p/googleearthtoolbox/'])
  end
end

keep_colorbar = 0;

% parse varargin
for k=1:length(varargin)
  if strcmp(varargin{k},'colorbar')
    keep_colorbar = 1;
  end
end

% turn off colorbar because it interferes with determining the axis 
% limits that define the edges of the image bitmap. If you want a 
% colorbar on the plot you need to place it inside the axis limits. 
% See the colorbar help for simple default placements of the colorbar
% inside the axes
if ~keep_colorbar
  colorbar off
end

% make the plot fill the frame. It won't matter if this gives a 
% dreadful aspect ratio because googleearth controls the georeferencing
% of the image
set(gca,'position',[0 0 1 1],'DataAspectRatioMode','auto');

% save the W, E, S, N limits
a = axis;

% make the NaN values and background transparent and remove axis lines
set(gca,'color','none')
set(gcf,'color','none')
axis off

% save the plot a png bitmap
pngfile = [filename '.png'];
export_fig(pngfile,'-png','-nocrop');

% googleearth toolbox
kmlStr = ge_groundoverlay(a(4),a(2),a(3),a(1),'imgURL',pngfile,'viewBoundScale',1e3);   
kmlFileName = [filename,'.kml'];
kmzFileName = [filename,'.kmz'];
ge_kml(kmlFileName,kmlStr);
ge_kmz(kmzFileName,'resourceURLs',{pngfile,kmlFileName});
           
% remove temporary kml file
system(['rm -f ' kmlFileName]);


function thegrid = roms_nget_grid(fnames,snames)
% thegrid = roms_nget_grid(gridfilelist,scoordfilelist)
%
% Like roms_get_grid but process a cell array of grid file names as would
% be defined for nested grid applications
%
% If scoordfilelist is not given the z-coordinates will not be computed
%
% John Wilkin - May 2018

get_z_coords = true;
if nargin < 2
  get_z_coords = false;
end

if iscell(fnames)
  nf = size(fnames,2);
else
  nf = 0;
end

if nf==0
  % treat as regular roms_get_grid call
  if get_z_coords
    thegrid = roms_get_grid(fnames,snames);
  else
    thegrid = roms_get_grid(fnames);
  end
else
  % preallocate outputs
  thegrid = cell([1 nf]);
  for i=1:nf
    if get_z_coords
      if iscell(snames)
        thegrid{i} = roms_get_grid(fnames{i},snames{i});
      else
        % assume snames was a 6-element s-coordinate vector
        thegrid{i} = roms_get_grid(fnames{i},snames);
      end
    else
      thegrid{i} = roms_get_grid(fnames{i});
    end
  end
end





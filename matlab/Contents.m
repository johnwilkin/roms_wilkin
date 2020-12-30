disp(' ')
disp('% This directory: ') 
disp(' ')
disp(pwd)
disp(' ')
echo on
% John Wilkin's m-files for working with the Regional Ocean Modeling 
% System (ROMS) (https://www.myroms.org)
%
% This collection of tools was originally created March 2008 (some are much
% older) for release to the ROMS User community.  Prior to December 2020 
% the code set was managed by svn and accessible with login credentials for
% the myroms.org portal ...
%    svn checkout https://www.myroms.org/svn/om/matlab/roms_wilkin
% 
% Presently the code is managed via github at: 
%    https://github.com/johnwilkin/roms_wilkin
%
% There is basic documentation for some of these files - especially for
% plotting - at https://romsmatlab.tiddlyspot.com/
%
% -------------------------------------------------------------------------
% A slightly incomplete listing of the Contents:
%
% For extracting and plotting ROMS model output:
%
% roms_get_grid.m    - load the coordinates (used by most functions)
% roms_get_ijgrid.m  - load coordinates in i,j space 
% roms_nget_grid.m   - for nested grids
% roms_zview.m       - plot a constant z slice
% roms_nzview.m      - nested grids
% roms_sview.m       - plot a constant s-level slice
% roms_nsview.m      - nested grids
% roms_iview.m	     - plot a constant i-index slice
% roms_jview.m       - plot a constant j-index slice
% roms_bview.m       - plot boundary conditions
% roms_bviews.m	     - stitch multiple boundaries together in one plot
% roms_genslice.m    - extract data at all depths along lon/lat/time path
% roms_zgenslice.m   - like genslice but 4-D extraction at specific depths
% roms_kmz.m         - old (might not work at all) 
%
% Supporting functions for the plots
%
% roms_islice.m
% roms_jslice.m
% roms_zslice.m
% roms_2dslice.m 
% roms_slice_var.m        - more versatile slicing
% roms_zslice_var.m       - z-slice of a 3-D ROMS variable in workspace
% roms_addvect.m          - add vectors (currents, winds, stresses, etc.)       
% roms_addvect_scale.m    - add vector scale
% roms_quivergrd.m        - plotting C-grid staggered vectors
% roms_curquivergrd.m     - "curvy" vectors (pseudo particle tracks)
% roms_get_date.m         - parse date/time (needs work)
% 
% Visualizing the grid
%
% roms_plot_bathy.m       roms_plot_coast.m          roms_plot_mesh.m
%
% Other utilities:
%
% roms_get_river_source_locations.m  	- interactively check river positions 
% roms_plot_river_source_locations.m  - plot point source locations 
%
% roms_uvrhotate.m        - rotate i,j vectors to east-north
% roms_cpplist.m          - parse CPP list from output file 
% roms_cgridpos.m         - determine C-grid position from the shape 
% roms_varlist.m          - cell array of common subsets of ROMS variables 
% roms_diff_grid.m        - get spatially variable viscosity
% roms_visc_grid.m        - get spatially variable diffusivity
% roms_rossbyradius.m     - compute baroclinic Rossby radius (slow)
% roms_lonlat2ij.m        - what it says
% roms_zint.m             - vertical integral between user defined layers         
%
% Other useful functions
%
% pcolorjw.m         - venerable function to plot centered pcolor patches
% amerc.m            - simple minded aspect ratio for lon/lat plots 
% callsfindstr.m     - reveal the path of all functions called 
% catstruct.m        - concatentate alike structures
% findstrinstruct.m  - find a particular string in a structure 
% closest.m          - find the i,j indicies for a lon/lat point
% rk4.m	             - Runge Kutta 4th order particle track integration
% strrep_.m          - escape _ to prevent interpretation as subscript
% griddata_lite.m    - old (don't use)
% parsetnc.m         - parse time units
%
% Reading data
%
% erddap_read.m                         - read an ERDDAP URL into workspace
% podaac_get_viirsL3U.m                 - SST from NASA PO-DAAC
% podaac_get_goes16L3C.m                - SST from NASA PO-DAAC
% roms_get_era5_NCARds633_bulkflux.m    - ERA5 from NCAR 
% roms_write_era5_NCARds633_frcfile.m   - write ERA5 forcing file
%
% For the doppio and pioneer ocean modeling projects
% 
% doppio_provenance_lookup.m
% pioneer_plot_gliders_nominal_tracks.m    - OOI Pioneer gliders
% pioneer_plot_mooring_sites.m             - OO Pioneer moorings
% pioneer_mplot_gliders_nominal_tracks.m   - with m_map toolbox
% pioneer_mplot_gliders_nominal_tracks.m   - with m_map toolbox
%
% Miscellaneous colormaps - I have many more
%
% rscolmap.m       - CSIRO ocean color
% coolavhrrmap.m   - RU COOL SST 
% zebra.m          
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id$

echo off

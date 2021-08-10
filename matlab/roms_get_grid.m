function grd = roms_get_grid(grd_file,scoord,tindex,~)
% grd = roms_get_grid(grd_file,outfile,zeta_input);
% grd = roms_get_grid(grd_file,scoord_info);
% grd = roms_get_grid(grd_file)
%
% Gets the lon,lat,mask,depth [and z coordinates] from roms netcdf
% grd_file or output file
%
% Input:
%     grd_file: The roms netcdf grid file name
%           or, an existing grd structure to which the vertical coordinates
%               are to be added or updated
%
% Optional inputs:
%     scoord_info: ROMS his/rst/avg output file from which the s-coord
%               params can be determined
%            or 4-element vector [theta_s theta_b Tcline N] where
%               Vtransform=1 and Vstretching=1 is assumed (compatibility)
%            or 6-element vector [theta_s theta_b Tcline N ...
%               Vtransform Vstretching] (recommended)
%
%     zeta_in:  How to obtain zeta information to use when including free
%                surface height in calculating z-coordinates
%            0 implies assume zeta=0
%            integer implies use this time index into his/rst/avg
%               file to read zeta
%               (In this case wet/dry masks will be read if they exist
%                and getpref('ROMS_WILKIN','USE_WETDRY_MASK') == true.
%            2-d array of zeta values
%
%     calc_zuv: If present, this argument (any value) activates computing
%               the depths z_u and z_v on the u and v points of the
%               ROMS C-grid
%
% Output is a structure containing all the grid information
%
% NOTE:     All arrays are stored in standard Matlab row-major order
% ******    (C-language) as they appear in netcdf ncdump and
%           NOT in column-major order (Fortran) as assumed by Arango's
%           ROMS Matlab utilities
%           For example,
%           grd.h  (lat,lon)          (j,i)   order
%           grd_z_r(N,lat,lon)        (k,j,i) order
%
% John Wilkin - 2000
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_get_grid.m 596 2020-12-29 16:46:14Z wilkin $

% This was confusing people so disable warning
warning('OFF','RomsGetGrid:NoVariable')

if exist('stretching','file')~=2
  disp('You need to add stretching.m to your matlab path')
  error('from e.g. https://www.myroms.org/svn/src/matlab/utility')
end
if exist('set_depth','file')~=2
  disp('You need to add set_depth.m to your matlab path')
  error('from e.g. https://www.myroms.org/svn/src/matlab/utility')
end

if isstruct(grd_file)
  % if the first input is already a grd structure
  % the intention is to add vertical coordinates below
  grd = grd_file;
  
else
  % get the horizontal grid information from a ROMS grid file (or
  % history, average etc file)
  grd.grd_file = grd_file;
  
  varlist = ...
    {'mask_rho','mask_psi','mask_u','mask_v','h','pm','pn','f','angle',...
    'dmde','dndx','visc_factor','diff_factor'};
  for v = varlist
    vname = char(v);
    try
      tmp = nc_varget(grd_file,vname);
      grd.(vname) = tmp; % grd = setfield(grd,vname,tmp);
    catch
      warning('RomsGetGrid:NoVariable',['Variable not found: ' vname])
      if strcmp(vname,'angle')
        grd.(vname) = zeros(size(grd.('h')));
      end
      if strcmp(vname,'h')
        warning('Using initial bath for h')
        grd.(vname) = squeeze(...
          nc_varget(grd_file,'bath',[0 0 0],[1 -1 -1]));
      end
    end
  end
  
  varlist = {'x_rho','y_rho','x_u','y_u','x_v','y_v','x_psi','y_psi'};
  if nc_isvar(grd_file,'x_rho')
    for v = varlist
      vname = char(v);
      try
        tmp = nc_varget(grd_file,vname);
        grd.(vname) = tmp; %replaces grd = setfield(grd,vname,tmp);
      catch
        warning('RomsGetGrid:NoVariable',['Variable not found: ' vname])
      end
    end
  end
  
  varlist = ...
    { 'lon_rho','lat_rho','lon_psi','lat_psi',...
    'lon_v','lat_v','lon_u','lat_u'};
  for v = varlist
    vname = char(v);
    try
      tmp = nc_varget(grd_file,vname);
      grd.(vname) = tmp; %replaces grd = setfield(grd,vname,tmp);
    catch
      %warning('RomsGetGrid:NoLonLat',...
      %  [vname ' not found. Substituting x/y coords instead'])
      if strcmp(vname(1:3),'lon')
        usevname = strrep(vname,'lon','x');
      else
        usevname = strrep(vname,'lat','y');
      end
      try
        grd.(vname) = grd.(usevname);
      catch
        % some files dont have psi coordinates
      end
      grd.nolatlon = true;
      grd.merc = false;
    end
  end
  try
    grd.bounding_box = [min(grd.lon_psi(:)) max(grd.lon_psi(:)) ...
      min(grd.lat_psi(:)) max(grd.lat_psi(:))];
  catch
  end
  try
    grd.perimeter(:,1) = ...
      [grd.lon_psi(1,:)' ; grd.lon_psi(:,end) ; ...
      grd.lon_psi(end,end:-1:1)'; grd.lon_psi(end:-1:1,1)];
    grd.perimeter(:,2) = ...
      [grd.lat_psi(1,:)' ; grd.lat_psi(:,end) ; ...
      grd.lat_psi(end,end:-1:1)'; grd.lat_psi(end:-1:1,1)];
  catch
  end
  try
    grd.corners = [grd.lon_rho(end,1) grd.lat_rho(end,1);
      grd.lon_rho(1,1) grd.lat_rho(1,1);
      grd.lon_rho(1,end) grd.lat_rho(1,end);
      grd.lon_rho(end,end) grd.lat_rho(end,end)];
  catch
  end
  
  varlist = {'rdrag','rdrag2','ZoBot'};
  for v = varlist
    vname = char(v);
    if nc_isvar(grd_file,vname)
      tmp = nc_varget(grd_file,vname);
      grd.(vname) = tmp;
    end
  end
  
  if isfield(grd,'mask_rho')
    grd.mask_rho_nan = grd.mask_rho;
    land = grd.mask_rho_nan==0;
    grd.mask_rho_nan(land) = NaN;
  else
    % if there is no mask information in the file so create unit masks
    % in case code tries to use them
    grd.mask_rho = ones(size(grd.h));
    grd.mask_rho_nan = grd.mask_rho;
    grd.mask_u = ones(size(grd.h(:,2:end)));
    grd.mask_v = ones(size(grd.h(2:end,:)));
    grd.mask_psi = ones(size(grd.h(2:end,2:end)));
    grd.nomask = 1;
  end
  
  % If the grid file includes coastline data, such as a file being used
  % with the Rutgers version of editmask.m, load this too
  try
    grd.lon_coast = ncread(grd_file,'lon_coast');
    grd.lat_coast = ncread(grd_file,'lat_coast');
  catch
  end
  
  if nargin > 2
    if length(tindex)==1 && tindex~=0
      % possibly get the appropriate land/sea or wet/dry mask for this time
      haswetdry = nc_isvar(grd_file,'wetdry_mask_rho');
      if haswetdry % wet dry mask exists in file
        try % override with preference
          usewetdry = getpref('ROMS_WILKIN','USE_WETDRY_MASK');
        catch % no preference - assume not
          usewetdry = false;
        end
        if usewetdry
          for v = {'wetdry_mask_rho','wetdry_mask_psi',...
              'wetdry_mask_u','wetdry_mask_v'}
            vname = char(v);
            tmp = nc_varget(grd_file,vname,[tindex-1 0 0],[1 -1 -1]);
            grd.(vname) = tmp; % grd = setfield(grd,vname,tmp);
          end
        end
      end
    end
  end
  
  % If the grid file is a refinement grid extract the info
  try grd.refine_factor = ncreadatt(grd_file,'/','refine_factor');
    varlist = {'parent_Imin','parent_Imax','parent_Jmin','parent_Jmax'};
    for v = varlist
      vname = char(v);
      grd.(vname) = ncreadatt(grd_file,'/',vname);
    end
  catch
  end
  
end % loading 2D fields independent of s-coordinate

if nargin > 1
  
  h = grd.h;
  [Mp,Lp] = size(h);
  L = Lp-1;
  M = Mp-1;
  
  % get z_r and z_w for the given s-coordinate parameters
  
  [theta_s,theta_b,Tcline,N,Vtransform,Vstretching] = ...
    roms_get_scoord(scoord);
  
  % This logic shifted to an embedded function at the end of this file
  %
  %   if ~ischar(scoord)
  %
  %     theta_s = scoord(1);
  %     theta_b = scoord(2);
  %     Tcline  = scoord(3);
  %     N       = scoord(4);
  %     if (length(scoord) < 5)
  %       Vtransform = 1;
  %       Vstretching = 1;
  %     else
  %       Vtransform = scoord(5);
  %       Vstretching = scoord(6);
  %     end
  %     hc = Tcline;
  %
  %   else
  %
  %     % input 'scoord' is a his/avg/rst file name or opendap url
  %     % attempt to get s-coord params from this file/url
  %
  %     theta_s = ncread(scoord,'theta_s');
  %     theta_b = ncread(scoord,'theta_b');
  %     hc      = ncread(scoord,'hc');
  %     try
  %       Tcline  = ncread(scoord,'Tcline');
  %     catch
  %       Tcline = hc;
  %     end
  %     N = length(ncread(scoord,'Cs_r'));
  %     if nc_isvar(scoord,'Vtransform')
  %       Vtransform = ncread(scoord,'Vtransform');
  %     else
  %       Vtransform = 1;
  %     end
  %     if nc_isvar(scoord,'Vstretching')
  %       Vstretching = ncread(scoord,'Vstretching');
  %     else
  %       Vstretching = 1;
  %     end
  %
  %   end
  
  hc = Tcline;
  if Vtransform == 1
    hc = min(Tcline,min(h(:)));
  end

  [s_rho,Cs_r] = stretching(Vstretching,theta_s,theta_b,hc,N,0,0);
  [s_w  ,Cs_w] = stretching(Vstretching,theta_s,theta_b,hc,N,1,0);
  sc_r = s_rho;
  sc_w = s_w;
  
  % zeta
  
  zeta = zeros(size(grd.h)); % default
  if nargin > 2 % option to include zeta in z calculation
    
    if tindex == 0 % if tindex==0 zeta defaults to zero
      % do nothing
      
    else % if tindex==0 zeta defaults to zero
      if length(tindex)==1
        % tindex is a single index to zeta in a roms output file
        if ischar(scoord)
          zeta = nc_varget(scoord,'zeta',[tindex-1 0 0],[1 -1 -1]);
          zeta = squeeze(zeta);
        else
          warning([ 'Cannot process zeta from file(2) in the case ' ...
            ' that scoord parameters are input as a vector or structure'])
          disp(['Reading zeta from ' grd_file ' for record ' int2str(tindex)])
          zeta = nc_varget(grd_file,'zeta',[tindex-1 0 0],[1 -1 -1]);
        end
        
        if isempty(zeta)
          warning([ 'zeta not found in ' scoord '. Assuming zeta=0.'])
          zeta = zeros(size(grd.h));
        end
      else
        % tindex should be a 2-d field of zeta values
        if any(size(tindex)-size(grd.h))
          % sizes of zeta and h don't match
          error('input tindex as zeta does not match grid dimensions')
        else
          zeta = tindex;
        end
      end
    end
    
  end
  grd.zeta = zeta;
  
  % rho-points
  h = grd.h;
  rgrid = 1;
  z_r = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, ...
    rgrid, h', zeta', 0);
  grd.z_r = permute(z_r,[3 2 1]);
  clear z_r
  
  % w-points
  wgrid = 5;
  z_w = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, ...
    wgrid, h', zeta', 0);
  grd.z_w = permute(z_w,[3 2 1]);
  
  % rho cell thicknesses
  grd.dz = diff(grd.z_w,1,1);
  
  % cell areas and volumes
  try
    grd.dA = 1./(grd.pm.*grd.pn);
    dA(1,:,:) = grd.dA;
    dV = bsxfun(@times,dA,grd.dz);
  catch
  end
  
  % compute the z depths on the velocity points as well
  % this used to be (before 2009/10/08) optional but there is little
  % reason not to do this and it's annoying when you've forgotten to
  
  % u-points (cell centres in vertical)
  ugrid = 3;
  z_u = set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,N, ...
    ugrid,h',zeta',0);
  grd.z_u = permute(z_u,[3 2 1]);
  
  % u-points (cell edges in vertical)
  z_uw = 0.5.*(z_w(1:L,1:Mp,:)+z_w(2:Lp,1:Mp,:));
  grd.z_uw = permute(z_uw,[3 2 1]);
  
  % v-points (cell centres in vertical)
  vgrid = 4;
  z_v = set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,N, ...
    vgrid,h',zeta',0);
  grd.z_v = permute(z_v,[3 2 1]);
  
  % v-points (cell edges in vertical)
  z_vw= 0.5.*(z_w(1:Lp,1:M,:)+z_w(1:Lp,2:Mp,:));
  grd.z_vw = permute(z_vw,[3 2 1]);
  
  clear z_u z_uw z_v z_vw
  
  grd.dV = dV;
  grd.Vtransform = Vtransform;
  grd.Vstretching = Vstretching;
  grd.theta_s = theta_s;
  grd.theta_b = theta_b;
  grd.Tcline = Tcline;
  grd.N = N;
  grd.hc = hc;
  grd.sc_w = sc_w;
  grd.Cs_w = Cs_w;
  grd.sc_r = sc_r;
  grd.Cs_r = Cs_r;
  grd.s_w = s_w;
  grd.s_rho = s_rho;
  
end

warning('OFF','RomsGetGrid:NoVariable')

function [theta_s,theta_b,Tcline,N,Vtransform,Vstretching] = ...
  roms_get_scoord(scoord)
% Parse the s-coordinate parameters

if ischar(scoord)
  
  % scoord is a his/avg/rst file name or opendap url
  
  theta_s = ncread(scoord,'theta_s');
  theta_b = ncread(scoord,'theta_b');
  hc = ncread(scoord,'hc');
  try
    Tcline  = ncread(scoord,'Tcline');
  catch
    Tcline = hc;
  end
  N = length(ncread(scoord,'Cs_r'));
  if nc_isvar(scoord,'Vtransform')
    Vtransform = ncread(scoord,'Vtransform');
  else
    Vtransform = 1;
  end
  if nc_isvar(scoord,'Vstretching')
    Vstretching = ncread(scoord,'Vstretching');
  else
    Vstretching = 1;
  end

elseif isstruct(scoord)
  
  % scoord is a previously loaded grd structure
  theta_s = scoord.theta_s;
  theta_b = scoord.theta_b;
  Tcline = scoord.Tcline;
  N = scoord.N;
  Vtransform = scoord.Vtransform;
  Vstretching = scoord.Vstretching;
  
else
  
  % scoord was a vector of parameters
  theta_s = scoord(1);
  theta_b = scoord(2);
  Tcline  = scoord(3);
  N       = scoord(4);
  if (length(scoord) < 5)
    Vtransform = 1;
    Vstretching = 1;
  else
    Vtransform = scoord(5);
    Vstretching = scoord(6);
  end
  % hc = Tcline;

end
  

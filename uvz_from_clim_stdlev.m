function clm = uvz_from_clim_stdlev(clm,pref,hmean)
% $Id: uvz_from_clim_stdlev.m 358 2008-04-07 14:15:03Z zhang $
% Compute u,v geostrophic from a standard level t,s,ssh climatology
%
% clm = uvz_from_clim_stdlev(clm,pref)
%
% Inputs:
%
%    clm = structure of 3D gridded (stdlev) temperature and salinity
%          clm.temp [time z lat lon]
%          clm.salt [time z lat lon]
%          clm.lon
%          clm.lat
%          clm.z [z(1) is assumed to be the deepest level]
%       optional:
%          clm.zeta [time lat lon]
%
%    pref = reference pressure (db) to compute dynamic height from
%    hmean = constant value subtracted from zeta before output (e.g. for
%          removing a large scale mean in dynamic height)
%
% Outputs:
%
%    Velocity components added to the climatology structure are 
%    defined on a Arakawa C grid
%          clm.u
%          clm.v
%          clm.h*** is the dynamic height (where *** denotes the 
%              reference pressure used)
%          clm.zeta (if it wasn't provided on input) is the surface 
%              height (h***) with the mean removed
%
% John Wilkin Dec 2000

% separate out the data -----------------------------------------------

lon  = clm.lon;
lat  = clm.lat;
z    = clm.z;

if ndims(clm.temp)<4
  ntimes = 1;
else
  ntimes = size(clm.temp,1);
end

tclm(1:ntimes,:,:,:) = clm.temp;
sclm(1:ntimes,:,:,:) = clm.salt;

if nargin < 2
  try
    hclm(1:ntimes,:,:) = clm.zeta;
    warning('Found clm.zeta input so using this to reference the flow')
  catch
    warning([lasterr])
    disp(' No pref requested or zeta variable in climatology')
    disp('  Using zeta = 0')
    hclm = zeros([ntimes size(tclm,3) size(tclm,4)]);
  end
end

% ---------------------------------------------------------------------

% allocate the u,v arrays (the geostrophic calculation gets v at 
% u-points and u at v-points of a C grid)

uclm = zeros([size(tclm)+[0 0 -1 0]]); 
vclm = zeros([size(tclm)+[0 0 0 -1]]);
GA = zeros(size(tclm));

disp('geopotential anomaly') % ----------------------------------------

for mon = 1:ntimes

  disp([' Doing time index ' int2str(mon)])

  for J = 1:size(vclm,3)

    % disp(['  Doing J index ' int2str(J)])
    
    x = lon(J,:);
    temp = squeeze(tclm(mon,:,J,:));
    temp = flipud(temp);
    salt = squeeze(sclm(mon,:,J,:));
    salt = flipud(salt);
    pres = flipud(abs(z));

    % geopotential anomaly
    ga = sw_gpan(salt,temp,pres);
    GA(mon,:,J,:) = ga;
    
  end

  if nargin > 1
    
    % compute dynamic height relative to pref    
    dha = GA/9.81; % to get dynamic metres divide by g
  
    % relative to input pref
    l_ref = min(find(pres >= abs(pref)));
    hclm(mon,:,:) = squeeze(dha(mon,l_ref,:,:)-dha(mon,1,:,:));

  end
  
end

disp('v-component of velocity') % -------------------------------------

for mon = 1:ntimes

  disp([' Doing time index ' int2str(mon)])

  for J = 1:size(vclm,3)

    % disp(['  Doing J index ' int2str(J)])

    x = lon(J,:);
    zeta = squeeze(hclm(mon,J,:));

    % geopotential anomaly
    ga = squeeze(GA(mon,:,J,:));
    
    % velocity relative to surface
    velz = sw_gvel(ga,lat(J,:),x);

    % g*zeta for surface
    ga0 = repmat(-9.81*zeta',length(z),1);

    % 'barotropic' velocity
    vel0 = sw_gvel(ga0,lat(J,:),x);

    % total velocity
    vel = vel0 + velz;

    % flip back to fill clm array
    vel = flipud(vel);

    % fill this slice
    vclm(mon,:,J,:) = vel;

  end

end

disp('u-component of velocity') % -------------------------------------

for mon = 1:size(uclm,1)

  disp([' Doing time index ' int2str(mon)])

  for I = 1:size(uclm,4)

    % disp(['  Doing I index ' int2str(I)])
    
    y = lat(:,I);
    zeta = squeeze(hclm(mon,:,I));

    % geopotential anomaly
    ga = squeeze(GA(mon,:,:,I));
     
    % velocity relative to surface
    velz = -sw_gvel(ga,y',lon(:,I));

    % g*zeta for surface
    ga0 = repmat(-9.81*zeta,length(z),1);

    % 'barotropic' velocity
    vel0 = -sw_gvel(ga0,y',lon(:,I));

    % total velocity
    vel = vel0 + velz;
	
    % flip back to fill clm array
    vel = flipud(vel);

    % fill this slice
    uclm(mon,:,:,I) = vel;

  end
end

ArakawaGrid = 'A';

switch ArakawaGrid

  case 'A'
    
    disp('averaging u to Arakawa-A u-points') % ---------------------------
    % and pad extra rows/columns

    M = size(uclm,3);
    Mm = M-1;    
    tmp = 0.5*(uclm(:,:,1:Mm,:)+uclm(:,:,2:M,:));
    uclm = cat(3,tmp(:,:,1,:),tmp,tmp(:,:,Mm,:));
    
    disp('averaging v to Arakawa-A v-points') % ---------------------------
    % and pad extra rows/columns
    
    L = size(vclm,4);
    Lm = L-1;    
    tmp = 0.5*(vclm(:,:,:,1:Lm)+vclm(:,:,:,2:L));
    vclm = cat(4,tmp(:,:,:,1),tmp,tmp(:,:,:,Lm));
    
  case 'C'

    disp('averaging u to Arakawa-C u-points') % ---------------------------
    % and pad extra rows/columns

    M = size(uclm,3);
    Mm = M-1;
    L = size(uclm,4);
    Lm = L-1;
    tmp = 0.25*(uclm(:,:,1:Mm,1:Lm)+uclm(:,:,1:Mm,2:L)+...
	uclm(:,:,2:M,1:Lm)+uclm(:,:,2:M,2:L));
    uclm = cat(3,tmp(:,:,1,:),tmp,tmp(:,:,Mm,:));
    
    disp('averaging v to Arakawa-C v-points') % ---------------------------
    % and pad extra rows/columns
    
    M = size(vclm,3);
    Mm = M-1;
    L = size(vclm,4);
    Lm = L-1;
    tmp = 0.25*(vclm(:,:,1:Mm,1:Lm)+vclm(:,:,1:Mm,2:L)+...
	vclm(:,:,2:M,1:Lm)+vclm(:,:,2:M,2:L));
    vclm = cat(4,tmp(:,:,:,1),tmp,tmp(:,:,:,Lm));
    
end

    
% append results to input clm structure % ---------------------------------

if ntimes == 1
  uclm = squeeze(uclm);
  vclm = squeeze(vclm);
  hclm = squeeze(hclm);
end
  
clm.u = uclm;
clm.v = vclm;
if nargin > 1
  eval(['clm.h' int2str(abs(pref)) ' = hclm;'])
end

if nargin < 3
  clm.zeta = hclm-nanmean(hclm(:));  
else  
  clm.zeta = hclm-hmean;
end


function [R,S] = roms_rossbyradius(T,S,g)
% [Ro,Uo] = roms_rossbyradius(temp,salt,grd)
% Compute 1st baroclinic Rossby radius (km) and wave speed (m/s) from ROMS
% 3-D temperature and salt fields and the GRD structure from roms_get_grid
%
% This is very slow because it has to step through each vertical profile
% and solve the eigenvalue problem for first baroclinic Rossby mode.
%
% The work has been placed in a parfor loop to accelerate the calculation
% if a pool of parallel workers is active (e.g. parpool(4))
%
% John Wilkin - June 2020
% $Id: roms_rossbyradius.m 587 2020-10-16 14:32:09Z wilkin $

p = gcp('nocreate');
if isempty(p)
  disp('This function is very slow on large data. It will run faster ',...
    'you create a parallel pool with the parpool command')
end

P = -g.z_r;
bfrq = sw_bfrq(S(:,:),T(:,:),P(:,:));
pres = -g.z_w(2:end-1,:,:);
bfrq = reshape(bfrq,size(pres));

% preallocate for speed
S = nan(size(g.lat_rho));
R = S;

% parallel loop
parfor k=1:numel(R)
  if g.mask_rho(k)==1
    [speed,radius] = rossby_modes(bfrq(:,k),pres(:,k),g.lat_rho(k),1,0);
    % sometimes there are complex values in shallow water with ill-defined
    % stratification
    S(k) = real(speed);
    R(k) = real(radius);
  end
end

% convert to km
R = R/1000;



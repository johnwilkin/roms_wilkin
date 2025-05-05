function [Radius,Speed] = roms_rossbyradius(T,S,g)
% [Ro,Uo] = roms_rossbyradius(temp,salt,grd)
% Compute 1st baroclinic Rossby radius (km) and wave speed (m/s) from ROMS
% 3-D temperature and salt fields and the GRD structure from roms_get_grid
%
% This is very slow because it has to step through each vertical profile.
%
% The approximate Chelton method is about twice as fast as the exact
% eigenvalue approach. They give similar results. 
%
% In previous versions this code could be acdelerated using a parfor loop 
% but presently parpool on Matlab 2023b and 2024a on Apple Silicon M1 is
% throwing memory probems so I have diabled the code. 
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_rossbyradius.m 587 2020-10-16 14:32:09Z wilkin $

p = gcp('nocreate');
if isempty(p)
  disp("This function is slow on large data. It will run faster "...
    +"if you create a parallel pool with the parpool command")
end

P = -g.z_r; % pressure at t,s observation depths
PW = -g.z_w(2:end-1,:,:); % pressure at depth where BFRQ will be defined

% preallocate for speed
Speed = nan(size(g.lat_rho));
Radius = Speed;

% only compute on water points
wet = find(g.mask_rho==1);

prog = 0;
waitbar(prog);

% Mehtod to use computing Rossby wave speed
opt = 'Chelton';

% parfor k=1:numel(g.mask_rho)
for k = wet' % only wet points

% if g.mask_rho(k)~=0 % Must do this with parfor
  
  % S, ptemp and pressure
  sp = S(:,k);
  tp = T(:,k);
  pr = P(:,k);  % rho points
  pw = PW(:,k); % w points

  % sw_bfrq assumes potential temperature input, which is what ROMS outputs
  bfrq = sw_bfrq(sp,tp,pr);
  bfrq = max(eps,bfrq); % eliminate negative NSQ due to unstable profiles

  switch lower(opt)
    case 'eigenvalue'
      [speed,radius] = rossby_modes(bfrq,pw,g.lat_rho(k),1);
    case 'chelton'
      % Chelton method:
      %      c_m = 1/(m*pi) integral_-h^0 (Nz) dz
      % Vertical integral is inner product of NSQ times DZ
      speed = 1/pi*transpose(sqrt(bfrq))*diff(-pr);
      radius = speed/g.f(k);
  end
  Speed(k) = real(speed);
  Radius(k) = real(radius);
  prog = k/wet(end);
  waitbar(prog)

% end % test of mask_rho

end

% convert to km
Radius = Radius/1000;

function [speed,radius,modes] = rossby_modes(bfrq,p,lat,keepmodes)
% [speed,radius,shape] = rossby_modes(bfrq,p,lat,[keepmodes],[plot])
%
% brfq = Brunt-Vaisala frequency vector such as computed by sw_bfrq at
% p    = depth values (<0) or pressures (dbar)
% lat  = latitude
% keepmodes (optional, default = 2)
%
% speed = mode speed in (m/s)
% radius = baroclinic Rossby radius (m)
% modes = matrix of mode structures of normalized vertical displacement
%         but with the first column being the depth coordinate (z)
%
% Calculates phase speeds, Rossby radii, and modal structures (of vertical
% displacement) for a stratification profile specified by Brunt-Vaisala 
% frequency profile BRFZ(P) under the Bousinnesq approximation, i.e. solves
% the Sturm-Liouville eigenvalue problem equation (6.11.18) of Gill (1982)
%
% Only the first  "keepmodes" modes are kept.
%
% John Wilkin, CSIRO Division of Oceanography, Australia
% wilkin@flood.ml.csiro.au
% April 24, 1991
%
% Converted to a function for a tabular dens(z) JLW 2019

if nargin < 4
  keepmodes = 2;
end

H = max(abs(p));  % allow for coordinate P being pressure (db)
f = sw_f(lat);
n    = 100;       % number of grid points in z in discretization
z = linspace(-abs(p(1)),-abs(p(end)),n-1);
delz = diff(z(1:2));

% interp tabulated values
% allow for coordinate P being pressure (db) - convert to z < 0 (m)
nsq = interp1(-abs(p),bfrq,z);

% simple 2nd-order finite differences
coef   = -1/delz^2  ./ (nsq+eps);
A      = diag(-2*ones(n-1,1),0)+diag(ones(n-2,1),1)+diag(ones(n-2,1),-1);
B      = diag(coef) * A;
[V,D]  = eig(B);
eigval = diag(D);
c      = 1.0./sqrt(eigval);
[Y,I]  = sort(-c);

% phase speed of each mode
speed  = -Y(1:keepmodes);

% Rossby radius of deformation of each mode
radius = speed(1:keepmodes)/f;

% vertical displacement profile each mode
shape  = V(:,I(1:keepmodes));

% add in points at z=-H,0 to make the plots complete
zplt   = [ -H ; z'; 0];
shape  = [zeros(1,keepmodes) ; shape; zeros(1,keepmodes) ];

% normalize the modes
shape = diag( 1.0 ./ max(abs(shape)))*shape';

modes = [zplt(:) shape'];

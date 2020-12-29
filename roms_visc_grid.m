function [visc2,uRe,vRe] = roms_visc_grid(g,uvnu2)
% [visc2,uRe,vRe] = roms_diff_grid(grd,UVNU2)
%
%    Inputs:   grd is a ROMS grid structure (see roms_get_grid)
%            uvn2 is the Laplacian viscosity coefficient given in ocean.in
%
% Calculate the spatially varying diffusivity when ROMS option VISC_GRID is
% used to scale the parameter value for UVNU2
%
% Outputs are on the psi points grid.
%
% Optionally returns the grid Reynolds number in each grid direction
% assuming U = 1 m/s
%
% From ini_hmixcoef.F (see metrics.F for grdscl and grdmax)
%
% #if defined VISC_GRID && defined SOLVE3D
% !
% !-----------------------------------------------------------------------
% !  Scale horizontal diffusion according to the grid size.
% !-----------------------------------------------------------------------
% !
% # ifdef UV_VIS2
%       cff=viscosity2/grdmax(ng)
%       DO j=JstrT,JendT
%         DO i=IstrT,IendT
%           visc2_r(i,j)=cff*grdscl(i,j)
%         END DO
%       END DO
%       cff=0.25_r8*cff
%       DO j=JstrP,JendT
%         DO i=IstrP,IendT
%           visc2_p(i,j)=cff*(grdscl(i-1,j-1)+grdscl(i,j-1)+              &
%      &                      grdscl(i-1,j  )+grdscl(i,j  ))
%         END DO
%       END DO
% # endif
%
% John Wilkin - during COVID-19 lockdown
%
% Copyright (c) 2020 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_visc_grid.m 596 2020-12-29 16:46:14Z wilkin $

dA = 1./(g.pm.*g.pn);
grdscl = sqrt(dA);
grdscl = 0.5*(grdscl(1:end-1,:)+grdscl(2:end,:));
grdscl = 0.5*(grdscl(:,1:end-1)+grdscl(:,2:end));
grdmax = max(grdscl(g.mask_psi==1));
visc2  = uvnu2*grdscl/grdmax;

U = 1;
if nargin > 1
  tmp = 0.5*(g.pm(1:end-1,:)+g.pm(2:end,:));
  pm = 0.5*(tmp(:,1:end-1)+tmp(:,2:end));
  uRe = U./pm./visc2;
end
if nargin > 1
  tmp = 0.5*(g.pn(1:end-1,:)+g.pn(2:end,:));
  pn = 0.5*(tmp(:,1:end-1)+tmp(:,2:end));
  vRe = U./pn./visc2;
end

function [diff2,uPec,vPec] = roms_diff_grid(g,tnu2)
% [diff2,Pecletu,Pecletv] = roms_diff_grid(grd,TNU2)
%
%    Inputs:   grd is a ROMS grid structure (see roms_get_grid)
%             tnu2 is the Laplacian diffusion coefficient given in ocean.in
% 
% Calculate the spatially varying diffusivity when ROMS option DIFF_GRID is
% used to scale the parameter value for TNU2
%
% Optionally returns the grid Peclet number in each grid direction
% assuming U = 1 m/s
% 
% From ini_hmixcoef.F (see metrics.F for grdscl and grdmax)
%
% #if defined DIFF_GRID && defined SOLVE3D
% !
% !-----------------------------------------------------------------------
% !  Scale horizontal diffusion according to the grid size.
% !-----------------------------------------------------------------------
% !
% # ifdef TS_DIF2
%       DO itrc=1,NT(ng)
%         cff=diffusion2(itrc)/grdmax(ng)
%         DO j=JstrT,JendT
%           DO i=IstrT,IendT
%             diff2(i,j,itrc)=cff*grdscl(i,j)
%           END DO
%         END DO
%       END DO
% # endif
%
% John Wilkin - January 2017
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_diff_grid.m 598 2020-12-29 16:50:22Z wilkin $

dA = 1./(g.pm.*g.pn);
grdscl = sqrt(dA);
grdmax = max(grdscl(g.mask_rho==1));
diff2  = tnu2*grdscl/grdmax;

U = 1;
if nargin > 1
    uPec = U./g.pm./diff2;
end
if nargin > 1
    vPec = U./g.pn./diff2;
end

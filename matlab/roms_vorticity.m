function out = roms_vorticity(u,v,g,vortype)
% vor = roms_vorticity(u,v,g,vortype)
% compute relative vorticity
%
% Options:
%   vortype = 'relative' 
%   vortype = 'okubo-weiss'
%   vortype = 'potential' is not supported at present
%
% Output is on psi points grid
%
% -------------------------------------------------------------------------
% !  From ROMS vorticity.F 
% !  This routine computes relative (s-1) and  potential (m-1 s-1)       !
% !  vorticity for an adiabatic Boussinesq fluid where the potential     !
% !  density is conserved:                                               !
% !                                                                      !
% !    pvor = 1/rho0 dot_product(avor, grad(pden))                       !
% !                                                                      !
% !  where "avor" is the absolute (relative plus planetary) vorticity    !
% !  and "pden" is the potential density (a conserved quantity).         !
% !                                                                      !
% !    avor = rvor + f                                                   !
% !                                                                      !
% !  In curvilinear coordinates, the vertical component of relative      !
% !  vorticity and potential vorticity are:                              !
% !                                                                      !
% !                                                                      !
% !    rvor = mn * [d(v/n)/d(xi) - d(u/m)/d(eta)]                        !
% !                                                                      !
% !    pvor = mn/rho0 * [f/mn +                                          !
% !                      d(v/n)/d(xi) - d(u/m)/d(eta)] * d(pden)/d(z) +  !
% !           1/rho0 * [1/n d(pden)/d(eta) d(u)/d(z) -                   !
% !                     1/m d(pden)/d(xi)  d(v)/d(z)]                    !
% !                                                                      !
% !  In addition, the vertically integrated (shallow water) relative     !
% !  and potential vorticity are computed.                               !
% !                                                                      !
% !  The relative and potential vorticity is discretized at horizontal   !
% !  PSI-points and vertical RHO-points.                                 !
% !                                                                      !
% !=======================================================================

if nargin < 4
  vortype = 'relative';
end

switch vortype(1)
  
  case 'r'
    
    % relative vorticity
    von = 2*v./(g.pn(1:(end-1),:)+g.pn(2:end,:));
    uom = 2*u./(g.pm(:,1:(end-1))+g.pm(:,2:end));
    mn_p = 0.0625*...
      (g.pm(1:(end-1),1:(end-1))+g.pm(1:(end-1),2:end)+...
       g.pm(2:end,1:(end-1))+g.pm(2:end,2:end)).*...
      (g.pn(1:(end-1),1:(end-1))+g.pn(1:(end-1),2:end)+...
       g.pn(2:end,1:(end-1))+g.pn(2:end,2:end));
    xi_p = mn_p.*(von(:,2:end)-von(:,1:(end-1))-uom(2:end,:)+uom(1:(end-1),:));
    out = xi_p;
    
  case 'o'
    % Okubo-Weiss parameter
    % positive is dominated by irrotational, straining deformation
    % negative is dominated by vorticity and rotation
    % (Samelson, 2013, Annu. Rev. Mar. Sci., 5:137?63)
    
    % 1/2(dv/dx-du/dy)
    von = 2*v./(g.pn(1:(end-1),:)+g.pn(2:end,:));
    uom = 2*u./(g.pm(:,1:(end-1))+g.pm(:,2:end));
    xi = von(:,2:end)-von(:,1:(end-1))-uom(2:end,:)+uom(1:(end-1),:);
    % average to rho points
    xi = av2(av2(xi')');
    % pad to correct dimension with NaNs at edges
    xi = 0.5*g.pm.*g.pn.*xi([1 1:end end],[1 1:end end]);
    
    % 1/2(du/dy+dvdx)
    gam = von(:,2:end)-von(:,1:(end-1))+uom(2:end,:)-uom(1:(end-1),:);
    % average to rho points
    gam = av2(av2(gam')');
    % pad to correct dimension with NaNs at edges
    gam = 0.5*g.pm.*g.pn.*gam([1 1:end end],[1 1:end end]);

    % delta terms
    vom = 2*v./(g.pm(1:(end-1),:)+g.pm(2:end,:));
    uon = 2*u./(g.pn(:,1:(end-1))+g.pn(:,2:end));
  
    % delta1 du/dx
    d1 = uon(:,2:end)-uon(:,1:end-1);
    d1 = g.pm.*g.pn.*d1(:,[1 1:end end]);
    
    % delta2 dv/dy
    d2 = vom(2:end,:)-vom(1:end-1,:);
    d2 = g.pm.*g.pn.*d2([1 1:end end],:);
  
    out = 0.25*(d1-d2).^2 + gam.^2 - xi.^2;
    
  case 'p'
    
    error('this function does not compute potential vorticity')
    
end

function [u,v] = uv2rho(u,v)

% average vector components to rho points
u = av2(u')';
v = av2(v);

% pad to correct dimension with NaNs at edges
u = u(:,[1 1:end end]);
u(:,[1 end]) = NaN;
v = v([1 1:end end],:);
v([1 end],:) = NaN;

function [xn,yn] = rk4(x,y,t,u,v,x0,y0,t0,tn,xmax,ymax);
% [xn,yn] = rk4(x,y,t,u,v,x0,y0,t0,tn,xmax,ymax)
% 4th ord Runge-Kutta integration of velocity to track Lagrangian particles
%
% Inputs:
%     x,y are 2-D arrays of horizontal position 
%     t is 1-D vector of times
%     u,v are 3-D arrays of velocity component
%     x0,y0 are 1-D vectors of initial position of the particles at time t0
%     xn,yn are 1-D vectors of final position of the particles at time tn
%
% The inputs must be defined such that d(x)/dt = u etc.  
%
% Curvilinear (e.g. lon/lat) grids can be handled by tracking the non-
% dimensional coordinates of particles and using the product of velocity and
% grid metrics. For example, if x,y are matrices of the non-dimensional
% horizontal positions of the grid points (i.e. grid indices), and dx,dy 
% (in metres) are matrices of the grid intervals on the x,y points, then 
% use u./dx and v./dy as the velocity data (in m/s):
%                [xn,yn] = rk4(x,y,t,u./dx,v./dy,x0,y0,t0,tn);
%
% INTERP2 can be used to switch between the x,y and lon/lat space.
% Lon0 = interp2(1:size(lon,2),1:size(lon,1),lon,x,y); % i.e. lon(y,x)
% Lat0 = interp2(1:size(lat,2),1:size(lat,1),lat,x,y); % i.e. lat(y,x)
% [X,Y] = meshgrid(1:size(lon,2),1:size(lon,1));
% x = interp2(lon,lat,X,Lon0,Lat0); % i.e. x(lon,lat) 
% y = interp2(lon,lat,Y,Lon0,Lat0); % i.e. y(lon,lat) 
%
% THIS ALGORITHM IS *WAY MORE* ACCURATE THAN SIMPLE FIRST-ORDER
% DIFFERENCES FOR PREDICTING PARTICLE TRAJECTORIES. 1ST-ORDER SCHEMES WILL
% HAVE SUBSTANTIAL ERRORS IN REGIONS Of STRONG FLOW CURVATURE.
%
% John Wilkin - 2006
% $Id: rk4.m 544 2019-11-23 17:08:13Z robertson $

x3d = x(:,:,ones([1 length(t)]));
y3d = y(:,:,ones([1 length(t)]));
t3d = reshape(t(ones(size(x)),:),size(x3d));

% dt is time step
dt = tn-t0;
t0 = t0(:,ones(size(x0)));
tn = tn(:,ones(size(x0)));

% 4-th order runge-kutta integration

k1x = dt*interp3(x3d,y3d,t3d,u,x0,y0,t0);
k1y = dt*interp3(x3d,y3d,t3d,v,x0,y0,t0);

x1 = x0+0.5*k1x;
y1 = y0+0.5*k1y;

if nargin > 9 
  x1 = max(1,min(x1,xmax));
  y1 = max(1,min(y1,ymax));
end
   
k2x = dt*interp3(x3d,y3d,t3d,u,x1,y1,t0+0.5*dt);
k2y = dt*interp3(x3d,y3d,t3d,v,x1,y1,t0+0.5*dt);

x2 = x0+0.5*k2x;
y2 = y0+0.5*k2y;

if nargin > 9 
  x2 = max(1,min(x2,xmax));
  y2 = max(1,min(y2,ymax));
end

k3x = dt*interp3(x3d,y3d,t3d,u,x2,y2,t0+0.5*dt);
k3y = dt*interp3(x3d,y3d,t3d,v,x2,y2,t0+0.5*dt);

x3 = x0+k3x;
y3 = y0+k3y;

if nargin > 9 
  x3 = max(1,min(x3,xmax));
  y3 = max(1,min(y3,ymax));
end

k4x = dt*interp3(x3d,y3d,t3d,u,x3,y3,tn);
k4y = dt*interp3(x3d,y3d,t3d,v,x3,y3,tn);

w = [1/6 1/3 1/3 1/6];
xn = x0 + w*[k1x; k2x; k3x; k4x];
yn = y0 + w*[k1y; k2y; k3y; k4y];

if nargin > 9 
  xn = max(1,min(xn,xmax));
  yn = max(1,min(yn,ymax));
end

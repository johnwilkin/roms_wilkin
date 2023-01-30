function [vel,coord] = roms_velslice(file,uname,vname,varargin)
% [vel,coord] = roms_velslice(file,uname,vname,lonTrk,latTrk,timTrk,opt)
%
% Get vertical slice of ROMS vector data at coordinates specified by inputs
% of lon, lat and time (e.g. glider track, estuary talweg, isobath). Data
% is returned on the ROMS native s-coordinate. The vector outputs are 
% rotated to along- and across-track components. 
% See below for the definition of the coordinate rotation convention. 
%
% This function works by making two calls to roms_genslice for the two
% vector components, and then reconciling the along/across-track coordinate
% projection in fractional i,j coordinates.
%
% Input:
%   file      netcdf or opendap url (aggregation across time)
%   uname     i-direction vector component variable name, e.g. 'u'
%   vname     j-direction vector component variable name, e.g. 'v'
%             (Other names might be e.g. the diagnostics of tracer fluxes 
%             on cell faces.)
%   lon, lat  coordinates along the track
%   timTrk    time values along the track
%             - default: assume Matlab datenum 
%             - can interpret as index with option 'index' (see below)
%             - must match dimension of lon/lat OR be a scalar, in
%               which case it is assumed the track is for a single time
%
% Optional inputs:
%   options   integer Ntrk - interpolate to this many points to refine
%               resolution of output along trajectory. Default is return
%               values strictly at the requested lon,lat coordinates.
%             string 'nearest' or 'natural' to change from linear interp
%             string 'date' to interpret timTrk as datenum (default)
%             string 'index' interpret timTrk as index into time
%             string 'verbose' to monitor progress
% Output:
%   vel (complex) interpolated vector data along the track
%               real(vel) along-track (>0 if in direction of travel)
%               imag(vel) cross-track (>0 to left of direction of travel)
%   coord     structure of data coordinates suited to plotting with pcolor
%               lon  2D lon along the track
%               lat  2D lat along the track
%               z    2D depth along the track in vertical stretched coords
%               zw   corresponding depths of the layer interfaces
%                    (if the requested variable is a 'w' points variable
%                    then z and zw are the same)
%               dz   thickness of centered layers (so vertical sum of dz
%                    is (h+zeta)
%               dis  2D distance along the track (in kilometers)
%               dislen 2D track segment lengths for alongtrack integral
%               en,ep unit vectors normal and parallel to the track for 
%                    use by roms_velslice when computing across-track and 
%                    along-track velocity
%
% John Wilkin 24-Oct-2013 / 13-Apr-2016 / 08-Jul-2016
%
% This code was re-written 08-Jul-2016 - there were errors in the 
% along/across-track rotation math. 
% Versions of roms_velslice older than 2016 SHOULD BE DELETED. 
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu 

% u,v in native ROMS coordinate directions along the track
[Tu,Coordu] = roms_genslice(file,uname,varargin{:});
[Tv,Coordv] = roms_genslice(file,vname,varargin{:});

% Unit vectors (complex) in cross (normal) and along-track (parallel) 
% directions 
en = Coordu.en;
ep = Coordu.ep;

% Use complex variable for vector data
vel = complex(Tu,Tv);

% take the Real part of multiplication by complex conjugate to get
% vector inner product

% the track parallel velocity 
Ep = repmat(ep,[size(Tu,1) 1]); 
Valong = real(conj(Ep).*vel);

% the track normal velocity 
En = repmat(en,[size(Tu,1) 1]); 
Vcross = real(conj(En).*vel);

% Yes, it's really that simple. No need to wrestle with the grid angle
% w.r.t. east or any of that because this code above works in the
% fractional i,j coordinates aligned with the grid. The magic is done in
% roms_genslice computing the unit vectors en,ep. 

% Output
vel      = complex(Valong,Vcross);
coord    = Coordu;
coord.z  = 0.5*(Coordu.z+Coordv.z);
coord.dz = 0.5*(Coordu.dz+Coordv.dz);
coord.h  = 0.5*(Coordu.h+Coordv.h);

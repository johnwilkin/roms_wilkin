function [xi,yi,zi] = griddata_lite(x,y,z,xi,yi,method,options)
% $Id: griddata_lite.m 395 2011-01-10 18:05:38Z wilkin $
% Stripped down version of GRIDDATA that skips a bunch of data checks and 
% saves the delaunay triangulation for use on subsequent calls with the
% exact same set of data coordinates
%
% Invoke griddata_lite the same way as griddata
% e.g. zi = griddata_lite(x,y,z,xi,yi,method)
% 
% Note that only one output is allowed (the interpolated data) and that 
% 'method' must be either 'linear' (the default) or 'nearest'.
%
% Global variable TRI is declared, and if it is empty then the Delaunay
% triangulation is performed as it would be by griddata. The triangulation
% information is saved in global variable TRI (actually a structure) for use 
% on a subsequent call to griddata_lite for input data at the same locations. 
%
% IMPORTANT: If the input data locations alter (e.g. by switching between
% ROMS Arakawa C-grid locations) then the Delaunay triangulation has to be
% redone. This is forced by clearing the global variable TRI and redeclaring 
% it (as empty):
%
%    clear global TRI (be sure to include option global)
%
%
% See help griddata for more info on the method
%
% John Wilkin, October 2002
%    Modified January 2008 to add the 'nearest' option
%    Modified January 2011
%      (1) to check if the input coordinates change in which 
%          case recalculation of the triangulation is forced
%      (2) to repeat the tsearch on every interpolation - this
%          is required for safety in case the data being gridded
%          has different different distributions of NaNs even if
%          the data coordinates have not altered

%GRIDDATA Data gridding and surface fitting.
%   ZI = GRIDDATA(X,Y,Z,XI,YI) fits a surface of the form Z = F(X,Y) to the
%   data in the (usually) nonuniformly-spaced vectors (X,Y,Z). GRIDDATA
%   interpolates this surface at the points specified by (XI,YI) to produce
%   ZI.  The surface always goes through the data points. XI and YI are
%   usually a uniform grid (as produced by MESHGRID) and is where GRIDDATA
%   gets its name.
%
%   XI can be a row vector, in which case it specifies a matrix with
%   constant columns. Similarly, YI can be a column vector and it specifies
%   a matrix with constant rows.
%
%   [XI,YI,ZI] = GRIDDATA(X,Y,Z,XI,YI) also returns the XI and YI formed
%   this way (the results of [XI,YI] = MESHGRID(XI,YI)).
%
%   [...] = GRIDDATA(X,Y,Z,XI,YI,METHOD) where METHOD is one of
%       'linear'    - Triangle-based linear interpolation (default)
%       'cubic'     - Triangle-based cubic interpolation
%       'nearest'   - Nearest neighbor interpolation
%       'v4'        - MATLAB 4 griddata method
%   defines the type of surface fit to the data. The 'cubic' and 'v4'
%   methods produce smooth surfaces while 'linear' and 'nearest' have
%   discontinuities in the first and zero-th derivative respectively.  All
%   the methods except 'v4' are based on a Delaunay triangulation of the
%   data.
%   If METHOD is [], then the default 'linear' method will be used.
%
%   [...] = GRIDDATA(X,Y,Z,XI,YI,METHOD,OPTIONS) specifies a cell array 
%   of strings OPTIONS to be used as options in Qhull via DELAUNAYN. 
%   If OPTIONS is [], the default DELAUNAYN options will be used.
%   If OPTIONS is {''}, no options will be used, not even the default.
%
%   Example:
%      x = rand(100,1)*4-2; y = rand(100,1)*4-2; z = x.*exp(-x.^2-y.^2);
%      ti = -2:.25:2; 
%      [xi,yi] = meshgrid(ti,ti);
%      zi = griddata(x,y,z,xi,yi);
%      mesh(xi,yi,zi), hold on, plot3(x,y,z,'o'), hold off
%
%   See also GRIDDATA3, GRIDDATAN, DELAUNAY, INTERP2, MESHGRID, DELAUNAYN.

%   Copyright 1984-2007 The MathWorks, Inc. 
%   $Revision: 5.33.4.9 $  $Date: 2008/06/20 08:00:53 $

if nargin < 6
  method = 'linear';
end
switch method
  case { 'linear','nearest'}
  otherwise
    error('only methods LINEAR and NEARSET are supported')
end

global TRI

% % error(nargchk(5,7,nargin,'struct'));
% % 
% % [msg,x,y,z,xi,yi] = xyzchk(x,y,z,xi,yi);
% % if ~isempty(msg), error(msg); end
% % if ndims(x) > 2 || ndims(y) > 2 || ndims(xi) > 2 || ndims(yi) > 2
% %     error('MATLAB:griddata:HigherDimArray',...
% %           'X,Y and XI,YI cannot be arrays of dimension greater than two.');
% % end
% % 
% % if ( issparse(x) || issparse(y) || issparse(z) || issparse(xi) || issparse(yi) )
% %     error('MATLAB:griddata:InvalidDataSparse',...
% %         'Input data cannot be sparse.');
% % end
% % 
% % if ( ~isreal(x) || ~isreal(y) || ~isreal(xi) || ~isreal(yi) )
% %     error('MATLAB:griddata:InvalidDataComplex',...
% %         'Input data cannot be complex.');
% % end
% % 
% % if ( nargin < 6 || isempty(method) ),  method = 'linear'; end
% % if ~ischar(method), 
% %   error('MATLAB:griddata:InvalidMethod',...
% %         'METHOD must be one of ''linear'',''cubic'',''nearest'', or ''v4''.');
% % end
% % 
% % if nargin == 7
% %     if ~iscellstr(options)
% %         error('MATLAB:OptsNotStringCell',...
% %               'OPTIONS should be cell array of strings.');           
% %     end
% %     opt = options;
% % else
% %     opt = [];
% % end
% % 
% % if numel(x) < 3 || numel(y) < 3
% %   error('MATLAB:griddata:NotEnoughSamplePts',...
% %         'Not enough unique sample-points specified.');
% % end

opt = [];

% Sort x and y so duplicate points can be averaged before passing to delaunay

%Need x,y and z to be column vectors
sz = numel(x);
x = reshape(x,sz,1);
y = reshape(y,sz,1);
z = reshape(z,sz,1);
sxyz = sortrows([x y z],[2 1]);
x = sxyz(:,1);
y = sxyz(:,2);
z = sxyz(:,3);

% JW: Disable duplicate data check for speed - this does not happen with model output

% % myepsx = eps(0.5 * (max(x) - min(x)))^(1/3);
% % myepsy = eps(0.5 * (max(y) - min(y)))^(1/3);
% % ind = [0; ((abs(diff(y)) < myepsy) & (abs(diff(x)) < myepsx)); 0];

% % if sum(ind) > 0
% %   warning('MATLAB:griddata:DuplicateDataPoints',['Duplicate x-y data points ' ...
% %             'detected: using average of the z values.']);
% %   fs = find(ind(1:end-1) == 0 & ind(2:end) == 1);
% %   fe = find(ind(1:end-1) == 1 & ind(2:end) == 0);
% %   for i = 1 : length(fs)
% %     % averaging z values
% %     z(fe(i)) = mean(z(fs(i):fe(i)));
% %   end
% %   x = x(~ind(2:end));
% %   y = y(~ind(2:end));
% %   z = z(~ind(2:end));
% % end

switch lower(method),
  case 'linear'
    zi = linear(x,y,z,xi,yi,opt);
  case 'cubic'
    if(isreal(z))
        zi = cubic(x,y,z,xi,yi,opt);
    else
        zi = complex(cubic(x,y,real(z),xi,yi,opt),cubic(x,y,imag(z),xi,yi,opt));
    end
  case 'nearest'
    zi = nearest(x,y,z,xi,yi,opt);
  case {'invdist','v4'}
    zi = gdatav4(x,y,z,xi,yi);
  otherwise
    error('MATLAB:griddata:UnknownMethod', 'Unknown method.');
end
  
if nargout<=1, xi = zi; end


%------------------------------------------------------------
function zi = linear(x,y,z,xi,yi,opt)
%LINEAR Triangle-based linear interpolation

%   Reference: David F. Watson, "Contouring: A guide
%   to the analysis and display of spacial data", Pergamon, 1994.

siz = size(xi);
xi = xi(:); yi = yi(:); % Treat these as columns
x = x(:); y = y(:); % Treat these as columns

global TRI

retri = 0;
if isempty(TRI)
  retri = 1;
else
  % check whether coordinates changed
  if ~isempty(setdiff(x,TRI.x)) || ~isempty(setdiff(y,TRI.y))
    disp('coordinates appear to have changed - must retriangulate')
    retri = 1;
  end
end

if retri
  disp([mfilename ' is triangulating data ...'])
  % triangularize the data
  tri = delaunayn([x y]);
  if isempty(tri),
    warning('Data cannot be triangulated.');
    zi = repmat(NaN,size(xi));
    return
  end
  
  % % % Triangularize the data
  % % if isempty(opt)
  % %     tri = delaunayn([x y]);
  % % else
  % %     tri = delaunayn([x y],opt);
  % % end
  % % if isempty(tri),
  % %   warning('MATLAB:griddata:CannotTriangulate','Data cannot be triangulated.');
  % %   zi = repmat(NaN,size(xi));
  % %   return
  % % end
  
  % save in global variable
  TRI.tri = tri;
  TRI.x = x;
  TRI.y = y;
  disp('      ... triangulation complete')
  
else
  
  % use the triangulation held in global from last
  % call to griddata_lite
  tri = TRI.tri;

end

% Find the nearest triangle (t)
t = tsearch(x,y,tri,xi,yi);

% Only keep the relevant triangles.
out = find(isnan(t));
if ~isempty(out), t(out) = ones(size(out)); end
tri = tri(t,:);
  
% Compute Barycentric coordinates (w).  P. 78 in Watson.
del = (x(tri(:,2))-x(tri(:,1))) .* (y(tri(:,3))-y(tri(:,1))) - ...
      (x(tri(:,3))-x(tri(:,1))) .* (y(tri(:,2))-y(tri(:,1)));
w(:,3) = ((x(tri(:,1))-xi).*(y(tri(:,2))-yi) - ...
          (x(tri(:,2))-xi).*(y(tri(:,1))-yi)) ./ del;
w(:,2) = ((x(tri(:,3))-xi).*(y(tri(:,1))-yi) - ...
          (x(tri(:,1))-xi).*(y(tri(:,3))-yi)) ./ del;
w(:,1) = ((x(tri(:,2))-xi).*(y(tri(:,3))-yi) - ...
          (x(tri(:,3))-xi).*(y(tri(:,2))-yi)) ./ del;
w(out,:) = zeros(length(out),3);

z = z(:).'; % Treat z as a row so that code below involving
            % z(tri) works even when tri is 1-by-3.
zi = sum(z(tri) .* w,2);

zi = reshape(zi,siz);

if ~isempty(out), zi(out) = NaN; end
%------------------------------------------------------------

%------------------------------------------------------------
function zi = cubic(x,y,z,xi,yi,opt)
%TRIANGLE Triangle-based cubic interpolation

%   Reference: T. Y. Yang, "Finite Element Structural Analysis",
%   Prentice Hall, 1986.  pp. 446-449.
%
%   Reference: David F. Watson, "Contouring: A guide
%   to the analysis and display of spacial data", Pergamon, 1994.

% Triangularize the data
if isempty(opt)
    tri = delaunayn([x(:) y(:)]);
else
    tri = delaunayn([x(:) y(:)],opt);
end
if isempty(tri), 
  warning('MATLAB:griddata:CannotTriangulate','Data cannot be triangulated.');
  zi = repmat(NaN,size(xi));
  return
end

% Find the nearest triangle (t)
t = tsearch(x,y,tri,xi,yi);

zi = cubicmx(x,y,z,xi,yi,tri,t);
%------------------------------------------------------------

%------------------------------------------------------------
function zi = nearest(x,y,z,xi,yi,opt)
%NEAREST Triangle-based nearest neightbor interpolation

%   Reference: David F. Watson, "Contouring: A guide
%   to the analysis and display of spacial data", Pergamon, 1994.

siz = size(xi);
xi = xi(:); yi = yi(:); % Treat these a columns
x = x(:); y = y(:); z = z(:); % Treat these as columns

global TRI

retri = 0;
if isempty(TRI)
  retri = 1;
else
  % check whether coordinates changed
  if ~isempty(setdiff(x,TRI.x)) || ~isempty(setdiff(y,TRI.y))
    disp('coordinates appear to have changed - must retriangulate')
    retri = 1;
  end
end

if retri
  disp([mfilename ' is triangulating data ...'])
  % triangularize the data
  tri = delaunayn([x y]);
  if isempty(tri),
    warning('Data cannot be triangulated.');
    zi = repmat(NaN,size(xi));
    return
  end
  
  % % % Triangularize the data
  % % if isempty(opt)
  % %     tri = delaunayn([x y]);
  % % else
  % %     tri = delaunayn([x y],opt);
  % % end
  % % if isempty(tri),
  % %   warning('MATLAB:griddata:CannotTriangulate','Data cannot be triangulated.');
  % %   zi = repmat(NaN,size(xi));
  % %   return
  % % end
  
  % save in global variable
  TRI.tri = tri;
  TRI.x = x;
  TRI.y = y;
  disp('      ... triangulation complete')
  
else
  
  % use the triangulation held in global from last
  % call to griddata_lite
  tri = TRI.tri;

end

% Find the nearest vertex
k = dsearch(x,y,tri,xi,yi);
zi = k;
d = find(isfinite(k));
zi(d) = z(k(d));
zi = reshape(zi,siz);
%----------------------------------------------------------


%----------------------------------------------------------
function [xi,yi,zi] = gdatav4(x,y,z,xi,yi)
%GDATAV4 MATLAB 4 GRIDDATA interpolation

%   Reference:  David T. Sandwell, Biharmonic spline
%   interpolation of GEOS-3 and SEASAT altimeter
%   data, Geophysical Research Letters, 2, 139-142,
%   1987.  Describes interpolation using value or
%   gradient of value in any dimension.

xy = x(:) + y(:)*sqrt(-1);

% Determine distances between points
d = xy(:,ones(1,length(xy)));
d = abs(d - d.');
n = size(d,1);
% Replace zeros along diagonal with ones (so these don't show up in the
% find below or in the Green's function calculation).
d(1:n+1:numel(d)) = ones(1,n);

non = find(d == 0);
if ~isempty(non),
  % If we've made it to here, then some points aren't distinct.  Remove
  % the non-distinct points by averaging.
  [r,c] = find(d == 0);
  k = find(r < c);
  r = r(k); c = c(k); % Extract unique (row,col) pairs
  v = (z(r) + z(c))/2; % Average non-distinct pairs
  
  rep = find(diff(c)==0);
  if ~isempty(rep), % More than two points need to be averaged.
    runs = find(diff(diff(c)==0)==1)+1;
    for i=1:length(runs),
      k = find(c==c(runs(i))); % All the points in a run
      v(runs(i)) = mean(z([r(k);c(runs(i))])); % Average (again)
    end
  end
  z(r) = v;
  if ~isempty(rep),
    z(r(runs)) = v(runs); % Make sure average is in the dataset
  end

  % Now remove the extra points.
  x(c) = [];
  y(c) = [];
  z(c) = [];
  xy(c,:) = [];
  xy(:,c) = [];
  d(c,:) = [];
  d(:,c) = [];
  
  % Determine the non distinct points
  ndp = sort([r;c]);
  ndp(find(ndp(1:length(ndp)-1)==ndp(2:length(ndp)))) = [];

  warning('MATLAB:griddata:NonDistinctPoints',['Averaged %d non-distinct ' ...
            'points.\n         Indices are: %s.'],length(ndp),num2str(ndp'))
end

% Determine weights for interpolation
g = (d.^2) .* (log(d)-1);   % Green's function.
% Fixup value of Green's function along diagonal
g(1:size(d,1)+1:numel(d)) = zeros(size(d,1),1);
weights = g \ z(:);

[m,n] = size(xi);
zi = zeros(size(xi));
jay = sqrt(-1);
xy = xy.';

% Evaluate at requested points (xi,yi).  Loop to save memory.
for i=1:m
  for j=1:n
    d = abs(xi(i,j)+jay*yi(i,j) - xy);
    mask = find(d == 0);
    if length(mask)>0, d(mask) = ones(length(mask),1); end
    g = (d.^2) .* (log(d)-1);   % Green's function.
    % Value of Green's function at zero
    if length(mask)>0, g(mask) = zeros(length(mask),1); end
    zi(i,j) = g * weights;
  end
end

if nargout<=1,
  xi = zi;
end
%----------------------------------------------------------


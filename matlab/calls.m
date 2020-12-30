function out = calls(theFcn, varargin)
% calls -- List of function calls.
%  calls('demo') demonstrates itself by listing the
%   routines that would be called by "bessel(2, 1)".
%  [...] = calls('theFcn', ...) executes 'theFcn' with
%	the given arguments, then displays the list of M-file
%   routines that were called, according to the Matlab
%   "which" function.  The list is also placed in the
%   caller's "ans".  If the list is empty, then all
%   the calls were made to "built-in" functions.
%
% NOTE: When using this routine to assist a bundling
%  procedure, be sure not to include files whose
%  respective copyrights would be violated.
%
% Copyright (C) 2000 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%
% Version of 18-Apr-2000 09:11:50.
% Updated    18-Apr-2000 09:46:21.
% 
% Modified J. Wilkin Jan 2020
% to work with more recent Matlab versions

if nargin > 1
  feval(theFcn, varargin{:});
else
  feval(theFcn);
end

theCalled = inmem;

for i = 1:length(theCalled)
  theCalled{i} = which(theCalled{i});
end

w = which(mfilename);
for i = length(theCalled):-1:1
  if isequal(w, theCalled{i})
    theCalled(i) = [];
    break
  end
end

if isempty(theCalled)
  theCalled = '    (none)';
end

assignin('caller', 'ans', theCalled)
if nargout == 0
disp(' ')
disp([' ## M-files called by "' theFcn '":'])
disp(' ')
disp(theCalled)
end

out = theCalled;



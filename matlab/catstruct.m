function B = catstruct(dim,sa,sb)
% CATSTRUCT(DIM,SA,SB) concatenates all fields in structures SA and SB
% along the dimension DIM
%
% Useful when a large data download task is broken into a series of smaller
% queries that return similar structures, for example several ERDDAP
% requests stepping through a sequence of time intervals. 
%
% A try-catch test is used to skip over fields in the structure for which
% the concatentation fails - typically string/character data where 
% one of the dimensions doesn't match
%
% John Wilkin - June 2019
%
% Copyright (c) 2020 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: catstruct.m 599 2020-12-29 16:53:00Z wilkin $

na = fieldnames(sa);

for name = na'
  if ~isfield(sb,name)
    warning(['No ' char(name) ' in SB - filling with NaNs to allow cat'])
    % hope first variable in SB is the same dimension as this
    sb.(char(name)) = NaN(size(sb.(char(na{1}))));
  end
  try
    B.(char(name)) = cat(dim,sa.(char(name)),sb.(char(name)));
  catch
    warning(['Unable to cat field ' char(name)])
  end
end

function index = findstrinstruct(S,field,string)
% index = findstrinstruct(S,field,string)
%
% find INDEX into a structure S for which S.FIELD matches STRING
%
% John Wilkin - Nov 2018
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: findstrinstruct.m 599 2020-12-29 16:53:00Z wilkin $

% original code: match string exactly
index = find(arrayfun(@(n) strcmp(S(n).(field),string),1:numel(S)));

% 2024/06/19 find results that START with STRING
% index = find(arrayfun(@(n) strncmp(S(n).(field),string,length(string)), 1:numel(S)));
  

function s = strrep_(s)
% $Id: strrep_.m 433 2016-02-09 14:10:36Z wilkin $
% Escape string characters with \ so that they are not interpretted as 
% TeX instructions.
% This function is used by several roms_*view routines so that ROMS
% filenames are not corrupted in the title string
%
% John Wilkin
s = strrep(s,'\','\\');
s = strrep(s,'_','\_');
s = strrep(s,'^','\^');


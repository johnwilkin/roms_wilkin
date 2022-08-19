function [cpp1,cpp2] = roms_cpplist(file,file2)
% [cpp1,cpp2] = roms_cpplist(file,file2)
% Extract the list of CPP options from a ROMS output file
% If two files are input the options are compared
%
% John Wilkin
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id$

cpp = ncreadatt(file,'/','CPP_options');
c = strfind(cpp,',');
c = [1 c];
clear s
for k = 2:length(c)
  tstr = cpp(c(k-1):c(k));
  tstr = strrep(tstr,',','');
  tstr = strrep(tstr,' ','');
  s{k-1} = tstr;
end
cpp1 = char(s);

if nargin==2
  cpp = nc_attget(file2,nc_global,'CPP_options');
  c = strfind(cpp,',');
  c = [1 c];
  clear s
  for k = 2:length(c)
    tstr = cpp(c(k-1):c(k));
    tstr = strrep(tstr,',','');
    tstr = strrep(tstr,' ','');
    s{k-1} = tstr;
  end
  cpp2 = char(s);
  for k=1:size(cpp1,1)
    str1 = cpp1(k,:);
    found = 0;
    for j=1:size(cpp2,1)
      if strcmp(deblank(str1),deblank(cpp2(j,:)))
        disp([str1 ' in common'])
        found = 1;
      end
    end
    if found == 0
      disp([str1 ' in file 1 but not file 2'])
    end
  end
  for k=1:size(cpp2,1)
    str1 = cpp2(k,:);
    found = 0;
    for j=1:size(cpp1,1)
      if strcmp(deblank(str1),deblank(cpp1(j,:)))
        found = 1;
      end
    end
    if found == 0
      disp([str1 ' in file 2 but not file 1'])
    end
  end
end


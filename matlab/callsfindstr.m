function callsfindstr(s,str,notstr)
% callsfindstr(s,str,notstr)
%
% Find entries in the structure S created by CALLS that contain
% string STR and, optionally, do not contain NOTSTR
%
% If STR is empty, the test is simply NOTSTR
%
% Use this to check the portability of a routine.
%
% For example: To find
%      S = calls('myfunction',args)
%      callsfindstr(S,'wilkin','roms_wilkin')
%
% John Wilkin
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu

if nargin == 2
  opt = 1;    % just find if str occurs
else
  if isempty(str)
    opt = 2;  % just find if notstr does not occur
  else
    opt = 3;  % find cases where str occurs and notstr does not
  end
end

for i=1:length(s)
  chkstr = char(s{i});
  switch opt
    case 1
      if contains(chkstr,str)
        disp(chkstr)
      end
    case 2
      if ~contains(chkstr,notstr)
        disp(chkstr)
      end
    case 3
      if contains(chkstr,str) && ~contains(chkstr,notstr)
        disp(chkstr)
      end
  end
end

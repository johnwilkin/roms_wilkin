function callsfindstr(s,str,notstr)
% callsfindstr(s,str,notstr)
%
% Find entries in the structure S - created by CALLS - that contain
% string STR and/or optionally do not contain NOTSTR
%
% If STR is empty, the test is simply not NOTSTR
%
% I use this to check the portability of a routine.
%
% Example usage:
%
%     Suppose you have a function that you want to share with a colleague. 
%     You know it calls many other utilities that are not part of the 
%     Matlab release and you want to know what those are so you can share a
%     complete package of all the necessary routines. 
%
%     Execute your function, 'myfunction.m', with the CALLS function, 
%     saving the result in a structure, S. Include any optional inputs 
%     your function requires in the argument list as shown below.  
%     Note: This must be a function, not a script. If it's a script, wrap 
%     in a function ensuring you pass along any necessary variables from 
%     your workspace as inputs. 
%
%          S = calls('myfunction',arg1,arg2,...)
%
%     CALLSFINDSTR will parse the list (in S) of ALL functions that your
%     MYFUNCTION called. You might search for your username to reveal
%     functions in your own paths (as opposed to those in the default 
%     Matlab path) that should be distributed with MYFUNCTION for 
%     it to be portable, e.g.:
%
%          callsfindstr(S,'wilkin','myroms')
%
%     But you might want to exclude routines you know a user would have 
%     from other contributed toolboxes, e.g.:
%     
%          callsfindstr(S,'wilkin','myroms')
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

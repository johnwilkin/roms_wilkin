function [file,time] = roms_filetime_fromctl(ctl,time)
% $Id: roms_filetime_fromctl.m 358 2008-04-07 14:15:03Z zhang $
% [file,time] = roms_filetime_fromctl(ctl,tinput)
% extracts the filename and time index from a control structure of multiple
% files created by roms_timectl
%
% John Wilkin

if isstr(time)
  fdnums = ctl.dnum;
  dnum = datenum(time);
  if dnum >= fdnums(1) & dnum <= fdnums(end)
    [tmp,time] = min(abs(dnum-fdnums));
    time = time(1);
  else
    warning(' ')
    disp(['Requested date ' time ' is not between the dates in '])
    disp(['ctl which are ' datestr(fdnums(1),0) ' to ' ])
    disp(datestr(fdnums(end),0))
    thedata = -1;
    return
  end
end
  
file = char(ctl.files{ctl.file_index(time)});
file = strrep(file,'\','/');
time = ctl.in_file_index(time)+1;


function ctl = roms_timectl( arg1, arg2, arg3)
% $Id: roms_timectl.m 426 2014-11-07 18:22:03Z wilkin $
% ROMS_TIMECTL: retrieves files in a directory matching a regular 
% expression and creates a structure noting filename/timevalue for managing 
% ROMS output split across multiple files. The ctl structure can be 
% passed to roms_zview etc in lieu of a single netcdf filename. 
%
% After retrieving matching files, they are sorted by time.  
%
% USAGE:  ctl = roms_timectl (direc,regular_expression, timevarname);
% USAGE:  ctl = roms_timectl (opendap_url);
%
% PARAMETERS:
% Input:
%     direc:  
%         directory path
%     regular_expression:
%         a regular expression suitable to be used with the REGEXP 
%         function.  Make sure you are fully versed with REGEXP!
%         This argument is optional, if not supplied, it is assumed to
%         be '.', which matches everything.
%     timevarname:
%         If not supplied, assume this to be 'ocean_time'
%
%     or
%
%     opendap_url:
%         
% Output:
%     ctl:
%         structure with fields
%
%         files:
%             Cell array of all files (full pathnames) that matched 
%             the regular expression.  If an opendap url is given,
%             then the cell array contains just this.
%         time:
%             Array of sorted time.
%         file_index:
%             Array of indices that tag each time value to each file.
%         in_file_index:
%             Array of indices that tag each time value to it's index
%             within the file specified by file_index
%
% Example:
%     direc = '.'
%     ctl = roms_timectl( '.','his');
%     ctl = roms_timectl( '.','his.*_008\d\.nc');
%
%     The first case would return all the history files in the current
%     directory.  The second example would return only those files
%     enumerated '0080', '0081', '0082', ... '0089'.
%
% John Evans / John Wilkin (originating with Manu DiLorenzo's rnt tools) 

files = [];

if nargin < 1
	msg = sprintf ( '%s:  must supply at least one input argument.\n', mfilename );
	error ( msg );
end

urlcheck = regexp ( arg1, 'http://.*/.*' );
if ~isempty(urlcheck)
	opendap_url = arg1;
	ctl.files = { opendap_url };
	ctl.time = nc_varget ( opendap_url, 'ocean_time' );
	ctl.file_index = ones(length(ctl.time),1);
	ctl.in_file_index = [1:length(ctl.time)]';
	return
end

direc = arg1;

%
% If not supplied, then match everything.
if nargin < 2
	rexp = '.';
	timevar = 'ocean_time';
else
	rexp = arg2;
end



if nargin < 3
	timevar = 'ocean_time';
else
  timevar = arg3;
end



if ~exist ( direc, 'dir' )
	msg = sprintf ( '%s:  directory ''%s'' does not exist.\n', mfilename, direc );
	error ( msg );
end

d_struct = dir ( direc );


%
% This keeps track of how many files we matched.
match_count = 0;

%
% go thru the directory listing and try to match them up
for j = 1:length(d_struct)

	%
	% skip directories, of course
	if ( d_struct(j).isdir )
		continue
	end

	start_match = regexp ( d_struct(j).name, rexp );
	if ~isempty(start_match)
		candidate = sprintf ( '%s%s%s', direc, filesep,  d_struct(j).name );

		[ncid, status] = mexnc ( 'open', candidate, nc_noclobber_mode );
		if ( status == 0 )
			match_count = match_count + 1;
			mexnc ( 'close', ncid );
			files{match_count,1} = candidate;
		else
			msg = sprintf ( 'Could not open %s...\n', candidate );
			error ( msg );
		end
	end

end

if isempty(files)
	return
end

%
% Now sort by time.  
ocean_time = [];
index = [];
in_file_index = [];
for j = 1:length(files)
	t = nc_varget ( files{j}, timevar );
	ocean_time = [ocean_time; t];
	index = [index; j*ones(length(t),1)];
	in_file_index = [in_file_index; [0:length(t)-1]'];
end

[ocean_time,I] = sort ( ocean_time );

index = index(I);
in_file_index = in_file_index(I);

ctl.files = files;
ctl.time = ocean_time;
ctl.file_index = index;
ctl.in_file_index = in_file_index;

% Wilkin added calculation of Matlab datenum equivalent for ocean_time
tunits = nc_attget(ctl.files{1},timevar,'units');
switch tunits(1:3)
  case 'day'
    fac = 1;
  case 'hou'
    fac = 1/24;
  case 'sec'
    fac = 1/86400;
  otherwise
    warning(' ')
    disp(['Got time units of ' tunits])
    disp(['but don''t have a rule for converting to datenum'])
    ctl.dnum = [];
    return
end
try
  ctl.dnum = datenum(parsetnc(tunits)) + ocean_time*fac;
catch
  ctl.dnum = ocean_time*fac;
end




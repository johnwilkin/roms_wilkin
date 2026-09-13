function [DATA,url] = erddap_read(varargin)
% DATA = erddap_read(URL)
% Load data from an ERDDAP Data Access Form query for filetype 'mat'
%      created with the "Just generate the URL" button.
%      The data are returned in the output structure named DATA.
%      Any time coordinate data is converted to a Matlab datenum.
%
% This function circumvents the inconvenience of needing to know the
% datasetID that ERDDAP uses to name the output structure.
%
% Some ERDDAP services require percent encoding of the URL, which this
% function does if the query at first throws an error.
%
% John Wilkin - August 2018
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu

% Download to temporary file
tmpfile = ['tmp_' randstr(10) '.mat'];

wopt = [];
encodeurl = false;

url = varargin{1};
for k = 2:nargin
	if isa(varargin{k},'weboptions')
		wopt = varargin{k};
	else
		encodeurl = true;
	end
end
if encodeurl
	disp('URL encoding the dataurl')
	url = [url(1:strfind(url,'.mat?')+4) ...
		urlencode(url(strfind(url,'.mat?')+5:end))];
end

% wopt = weboptions('Timeout',600);

try
	if isempty(wopt)
		tmp = load(websave(tmpfile,url));
	else
		tmp = load(websave(tmpfile,url,wopt));
	end
catch
	warning("EDRRAP read failed"+newline+ ...
		"Try extending timeout "+newline+ ...
		"erddap_read(url,weboptions('Timeout',600)); "+newline+ ...
		"and/or encoding special characters in the url "+newline+ ...
		"erddap_read(url,true)")
	DATA = NaN;
	return
end

% There will be one field in the structure - get it, and call it DATA
names = fieldnames(tmp);
DATA = tmp.(names{1});

if isfield(DATA,'time') % convert time to a datenum
	% ERDDAP always returns time data "seconds since 1970-01-01"
	DATA.time = DATA.time/86400+datenum(1970,1,1);
	DATA.datetime = datetime(DATA.time,'ConvertFrom','datenum');
end

delete(tmpfile) % remove the temporary file

function str = randstr(len)
% str = randstr(len) creates a random lowercase string of length LEN
% John Wilkin

% ASCII codes 97-122 are a-z
str = char(floor(25*rand(1,len)+97));


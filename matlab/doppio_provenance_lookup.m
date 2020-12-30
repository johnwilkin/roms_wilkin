function [short,long,num] = doppio_provenance_lookup(prov)
% [Shortname,Description,Prov_number] = doppio_provenance_lookup(PROV)
%
% Query the Doppio provenance table in the ERDDAP service at
% http://tds.marine.rutgers.edu/erddap/tabledap/DOPPIO_PROVENANCE.html
% to translate a provenance number into its short name (short_CODE) and
% longer name (long_Description), or for an input code return its number.
%
% Input PROV can be a scalar, or a 2-element vector defining a range
% request for all valid provenances in that interval. In the latter case
% the outputs are a vector of valid provenance numbers, and string
% matrices of corresponding short and long names.
%
% In the range vector input case outputs are only given for valid defined
% provenances in the table.
%
% When the input provenance is a scalar, the num output just echoes this,
% and if no entry exists in the table for that value the function returns
% empty strings, e.g.
%     [short,~] = doppio_provenance_lookup(560)
%     returns short = 'AMAG'
%
% When the input is a string, the num output is the numeric code, e.g.:
%     [~,~,num] = doppio_provenance_lookup('AMAG')
%     returns num = 560
%
% The ERDDAP catalog entry is maintained by Eli Hunter. It reads the file
% /home/om/roms/doppio/Provenancetable2018.txt
% edited by Julia Levin and John Wilkin
%
% John Wilkin - 18 April 2018
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: doppio_provenance_lookup.m 545 2020-01-09 19:41:09Z wilkin $

server =  'http://tds.marine.rutgers.edu/erddap/tabledap/';

if isnumeric(prov)
  url =  [server 'DOPPIO_PROVENANCE.mat?' ...
    'Provenance%2CCODE%2CDescription&Provenance=' int2str(prov(1))];
  
  try
    load(urlwrite(url,'doppio_provenance_lookup_tmpfile.mat'))
    num = prov(1);
    short = DOPPIO_PROVENANCE.CODE;
    long = DOPPIO_PROVENANCE.Description;
    unix('rm doppio_provenance_lookup_tmpfile.mat');
  catch
    num = prov;
    short = '';
    long = '';
    warning(['There is no entry for number ' int2str(prov)])
  end
  
  if length(prov)>1
    for k=2:length(prov)
      url =  [server 'DOPPIO_PROVENANCE.mat?' ...
        'Provenance%2CCODE%2CDescription&Provenance=' int2str(prov(k))];
      try
        load(urlwrite(url,'doppio_provenance_lookup_tmpfile.mat'))
        num(k) = DOPPIO_PROVENANCE.Provenance;
        short = char(short,DOPPIO_PROVENANCE.CODE);
        long = char(long,DOPPIO_PROVENANCE.Description);
        unix('rm doppio_provenance_lookup_tmpfile.mat');
      catch
        warning(['There is no entry for number ' int2str(prov)])
      end
    end
  end
  
else
  url =  [server 'DOPPIO_PROVENANCE.mat?' ...
    'Provenance%2CCODE%2CDescription&CODE=%22' prov '%22'];
  try
    load(urlwrite(url,'doppio_provenance_lookup_tmpfile.mat'))
    num = DOPPIO_PROVENANCE.Provenance;
    short = DOPPIO_PROVENANCE.CODE;
    long = DOPPIO_PROVENANCE.Description;
    unix('rm doppio_provenance_lookup_tmpfile.mat');
  catch
    num = NaN;
    short = '';
    long = '';
    warning(['There is no provenance code for ' prov])
  end
end

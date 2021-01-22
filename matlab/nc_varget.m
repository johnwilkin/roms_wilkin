function outf=nc_vargetj(f,u,varargin)
% Subvert a call to nc_varget (from snctools) and instead use Matlab 
% native ncread
% John Wilkin and Jack McSweeney - May 2016
%
% nc_varget  Retrieve data from netCDF variable or HDF4 data set.
%     DATA = nc_varget(NCFILE,VARNAME) retrieves all the data from the 
%     variable VARNAME in the netCDF file NCFILE.  
%  
%     DATA = nc_varget(NCFILE,VARNAME,START,COUNT) retrieves the contiguous
%     portion of the variable specified by the index vectors START and 
%     COUNT.  Remember that SNCTOOLS indexing is zero-based, not 
%     one-based.  Specifying Inf in COUNT means to retrieve everything 
%     along that dimension from the START coordinate.
%  
%     DATA = nc_varget(NCFILE,VARNAME,START,COUNT,STRIDE) retrieves 
%     a non-contiguous portion of the dataset.  The amount of
%     skipping along each dimension is given through the STRIDE vector.
%  
%     DATA is returned as double precision when this is reasonable. 
%     Consequently, '_FillValue' and 'missing_value' attributes are honored 
%     by flagging those datums as NaN.  Any 'scale_factor' and 'add_offset' 
%     attributes are honored by applying the linear transformation.
%  
%     EXAMPLE:  This example file is shipped with R2008b.
%         data = nc_varget('example.nc','peaks',[0 0], [20 30]);
%  
%     EXAMPLE: Retrieve data from a URL. This requires the netcdf-java 
%     backend.
%         url = 'http://coast-enviro.er.usgs.gov/models/share/balop.nc';
%         data = nc_varget(url,'lat_rho');
%  
%     Example:  Retrieve data from the example HDF4 file.
%         data = nc_varget('example.hdf','Example SDS');

nargs=length(varargin);

% For 3rd and subsequent arguments reverse the order of indices
% Replace 0-based with 1-based index counting
% Replace -1 with Inf to request all data for a dimension
if nargs >0
  start = fliplr(1+varargin{1});
  count = fliplr(varargin{2});
  count(count==-1)=Inf;
end

% We noticed snctools nc_varget can give different (wrong) results here
% by getting all values for stride ~=1
if nargs == 3
  stride=fliplr(varargin{3});
end

switch nargs
  case 0
    out=ncread(f,u);
  case 1
    error('Invalid numberof inputs')
  case 2
    out=ncread(f,u,start,count);
  case 3
    out=ncread(f,u,start,count,stride);
end

% apply double in case the code using output from this function anticipates
% double instead single
outf = double(permute(out,ndims(out):-1:1));
if isrow(outf)
  outf = outf'; % to be consistent with snctools nc_varget
end






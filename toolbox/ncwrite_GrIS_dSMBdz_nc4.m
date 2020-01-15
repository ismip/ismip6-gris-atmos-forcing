function ncwrite_GrIS(ancfile,avar,varname,dimnames,res,time, time_bounds)
% write variables to a netcdf file and add ISMIP6 specific information
% function ncwrite_GrIS(ancfile,avar,varname, dimnames, res, time, time_bounds)

% TODO: check input arguments and sizes

% Always overwrite, remove destination if exists
if exist(ancfile, 'file') ~= 0;
    delete(ancfile)
end

% Overwrite coordinate dim names assuming standard ordering (x, y, ...)
dimnames{1}='x';
dimnames{2}='y';
dimnames{3}='time';
dimnames{4}='nv';

% Open netCDF file.
ncid = netcdf.create(ancfile, bitor(netcdf.getConstant('NETCDF4'),netcdf.getConstant('CLOBBER')));

% Define the dimensions of the variable.
dimids=[];
for i=1:3
  dimids = [dimids,netcdf.defDim(ncid,dimnames{i},size(avar,i))];
end
netcdf.defDim(ncid,dimnames{4},2);

% Define a new variable in the file.
my_varID = netcdf.defVar(ncid,varname,'double',dimids);

% define compression level
netcdf.defVarDeflate(ncid,my_varID,false,true,1);

ncwriteatt(ancfile,varname, '_FillValue', 9.96921e36) ; 

% Leave define mode and enter data mode to write data.
netcdf.endDef(ncid);

% Write data to variable.
netcdf.putVar(ncid,my_varID,avar);

netcdf.close(ncid);

% Standard ISMIP6 coordinates
dx = -720000;
dy = -3450000;
nx1=1681;
ny1=2881;
nx=(nx1-1)/res+1;
ny=(ny1-1)/res+1;
xd=single(zeros(nx,1));
yd=single(zeros(ny,1));
for ip=1:nx
    xd(ip) = (dx + (ip-1) * res*1000);
end
for jp=1:ny
    yd(jp) = (dy + (jp-1) * res*1000);
end

%% Coordinates
nccreate(ancfile,'x','Dimensions',{'x',nx}, 'Datatype','single', 'Format','netcdf4');
nccreate(ancfile,'y','Dimensions',{'y',ny}, 'Datatype','single', 'Format','netcdf4');
nccreate(ancfile,'time','Dimensions',{'time',1}, 'Datatype','single', 'Format','netcdf4');
nccreate(ancfile,'time_bounds','Dimensions',{'nv',2,'time',1}, 'Datatype','single', 'Format','netcdf4');
ncwrite(ancfile,'x',xd);
ncwrite(ancfile,'y',yd);
ncwrite(ancfile,'time',time);
ncwrite(ancfile,'time_bounds',time_bounds(:));

ncwriteatt(ancfile,'x', 'units', 'm') ;
ncwriteatt(ancfile,'x', 'standard_name', 'projection_x_coordinate') ;
ncwriteatt(ancfile,'x', 'axis', 'x') ;
ncwriteatt(ancfile,'y', 'units', 'm') ;
ncwriteatt(ancfile,'y', 'standard_name', 'projection_y_coordinate') ;
ncwriteatt(ancfile,'y', 'axis', 'y') ;
ncwriteatt(ancfile,'time', 'units', 'days since 1900-1-1 00:00:00') ;
ncwriteatt(ancfile,'time', 'standard_name', 'time') ;
ncwriteatt(ancfile,'time', 'calendar', 'gregorian') ;
ncwriteatt(ancfile,'time', 'bounds', 'time_bounds') ;
ncwriteatt(ancfile,'time', 'axis', 'time') ;

%% Attributes
ncwriteatt(ancfile,'/','proj4','+init=epsg:3413')
ncwriteatt(ancfile,'/','Description',['MAR SMB change with surface elevation. Prepared for ISMIP6 by Heiko Goelzer (h.goelzer@uu.nl), IMAU, ', date ])
ncwriteatt(ancfile,'/','MAR-contact','xavierfettweis@uliege.be')
ncwriteatt(ancfile,'/','MAR-institute','University of Liege (Belgium)')

% Add attributes to variables
ncwriteatt(ancfile,varname, 'units', 'kg m-2 s-1 m-1') ;
ncwriteatt(ancfile,varname, 'long_name', 'Surface mass balance change with surface elevation') ;
ncwriteatt(ancfile,varname, 'grid_mapping', 'mapping') ;

% Mapping
ncwrite_mapping_GrIS(ancfile);

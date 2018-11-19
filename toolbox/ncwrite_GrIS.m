function ncwrite_GrIS(ancfile,avar,varname,dimnames,res)
% write variables to a netcdf file and add ISMIP6 specific information
% function ncwrite_GrIS(ancfile,avar,varname, dimnames, res)

% TODO: check input arguments and sizes

% Always overwrite, remove destination if exists
if exist(ancfile, 'file') ~= 0;
    delete(ancfile)
end

% Overwrite coordinate dim names assuming standard ordering (x, y, ...)
dimnames{1}='x';
dimnames{2}='y';

% Open netCDF file.
ncid = netcdf.create(ancfile,'CLOBBER');

% Define the dimensions of the variable.
dimids=[];
for i=1:ndims(avar)
  dimids = [dimids,netcdf.defDim(ncid,dimnames{i},size(avar,i))];
end
% Define a new variable in the file.
my_varID = netcdf.defVar(ncid,varname,'double',dimids);

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
nccreate(ancfile,'x','Dimensions',{'x',nx}, 'Datatype','single', 'Format','classic');
nccreate(ancfile,'y','Dimensions',{'y',ny}, 'Datatype','single', 'Format','classic');
ncwrite(ancfile,'x',xd);
ncwrite(ancfile,'y',yd);
ncwriteatt(ancfile,'x', 'units', 'm') ;
ncwriteatt(ancfile,'y', 'units', 'm') ;
ncwriteatt(ancfile,'x', 'standard_name', 'projection_x_coordinate') ;
ncwriteatt(ancfile,'y', 'standard_name', 'projection_y_coordinate') ;
ncwriteatt(ancfile,'x', 'axis', 'x') ;
ncwriteatt(ancfile,'y', 'axis', 'y') ;

%% Attributes
ncwriteatt(ancfile,'/','proj4','+init=epsg:3413')
ncwriteatt(ancfile,'/','Description',['Created for ISMIP6 by Heiko Goelzer (h.goelzer@uu.nl), IMAU, ', date ])

% Mapping information
mapping = 'mapping';
nccreate(ancfile,'mapping','Datatype','char');
ncwriteatt(ancfile,'mapping', 'ellipsoid', 'WGS84') ;
ncwriteatt(ancfile,'mapping', 'false_easting', '0.') ;
ncwriteatt(ancfile,'mapping', 'false_northing', '0.') ;
ncwriteatt(ancfile,'mapping', 'grid_mapping_name', 'polar_stereographic') ;
ncwriteatt(ancfile,'mapping', 'latitude_of_projection_origin', '90.') ;
ncwriteatt(ancfile,'mapping', 'standard_parallel', '70.') ;
ncwriteatt(ancfile,'mapping', 'straight_vertical_longitude_from_pole', '-45') ;

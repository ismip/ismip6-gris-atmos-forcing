function ncwrite_TDSMBTable(ancfile,dimz,dimt,avar,varname,dimnames)
% write variables to a netcdf file and add ISMIP6 specific information
% function ncwrite_TDSMBTable(ancfile,dimz,dimt,avar,varname, dimnames)

% TODO: check input arguments and sizes

nz = length(dimz);
zd = dimz;
nt = length(dimt);
td = dimt;


% Always overwrite, remove destination if exists
if exist(ancfile, 'file') ~= 0;
    delete(ancfile)
end

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

%% Attributes
ncwriteatt(ancfile,'/','Description',['Created for ISMIP6 by Heiko Goelzer (h.goelzer@uu.nl), IMAU, ', date ])

%% Coordinates
nccreate(ancfile,'z','Dimensions',{'z',nz}, 'Datatype','single', 'Format','classic');
ncwrite(ancfile,'z',zd);
ncwriteatt(ancfile,'z', 'units', 'm') ;
ncwriteatt(ancfile,'z', 'axis', 'z') ;
ncwriteatt(ancfile,'z', 'positive', 'up') ;

nccreate(ancfile,'time','Dimensions',{'time',nt}, 'Datatype','single', 'Format','classic');
ncwrite(ancfile,'time',td);
ncwriteatt(ancfile,'time', 'units', 'seconds') ;
ncwriteatt(ancfile,'time', 'axis', 'time') ;


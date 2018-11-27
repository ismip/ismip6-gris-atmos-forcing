% Extract DSMB data
% Here schematic interpolation in time based on end of 19th and 20th century results 
clear

addpath('../toolbox')

%%% new MAR data
pdata='../Data/MAR/';

d19 = ncload([pdata 'MARv3.7-yearly-MIROC5-19xx.nc']);
d20 = ncload([pdata 'MARv3.7-yearly-MIROC5-20xx.nc']);
DSMB = (d20.SMB(1:5:end,1:5:end)-d19.SMB(1:5:end,1:5:end))/1000;
LAT=d19.LAT(1:5:end,1:5:end);
SH=d19.SRF(1:5:end,1:5:end);
dSMBdz19=-d19.dRU(1:5:end,1:5:end)/1000;
dSMBdz20=-d20.dRU(1:5:end,1:5:end)/1000;

% Time and amplitude
TIME=2000:1:2100;
nt=length(TIME);
ts=1:nt;
% linear
%amp = min(1, 0 + 1 * (floor(ts-1) /100)); 
% quadratic
amp = min(1, 0 + 1 * (floor(ts-1) /100).^2); 

% SMB anomaly
TDSMB = repmat(DSMB,[1,1,nt]);
for t=ts
    TDSMB(:,:,t) = TDSMB(:,:,t)*amp(t);
end
TDSMB(isnan(TDSMB))=0.0;

% dSMB/dz
TdSMBdz = repmat(dSMBdz19,[1,1,nt]);
for t=ts
    TdSMBdz(:,:,t) = TdSMBdz(:,:,t) + (dSMBdz20-dSMBdz19) * amp(t);
end

% write out
%save ../Data/MAR/trans_DSMB_MAR37.mat TDSMB SH TIME LAT TdSMBdz


% write out as netcdf
td = (0:100)*31556926;
nt = length(td);
nx = size(SH,1);
x = 1:nx;
ny = size(SH,2);
y = 1:ny;

%% dSMB/dz
ancfile = ['../Data/MAR/TdSMBdz_MAR37_MIROC5_rcp85_05000m.nc'];
ncwrite_GrIS(ancfile,TdSMBdz,'dSMBdz',{'x','y','time'},5);
% time
nccreate(ancfile,'time','Dimensions',{'time',nt}, 'Datatype','single', 'Format','classic');
ncwrite(ancfile,'time',td);
ncwriteatt(ancfile,'time', 'units', 'seconds') ;
ncwriteatt(ancfile,'time', 'axis', 'time') ;
% topg
nccreate(ancfile,'topg','Dimensions',{'x',nx,'y',ny}, 'Datatype','single', 'Format','classic');
ncwrite(ancfile,'topg',SH);
ncwriteatt(ancfile,'topg', 'units', 'm') ;
% lat
nccreate(ancfile,'lat','Dimensions',{'x',nx,'y',ny}, 'Datatype','single', 'Format','classic');
ncwrite(ancfile,'lat',LAT);


%% DSMB
ancfile = ['../Data/MAR/TDSMB_MAR37_MIROC5_rcp85_05000m.nc'];
ncwrite_GrIS(ancfile,TDSMB,'aSMB',{'x','y','time'},5);
% time
nccreate(ancfile,'time','Dimensions',{'time',nt}, 'Datatype','single', 'Format','classic');
ncwrite(ancfile,'time',td);
ncwriteatt(ancfile,'time', 'units', 'seconds') ;
ncwriteatt(ancfile,'time', 'axis', 'time') ;
% topg
nccreate(ancfile,'topg','Dimensions',{'x',nx,'y',ny}, 'Datatype','single', 'Format','classic');
ncwrite(ancfile,'topg',SH);
ncwriteatt(ancfile,'topg', 'units', 'm') ;
% lat
nccreate(ancfile,'lat','Dimensions',{'x',nx,'y',ny}, 'Datatype','single', 'Format','classic');
ncwrite(ancfile,'lat',LAT);
ncwriteatt(ancfile,'topg', 'units', 'deg') ;

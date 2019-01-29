% Extract DSMB data
% Here schematic interpolation in time based on end of 19th and 20th century results 
clear

% resolution of output files in km
res = 1;

datapath = '/Volumes/ISMIP6/Data/Raw/SMB/MAR3.9';

gcm = 'MIROC5';
scen = 'rcp85';

%gcm = 'NorESM1';
%scen = 'rcp85';

% timer
%time = 1995:2005; % needs histo srcscen below

%time = 2006:2014;
time = 2015:2100;
nt = length(time);


%%%%%%%
%srcscen = 'histo';
%scenpath = [ datapath '/' gcm '-' srcscen '_1950_2005'];
%file_root = ['MARv3.9-yearly-' gcm '-' srcscen '-'];

scenpath = [ datapath '/' gcm '-' scen '_2006_2100'];
file_root = ['MARv3.9-yearly-' gcm '-' scen '-'];

outpath = '../Data/dSMB';

outfile_root_a = [ 'aSMB_MARv3.9-yearly-' gcm '-' scen ];
outfile_root_d = [ 'dSMBdz_MARv3.9-yearly-' gcm '-' scen ];
outfile_root_r = [ 'dRUNdz_MARv3.9-yearly-' gcm '-' scen ];

addpath('../toolbox')

% Load reference SMB
d0 = ncload(['../Data/MAR/MARv3.9-yearly-' gcm '-' scen '-ltm1995-2014.nc']);

%% Time loop, scenario from 2015-2100, hist from 1950-2005, present 2006-2014
%for t = 1:5 
for t = 1:nt 
    time(t)
    
    %% Load forcing file
    d1 = ncload([scenpath '/' file_root num2str(time(t)) '.nc']);

    %% aSMB [mmWE/yr] 
    aSMB = (d1.SMB(1:res:end,1:res:end)-d0.SMB(1:res:end,1:res:end));
    %% write out aSMB
    ancfile = [outpath '/' outfile_root_a  '-' num2str(time(t)) '.nc'];
    ncwrite_GrIS_aSMB(ancfile, aSMB, 'aSMB', {'x','y','t'}, res, time(t));

    %% dSMB/dz [mmWE/yr /m] 
    dSMBdz = (d1.dSMB(1:res:end,1:res:end));
    %% write out dSMBdz
    ancfile = [outpath '/' outfile_root_d  '-' num2str(time(t)) '.nc'];
    ncwrite_GrIS_dSMBdz(ancfile, dSMBdz, 'dSMBdz', {'x','y','t'}, res, time(t));

    %% dRUN/dz [mmWE/yr /m] 
    dRUNdz = (d1.dRU(1:res:end,1:res:end));
    %% write out dRUNdz
    ancfile = [outpath '/' outfile_root_r  '-' num2str(time(t)) '.nc'];
    ncwrite_GrIS_dRUNdz(ancfile, dRUNdz, 'dRUNdz', {'x','y','t'}, res, time(t));
    

end
%% end time loop


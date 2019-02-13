% Extract DSMB data

clear

% resolution of output files in km
res = 1;
secpyear = 31556926;

%datapath = '/Volumes/ISMIP6/Data/Raw/SMB/MAR3.9';
datapath = '/work/hgoelzer/Processing/RCM/MAR3.9';

%gcm = 'MIROC5';
%scen = 'rcp85';

gcm = 'NorESM1';
scen = 'rcp85';

% timer
time = 1950:2005; % needs histo srcscen below
%time = 1995:2005; % needs histo srcscen below

%time = 2006:2014;
%time = 2015:2100;
%time = 2006:2100;

nt = length(time);

% read days for time axis 
caldays = load('../Data/Grid/days_1900-2300.txt');

%%%%%%%
srcscen = 'histo';
scenpath = [ datapath '/' gcm '-' srcscen '_1950_2005'];
file_root = ['MARv3.9-yearly-' gcm '-' srcscen '-'];

%scenpath = [ datapath '/' gcm '-' scen '_2006_2100'];
%file_root = ['MARv3.9-yearly-' gcm '-' scen '-'];

outpath = ['../Data/dSMB/' gcm '-' scen ];
mkdir(outpath);
mkdir([outpath '/aSMB' ]);
mkdir([outpath '/dSMBdz' ]);
mkdir([outpath '/dRUNdz' ]);


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
    timestamp = caldays(time(t)-1900+1,3)
    time_bounds = [caldays(time(t)-1900+1,2), caldays(time(t)-1900+2,2)]
    %% Load forcing file
    d1 = ncload([scenpath '/' file_root num2str(time(t)) '.nc']);

    %% aSMB, convert [mmWE/yr] to [kg m-2 s-1], devinde by seconds-per-year
    aSMB = (d1.SMB(1:res:end,1:res:end)-d0.SMB(1:res:end,1:res:end))/secpyear;
    %% write out aSMB
    ancfile = [outpath '/aSMB/' outfile_root_a  '-' num2str(time(t)) '.nc'];
    ncwrite_GrIS_aSMB(ancfile, aSMB, 'aSMB', {'x','y','t'}, res, timestamp, time_bounds);

    %% dSMB/dz convert [mmWE/yr /m] to [kg m-2 s-1 m-1], devinde by seconds-per-year
    dSMBdz = (d1.dSMB(1:res:end,1:res:end))/secpyear;
    %% write out dSMBdz
    ancfile = [outpath '/dSMBdz/' outfile_root_d  '-' num2str(time(t)) '.nc'];
    ncwrite_GrIS_dSMBdz(ancfile, dSMBdz, 'dSMBdz', {'x','y','t'}, res, timestamp, time_bounds);

    %% dRUN/dz convert [mmWE/yr /m] to [kg m-2 s-1 m-1], devinde by seconds-per-year
    dRUNdz = (d1.dRU(1:res:end,1:res:end))/secpyear;
    %% write out dRUNdz
    ancfile = [outpath '/dRUNdz/' outfile_root_r  '-' num2str(time(t)) '.nc'];
    ncwrite_GrIS_dRUNdz(ancfile, dRUNdz, 'dRUNdz', {'x','y','t'}, res, timestamp, time_bounds);
    

end
%% end time loop


% Extract DSMB data

clear

% resolution of output files in km
res = 1;
secpyear = 31556926;

%datapath = '/Volumes/ISMIP6/Data/Raw/SMB/MAR3.9';
datapath = '/work/hgoelzer/Processing/RCM/MAR3.9';

gcm = 'CNRM-CM6';
scen = 'ssp585';
%scen = 'ssp126';

%gcm = 'UKESM1-CM6';
%scen = 'ssp585';

%gcm = 'CNRM-ESM2';
%scen = 'ssp585';

%gcm = 'CESM2';
%scen = 'ssp585';


% timer
time = 1950:2100; 

nt = length(time);

% read days for time axis 
caldays = load('../Data/Grid/days_1900-2300.txt');

%%%%%%%
outpath = ['../Data/dSMB/' gcm '-' scen ];
mkdir(outpath);
mkdir([outpath '/aSMB' ]);
mkdir([outpath '/dSMBdz' ]);


outfile_root_a = [ 'aSMB_MARv3.9-yearly-' gcm '-' scen ];
outfile_root_d = [ 'dSMBdz_MARv3.9-yearly-' gcm '-' scen ];

addpath('../toolbox')

% Load reference SMB
d0 = ncload(['../Data/MAR/MARv3.9-yearly-' gcm '-ltm1960-1989.nc']);

%% Time loop, scenario from 2015-2100, hist from 1950-2005, present 2006-2014
%for t = 1:5 
for t = 1:nt 
    time(t)
    if ( time(t) < 2015 )
      scenpath = [ datapath '/' gcm '-histo_1950_2014'];
      file_root = ['MARv3.9-yearly-' gcm '-histo-'];
    else
      scenpath = [ datapath '/' gcm '-' scen '_2015_2100'];
      file_root = ['MARv3.9-yearly-' gcm '-' scen '-'];
    end

    timestamp = caldays(time(t)-1900+1,3)
    time_bounds = [caldays(time(t)-1900+1,2), caldays(time(t)-1900+2,2)]
    %% Load forcing file
    d1 = ncload([scenpath '/' file_root num2str(time(t)) '.nc']);

    %% aSMB, convert [mmWE/yr] to [kg m-2 s-1], devinde by seconds-per-year
    aSMB = (d1.SMB(1:res:end,1:res:end)-d0.SMB(1:res:end,1:res:end))/secpyear;
    %% write out aSMB
    ancfile = [outpath '/aSMB/' outfile_root_a  '-' num2str(time(t)) '.nc'];
    ncwrite_GrIS_aSMB(ancfile, aSMB, 'aSMB', {'x','y','t'}, res, timestamp, time_bounds);

    %% dSMB/dz based on runoff convert [mmWE/yr /m] to [kg m-2 s-1 m-1], devinde by seconds-per-year
    dSMBdz = (-d1.dRU(1:res:end,1:res:end))/secpyear;
    %% write out dSMBdz
    ancfile = [outpath '/dSMBdz/' outfile_root_d  '-' num2str(time(t)) '.nc'];
    ncwrite_GrIS_dSMBdz(ancfile, dSMBdz, 'dSMBdz', {'x','y','t'}, res, timestamp, time_bounds);
    

end
%% end time loop


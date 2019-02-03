% Sum time evolving aSMB over the scenario

clear

addpath('../toolbox')

%% Settings

% Scenario
rcm = 'MARv3.9';
gcm = 'MIROC5';
scen = 'rcp85';

%gcm = 'NorESM1';
%scen = 'rcp85';

%%%%%%%

% scenario specific
aSMBpath = ['../Data/dSMB/' gcm '-' scen ];
aSMBfile_root = [ 'aSMB_MARv3.9-yearly-' gcm '-' scen ];

outpath = ['../Data/dSMB/' gcm '-' scen ];
outfile_root = [ 'aSMB_MARv3.9-yearly-' gcm '-' scen ];

mkdir(outpath, 'aSMB_mean'); 

secpyear = 31556926;

% read days for time axis 
caldays = load('../Data/Grid/days_1900-2300.txt');

% timer
time = 2015:2100;
nt = length(time);

aSMBsum = zeros(1681,2881);
count = 0;

msg = (['running year, 0000']);
fprintf(msg);
%for t=1:5 % year loop
for t=1:nt % year loop

    fprintf(['\b\b\b\b']);
    fprintf([sprintf('%02d',time(t))]);
    %% original aSMB
    d1 = ncload([aSMBpath '/aSMB/' aSMBfile_root  '-' num2str(time(t)) '.nc']);

    aSMBsum = aSMBsum + d1.aSMB(:,:);

end
%% end time loop
fprintf('\n')

% Unit in [kg m-2 s-1] 
aSMBmean = aSMBsum/nt;

%% write out aSMB
timestamp = caldays(floor(mean(time))-1900+1,3);
time_bounds = [caldays(time(1)-1900+1,2), caldays(time(end)-1900+2,2)];
ancfile = [outpath '/aSMB_mean/' outfile_root  '-mean.nc'];
ncwrite_GrIS_aSMB(ancfile, aSMBmean, 'aSMB', {'x','y','t'}, 1, timestamp, time_bounds);

% time dependent approximation for dsmb

clear

addpath('../toolbox')

%% Settings
rcm = 'MARv3.9';
gcm = 'MIROC5';
scen = 'rcp85';

%gcm = 'NorESM1';
%scen = 'rcp85';

%%%%%%%

% replace nan for high elevation with last finite value?
flg_nanfill = 1;

colors=get(0,'DefaultAxesColorOrder');

% basin definition
load ../Data/Basins/ExtBasinMasks25.mat
nb=length(bas.ids);

% area factors
da = ncload('../Data/Grid/af2_ISMIP6_GrIS_01000m.nc');
af = double(da.af2);
% res
dx=1000;dy=1000;

outpath = '../Data/dSMB';

% MAR surface 
d0 = ncload(['../Data/MAR/MARv3.9_topg_01000m.nc']);
sur = d0.topg;

% scenario specific 
outfile_root_d = [ 'dSMBdz_MARv3.9-yearly-' gcm '-' scen ];
lookup_file = ['trans_lookup_b25_MARv3.9-' gcm '-' scen ];

% timer
time = 2015:2100;
nt = length(time);

ds = ncload('../Models/OBS/sftgif_01000m.nc');
mask = double(ds.sftgif);

%% fit a lookup table to the data
%% centers at sstep intervals with ds range
sstep = 100;
ds = 100;
ss=0:sstep:3500;
ns=length(ss);

table=zeros([nb,ns,nt]);
bint=zeros(nb,nt);

msg = (['running year, basin: 00,00']);
fprintf(msg);
%for t = 1 % year loop
%for t = 1:5 % year loop
for t = 1:nt % year loop

%    t
    fprintf(['\b\b\b\b\b']);
    fprintf([sprintf('%02d',t), ',00']);
    d1 = ncload([outpath '/' outfile_root_d  '-' num2str(time(t)) '.nc']);
    dSMBdz = d1.dSMBdz;
    
%    figure
    for b = 1:nb
        
        fprintf(['\b\b\b']);
        fprintf([',' sprintf('%02d',b)]);
        eval(['dSMBdz_b=dSMBdz.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
        eval(['sur_b=sur.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
        
%        subplot(5,4,b)
%        hold on; box on;
%        plot(sur_b(:),dSMBdz_b(:),'.');

        %% integral dsmb for this basin
        bint(b,t)=nansum(nansum(dSMBdz_b.*af.*mask))*dx*dy;
        
        %% fit a lookup table to the data
        %% centers at sstep intervals with ds range
        look = zeros(2,ns); 
        n=1;
        yold=0.;
        for s0=ss(2:end)
            n=n+1;
            %% find local average DSMB for given elevation range
%            [y1,ysel,s1,ssel]= find_local_average(sur_b,dSMBdz_b,s0,ds);
            [y1,ysel,s1,ssel]= find_local_median(sur_b,dSMBdz_b,s0,ds);
            %% fill nans with last value
            if(isnan(y1) && flg_nanfill)
                look(:,n)=[s0,yold];
%                plot(s0,yold,'o','Color',colors(n,:))
            else
                look(:,n)=[s0,y1];
                yold=y1;
%                plot(s0,y1,'o','Color',colors(n,:))
            end
%            plot(ssel,ysel,'.','Color',colors(n,:))
        end
        %% fill first value for h=0 from second
        look(2,1) = look(2,2);
        if(b==12)
            look(2,29:end) = look(2,28);
        end

        table(b,:,t)=look(2,:);

%        plot(look(1,:),look(2,:),'-k')
%        axis([0 3300 -5 2])
%        axis([0 3300 -2 1])
%        title(['B' num2str(bas.ids{b}) ' t:' num2str(t)])
        
    end % end basin loop
    %% manual correction for some basins needed?
    if(t==nt)
        warning('Manual corrections active, check !!!'); 
    end
    table(9,32:end,t) = table(9,31,t);
    table(7,35:end,t) = table(7,34,t);

end % end year loop

%% Write netcdf
nz = length(ss);
td = (1:86)*31556926;
nt = length(td);

% permute to get dsmb_table(h,b,t)
table_out = permute(table,[2,1,3]);

% Write netcdf
ancfile = ['../Data/lookup/TdSMBdz_' lookup_file '.nc'];
ncwrite_TDSMBTable(ancfile,ss,td,table_out,'dSMBdz_ltbl',{'z','b','time'});
nccreate(ancfile,'bint','Dimensions',{'b',nb,'time',nt}, 'Datatype','single', 'Format','classic');
ncwrite(ancfile,'bint',bint);

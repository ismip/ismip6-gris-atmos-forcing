% time dependent approximation for dsmb

clear

addpath('../toolbox')

%%%% ! Note special treatment of some basins below !!!

% fill with nan for high elevation?
flg_nanfill = 1;

colors=get(0,'DefaultAxesColorOrder');

% basin definition
load ../Data/Basins/ExtBasinMasks25.mat
nb=length(bas.ids);

% area factors
load ../Data/Grid/af_e05000m.mat af2
% dim
dx=5000;dy=5000;

% masks
obs=ncload('../Data/Grid/obs1_05.nc');

iscen = 5;

if (iscen ==5)
d0 = ncload('../Data/MAR/TDSMB_MAR37_MIROC5_rcp85_05000m.nc');
lookup_file='trans_lookup_MAR37_b25';
end

sur=d0.topg;
mask=obs.zmask;
lat=d0.lat;
time=d0.time;
nt=length(time);

%% fit a lookup table to the data
%% centers at sstep intervals with ds range
sstep = 100;
ds = 100;
ss=0:sstep:3500;
ns=length(ss);

table=zeros([nb,ns,nt]);
bint=zeros(nb,nt);

%for t=1 % year loop
%for t=1:5 % year loop
for t=1:nt % year loop
           

    dsd=d0.aSMB(:,:,t).*(mask./double(mask));

%    figure
    for b=1:nb

        eval(['dsd_b=dsd.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
        eval(['sur_b=sur.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
        
%        subplot(5,4,b)
%        hold on; box on;
%        plot(sur_b(:),dsd_b(:),'.');

        %% integral dsmb for this basin
        bint(b,t)=nansum(nansum(dsd_b.*af2))*dx*dy;
        
        %% fit a lookup table to the data
        %% centers at sstep intervals with ds range
        look = zeros(2,ns); 
        n=1;
        yold=0.;
        for s0=ss(2:end)
            n=n+1;
            %% find local average DSMB for given elevation range
%            [y1,ysel,s1,ssel]= find_local_average(sur_b,dsd_b,s0,ds);
            [y1,ysel,s1,ssel]= find_local_median(sur_b,dsd_b,s0,ds);
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
%    table(12,2,29:end) = table(12,2,28);

end % end year loop

lookup.table = table;

%% add basin totals to lookup
lookup.bint = bint;
lookup.ss = ss;
lookup.ds = ds;
lookup.time = time;

%% write lookup table for model/scenario
%save(lookup_file,'lookup')

%% Write netcdf
nz = length(ss);
td = (1:101)*31556926;
nt = length(td);

% permute to get dsmb_table(h,b,t)
table_out = permute(table,[2,1,3]);

ncwrite_TDSMBTable(['../Data/lookup/TDSMB_' lookup_file '.nc'],ss,td,table_out,'aSMB_ltbl',{'z','b','time'});

% Calculate time evolving DSMB at a given initial geometry

clear

addpath('../toolbox')

%% Settings

% Scenario
rcm = 'MARv3.9';
gcm = 'MIROC5';
scen = 'rcp85';

%gcm = 'NorESM1';
%scen = 'rcp85';

% Model
%amod = 'OBS';
amod = 'BISICLES1';
%amod = 'IMAUICE08';

%%%%%%%

% Parameters
flg_weigh = 1;

% flag for plotting 
flg_plot=0;
load cmap_dsmb

% scenario specific
aSMBpath = ['../Data/dSMB/' gcm '-' scen ];
aSMBfile_root = [ 'aSMB_MARv3.9-yearly-' gcm '-' scen ];
outpath = ['../Models/' amod '/' gcm '-' scen ];
outfile_root = [ 'aSMB_MARv3.9-yearly-' gcm '-' scen '-' amod ];

mkdir(outpath, 'aSMB'); 

secpyear = 31556926;

% basin definition
load(['../Data/Basins/ExtBasinMasks25.mat']);
x1 = 1:size(bas.basinIDs,1);
y1 = 1:size(bas.basinIDs,2);
nb = length(bas.ids);
[y,x] = meshgrid(y1,x1);

% area factors
da = ncload('../Data/Grid/af2_ISMIP6_GrIS_01000m.nc');
af = double(da.af2);
% dim
dx=1000;dy=1000;

% basin weights
load(['../Data/Basins/ExtBasinScale25_nn7_50.mat'], 'wbas');

% original forcing
lookup = ncload(['../Data/lookup/TaSMB_trans_lookup_b25_MARv3.9-' gcm '-' scen '.nc']);
modscen='MAR39';

% dummy lookup for zero
dummy0 = lookup.aSMB_ltbl(:,1,1);

% Load a modelled geometry for reconstruction
nc=ncload(['../Models/' amod '/orog_01000m.nc']);
nc1=ncload(['../Models/' amod '/sftgif_01000m.nc']);

% Load observed geometry 
nco1=ncload(['../Models/OBS/sftgif_01000m.nc']);

sur = max(0,double(nc.orog));

% Masks
ima_mod = double(nc1.sftgif);
ima_obs = double(nco1.sftgif);

nt=length(lookup.time);
time = lookup.time;

bint_obs=zeros(size(lookup.bint));
bint_ext=zeros(size(lookup.bint));
bint_map=zeros(size(lookup.bint));

msg = (['running year, basin: 00,00']);
fprintf(msg);
%for t=1:5 % year loop
for t=1:nt % year loop

    timestamp = (time(t)-1900)*secpyear;

    fprintf(['\b\b\b\b\b']);
    fprintf([sprintf('%02d',t), ',00']);
    %% original aSMB
    d1 = ncload([aSMBpath '/aSMB/' aSMBfile_root  '-' num2str(time(t)) '.nc']);

    aSMB = d1.aSMB(:,:);
    aSMB_re = zeros(size(aSMB));
    bint_o = zeros(1,nb);
    bint_e = zeros(1,nb);
    bint_m = zeros(1,nb);

    %% loop through basins
    for b=1:nb

        fprintf(['\b\b\b']);
        fprintf([',' sprintf('%02d',b)]);
        %% set current basin and lookup
        eval(['sur_b=sur.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
        eval(['mask_b =       (bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
        eval(['ima_b=ima_mod.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);

        %% set neighbor basin and lookup
        look0 = dummy0;
        if (wbas.n0(b)>0)
            look0=lookup.aSMB_ltbl(:,wbas.n0(b),t);
        end
        look1 = dummy0;
        if (wbas.n1(b)>0)
            look1=lookup.aSMB_ltbl(:,wbas.n1(b),t);
        end
        look2 = dummy0;
        if (wbas.n2(b)>0)
            look2=lookup.aSMB_ltbl(:,wbas.n2(b),t);
        end
        look3 = dummy0;
        if (wbas.n3(b)>0)
            look3=lookup.aSMB_ltbl(:,wbas.n3(b),t);
        end
        look4 = dummy0;
        if (wbas.n4(b)>0)
            look4=lookup.aSMB_ltbl(:,wbas.n4(b),t);
        end
        look5 = dummy0;
        if (wbas.n5(b)>0)
            look5=lookup.aSMB_ltbl(:,wbas.n5(b),t);
        end
        look6 = dummy0;
        if (wbas.n6(b)>0)
            look6=lookup.aSMB_ltbl(:,wbas.n6(b),t);
        end
        
        %% use lookup table to determine DSMB
        aSMB_b0 = interp1(lookup.z,look0(:),sur_b);
        aSMB_b1 = interp1(lookup.z,look1(:),sur_b);
        aSMB_b2 = interp1(lookup.z,look2(:),sur_b);
        aSMB_b3 = interp1(lookup.z,look3(:),sur_b);
        aSMB_b4 = interp1(lookup.z,look4(:),sur_b);
        aSMB_b5 = interp1(lookup.z,look5(:),sur_b);
        aSMB_b6 = interp1(lookup.z,look6(:),sur_b);

        if (flg_weigh == 0)
            %% combine according to weights
            aSMB_b = aSMB_b0.*wbas.wg;
        else
            aSMB_b = aSMB_b0.*wbas.wgc0 + aSMB_b1.*wbas.wgc1 + aSMB_b2.*wbas.wgc2 + aSMB_b3.*wbas.wgc3 + aSMB_b4.*wbas.wgc4 + aSMB_b5.*wbas.wgc5 + aSMB_b6.*wbas.wgc6;
        end
%    shade(aSMB_b)

        %% replace nan by zeros to add all basins together
        aSMB_b(isnan(aSMB_b)) = 0;
        aSMB_re = aSMB_re+aSMB_b;
        %% integral remapped aSMB for this basin
        bint_m(b) = nansum(nansum(aSMB_b.*ima_b.*af))*dx*dy;

        %% integral extended aSMB for this basin
        bint_e(b) = nansum(nansum(aSMB.*ima_b.*af))*dx*dy;

        %% integral observed aSMB for this basin
        bint_o(b) = nansum(nansum(aSMB.*mask_b.*ima_obs.*af))*dx*dy;
        
    end
    %% end basin loop

    %% collect results
    bint_obs(:,t) = bint_o(:);
    bint_ext(:,t) = bint_e(:);
    bint_map(:,t) = bint_m(:);

    %% aSMB [kg m-2 s-1] 
    %% write out aSMB
    ancfile = [outpath '/aSMB/' outfile_root  '-' num2str(time(t)) '.nc'];
    ncwrite_GrIS_aSMB(ancfile, aSMB_re, 'aSMB', {'x','y','t'}, 1, timestamp);

    if (flg_plot) 
        shade_bg(aSMB_re)
        colormap(cmap)
        caxis([-4,1])
        print('-dpng', '-r300', ['dsmb_' modscen '_re' sprintf('%02d',t)]) 
        close
        shade_bg(aSMB)
        colormap(cmap)
        caxis([-4,1])
        print('-dpng', '-r300', ['dsmb_' modscen '_or' sprintf('%02d',t)]) 
        close
    end


end
%% end time loop

save(['../Models/' amod '/biastest_' gcm '-' scen '-' amod ], 'bint_obs', 'bint_ext', 'bint_map');

% Plot
if (flg_plot) 
figure
bar([mean(bint_obs,2), mean(bint_ext,2), mean(bint_map,2)]*31556926/1e12)
axis tight
ylabel('Integrated aSMB [Gt yr-1]')
legend({'observed', 'extended', 'remapped'},'Location','southeast')
xlabel('Basin Id')
print('-dpng', '-r300', ['../Models/' amod '/dsmb_basinint_' gcm '-' scen '_sum']) 

figure
bar([mean(bint_ext,2)-mean(bint_obs,2), mean(bint_map,2)-mean(bint_obs,2)]*31556926/1e12)
axis tight
ylabel('Integrated aSMB biases [Gt yr-1]')
legend({'extended', 'remapped'},'Location','southeast')
xlabel('Basin Id')
print('-dpng', '-r300', ['../Models/' amod '/dsmb_basinint_' gcm '-' scen '_diff']) 
end

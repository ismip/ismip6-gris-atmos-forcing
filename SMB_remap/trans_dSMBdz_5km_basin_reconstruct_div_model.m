% Calculate time evolving dSMBdz at a given initial geometry

clear

addpath('../toolbox')

% Parameters
flg_weigh = 1;

% 0=initMIP; 1=MIROC8.5; 2=NorESM8.5; 3=CANSM8.5; 4=MIROC4.5; 5=new MIROC8.5
iscen = 5;

% flag for plotting 
flg_plot=0;
load cmap_dsmb

colors=get(0,'DefaultAxesColorOrder');

% basin definition
load ../Data/Basins/ExtBasinMasks25.mat
x1 = 1:size(bas.basinIDs,1);
y1 = 1:size(bas.basinIDs,2);
nb = length(bas.ids);
[y,x] = meshgrid(y1,x1);

% area factors
load ../Data/Grid/af_e05000m.mat af2
% dim
dx=5000;dy=5000;

% masks
obs=ncload('../Data/Grid/obs1_05.nc');

% basin weights
load ../Data/Basins/ExtBasinScale25_div6_50 wbas

% Forcing scenario
if (iscen ==5)
    % original forcing
    d0 = ncload('../Data/MAR/TdSMBdz_MAR37_MIROC5_rcp85_05000m.nc');
    lookup = ncload('../Data/lookup/TdSMBdz_trans_lookup_MAR37_b25.nc');
    modscen='MAR37';
end

% Load a modelled geometry for reconstruction
amod = 'VUBGISM';
nc=ncload(['../Models/' amod '/lithk_05000m.nc']);
nc1=ncload(['../Models/' amod '/topg_05000m.nc']);
nc2=ncload(['../Models/' amod '/sftgif_05000m.nc']);

% Operate on ice thickness
%thi = nc.lithk.*nc2.sftgif;
thi = nc.lithk; thi(isnan(thi))=0;
bed = nc1.topg;
sur = max(0,bed+thi);

ima=nc2.sftgif;
mask=nc2.sftgif;

lat=d0.lat;
nt=length(lookup.time);

bint_re=zeros(nb,nt);

% output array
tdsmb_re=zeros(size(d0.dSMBdz));

%for t=1:5 % year loop
for t=1:nt % year loop

    dsd=d0.dSMBdz(:,:,t);
    dsd_re=zeros(size(dsd));

    %% loop through basins
    for b=1:nb
        %% set current basin and lookup
        eval(['sur_b=sur.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
        eval(['ima_b=ima.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);

        look0=lookup.dSMBdz_ltbl(:,wbas.n0(b),t);
        %% set neighbor basin and lookup
        look1=lookup.dSMBdz_ltbl(:,wbas.n1(b),t);
        look2=lookup.dSMBdz_ltbl(:,wbas.n2(b),t);
        look3=lookup.dSMBdz_ltbl(:,wbas.n3(b),t);
        look4=lookup.dSMBdz_ltbl(:,wbas.n4(b),t);
        look5=lookup.dSMBdz_ltbl(:,wbas.n5(b),t);
        look6=lookup.dSMBdz_ltbl(:,wbas.n6(b),t);
        
        %% use lookup table to determine DSMB
        dsd_b0 = interp1(lookup.z,look0(:),sur_b);
        dsd_b1 = interp1(lookup.z,look1(:),sur_b);
        dsd_b2 = interp1(lookup.z,look2(:),sur_b);
        dsd_b3 = interp1(lookup.z,look3(:),sur_b);
        dsd_b4 = interp1(lookup.z,look4(:),sur_b);
        dsd_b5 = interp1(lookup.z,look5(:),sur_b);
        dsd_b6 = interp1(lookup.z,look6(:),sur_b);

        if (flg_weigh == 0)
            %% combine according to weights
            dsd_b = dsd_b0.*wbas.wg;
        else
            dsd_b = dsd_b0.*wbas.wgc0 + dsd_b1.*wbas.wgc1 + dsd_b2.*wbas.wgc2 + dsd_b3.*wbas.wgc3 + dsd_b4.*wbas.wgc4 + dsd_b5.*wbas.wgc5 + dsd_b6.*wbas.wgc6;
        end
%    shade(dsd_b)

        %% extended integral dsmb for this basin
        dsd_ex = dsd.*ima_b;
        bint_ex(b)=nansum(nansum(dsd_ex.*af2))*dx*dy;
        
        %% mask
        %dsd_b = dsd_b.*mask;

        %% integral dsmb for this basin
        bint_re(b)=nansum(nansum(dsd_b.*af2))*dx*dy;

        %% check integral again 
        bint_out(b)=nansum(nansum(dsd_b.*af2))*dx*dy;

        %% replace nan by zeros to add all basins together
        dsd_b(isnan(dsd_b))=0;
        dsd_re = dsd_re+dsd_b;

    end
    %% end basin loop

    %% collect results
    tdsmb_re(:,:,t) = dsd_re;    


    if (flg_plot) 
        shade_bg(dsd_re)
        colormap(cmap)
        caxis([-4,1])
        print('-dpng', '-r300', ['dsmb_' modscen '_re' sprintf('%02d',t)]) 
        close
%        shade_bg(dsd.*mask)
        shade_bg(dsd)
        colormap(cmap)
        caxis([-4,1])
        print('-dpng', '-r300', ['dsmb_' modscen '_or' sprintf('%02d',t)]) 
        close
    end


end
%% end time loop

ncwrite_GrIS(['../Models/' amod '/TdSMBdz_' modscen '_' amod '.nc'], tdsmb_re, 'TdSMBdz',{'x','y','time'},5)

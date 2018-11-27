% plot dRUN/dz
clear

addpath('../toolbox')

%%% new MAR data
pdata='../Data/MAR/';

d19 = ncload([pdata 'MARv3.7-yearly-MIROC5-19xx.nc']);
d20 = ncload([pdata 'MARv3.7-yearly-MIROC5-20xx.nc']);
dSMBdz19=-d19.dRU(1:5:end,1:5:end)/1000;
dSMBdz20=-d20.dRU(1:5:end,1:5:end)/1000;

obs2 = ncload('../Data/obs1_05.nc');
zma = obs2.zmask;
cma = (obs2.topg>0 | obs2.lithk>10);


shade_nt(dSMBdz19*1000)         
hold on
load cmap_polar.mat
caxis([-5,5])
colormap(cmap)
title('dRUN/dz 1990')
contour(cma',[0.5,0.5],'Color',[0.7,0.7,0.7],'Linewidth',1.5)
contour(zma',[0.5,0.5],'Linewidth',1,'Color','k')
print -dpng -r300 dRUNdz_1990


shade_nt(dSMBdz20*1000)      
hold on
print -dpng -r300 dRUNdz_1990
caxis([-5,5])
colormap(cmap)
title('dRUN/dz 2090')
contour(cma',[0.5,0.5],'Color',[0.7,0.7,0.7],'Linewidth',1.5)
contour(zma',[0.5,0.5],'Linewidth',1,'Color','k')
print -dpng -r300 dRUNdz_2090



shade_nt(dSMBdz20*1000-dSMBdz19*1000)
hold on
caxis([-5,5])
colormap(cmap)
title('dRUN/dz difference 2090-1990')
contour(cma',[0.5,0.5],'Color',[0.7,0.7,0.7],'Linewidth',1.5)
contour(zma',[0.5,0.5],'Linewidth',1,'Color','k')
print -dpng -r300 dRUNdz_2090-1990

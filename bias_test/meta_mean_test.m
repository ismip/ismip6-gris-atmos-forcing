% Run bias test for a number of models

clear

% initMIP labs and models  
labs={'ARC','AWI','AWI','BGC','BGC','BGC','DMI','DMI','DMI','DMI','DMI','IGE','IGE','ILTS','ILTSPIK','IMAU','IMAU','IMAU','JPL','LANL','LSCE','MIROC','MIROC','MPIM','UAF','UAF','UAF','UAF','UAF','UAF','UCIJPL','ULB','ULB','VUB','VUB'};
models={'PISM','ISSM1','ISSM2','BISICLES1','BISICLES2','BISICLES3','PISM1','PISM2','PISM3','PISM4','PISM5','ELMER1','ELMER2','SICOPOLIS','SICOPOLIS','IMAUICE1','IMAUICE2','IMAUICE3','ISSM','CISM','GRISLI','ICIES1','ICIES2','PISM','PISM1','PISM2','PISM3','PISM4','PISM5','PISM6','ISSM','FETISH1','FETISH2','GISM1','GISM2'};

% basic testing
ll = length(labs);
if (ll ~= length(models))
    error('labs and models vectors must have the same length')
end

% Results container
bints = zeros(length(labs),3,25);

% model loop
for m = 1:length(labs)
    amod = [labs{m}, '_' models{m}];
    mean_test_aSMB1;
    bints(m,1,:) = bint_obs;
    bints(m,2,:) = bint_ext;
    bints(m,3,:) = bint_map;
end

save('meta_bints','bints','labs','models');

% Plotting
bar(mean(bints(:,:,:),3))
set(gca,'Xtick',1:length(models))
set(gca,'xticklabels',models,'XTickLabelRotation',90)
axis tight
ax = axis;
axis([ax(1),ax(2),ax(3),0])
print('-dpng' ,'-r300','meta_mean_bias')

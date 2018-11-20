% Save basins in useful mask format 
clear

bas=ncload('../Data/Basins/ISMIP6Masks25_05000m.nc');

bas.ids = 1:25;
bas.basinIDs=zeros(size(bas.basin1));
check=zeros(size(bas.basin1));
for i=1:25
    eval(['bas.basin' num2str(i) '(isnan(bas.basin' num2str(i) '))=0.;']);
    eval(['bas.basinIDs=bas.basinIDs+bas.basin' num2str(i) '*i;']);
    eval(['check = check + bas.basin' num2str(i) ';']);
end
%shade(check)

x1=1:size(bas.basinIDs,1);
y1=1:size(bas.basinIDs,2);
[y,x] = meshgrid(y1,x1);

for i=1:25
    %% mean basin position for label
    eval(['xcb=x.*(bas.basin' num2str(i) './bas.basin' num2str(i) ');']);
    eval(['ycb=y.*(bas.basin' num2str(i) './bas.basin' num2str(i) ');']);
    bc(2,i)=nanmedian(ycb(:));
    bc(1,i)=nanmedian(xcb(:));
end

% make a zero basin id
bas.basinIDs(1,1)=0;
bas.basin0=bas.basin1*0;
bas.basin0(1,1) = 1;

% Plotting
shade(bas.basinIDs)
colormap(lines(25))
caxis([1,25])
text(bc(1,:),bc(2,:),num2str([1:25]'),'Color',[1,1,1],'Fontsize',12)

save ../Data/Basins/ExtBasinMasks25 bas

save ../Data/Basins/bc bc

print -dpng -r300 extbasins

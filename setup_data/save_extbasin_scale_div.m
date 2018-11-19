% save basin connection info
% For given basin, give neighbors and scaling factors according to distance
% include weighting across divide, other basins
clear

% Params
% falloff distance and resolution in km
foff = 50;
res = 5;

% load basin division
load ../Data/Basins/ExtBasinMasks25 
[y,x]= meshgrid(1:size(bas.basin1,2),1:size(bas.basin1,1));

% basin id
bas.ids = 1:25;
bas.n0 =  [ 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25];
% locking to the coast left neighbor
bas.n1 =  [25, 1:24];
% locking to the coast right neighbor
bas.n2 =  [ 2:25, 1];
% Any other basins of interest? 
bas.n3 =  [24  4 13  2 22 13 19 18 18 17 16 15 25 12 12 11  9  8  7  7  5  5  4  1 13];
% Any other basins of interest? 
bas.n4 =  [ 4 24 13 24 23 13  5 17 17 16 15 14 25 25 11 10 10  9  5  5  7  7  5  2 13];
% Any other basins of interest? 
bas.n5 =  [23 23 13 23 21 13 18 13 13 13 13 25 25 25 25 25  8  7 21 22 19 19  2  4 13];
% Any other basins of interest? 
bas.n6 =  [13 13 13  1  7 13 20 13 13 13 13 25 25 25 25 25 25 25 25 25 25  4  1 13 13];

bas.basn = [bas.n0; bas.n1; bas.n2; bas.n3; bas.n4; bas.n5; bas.n6];
% number of neighbors including basin itself
bas.nn = 7;

%% default weights: full weight on own basin
bas.wgn0 = ones(size(bas.basin1));
bas.wgn1 = zeros(size(bas.basin1));
bas.wgn2 = zeros(size(bas.basin1));
bas.wgn3 = zeros(size(bas.basin1));
bas.wgn4 = zeros(size(bas.basin1));
bas.wgn5 = zeros(size(bas.basin1));
bas.wgn6 = zeros(size(bas.basin1));

bas.wgc0 = ones(size(bas.basin1));
bas.wgc1 = zeros(size(bas.basin1));
bas.wgc2 = zeros(size(bas.basin1));
bas.wgc3 = zeros(size(bas.basin1));
bas.wgc4 = zeros(size(bas.basin1));
bas.wgc5 = zeros(size(bas.basin1));
bas.wgc6 = zeros(size(bas.basin1));

% determine weights from distance to closest neighbor for each point in basin
% calculate shortest distance to n1, n2, n3
%for i=[1:10]
%for i=[1:1]
for i=[1:25]
i
    eval(['bin0 = find(bas.basin' num2str(bas.n0(i)) '>0);']);
    eval(['bin1 = find(bas.basin' num2str(bas.n1(i)) '>0);']);
    eval(['bin2 = find(bas.basin' num2str(bas.n2(i)) '>0);']);
    eval(['bin3 = find(bas.basin' num2str(bas.n3(i)) '>0);']);
    eval(['bin4 = find(bas.basin' num2str(bas.n4(i)) '>0);']);
    eval(['bin5 = find(bas.basin' num2str(bas.n5(i)) '>0);']);
    eval(['bin6 = find(bas.basin' num2str(bas.n6(i)) '>0);']);

    %% big matrix of differences    
    xd = repmat(x(bin1),[1,length(bin0)])'-repmat(x(bin0),[1,length(bin1)]);
    yd = repmat(y(bin1),[1,length(bin0)])'-repmat(y(bin0),[1,length(bin1)]);
    vd = sqrt(xd.^2+yd.^2);
    sd1 = min(vd,[],2)*res;

    xd = repmat(x(bin2),[1,length(bin0)])'-repmat(x(bin0),[1,length(bin2)]);
    yd = repmat(y(bin2),[1,length(bin0)])'-repmat(y(bin0),[1,length(bin2)]);
    vd = sqrt(xd.^2+yd.^2);
    sd2 = min(vd,[],2)*res;

    xd = repmat(x(bin3),[1,length(bin0)])'-repmat(x(bin0),[1,length(bin3)]);
    yd = repmat(y(bin3),[1,length(bin0)])'-repmat(y(bin0),[1,length(bin3)]);
    vd = sqrt(xd.^2+yd.^2);
    sd3 = min(vd,[],2)*res;

    xd = repmat(x(bin4),[1,length(bin0)])'-repmat(x(bin0),[1,length(bin4)]);
    yd = repmat(y(bin4),[1,length(bin0)])'-repmat(y(bin0),[1,length(bin4)]);
    vd = sqrt(xd.^2+yd.^2);
    sd4 = min(vd,[],2)*res;

    xd = repmat(x(bin5),[1,length(bin0)])'-repmat(x(bin0),[1,length(bin5)]);
    yd = repmat(y(bin5),[1,length(bin0)])'-repmat(y(bin0),[1,length(bin5)]);
    vd = sqrt(xd.^2+yd.^2);
    sd5 = min(vd,[],2)*res;

    xd = repmat(x(bin6),[1,length(bin0)])'-repmat(x(bin0),[1,length(bin6)]);
    yd = repmat(y(bin6),[1,length(bin0)])'-repmat(y(bin0),[1,length(bin6)]);
    vd = sqrt(xd.^2+yd.^2);
    sd6 = min(vd,[],2)*res;

    %% weight decreases with increasing distance from region 
    bas.wgn1(bin0)=1-min((sd1)/foff,1);
    bas.wgn2(bin0)=1-min((sd2)/foff,1);
    bas.wgn3(bin0)=1-min((sd3)/foff,1);
    bas.wgn4(bin0)=1-min((sd4)/foff,1);
    bas.wgn5(bin0)=1-min((sd5)/foff,1);
    bas.wgn6(bin0)=1-min((sd6)/foff,1);
    %% combined weights
    bas.wgc0(bin0)=1-(bas.wgn1(bin0)+bas.wgn2(bin0)+bas.wgn3(bin0)+bas.wgn4(bin0)+bas.wgn5(bin0)+bas.wgn6(bin0))./(bas.wgn0(bin0)+bas.wgn1(bin0)+bas.wgn2(bin0)+bas.wgn3(bin0)+bas.wgn4(bin0)+bas.wgn5(bin0)+bas.wgn6(bin0));
    bas.wgc1(bin0)=(bas.wgn1(bin0))./(bas.wgn0(bin0)+bas.wgn1(bin0)+bas.wgn2(bin0)+bas.wgn3(bin0)+bas.wgn4(bin0)+bas.wgn5(bin0)+bas.wgn6(bin0));
    bas.wgc2(bin0)=(bas.wgn2(bin0))./(bas.wgn0(bin0)+bas.wgn1(bin0)+bas.wgn2(bin0)+bas.wgn3(bin0)+bas.wgn4(bin0)+bas.wgn5(bin0)+bas.wgn6(bin0));
    bas.wgc3(bin0)=(bas.wgn3(bin0))./(bas.wgn0(bin0)+bas.wgn1(bin0)+bas.wgn2(bin0)+bas.wgn3(bin0)+bas.wgn4(bin0)+bas.wgn5(bin0)+bas.wgn6(bin0));
    bas.wgc4(bin0)=(bas.wgn4(bin0))./(bas.wgn0(bin0)+bas.wgn1(bin0)+bas.wgn2(bin0)+bas.wgn3(bin0)+bas.wgn4(bin0)+bas.wgn5(bin0)+bas.wgn6(bin0));
    bas.wgc5(bin0)=(bas.wgn5(bin0))./(bas.wgn0(bin0)+bas.wgn1(bin0)+bas.wgn2(bin0)+bas.wgn3(bin0)+bas.wgn4(bin0)+bas.wgn5(bin0)+bas.wgn6(bin0));
    bas.wgc6(bin0)=(bas.wgn6(bin0))./(bas.wgn0(bin0)+bas.wgn1(bin0)+bas.wgn2(bin0)+bas.wgn3(bin0)+bas.wgn4(bin0)+bas.wgn5(bin0)+bas.wgn6(bin0));

end

% Group in one matrix
bas.wbasWGTs = zeros(size(bas.basinIDs,1),size(bas.basinIDs,2),bas.nn); 
bas.wbasWGTs(:,:,1) = bas.wgc0;
bas.wbasWGTs(:,:,2) = bas.wgc1;
bas.wbasWGTs(:,:,3) = bas.wgc2;
bas.wbasWGTs(:,:,4) = bas.wgc3;
bas.wbasWGTs(:,:,5) = bas.wgc4;
bas.wbasWGTs(:,:,6) = bas.wgc5;
bas.wbasWGTs(:,:,7) = bas.wgc6;

% Write out 
wbas = bas;
save(['../Data/Basins/ExtBasinScale25_div6_' num2str(foff)], 'wbas');

%%% Write in netcdf
% BasinIDs
wbas.wbasIDs = zeros(size(wbas.basinIDs,1),size(wbas.basinIDs,2),wbas.nn); 
%wbasIDs(:,:,1) = basinIDs;
for i=1:size(wbas.basinIDs,1)
    for j=1:size(wbas.basinIDs,2)
        in = wbas.basinIDs(i,j);
        for n = 1:wbas.nn
            wbas.wbasIDs(i,j,n) = wbas.basn(n,in);
        end
    end
end
ncwrite_GrIS('../Data/Basins/BasinIDs_nb25_nn7_05000m.nc',wbas.wbasIDs,'BasinIDs',{'x','y','n'},5);
ncwrite_GrIS('../Data/Basins/BasinWGTs_nb25_nn7_05000m.nc',wbas.wbasWGTs,'BasinWGTs',{'x','y','n'},5);


% Plotting
shade(bas.wgc0)
caxis([0 1])
print -dpng -r300 scale_div_wgc0

shade(bas.wgc1+bas.wgc2+bas.wgc3+bas.wgc4+bas.wgc5+bas.wgc6)
caxis([0 1])
print -dpng -r300 scale_div_wgc1-6

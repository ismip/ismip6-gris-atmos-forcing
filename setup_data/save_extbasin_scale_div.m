% Calculate basin weights
% For given basin, give neighbors and scaling factors according to distance
% include weighting across divide, other basins
clear

addpath('../toolbox')

% Params
% falloff distance and resolution in km
foff = 50;
res = 5;

% load basin division
load(['../Data/Basins/ExtBasinNeighbours25_nn7'], 'bas');
[y01,x01]= meshgrid(1:size(bas.basin1,2),1:size(bas.basin1,1));
[y05,x05]= meshgrid(1:5:size(bas.basin1,2),1:5:size(bas.basin1,1));

%% subsample to 5 km, too big operation on 1 km
bas05.ids = bas.ids; 
bas05.basn = bas.basn; 
for n = 0:6
    eval(['bas05.n' num2str(n) '= bas.n' num2str(n) ';']);
end
bas05.nn = bas.nn; 
bas05.basinIDs = bas.basinIDs(1:5:end,1:5:end);
for b = 0:25
    eval(['bas05.basin' num2str(b) '= bas.basin' num2str(b) '(1:5:end,1:5:end);']);
end
[y,x]= meshgrid(1:size(bas05.basin1,2),1:size(bas05.basin1,1));

%% default weights: full weight on own basin
bas05.wgn0 =  ones(size(bas05.basin1));
bas05.wgn1 = zeros(size(bas05.basin1));
bas05.wgn2 = zeros(size(bas05.basin1));
bas05.wgn3 = zeros(size(bas05.basin1));
bas05.wgn4 = zeros(size(bas05.basin1));
bas05.wgn5 = zeros(size(bas05.basin1));
bas05.wgn6 = zeros(size(bas05.basin1));

bas05.wgc0 =  ones(size(bas05.basin1));
bas05.wgc1 = zeros(size(bas05.basin1));
bas05.wgc2 = zeros(size(bas05.basin1));
bas05.wgc3 = zeros(size(bas05.basin1));
bas05.wgc4 = zeros(size(bas05.basin1));
bas05.wgc5 = zeros(size(bas05.basin1));
bas05.wgc6 = zeros(size(bas05.basin1));

% determine weights from distance to closest neighbor for each point in basin
% calculate shortest distance to n1, n2, n3
%for b=[1:10]
%for b=[4]
for b=[1:25]
    b
    tic
    % neighbour loop
%    for n = 1:7
%        eval(['bin' num2str(n-1) '= find(bas05.basin' num2str(bas05.basn(n,b)) '>0);']);
%    end
    eval(['bin0 = find(bas05.basin' num2str(bas05.n0(b)) '>0);']);
    eval(['bin1 = find(bas05.basin' num2str(bas05.n1(b)) '>0);']);
    eval(['bin2 = find(bas05.basin' num2str(bas05.n2(b)) '>0);']);
    eval(['bin3 = find(bas05.basin' num2str(bas05.n3(b)) '>0);']);
    eval(['bin4 = find(bas05.basin' num2str(bas05.n4(b)) '>0);']);
    eval(['bin5 = find(bas05.basin' num2str(bas05.n5(b)) '>0);']);
    eval(['bin6 = find(bas05.basin' num2str(bas05.n6(b)) '>0);']);
    
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
    bas05.wgn1(bin0)=1-min((sd1)/foff,1);
    bas05.wgn2(bin0)=1-min((sd2)/foff,1);
    bas05.wgn3(bin0)=1-min((sd3)/foff,1);
    bas05.wgn4(bin0)=1-min((sd4)/foff,1);
    bas05.wgn5(bin0)=1-min((sd5)/foff,1);
    bas05.wgn6(bin0)=1-min((sd6)/foff,1);

    %% combined weights
    bas05.wgc0(bin0)=1-(bas05.wgn1(bin0)+bas05.wgn2(bin0)+bas05.wgn3(bin0)+bas05.wgn4(bin0)+bas05.wgn5(bin0)+bas05.wgn6(bin0))./(bas05.wgn0(bin0)+bas05.wgn1(bin0)+bas05.wgn2(bin0)+bas05.wgn3(bin0)+bas05.wgn4(bin0)+bas05.wgn5(bin0)+bas05.wgn6(bin0));
    bas05.wgc1(bin0)=(bas05.wgn1(bin0))./(bas05.wgn0(bin0)+bas05.wgn1(bin0)+bas05.wgn2(bin0)+bas05.wgn3(bin0)+bas05.wgn4(bin0)+bas05.wgn5(bin0)+bas05.wgn6(bin0));
    bas05.wgc2(bin0)=(bas05.wgn2(bin0))./(bas05.wgn0(bin0)+bas05.wgn1(bin0)+bas05.wgn2(bin0)+bas05.wgn3(bin0)+bas05.wgn4(bin0)+bas05.wgn5(bin0)+bas05.wgn6(bin0));
    bas05.wgc3(bin0)=(bas05.wgn3(bin0))./(bas05.wgn0(bin0)+bas05.wgn1(bin0)+bas05.wgn2(bin0)+bas05.wgn3(bin0)+bas05.wgn4(bin0)+bas05.wgn5(bin0)+bas05.wgn6(bin0));
    bas05.wgc4(bin0)=(bas05.wgn4(bin0))./(bas05.wgn0(bin0)+bas05.wgn1(bin0)+bas05.wgn2(bin0)+bas05.wgn3(bin0)+bas05.wgn4(bin0)+bas05.wgn5(bin0)+bas05.wgn6(bin0));
    bas05.wgc5(bin0)=(bas05.wgn5(bin0))./(bas05.wgn0(bin0)+bas05.wgn1(bin0)+bas05.wgn2(bin0)+bas05.wgn3(bin0)+bas05.wgn4(bin0)+bas05.wgn5(bin0)+bas05.wgn6(bin0));
    bas05.wgc6(bin0)=(bas05.wgn6(bin0))./(bas05.wgn0(bin0)+bas05.wgn1(bin0)+bas05.wgn2(bin0)+bas05.wgn3(bin0)+bas05.wgn4(bin0)+bas05.wgn5(bin0)+bas05.wgn6(bin0));
    toc
end

% Interpolate to 1 km
bas.wgc0 = interp2(y05,x05,bas05.wgc0,y01,x01);
bas.wgc1 = interp2(y05,x05,bas05.wgc1,y01,x01);
bas.wgc2 = interp2(y05,x05,bas05.wgc2,y01,x01);
bas.wgc3 = interp2(y05,x05,bas05.wgc3,y01,x01);
bas.wgc4 = interp2(y05,x05,bas05.wgc4,y01,x01);
bas.wgc5 = interp2(y05,x05,bas05.wgc5,y01,x01);
bas.wgc6 = interp2(y05,x05,bas05.wgc6,y01,x01);

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
save(['../Data/Basins/ExtBasinScale25_nn7_' num2str(foff)], 'wbas');

%%% Write in netcdf
% BasinIDs
wbas.wbasIDs = zeros(size(wbas.basinIDs,1),size(wbas.basinIDs,2),wbas.nn); 
for i=1:size(wbas.basinIDs,1)
    for j=1:size(wbas.basinIDs,2)
        in = wbas.basinIDs(i,j);
        if (in > 0)
            for n = 1:wbas.nn
                wbas.wbasIDs(i,j,n) = wbas.basn(n,in);
            end
        end
    end
end
ncwrite_GrIS('../Data/Basins/BasinIDs_nb25_nn7.nc',wbas.wbasIDs,'BasinIDs',{'x','y','n'},1);
ncwrite_GrIS('../Data/Basins/BasinWGTs_nb25_nn7.nc',wbas.wbasWGTs,'BasinWGTs',{'x','y','n'},1);


% Plotting
shade(bas.wgc0)
caxis([0 1])
print -dpng -r300 scale_div_wgc0

shade(bas.wgc1+bas.wgc2+bas.wgc3+bas.wgc4+bas.wgc5+bas.wgc6)
caxis([0 1])
print -dpng -r300 scale_div_wgc1-6

% Identify basin neighbors

clear

addpath('../toolbox')

% Params
% falloff distance and resolution in km
res = 5;

% load basin division
load ../Data/Basins/ExtBasinMasks25 
[y,x]= meshgrid(1:size(bas.basin1,2),1:size(bas.basin1,1));

nx = size(y,1);
ny = size(y,2);

% find neighbors
all_nids = zeros(25,25);
for k=1:25
    all_nids_bas = [];
    for i=2:(nx-1)
        for j=2:(ny-1)
            if (bas.basinIDs(i,j) == k)
                nids = [bas.basinIDs(i+1,j),bas.basinIDs(i-1,j),bas.basinIDs(i,j+1),bas.basinIDs(i,j-1),bas.basinIDs(i+1,j+1),bas.basinIDs(i+1,j-1),bas.basinIDs(i-1,j-1),bas.basinIDs(i-1,j-1)];
                all_nids_bas= unique([all_nids_bas, nids]);
            end
        end
    end
    all_nids_bas = all_nids_bas(all_nids_bas~=k);
    all_nids(k,1:length(all_nids_bas)) = all_nids_bas;
end
%all_nids

% basin id
bas.ids = 1:25;
bas.n0 =  1:25;
% locking to the coast left neighbor
bas.n1 =  all_nids(:,1)';
% locking to the coast right neighbor
bas.n2 =  all_nids(:,2)';
% Any other basins of interest? 
bas.n3 =  all_nids(:,3)';
% Any other basins of interest? 
bas.n4 =  all_nids(:,4)';
% Any other basins of interest? 
bas.n5 =  all_nids(:,5)';
% Any other basins of interest? 
bas.n6 =  all_nids(:,6)';

bas.basn = [bas.n0; bas.n1; bas.n2; bas.n3; bas.n4; bas.n5; bas.n6];
% number of neighbors including basin itself
bas.nn = 7;

save(['../Data/Basins/ExtBasinNeighbours25_nn7'], 'bas');


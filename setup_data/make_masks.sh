#!/bin/bash
# Save a series of basin masks

# basin identifiers
# 1:25

# subsample
ncks -F -d x1,1,-1,5 -d y1,1,-1,5 ../Data/Basins/ISMIP6_Extensions.nc ../Data/Basins/ISMIP6_Extensions_e05000m.nc

for id in `seq 1 25`; do
    ncap2 -A -s "basin${id}=IDs*0.; where(IDs==${id}) basin${id}=1" -v ../Data/Basins/ISMIP6_Extensions_e05000m.nc ../Data/Basins/ISMIP6Masks25_05000m.nc
done


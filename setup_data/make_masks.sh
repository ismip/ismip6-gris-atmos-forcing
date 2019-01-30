#!/bin/bash
# Save a series of basin masks

# basin identifiers
# 1:25

for id in `seq 1 25`; do
    echo ${id}
    ncap2 -A -s "basin${id}=IDs*0.; where(IDs==${id}) basin${id}=1" -v ../Data/Basins/ISMIP6_Extensions.nc ../Data/Basins/ISMIP6Masks25.nc
done

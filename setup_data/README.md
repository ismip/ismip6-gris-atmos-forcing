# Prepare input data 

### Copy input files from external archive (search ISMIP6 ftp)
```cp <ExtArchive>/ISMIP6_Extensions.nc ../Data/Basins/```

```cp <ExtArchive>/af2_ISMIP6_GrIS_01000m.nc ../Data/Grid/```

```cp <ExtArchive>/orog_01000m.nc ../Models/OBS/```

```cp <ExtArchive>/sftgif_01000m.nc ../Models/OBS/```

### Prepare Basins (done only once)

`./make_masks.sh`

`matlab`

% Save basins in useful mask format 

`save_extbasins.m`

% Calculate basin weights

`save_extbasin_scale_div.m`


### Prepare MAR data (done once per scenario) 
`matlab`

% Build a forcing time series 

`save_trans_DSMB.m`


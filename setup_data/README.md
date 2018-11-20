# Prepare input data 

### Copy input files from external archive (search ISMIP6 ftp)
```cp <ExtArchive>/MARv3.7-yearly-MIROC5-19xx.nc ../Data/MAR/```

```cp <ExtArchive>/MARv3.7-yearly-MIROC5-20xx.nc ../Data/MAR/```

```cp <ExtArchive>/ISMIP6_Extensions.nc ../Data/Basins/```

```cp <ExtArchive>/af_e05000m.mat ../Data/Grid/```

```cp <ExtArchive>/orog_05000m.nc ../Model/OBS/```

```cp <ExtArchive>/sftgif_05000m.nc ../Model/OBS/```

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
`save_trans_DSMB_MAR37.m`


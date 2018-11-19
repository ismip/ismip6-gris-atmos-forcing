# ismip6-gris-atmos-forcing
Matlab/Shell workflow for generating Greenland atmospheric forcing for the ISMIP6 activity

# Workflow
### create input data (once and once per scenario) 
```setup_data/```

### Get modelled ice geometry and mask (once per initial state)
```Models/```



### Build lookup table (once per scenario)
```make_lookup/```

### Apply SMB scenario using lookup table 
```SMB_remap/```



# Input data
```Data/```
# Model specific 
```Models/```


# matlab toolbox
`toolbox/`

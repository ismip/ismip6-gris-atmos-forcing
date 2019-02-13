#!/bin/bash
# Make SMB forcing from original MAR files 
# Heiko Goelzer 2019 (h.goelzer@uu.nl)

#### Present-day reference is specific per scenario !

# Expected input data on archive
# MIROC5-histo_1950_2005/
#	MARv3.9-yearly-MIROC5-histo-1950.nc
# MIROC5-rcp85_2006_2100/
#	MARv3.9-yearly-MIROC5-rcp85-2006.nc

#datapath=/Volumes/ISMIP6/Data/Raw/SMB/MAR3.9/
datapath=/work/hgoelzer/Processing/RCM/MAR3.9/

gcm=MIROC5
scen=rcp26

#gcm=NorESM1
#scen=rcp85

mkdir -p proc
cd proc
# Collect files from hist and rcp 
for i in `seq 1995 2005`; do
	cp ${datapath}/${gcm}-histo_1950_2005/MARv3.9-yearly-${gcm}-histo-${i}.nc ./MAR_${i}.nc
done
for i in `seq 2006 2014`; do
	cp ${datapath}/${gcm}-${scen}_2006_2100/MARv3.9-yearly-${gcm}-${scen}-${i}.nc ./MAR_${i}.nc
done

# add time information 
for i in `seq 1995 2014`; do
    filename=MAR_${i}.nc
    ncks -O --mk_rec_dmn time $filename $filename
    ncap2 -O -s "time=time*0+$i-1900" $filename $filename
done

# concat them to one time series
ncrcat -O  MAR_1995.nc MAR_1996.nc MAR_1997.nc MAR_1998.nc MAR_1999.nc MAR_2000.nc MAR_2001.nc MAR_2002.nc MAR_2003.nc MAR_2004.nc MAR_2005.nc MAR_2006.nc MAR_2007.nc MAR_2008.nc MAR_2009.nc MAR_2010.nc MAR_2011.nc MAR_2012.nc MAR_2013.nc MAR_2014.nc  MARv3.9-yearly-${gcm}-${scen}-1995-2014.nc

# Long-term average
ncra -O  MARv3.9-yearly-${gcm}-${scen}-1995-2014.nc MARv3.9-yearly-${gcm}-${scen}-ltm1995-2014.nc

# Copy to destination
/bin/cp MARv3.9-yearly-${gcm}-${scen}-ltm1995-2014.nc ../../Data/MAR/

# Extract topg
ncks -v SRF MARv3.9-yearly-${gcm}-${scen}-ltm1995-2014.nc ../../Data/MAR/MARv3.9_topg_01000m.nc
ncrename -v SRF,topg ../../Data/MAR/MARv3.9_topg_01000m.nc


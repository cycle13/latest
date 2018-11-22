
### process observation data
#--------------------------------Settings Start-------------------------------------------------------------
project=cfmip

vars=(SWCF LWCF PRECT)
cases=(CERES-EBAF CERES-EBAF GPCP)

append=ANN_climo.nc
#--------------------------------Settings End-------------------------------------------------------------

nvars=${#vars[@]}

datadir=/home/lyl/WORK1/cesm1_2_1/amwg_diag/obs_data/
workdir=/home/share3/lyl/work3/qinyi/mid-data/${project}/

cd $workdir

for ivar in `seq 0 $[$nvars-1]`
do 
	echo ${vars[ivar]}
	echo ${cases[ivar]}
	ncks -O -v ${vars[ivar]} $datadir/${cases[ivar]}_$append ${vars[ivar]}_$append
done

# form interpolation grid from SWCF file: 1 degree x 1 degree
if [ ! -s "atmos_grid.txt" ]; then
ncks -O -v SWCF $datadir/CERES-EBAF_$append SWCF_$append
cdo griddes SWCF_ANN_climo.nc > atmos_grid.txt
rm SWCF_ANN_climo.nc
fi

for ivar in `seq 0 $[$nvars-1]`
do 
	# interpolate each file to the standard file grid
	cdo remapbil,atmos_grid.txt ${vars[ivar]}_$append tests.nc
	# remove time dimension
	ncwa -O -a time tests.nc tests.nc
	# append variable into one new file, which stores all variables with the same grid resolution
	ncks -A tests.nc lat-lon-atmos-obs.nc
	rm tests.nc
	rm ${vars[ivar]}_$append
done


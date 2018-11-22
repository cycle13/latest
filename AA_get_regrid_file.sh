
### product standard grid file for interpolation
#--------------------------------Settings Start-------------------------------------------------------------

project=cfmip


#--------------------------------Settings End-------------------------------------------------------------
workdir=/home/share3/lyl/work3/qinyi/mid-data/${project}/

# from observation
obsdir=/home/lyl/WORK1/cesm1_2_1/amwg_diag/obs_data/
# from model
moddir=/home/share3/lyl/work3/qinyi/mid-data/${project}/

cd $workdir

# form interpolation grid from SWCF file: 1 degree x 1 degree
if [ ! -s "atmos_grid_obs.txt" ]; then
ncks -O -v SWCF $obsdir/CERES-EBAF_ANN_climo.nc SWCF_ANN_climo_obs.nc
cdo griddes SWCF_ANN_climo_obs.nc > atmos_grid.txt
rm SWCF_ANN_climo_obs.nc
fi

# from atmo model result
if [ ! -s "atmos_grid_mod.txt" ]; then
ncks -O -v SWCF $moddir/FAMIPC5_f09f09_MG15_amip-1980.1984.nc SWCF_ANN_climo_mod.nc
cdo griddes SWCF_ANN_climo_mod.nc > atmos_grid_mod.txt
rm SWCF_ANN_climo_mod.nc
fi



## process model data
#--------------------------------Settings Start-------------------------------------------------------------
#casename=(FAMIPC5_f19f19_MG15_amip FAMIPC5_f19f19_MG15_amip-p4K)
casename=(FAMIPC5_f19f19_mac2_amip_outSGM FAMIPC5_f19f19_mac2_amip-p4K_debug_outSGM)

#casename=(FAMIPC5_f09f09_MG15_amip FAMIPC5_f09f09_MG15_amip-p4K FAMIPC5_f09f09_MG15_amip-4xCO2 \
#	  FAMIPC5_f09f09_mac2_amip FAMIPC5_f09f09_mac2_amip-p4K_debug FAMIPC5_f09f09_mac2_amip-4xCO2)

#casename=(FAMIPC5_f09f09_MG15_amip FAMIPC5_f09f09_MG15_amip-p4K \
#	  FAMIPC5_f09f09_mac2_amip_outSGM FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM)

#casename=(FAMIPC5_f09f09_MG15_amip FAMIPC5_f09f09_MG15_amip-p4K)
#casename=(FAMIPC5_f09f09_mac2_amip_outSGM FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM)

#casename=(FAMIPC5_f09f09_Park-RH_add_Gauss-PDF-output FAMIPC5_f09f09_Park-RH_add_Gauss-PDF-output-p4K)

#casename=(FAMIPC5_f09f09_MG10_amip FAMIPC5_f09f09_MG10_amip-p4K FAMIPC5_f09f09_mac2_amip_outSGM_MG10 FAMIPC5_f09f09_mac2_amip-p4K_debugoff_outSGM_MG10 \
#		  FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM_MG10 \
#		  FAMIPC5_f09f09_MG15_amip FAMIPC5_f09f09_MG15_amip-p4K \
#	      FAMIPC5_f09f09_mac2_amip_outSGM FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM)

#casename=(FAMIPC5_f09f09_MG10_amip FAMIPC5_f09f09_MG10_amip-p4K FAMIPC5_f09f09_mac2_amip_outSGM_MG10 FAMIPC5_f09f09_mac2_amip-p4K_debugoff_outSGM_MG10)
#casename=(FAMIPC5_f09f09_mac2_amip_outSGM_MG10 FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM_MG10 FAMIPC5_f09f09_mac2_amip-p4K_debugoff_outSGM_MG10)

#casename=(FAMIPC5_f09f09_mac2_amip)

ann_mean=true
climatology=true
refine=false
regrid=false
glb_mean=false
db_check_TS=false

int_year=1980
#end_year=1984	# for f09f09
end_year=1988	# for f19f19

int_year_4d=`printf %04d $int_year`
end_year_4d=`printf %04d $end_year`

#vars=(gw,FSNT,FSNTC,FLNT,FLNTC,FSNS,FSNSC,FLNS,FLNSC,LHFLX,SHFLX,PRECC,PRECL,SWCF,LWCF,U,V,OMEGA,CLDLOW,CLDHGH,TGCLDLWP,T,Q,AREL,sgm_tota,sgm_shal,sgm_turb,delta_q,deltaq_uns,deltaq_sat,N1,DCCLDLIQ,ZMDLIQ,DPDLFLIQ,CMFDLIQ,SHDLFLIQ,MACPDLIQ,MPDLIQ,VDCLDLIQ,DCQ,ZMDQ,CMFDQ,MACPDQ,MPDQ,VD01,DTCOND,ZMDT,DPDLFT,CMFDT,SHDLFT,MACPDT,MPDT,DTV,RELHUM,OCNFRAC,alst_def,alst_gp,CLOUD,TS,TREFHT,aa_gp,delta_q,deltaq_uns,deltaq_sat,qtu_shal,thlu_shal,umf_shal,clddep2,wstarPBL,lengi,shi,FREQSH,CMFMC)
vars=(gw,FSNT,FSNTC,FLNT,FLNTC,FSNS,FSNSC,FLNS,FLNSC,LHFLX,SHFLX,PRECC,PRECL,SWCF,LWCF,U,V,OMEGA,CLDLOW,CLDHGH,TGCLDLWP,T,Q,AREL,sgm_tota,sgm_shal,sgm_turb,DCCLDLIQ,ZMDLIQ,DPDLFLIQ,CMFDLIQ,SHDLFLIQ,MACPDLIQ,MPDLIQ,VDCLDLIQ,DCQ,ZMDQ,CMFDQ,MACPDQ,MPDQ,VD01,DTCOND,ZMDT,DPDLFT,CMFDT,SHDLFT,MACPDT,MPDT,DTV,RELHUM,OCNFRAC,alst_def,alst_gp,CLOUD,TS,TREFHT,aa_gp,bb_gp,beta_gp,delta_q,deltaq_uns,deltaq_sat,N1,qtu_shal,thlu_shal,umf_shal,clddep2,wstarPBL,lengi,shi,FREQSH,CMFMC,dqwdz,dthldz,qw2,thl2,qwthl,qwu_qw,thlu_thl,qwu_thlu,FISCCP1_COSP,CLMODIS,CLD_MISR) # for f19f19

#vars=(gw,FSNT,FSNTC,FLNT,FLNTC,FSNS,FSNSC,FLNS,FLNSC,LHFLX,SHFLX,PRECC,PRECL,SWCF,LWCF,U,V,OMEGA,CLDLOW,CLDHGH,TGCLDLWP,T,Q,AREL,DCCLDLIQ,ZMDLIQ,DPDLFLIQ,CMFDLIQ,SHDLFLIQ,MACPDLIQ,MPDLIQ,VDCLDLIQ,DCQ,ZMDQ,CMFDQ,MACPDQ,MPDQ,VD01,DTCOND,ZMDT,DPDLFT,CMFDT,SHDLFT,MACPDT,MPDT,DTV,RELHUM,OCNFRAC,CLOUD,TS,TREFHT,FREQSH,CMFMC,FISCCP1_COSP,CLMODIS,CLD_MISR)
#vars=(gw,FSNT,FSNTC,FLNT,FLNTC,FSNS,FSNSC,FLNS,FLNSC,LHFLX,SHFLX,PRECC,PRECL,SWCF,LWCF,U,V,OMEGA,CLDLOW,CLDHGH,TGCLDLWP,T,Q,AREL,DCCLDLIQ,ZMDLIQ,DPDLFLIQ,CMFDLIQ,SHDLFLIQ,MACPDLIQ,MPDLIQ,VDCLDLIQ,DCQ,ZMDQ,CMFDQ,MACPDQ,MPDQ,VD01,DTCOND,ZMDT,DPDLFT,CMFDT,SHDLFT,MACPDT,MPDT,DTV,RELHUM,OCNFRAC,CLOUD,al_st_prh,al_st_gpd,al_st_nc_prh,al_st_nc_gpd,conden_st_nc_gpd,sgm_gpd,FISCCP1_COSP,CLMODIS,CLD_MISR)

echo $vars
#--------------------------------Settings End-------------------------------------------------------------


ncase=${#casename[@]}

for icase in `seq 0 $[$ncase-1]`
do

casedir=/home/lyl/WORK4/cesm1_2_1/archive/${casename[icase]}/atm/hist/
echo $casedir

workdir=/home/share3/lyl/work3/qinyi/data/cfmip/
if [ ! -d "$workdir" ]; then
mkdir $workdir
fi

outdir=/home/share3/lyl/work3/qinyi/mid-data/cfmip/
if [ ! -d "$outdir" ]; then
mkdir $outdir
fi

#----------------------------------------------------------

cd $workdir

if [ ${ann_mean} == "true" ] ; then
echo "qinyi"

	for i in `seq $int_year $end_year`;
	do
	year=`printf %04d $i`
	echo $year
	
	# ann mean
	ncra -O  -v $vars $casedir/${casename[icase]}.cam.h0.${year}-??.nc $workdir/${casename[icase]}.${year}.ann.mean.nc
#	ncra -O  $casedir/${casename[icase]}.cam.h0.${year}-??.nc $workdir/${casename[icase]}.${year}.ann.mean.nc

	done  # i
fi


if [ ${climatology} == "true" ]; then

# find all files in the request time range.
# caution: the list_tmp is an arry instead of one string.

all_file_list=''
for iyr0 in `seq $int_year $end_year`
do
        iyr0_4d=`printf %04d $iyr0`
    
        list_tmp=`ls $workdir/${casename[icase]}.${iyr0_4d}.ann.mean.nc`
        all_file_list="${all_file_list} ${list_tmp}"

done

file_list=${all_file_list}

echo `ls ${file_list}`

#ncra -O -v $vars ${file_list[*]} $workdir/${casename[icase]}-${int_year_4d}.${end_year_4d}.nc
ncra -O ${file_list[*]} $workdir/${casename[icase]}.ANN.climo.nc

#cp $workdir/${casename[icase]}-${int_year_4d}.${end_year_4d}.nc $outdir
cp $workdir/${casename[icase]}.ANN.climo.nc $outdir

fi

if [ ${glb_mean} == "true" ]; then
	ncwa -a lat,lon -w gw $outdir/${casename[icase]}.ANN.climo.nc $outdir/${casename[icase]}.ANN.climo.globmean.nc
fi

#-------------------------------------------------------------------------
cd $outdir


if [ ${refine} == "true" ]; then
# eliminate time dimension
ncwa -O -a time ${casename[icase]}-${int_year_4d}.${end_year_4d}.nc ${casename[icase]}-${int_year_4d}.${end_year_4d}.nc
# get total precipitation rate
ncap2 -O -s 'PRECT=PRECC+PRECL' ${casename[icase]}-${int_year_4d}.${end_year_4d}.nc PRECT-v1.nc
ncks -A PRECT-v1.nc ${casename[icase]}-${int_year_4d}.${end_year_4d}.nc
rm PRECT-v1.nc
# change unit of PRECT,PRECC,PRECL
ncap2 -O -s 'PRECT=PRECT*86400000' ${casename[icase]}-${int_year_4d}.${end_year_4d}.nc ${casename[icase]}-${int_year_4d}.${end_year_4d}.nc
ncap2 -O -s 'PRECC=PRECC*86400000' ${casename[icase]}-${int_year_4d}.${end_year_4d}.nc ${casename[icase]}-${int_year_4d}.${end_year_4d}.nc
ncap2 -O -s 'PRECL=PRECL*86400000' ${casename[icase]}-${int_year_4d}.${end_year_4d}.nc ${casename[icase]}-${int_year_4d}.${end_year_4d}.nc
ncatted -O -a units,PRECT,o,c,"mm/day" ${casename[icase]}-${int_year_4d}.${end_year_4d}.nc
ncatted -O -a units,PRECC,o,c,"mm/day" ${casename[icase]}-${int_year_4d}.${end_year_4d}.nc
ncatted -O -a units,PRECL,o,c,"mm/day" ${casename[icase]}-${int_year_4d}.${end_year_4d}.nc

fi


if [ ${regrid} == "true" ]; then
# interpolate from model grid to defined standard grid
# because gw is only one-dimension, it cannot be interpolated in CDO. so I eliminate it from the file first.
ncks -O -x -v gw ${casename[icase]}-${int_year_4d}.${end_year_4d}.nc tmp.nc
cdo remapbil,atmos_grid.txt tmp.nc ${casename[icase]}-${int_year_4d}.${end_year_4d}_regrid.nc
rm tmp.nc
fi

done # icase

if [ ${db_check_TS} == "true" ]; then
	cd $outdir
	appends=(ANN.climo.globmean)
	ncdiff -O -v TS FAMIPC5_f09f09_MG10_amip-p4K.$appends.nc FAMIPC5_f09f09_MG10_amip.$appends.nc MG10_diff_TS_globmean.nc
	ncdiff -O -v TS FAMIPC5_f09f09_MG15_amip-p4K.$appends.nc FAMIPC5_f09f09_MG15_amip.$appends.nc MG15_diff_TS_globmean.nc
	ncdiff -O -v TS FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM.$appends.nc FAMIPC5_f09f09_mac2_amip_outSGM.$appends.nc Mac2_MG15_diff_TS_globmean.nc
	ncdiff -O -v TS FAMIPC5_f09f09_mac2_amip-p4K_debugoff_outSGM_MG10.$appends.nc FAMIPC5_f09f09_mac2_amip_outSGM_MG10.$appends.nc Mac2_MG10_debugoff_diff_TS_globmean.nc
	ncdiff -O -v TS FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM_MG10.$appends.nc FAMIPC5_f09f09_mac2_amip_outSGM_MG10.$appends.nc Mac2_MG10_debugon_diff_TS_globmean.nc

	
fi 

echo "Well done!"

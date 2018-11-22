load "./functions_read_file.ncl"
load "./functions_calculate.ncl"
load "./functions_plot.ncl"
load "./functions_discretize.ncl"
load "./functions_regrid.ncl"
load "./functions_bony_method.ncl"

;;;;;processing steps
; 1. read data of model results: function read_data_monthly()
; 2. read data of observation
; 3. regrid all data into the same grid: function time_lat_lon_regrid_to_specific()
; 4. do bony-type binning: function bony_method():30S~30N
; 5. organize data into specific format: function do_pdf_var()


begin

fdir = "/home/lyl/WORK4/cesm1_2_1/archive/"
fappend = "/atm/hist/"

case=(/(/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_mac2_amip"/),\
       (/"FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_mac2_amip-p4K_debug"/),\
       (/"FAMIPC5_f09f09_MG15_amip-4xCO2","FAMIPC5_f09f09_mac2_amip-4xCO2"/) \
     /)

     printVarSummary(case)

var = (/"OMEGA","SWCF","LWCF","CLDLOW","CLDMED","CLDHGH","CLDTOT"/)

tag = (/"nor","p4K","4xCO2"/)

modeldata 				= "T"
obsdata   				= "T"
regrid    				= "T"

do icase=0,2

if (modeldata.eq."T")then

OMEGA_def = read_data_monthly(fdir,fappend,case(icase,0),var(0))
OMEGA_new = read_data_monthly(fdir,fappend,case(icase,1),var(0))

OMEGA_def_500 = OMEGA_def(:,{500},:,:)
OMEGA_new_500 = OMEGA_new(:,{500},:,:)

;=======change unit from Pa/s to hPa/day
OMEGA_def_500 = OMEGA_def_500*24*3600/100.
OMEGA_new_500 = OMEGA_new_500*24*3600/100.

SWCF_def = read_data_monthly(fdir,fappend,case(icase,0),var(1))
SWCF_new = read_data_monthly(fdir,fappend,case(icase,1),var(1))

LWCF_def = read_data_monthly(fdir,fappend,case(icase,0),var(2))
LWCF_new = read_data_monthly(fdir,fappend,case(icase,1),var(2))

CLDLOW_def = read_data_monthly(fdir,fappend,case(icase,0),var(3))
CLDLOW_new = read_data_monthly(fdir,fappend,case(icase,1),var(3))
CLDLOW_def = CLDLOW_def*100.
CLDLOW_new = CLDLOW_new*100.

CLDMED_def = read_data_monthly(fdir,fappend,case(icase,0),var(4))
CLDMED_new = read_data_monthly(fdir,fappend,case(icase,1),var(4))
CLDMED_def = CLDMED_def*100.
CLDMED_new = CLDMED_new*100.

CLDHGH_def = read_data_monthly(fdir,fappend,case(icase,0),var(5))
CLDHGH_new = read_data_monthly(fdir,fappend,case(icase,1),var(5))
CLDHGH_def = CLDHGH_def*100.
CLDHGH_new = CLDHGH_new*100.

CLDTOT_def = read_data_monthly(fdir,fappend,case(icase,0),var(6))
CLDTOT_new = read_data_monthly(fdir,fappend,case(icase,1),var(6))
CLDTOT_def = CLDTOT_def*100.
CLDTOT_new = CLDTOT_new*100.


end if


init = addfile(fdir+case(icase,0)+fappend+case(icase,0)+".cam.h0.1980-01.nc","r")

;======read observation data
fdir1 = "/home/lyl/WORK3/qinyi/data/ERAI/"
fdir2 = "/home/lyl/WORK3/qinyi/data/ISCCP/ftp.climserv.ipsl.polytechnique.fr/cfmip/ISCCP/old_ISCCP/"
fdir3 = "/home/lyl/WORK3/qinyi/data/CERES-EBAF/"


fn1 = "Omega_200001-201012_1deg.nc"
fn2 = "ISCCP_totlowmidhgh_cldcover_198307-200806-unit.nc"
fn3 = "CERES_EBAF-TOA_Ed2.8_Subset_200101-201001.nc"


f1 = addfile(fdir1+fn1,"r")
f2 = addfile(fdir2+fn2,"r")
f3 = addfile(fdir3+fn3,"r")

if (obsdata.eq."T")then

; CERES data starts from 2001
; so I also make other observation data start from 2001
; 2001-2008 = 8*12 = 96 months
OMEGA_obs = short2flt(f1->w(12:71,:,:,:))
;=======change unit from Pa/s to hPa/day
OMEGA_obs = OMEGA_obs*24*3600/100.
;======get 500hPa variables
OMEGA_obs_500 = OMEGA_obs(:,{500},:,:)

CLDLOW_obs = f2->lowcld(198:257,:,:)
CLDMED_obs = f2->midcld(198:257,:,:)
CLDHGH_obs = f2->highcld(198:257,:,:)
CLDTOT_obs = f2->cldtot(198:257,:,:)

SWCF_obs = f3->toa_cre_sw_mon(:59,:,:)
LWCF_obs = f3->toa_cre_lw_mon(:59,:,:)

printVarSummary(OMEGA_obs)
printVarSummary(CLDLOW_obs)
printVarSummary(SWCF_obs)

; there is no problem about regrid and data reading.
; because of the missing value in different time focuses on different regions: north pole or south pole.
; be careful about this.
CLDLOW_obs = lonPivot(CLDLOW_obs,1.25)
CLDMED_obs = lonPivot(CLDMED_obs,1.25)
CLDHGH_obs = lonPivot(CLDHGH_obs,1.25)
CLDTOT_obs = lonPivot(CLDTOT_obs,1.25)

end if


if(regrid.eq."T")then

;=====pay attention to here! the raw data latitude is from positive
;=====to negative. Using (::-1) can flip this.
lat1 = f1->latitude(::-1)
lon1 = f1->longitude

;=====ISCCP: lon2 from -180 to 180, which is different from other dataset.
lat2 = f2->lat
;lon2 = f2->lon
lon2 = CLDLOW_obs&lon

;=====CERES-EBAF
lat3 = f3->lat
lon3 = f3->lon

; this is the model resolution.
lat0 = init->lat
lon0 = init->lon

; regrid grid information
lat_regrid = fspan(-90,90,72)
lon_regrid = fspan(0,360,144)
;print(lat_regrid)
;print(lon_regrid)

;lat_regrid = lat2
;lon_regrid = lon2


dsize = dimsizes(OMEGA_obs)
ntime = dsize(0)

;======regrid to prescribed grid information 2.5deg * 2.5 deg

OMEGA_500_regrid = time_lat_lon_regrid_to_specific(ntime,\
				       OMEGA_obs_500(:,::-1,:),lat1,lon1,\
				       OMEGA_def_500,lat0,lon0,\
				       OMEGA_new_500,lat0,lon0,\
				       lat_regrid,lon_regrid)
printVarSummary(OMEGA_500_regrid)

SWCF_regrid = time_lat_lon_regrid_to_specific(ntime,\
				       SWCF_obs,lat3,lon3,\
				       SWCF_def,lat0,lon0,\
				       SWCF_new,lat0,lon0,\
				       lat_regrid,lon_regrid)
printVarSummary(SWCF_regrid)

LWCF_regrid = time_lat_lon_regrid_to_specific(ntime,\
				       LWCF_obs,lat3,lon3,\
				       LWCF_def,lat0,lon0,\
				       LWCF_new,lat0,lon0,\
				       lat_regrid,lon_regrid)
printVarSummary(LWCF_regrid)


CLDLOW_regrid = time_lat_lon_regrid_to_specific(ntime,\
				       CLDLOW_obs,lat2,lon2,\
				       CLDLOW_def,lat0,lon0,\
				       CLDLOW_new,lat0,lon0,\
				       lat_regrid,lon_regrid)
printVarSummary(CLDLOW_regrid)

CLDMED_regrid = time_lat_lon_regrid_to_specific(ntime,\
				       CLDMED_obs,lat2,lon2,\
				       CLDMED_def,lat0,lon0,\
				       CLDMED_new,lat0,lon0,\
				       lat_regrid,lon_regrid)
printVarSummary(CLDMED_regrid)

CLDHGH_regrid = time_lat_lon_regrid_to_specific(ntime,\
				       CLDHGH_obs,lat2,lon2,\
				       CLDHGH_def,lat0,lon0,\
				       CLDHGH_new,lat0,lon0,\
				       lat_regrid,lon_regrid)

printVarSummary(CLDHGH_regrid)

CLDTOT_regrid = time_lat_lon_regrid_to_specific(ntime,\
				       CLDTOT_obs,lat2,lon2,\
				       CLDTOT_def,lat0,lon0,\
				       CLDTOT_new,lat0,lon0,\
				       lat_regrid,lon_regrid)

printVarSummary(CLDTOT_regrid)



;;;;; output files
system("/bin/rm -f ./tmp-data/data-after-regrid-"+tag(icase)+"-2.5deg.nc")
outfile = addfile("./tmp-data/data-after-regrid-"+tag(icase)+"-2.5deg.nc","c")


outfile->OMEGA_500_regrid = OMEGA_500_regrid
outfile->SWCF_regrid = SWCF_regrid
outfile->LWCF_regrid = LWCF_regrid
outfile->CLDLOW_regrid = CLDLOW_regrid
outfile->CLDMED_regrid = CLDMED_regrid
outfile->CLDHGH_regrid = CLDHGH_regrid
outfile->CLDTOT_regrid = CLDTOT_regrid 

delete(outfile)

end if ; regrid

end do ; do icase=

end


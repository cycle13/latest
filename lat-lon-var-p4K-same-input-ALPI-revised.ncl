
load "./functions_regrid.ncl"
load "./functions_getvar.ncl"
load "./functions_calculate.ncl"

begin

;***************file 1********************

;fdir 				= "/home/lyl/WORK4/cesm1_2_1/archive/"
;fdir_append			= "/atm/hist/"

horizon_interp			= False ; used for different input data resolution


pattern_correlation	= False
regime_contribution = False 
lat_lon_plot		= False

ALPI_revised			= True
ALPI_readout_data		= True
ALPI_mask_ocn			= True
ALPIsub_revised			= True

ALPIsub_mask_plot		= False
ALPIsub_bin_1d_plot		= True
ALPIsub_bin_2d_plot		= True



change_Yaxis		= True

latS_trop			= -30.
latE_trop			= 30.
lonS_scatter		= 220.
lonE_scatter		= 270.

spec_lev			= 887


nbins				= 10
;ALPI_mask_val		= fspan(0,nbins-1,nbins)
ALPI_mask_val		= fspan(0,nbins-1,nbins)+0.5 ;-- the addition of 0.5 is to make sure that 
												 ;-- the value is in the center of cnLevels rather 
												 ;-- than the edge. After this change, RasterFill
												 ;-- is not needed and masking plot would be right as well.

plot_type			= "pdf"


select_idx					= 0	; choose how many regimes would be plotted 
scale_factor				= 2.	; use to amplify diff to make them able to plot in non-diff plot
lev_max						= 887.	; the level where the cloud fraction in Gauss-PDF is largest than that in Park-RH

fdir				= "/home/share3/lyl/work3/qinyi/mid-data/cfmip/"

;case_nor 				= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_MG15_amip-p4K"/)
;case_p4K		 		= (/"FAMIPC5_f09f09_mac2_amip_outSGM","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM"/)

; !!!!!! pay attention: the case_nor(0) and case_nor(1) are the same case. 
; therefore, even though two or more plots are obtained in one figure, they are the same thing. 
;case_nor 				= (/"FAMIPC5_f09f09_mac2_amip_outSGM","FAMIPC5_f09f09_mac2_amip_outSGM"/)
;case_p4K		 		= (/"FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM"/)
; --- 2019-04-29 22:40:55
; --- it is pity that I did not output the two cloud fractions from new 2deg simulation
case_nor 				= (/"FAMIPC5_f09f09_Park-RH_add_Gauss-PDF-output","FAMIPC5_f09f09_Park-RH_add_Gauss-PDF-output"/)
case_p4K		 		= (/"FAMIPC5_f09f09_Park-RH_add_Gauss-PDF-output-p4K","FAMIPC5_f09f09_Park-RH_add_Gauss-PDF-output-p4K"/)


case_name				= (/"RH","PDF"/)

case					= (/case_nor,case_p4K/)
tag					= (/"nor", "p4K"/)

append 				= ".ANN.climo.nc"

SWCF_levs	= (/-170,-150,-135,-120,-105,-90,-75,-60,-45,-30,-15,0,15,30,45/)
;SWCF_dlevs	= (/-80,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,80/)
SWCF_dlevs_1	= (/-50,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,50/)/2.
SWCF_dlevs_2	= (/-50,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,50/)/2.
SWCF_dlevs_sp	= (/-10,-5,5,10/)

LWCF_levs	= (/-45,-30,-20,-10,0,10,20,30,40,50,60,70,85,100,115/) 
LWCF_dlevs_1	= (/-35,-30,-15,-10,-6,-4,-2, 0, 2,4,6,10,15,30,35/)/2.
LWCF_dlevs_2	= (/-35,-30,-15,-10,-6,-4,-2, 0, 2,4,6,10,15,30,35/)/2.
LWCF_dlevs_sp	= (/-10,-5,5,10/)

PRECT_levs	= (/.2,.5,1,2,3,4,5,6,7,8,9,10,12,14,17/)
;PRECT_dlevs	= (/-6,-5,-4,-3,-2,-1,-.5,0,.5,1,2,3,4,5,6/)
PRECT_dlevs_1	= (/-4,-3,-2.5,-2,-1.5,-1,-.5,0,.5,1,1.5,2,2.5,3,4/)
PRECT_dlevs_2	= (/-4,-3,-2.5,-2,-1.5,-1,-.5,0,.5,1,1.5,2,2.5,3,4/)
PRECT_dlevs_sp	= (/-2,-1,1,2/)

CLDLOW_levs	= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
CLDLOW_dlevs_1	= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
CLDLOW_dlevs_2	= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
CLDLOW_dlevs_sp	= (/-0.1,-0.05,0.05,0.1/)

CLDHGH_levs	= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
CLDHGH_dlevs_1	= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
CLDHGH_dlevs_2	= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
CLDHGH_dlevs_sp	= (/-0.1,-0.05,0.05,0.1/)

LTS_levs			= (/5.,6,7,8,9,10,12,14,16,18,20,22,24,26,28/)
LTS_dlevs_1			= (/-1.2,-1.,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1,1.2/)
LTS_dlevs_2			= (/-1.2,-1.,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1,1.2/)
LTS_dlevs_sp		= (/-1.0, -0.5, 0.5, 1.0/)

EIS_levs			= (/-3.,-2,-1,0,1,2,3,4,5,6,8,10,12,13,14/)
EIS_dlevs_1			= (/-1.2,-1.,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1,1.2/)
EIS_dlevs_2			= (/-1.2,-1.,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1,1.2/)
EIS_dlevs_sp		= (/-1.0, -0.5, 0.5, 1.0/)

ECTEI_levs			= (/-3.,-2,-1,0,1,2,3,4,5,6,8,10,12,13,14/)
ECTEI_dlevs_1		= (/-1.2,-1.,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1,1.2/)
ECTEI_dlevs_2		= (/-1.2,-1.,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1,1.2/)
ECTEI_dlevs_sp		= (/-1.0, -0.5, 0.5, 1.0/)

TGCLDLWP_levs 		= (/10,25,50,75,100,150,175,200,250,300,350,400,450,500,550/)
TGCLDLWP_dlevs_1	 	= (/-120,-100,-80,-60,-40,-20,-10,0,10,20,40,60,80,100,120/)
TGCLDLWP_dlevs_2	 	= (/-120,-100,-80,-60,-40,-20,-10,0,10,20,40,60,80,100,120/)
TGCLDLWP_dlevs_sp	= (/-100.,-50., 50., 100. /)

alst_def_levs		= CLDLOW_levs
alst_def_dlevs_1		= CLDLOW_dlevs_1
alst_def_dlevs_2		= CLDLOW_dlevs_2
alst_def_dlevs_sp	= CLDLOW_dlevs_sp

alst_gp_levs		= CLDLOW_levs
alst_gp_dlevs_1		= CLDLOW_dlevs_1
alst_gp_dlevs_2		= CLDLOW_dlevs_2
alst_gp_dlevs_sp	= CLDLOW_dlevs_sp

al_st_gpd_levs		= CLDLOW_levs*0.5
al_st_gpd_dlevs_1		= CLDLOW_dlevs_1
al_st_gpd_dlevs_2		= CLDLOW_dlevs_2
al_st_gpd_dlevs_sp	= CLDLOW_dlevs_sp

al_st_prh_levs		= CLDLOW_levs*0.5
al_st_prh_dlevs_1		= CLDLOW_dlevs_1
al_st_prh_dlevs_2	= CLDLOW_dlevs_2
al_st_prh_dlevs_sp	= CLDLOW_dlevs_sp

Vars_mid	= (/"T", "Q", "OMEGA"/)

;				   0            1          2           3           4            5          6	    7      8       9		10          11           12        13          14
varsA		= (/"PRECT",	"PRECC",	"PRECL",	"SWCF",		"LWCF",		"CLDLOW",	"CLDHGH", "LTS", "EIS", "ECTEI","TGCLDLWP",	"alst_def", "alst_gp", "al_st_prh", "al_st_gpd"/)
varsA_levs	= (/PRECT_levs, PRECT_levs,     PRECT_levs,     SWCF_levs,      LWCF_levs, 	CLDLOW_levs,	CLDHGH_levs, LTS_levs, EIS_levs, ECTEI_levs,TGCLDLWP_levs,alst_def_levs, alst_gp_levs, al_st_prh_levs, al_st_gpd_levs/)
varsA_dlevs_1	= (/PRECT_dlevs_1,PRECT_dlevs_1,    PRECT_dlevs_1,    SWCF_dlevs_1,     LWCF_dlevs_1,	CLDLOW_dlevs_1,	CLDHGH_dlevs_1, LTS_dlevs_1, EIS_dlevs_1,	ECTEI_dlevs_1,TGCLDLWP_dlevs_1, alst_def_dlevs_1, alst_gp_dlevs_1, al_st_prh_dlevs_1, al_st_gpd_dlevs_1/)
varsA_dlevs_2	= (/PRECT_dlevs_2,PRECT_dlevs_2,    PRECT_dlevs_2,    SWCF_dlevs_2,     LWCF_dlevs_2,	CLDLOW_dlevs_2,	CLDHGH_dlevs_2, LTS_dlevs_2, EIS_dlevs_2,	ECTEI_dlevs_2,TGCLDLWP_dlevs_2, alst_def_dlevs_2, alst_gp_dlevs_2, al_st_prh_dlevs_2, al_st_gpd_dlevs_2/)
varsA_dlevs_sp	= (/PRECT_dlevs_sp,PRECT_dlevs_sp,    PRECT_dlevs_sp,    SWCF_dlevs_sp,     LWCF_dlevs_sp,	CLDLOW_dlevs_sp,	CLDHGH_dlevs_sp, LTS_dlevs_sp, EIS_dlevs_sp, ECTEI_dlevs_sp,TGCLDLWP_dlevs_sp, alst_def_dlevs_sp, alst_gp_dlevs_sp, al_st_prh_dlevs_sp, al_st_gpd_dlevs_sp/)

varsA_units	= (/"mm/day",   "mm/day",       "mm/day",       "W/m~S~2~N~",   "W/m~S~2~N~",	"fraction",	"fraction", "K",	"K",	"K", "g/m2", "fraction", "fraction", "fraction" , "fraction"/)

plotsA_strings		= (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"/)

Vars_mid_4d		= (/"T", "Q", "OMEGA", "RELHUM"/)
Vars_mid_3d		= (/"PRECT", "PRECL", "PRECC"/)


;use_idx		= (/11,12/)
use_idx		= (/0,13,14/)

vars				= varsA(use_idx)
vars_levs			= varsA_levs(use_idx,:)
vars_dlevs_1		= varsA_dlevs_1(use_idx,:)
vars_dlevs_2		= varsA_dlevs_2(use_idx,:)
vars_dlevs_sp		= varsA_dlevs_sp(use_idx,:)
vars_units			= varsA_units(use_idx)

nlat				= 192
nlon				= 288

nlev				= 30

vars_mid			= new((/dimsizes(Vars_mid),dimsizes(tag),dimsizes(case_nor),nlev,nlat,nlon/),float)
vars_mid!0			= "var"
vars_mid&var		= fspan(0,dimsizes(Vars_mid)-1,dimsizes(Vars_mid))
vars_mid!1		 	= "tag"
vars_mid&tag		= fspan(0,dimsizes(tag)-1,dimsizes(tag))
if(dimsizes(case_nor).eq.1)then
vars_mid!2		 	= "case"
vars_mid&case		= (/0/)
else
vars_mid!2		 	= "case"
vars_mid&case		= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))
end if

vars_all			= new((/dimsizes(vars),dimsizes(tag),dimsizes(case_nor),nlat,nlon/),float)
vars_all!0			= "var"
vars_all&var		= fspan(0,dimsizes(vars)-1,dimsizes(vars))
vars_all!1		 	= "tag"
vars_all&tag		= fspan(0,dimsizes(tag)-1,dimsizes(tag))
if(dimsizes(case_nor).eq.1)then
vars_all!2		 	= "case"
vars_all&case		= (/0/)
else
vars_all!2		 	= "case"
vars_all&case		= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))
end if


do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1

	ftmp1								= addfile(fdir+case(itag,icase)+append,"r")
	lev									= ftmp1->lev
	OCNFRAC								= ftmp1->OCNFRAC(0,:,:)
	printVarSummary(OCNFRAC)

	do ivar=0,dimsizes(Vars_mid)-1
		vars_mid(ivar,itag,icase,:,:,:)		= ftmp1->$Vars_mid(ivar)$(0,:,:,:)
	end do ; do ivar=

	do ivar=0,dimsizes(vars)-1
		if(isfilevar(ftmp1,vars(ivar))) then
			if(vars(ivar).eq."alst_def".or.vars(ivar).eq."alst_gp" .or. vars(ivar).eq."al_st_prh" .or. vars(ivar).eq."al_st_gpd")then
				tmp								= ftmp1->$vars(ivar)$(0,{lev_max},:,:)
;				tmp								= tmp*100.	; change units
				vars_all(ivar,itag,icase,:,:)	= tmp
			else
				vars_all(ivar,itag,icase,:,:)	 	= ftmp1->$vars(ivar)$(0,:,:)
			end if
			delete(tmp)

			if(vars(ivar).eq."TGCLDLWP")then
				vars_all(ivar,itag,icase,:,:)	= vars_all(ivar,itag,icase,:,:)*1e3
			end if
		else
			if(vars(ivar).eq."PRECT")then
				precl_tmp				= ftmp1->PRECL(0,:,:)
				precc_tmp				= ftmp1->PRECC(0,:,:)
				prect_tmp				= precl_tmp+precc_tmp
				copy_VarCoords(precl_tmp,prect_tmp)
				vars_all(ivar,itag,icase,:,:)	= prect_tmp
				vars_all(ivar,itag,icase,:,:)	= vars_all(ivar,itag,icase,:,:)*8.64e7
				delete(precl_tmp)
				delete(precc_tmp)
				delete(prect_tmp)
			end if
			if(vars(ivar).eq."LTS")then
				result					= cal_ECTEI(vars_mid(0,itag,icase,:,:,:),vars_mid(1,itag,icase,:,:,:),lev)
				LTS_tmp					= result[0]
				vars_all(ivar,itag,icase,:,:)	= LTS_tmp
				delete(LTS_tmp)
			end if
			if(vars(ivar).eq."EIS")then
				result					= cal_ECTEI(vars_mid(0,itag,icase,:,:,:),vars_mid(1,itag,icase,:,:,:),lev)
				EIS_tmp					= tofloat(result[1])
				vars_all(ivar,itag,icase,:,:)	= (/EIS_tmp/)
				delete(EIS_tmp)
			end if
			if(vars(ivar).eq."ECTEI")then
				result					= cal_ECTEI(vars_mid(0,itag,icase,:,:,:),vars_mid(1,itag,icase,:,:,:),lev)
				ECTEI_tmp						= tofloat(result[2])
				vars_all(ivar,itag,icase,:,:)	= (/ECTEI_tmp/)
				delete(ECTEI_tmp)
			end if
		end if
		printMinMax(vars_all(ivar,itag,icase,:,:),True)
	end do ; do ivar=

end do ; do icase=
end do ; do itag=
printVarSummary(vars_all)
printVarSummary(vars_mid)

; ===================== difference between P4K and CTL ==========================
vars_all_tagdiff		= new((/dimsizes(vars),dimsizes(case_nor),nlat,nlon/),float)
vars_all_tagdiff		= (/vars_all(:,1,:,:,:)-vars_all(:,0,:,:,:)/)
copy_VarCoords_1(vars_all(var|:,case|:,lat|:,lon|:,tag|:),vars_all_tagdiff)
printVarSummary(vars_all_tagdiff)

; ==================== difference between PDF and RH ===============================
vars_all_casediff		= new((/dimsizes(vars),dimsizes(tag),nlat,nlon/),float)
vars_all_casediff(:,:,:,:)		= (/vars_all(:,:,1,:,:)-vars_all(:,:,0,:,:)/)
copy_VarCoords_1(vars_all(var|:,tag|:,lat|:,lon|:,case|:),vars_all_casediff)
printVarSummary(vars_all_casediff)

; ==================== difference between PDF_4K-PDF and RH_4K-RH ===============================
vars_all_diff2		= new((/dimsizes(vars),nlat,nlon/),float)
vars_all_diff2		= (/vars_all_tagdiff(:,1,:,:)-vars_all_tagdiff(:,0,:,:)/)
copy_VarCoords_1(vars_all_tagdiff(var|:,lat|:,lon|:,case|:),vars_all_diff2)
printVarSummary(vars_all_diff2)

; ===================== var difference ==========================
vars_all_vardiff		= new((/dimsizes(tag),dimsizes(case_nor),nlat,nlon/),float)
vars_all_vardiff		= (/vars_all(1,:,:,:,:)-vars_all(0,:,:,:,:)/)
copy_VarCoords_1(vars_all(tag|:,case|:,lat|:,lon|:,var|:),vars_all_vardiff)
printVarSummary(vars_all_vardiff)

; ===================== difference of var difference ==========================
vars_all_diff3		= new((/dimsizes(case_nor),nlat,nlon/),float)
vars_all_diff3		= (/vars_all_vardiff(1,:,:,:)-vars_all_vardiff(0,:,:,:)/)
copy_VarCoords_1(vars_all_vardiff(case|:,lat|:,lon|:,tag|:),vars_all_diff3)
printVarSummary(vars_all_diff3)


if (regime_contribution)then
; ================= broken down by regions and regimes =====================================
; outside of the tropics --- broken down into: 90-60S, 60-30S, 30-60N, 60-90N
; tropics --- broken down by 500hPa vertical velocity and LTS
; 1. deep convection and land (land .or. omega < 0)
; 2. trade cumulus (ocean .and. omega > 0 .and. LTS < 17K)
; 3. stratocumulus (ocean .and. omega > 0 .and. LTS > 17K)

; 1.1 calculate LTS
LTS					= new((/dimsizes(tag),dimsizes(case_nor),nlat,nlon/),float) ; unit is degreeC
LTS!0		 		= "tag"
LTS&tag				= fspan(0,dimsizes(tag)-1,dimsizes(tag))
LTS!1		 		= "case"
LTS&case			= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))

do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1
result					= cal_eis(vars_mid(0,itag,icase,:,:,:),vars_mid(1,itag,icase,:,:,:),lev)
LTS(itag,icase,:,:)		= tofloat(result[0])
end do ; do icase=
end do ; do itag=
copy_VarCoords(vars_all(0,:,:,:,:),LTS)
printVarSummary(LTS)
printMinMax(LTS,True)

; 1.2 regime separation
regimes_name		= (/"SH-high", "SH-mid", "NH-mid", "NH-high", "Trop_deep_lnd", "Trop_Shcu", "Trop_StCu"/)
nregimes			= dimsizes(regimes_name)
regimes				= new((/nregimes,dimsizes(vars),dimsizes(tag),dimsizes(case_nor)/),float) 
data2regimes		= new((/nregimes,dimsizes(vars),dimsizes(tag),dimsizes(case_nor),nlat,nlon/),float) 
data2regimes(0,:,:,:,:,:)		= vars_all

lon					= vars_all&lon
lat                 = vars_all&lat
clat                = tofloat(cos(0.01745329*lat))
clat!0              = "lat"
clat&lat            = vars_all&lat

latSs				= (/-90,-60,30,60/)
latEs				= (/-60,-30,60,90/)
lonSs				= (/0,0,0,0/)
lonEs				= (/360,360,360,360/)

do ivar=0,dimsizes(vars)-1
do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1

	omega500				= vars_mid(2,itag,icase,{500},:,:)
	lat2d					= conform_dims(dimsizes(omega500),lat,0)
	lon2d					= conform_dims(dimsizes(omega500),lon,1)

	; mid- and high- latitudes: SH-high, SH-mid, NH-mid, NH-high
	do ilat=0,dimsizes(latSs)-1
		tmp1			= where(lat2d.gt.latSs(ilat) .and. lat2d.lt.latEs(ilat) .and. lon2d.gt.lonSs(ilat) .and. lon2d.lt.lonEs(ilat),vars_all(ivar,itag,icase,:,:),vars_all@_FillValue)
		data2regimes(ilat,ivar,itag,icase,:,:)		= (/tmp1/)
;		regimes(ilat,:,:,:)		= wgt_areaave_Wrap(vars_all(:,:,:,{latSs(ilat):latEs(ilat)},{lonSs(ilat):lonEs(ilat)}),clat({latSs(ilat):latEs(ilat)}),1,0)
		regimes(ilat,ivar,itag,icase)		= wgt_areaave_Wrap(tmp1,clat,1,0)
		delete(tmp1)
	end do ; do ilat=
	
	; tropical deep convection and land regions
	crit1					= (OCNFRAC.lt.0.2).and.(lat2d.gt.-30.and.lat2d.lt.30.)
	crit2					= (omega500.lt.0).and.(lat2d.gt.-30.and.lat2d.lt.30.)
	tmp1					= where(crit1 .or. crit2, vars_all(ivar,itag,icase,:,:),vars_all@_FillValue)
	data2regimes(4,ivar,itag,icase,:,:)		= (/tmp1/)
	regimes(4,ivar,itag,icase)		= wgt_areaave_Wrap(tmp1,clat,1,0)
	; --- 2019-07-24 23:35:01
	delete([/crit1, crit2, tmp1/])
	
	; tropical trade cumulus
	crit1					= (OCNFRAC.gt.0.8) .and. (omega500.gt.0).and.(lat2d.gt.-30.and.lat2d.lt.30.).and.(LTS(itag,icase,:,:).lt.17)
	tmp1					= where(crit1, vars_all(ivar,itag,icase,:,:),vars_all@_FillValue)
	data2regimes(5,ivar,itag,icase,:,:)		= (/tmp1/)
	regimes(5,ivar,itag,icase)		= wgt_areaave_Wrap(tmp1,clat,1,0)
	delete([/crit1, tmp1/])
	
	; tropical stratocumulus
	crit1					= (OCNFRAC.gt.0.8) .and. (omega500.gt.0).and.(lat2d.gt.-30.and.lat2d.lt.30.).and.(LTS(itag,icase,:,:).gt.17)
	tmp1					= where(crit1, vars_all(ivar,itag,icase,:,:),vars_all@_FillValue)
	data2regimes(6,ivar,itag,icase,:,:)		= (/tmp1/)
	regimes(6,ivar,itag,icase)		= wgt_areaave_Wrap(tmp1,clat,1,0)
	delete([/crit1, tmp1/])
end do ; do itag=
end do ; do icase=
end do ; do ivar=

printVarSummary(data2regimes)

;;; ======== the following part is used to verify whether the regime-sorting is right.
; ---- yes, plots indicate that all regimes are right.
;;wks			= gsn_open_wks(plot_type,"test")
;;
;;res			= True
;;res@cnFillOn		= True
;;res@cnLinesOn		= False
;;res@cnLineLabelsOn	= False
;;
;;res@mpCenterLonF	= 180
;;
;;do ir=0,dimsizes(regimes_name)-1
;;do ivar=0,0;dimsizes(vars)-1
;;do itag=0,dimsizes(tag)-1
;;do icase=0,0;dimsizes(case_nor)-1
;;plot		= gsn_csm_contour_map_ce(wks,data2regimes(ir,ivar,itag,icase,:,:),res)
;;end do 
;;end do
;;end do
;;end do
;;; ======================================================

printVarSummary(regimes)
;print(regimes_name+"		"+regimes(:,0,0,0))

; -================ area-weighted regime's contribution
re   = 6.37122e06
rad  = 4.0 * atan(1.0) / 180.
con  = re * rad                 
clat = tofloat(cos(lat * rad))           ; cosine of latitude

dlon = (lon(2) - lon(1))        ; assume dlon is constant
dlat = (lat(2) - lat(1))        ; assume dlat is constant

dx   = con * dlon * clat        ; dx at each latitude
dy   = con * dlat               ; dy is constant
dydx = dy * dx                  ; dydx(nlat)

wgt  = new((/dimsizes(lat),dimsizes(lon)/), float)
wgt  = tofloat(conform(wgt, dydx, 0))
opt	= 0

; change the non-missing region value into 1.0
data2regimes_1				= where(.not.ismissing(data2regimes),1.0,vars_all@_FillValue)
; change all grids into 1. to get the total earth area.
data2regimes_full			= data2regimes_1
data2regimes_full			= where(ismissing(data2regimes_full),1.0,1.0)
printVarSummary(data2regimes_full)
; get the area weigting for each regimes
data2regimes_1_percent = wgt_areasum2(data2regimes_1, wgt, opt)/wgt_areasum2(data2regimes_full,wgt,opt)
printVarSummary(data2regimes_1_percent)

; ================= difference between P4K and CNTL

regimes_tagdiff				= (/regimes(:,:,1,:)-regimes(:,:,0,:)/) ; (regime, var, case_nor)
; 2018.10.04 difference of tagdiff between alst_gp and alst_def
regimes_diff2				= (/regimes_tagdiff(:,1,:)-regimes_tagdiff(:,0,:)/) ; (regime, case_nor)

regimes_areawgt			 	= (/data2regimes_1_percent*regimes/)
regimes_tagdiff_areawgt		= (/regimes_areawgt(:,:,1,:)-regimes_areawgt(:,:,0,:)/) ; (regime, var, case_nor)
; 2018.10.04
regimes_diff2_areawgt		= (/regimes_tagdiff_areawgt(:,1,:)-regimes_tagdiff_areawgt(:,0,:)/) ; (regime, case_nor)
print(regimes_diff2_areawgt)

Vars_out                    = regimes_name
Vars_out_fix                = Vars_out
do ivar = 0,dimsizes(Vars_out)-1
    lens                    = strlen(Vars_out(ivar))
    Vars_out_fix(ivar)      = str_concat((/conform_dims(20-lens, " ",-1)/))
    Vars_out_fix(ivar)      = str_concat((/Vars_out(ivar),Vars_out_fix(ivar)/))
    delete(lens)
end do ; do ivar=

do icase = 0,dimsizes(case_nor)-1
; absolute change
;print("For "+case_nor(icase))
;print("regimes	area-weight_Var1	area-weight_Var2	tagdiff_Var1	tagdiff_Var2 tagdiff_Var2-Var1	tagdiff_areawgt_Var1	tagdiff_areawgt_Var2	tagdiff_areawgt_Var2-Var1")
;print(Vars_out_fix+" "+sprintf("%9.5f",data2regimes_1_percent(:,0,0,icase))+"	"+sprintf("%9.5f",data2regimes_1_percent(:,1,0,icase))+"	"+sprintf("%9.5f",regimes_tagdiff(:,0,icase))+"	"+sprintf("%9.5f",regimes_tagdiff(:,1,icase))+"	"+sprintf("%9.5f",regimes_diff2(:,0))+"	"+sprintf("%9.5f",regimes_tagdiff_areawgt(:,0,icase))+"	"+sprintf("%9.5f",regimes_tagdiff_areawgt(:,1,icase))+"	"+sprintf("%9.5f",regimes_diff2_areawgt(:,1)))

; relative change
print("For "+case_nor(icase))
print("regimes	area-weight_Var1	area-weight_Var2	tagdiff_Var1	tagdiff_Var2 tagdiff_Var2-Var1	tagdiff_areawgt_Var1	tagdiff_areawgt_Var2	tagdiff_areawgt_Var2-Var1")
print(Vars_out_fix+" "+sprintf("%9.5f",data2regimes_1_percent(:,0,0,icase))+"	"+sprintf("%9.5f",data2regimes_1_percent(:,1,0,icase))+"	"+sprintf("%9.5f",regimes_tagdiff(:,0,icase))+"	"+sprintf("%9.5f",regimes_tagdiff(:,1,icase))+"	"+sprintf("%9.5f",regimes_diff2(:,0)*100./dim_sum_n(regimes_tagdiff(:,0,icase),0))+"	"+sprintf("%9.5f",regimes_tagdiff_areawgt(:,0,icase))+"	"+sprintf("%9.5f",regimes_tagdiff_areawgt(:,1,icase))+"	"+sprintf("%9.5f",regimes_diff2_areawgt(:,1)*100./dim_sum_n(regimes_tagdiff_areawgt(:,0,icase),0)))

; ---- check whether adding the weigting could equal to the global mean value got from direct calculation. -- the answer is yes.
;print("verify the correctness of weighting:		"+dim_sum_n(regimes_tagdiff_areawgt(:,0,icase),0)+"	"+dim_sum_n(regimes_tagdiff_areawgt(:,1,icase),0))
; -----
end do ; do icase=

nregimes_used				= dimsizes(regimes_name(select_idx:))

wks					 	= gsn_open_wks(plot_type, "figure/Global-regime-contribution-p4K-"+vars(dimsizes(vars)-1)+"_"+case_nor(1))

sres 						= True
sres@gsnDraw               = False
sres@gsnFrame              = False

sres@vpWidthF              = 0.7
sres@vpHeightF             = 0.5

sres@trXMinF               = 0.
sres@trXMaxF               = nregimes_used+1
;sres@trYMinF               = -1.
;sres@trYMaxF               = 1.

barwidth					= 0.20
sres@gsnXYBarChart         = True
sres@gsnXYBarChartBarWidth = barwidth           ; change bar widths

sres@gsnYRefLine 			= 0

sres@tmXBMode              = "Explicit"         ; explicit labels
sres@tmXBValues            = ispan(1,nregimes_used,1)
sres@tmXBLabels            = regimes_name(select_idx:)

sres@tmXBLabelAngleF		= 45.

sres@tmXBLabelFontHeightF  = 0.0205
sres@tmXTLabelFontHeightF  = 0.0205
sres@tmYLLabelFontHeightF  = 0.0225

sres@gsnStringFontHeightF	= 0.02
sres@gsnRightString        = ""
;sres@tiYAxisString         = "W/m2"
sres@tiYAxisString			= "CLOUD_887hPa (fraction)"
sres@tiMainFontHeightF		= 0.020
sres2						= sres

plots 						= new((/dimsizes(case_nor),dimsizes(case_nor)+1/),graphic)
plots_wgt					= new((/dimsizes(case_nor),dimsizes(case_nor)+1/),graphic)

;colors 						= (/"lightpink",   "lightblue"/)
colors 						= (/"grey",   "black",	"blue"/)

; The X values for each plot is different, ensuring that the bars don't overlap.
do icase=0,dimsizes(case_nor)-1
	sres@tiMainString		= case_nor(icase)
	sres@trYMinF			= -1.*max((/max(regimes_tagdiff(select_idx:,:,icase)),abs(min(regimes_tagdiff(select_idx:,:,icase))),max(scale_factor*regimes_diff2(select_idx:,icase)),abs(min(scale_factor*regimes_diff2(select_idx:,icase)))/))*1.1
	sres@trYMaxF			= 1.*max((/max(regimes_tagdiff(select_idx:,:,icase)),abs(min(regimes_tagdiff(select_idx:,:,icase))),max(scale_factor*regimes_diff2(select_idx:,icase)),abs(min(scale_factor*regimes_diff2(select_idx:,icase)))/))*1.1
	sres@gsnLeftString		= "un-weighted"

	sres2@tiMainString		= case_nor(icase)
	sres2@trYMinF			= -1.*max((/max(regimes_tagdiff_areawgt(select_idx:,:,icase)),abs(min(regimes_tagdiff_areawgt(select_idx:,:,icase))),max(scale_factor*regimes_diff2_areawgt(select_idx:,icase)),abs(min(scale_factor*regimes_diff2_areawgt(select_idx:,icase)))/))*1.1
	sres2@trYMaxF			= 1.*max((/max(regimes_tagdiff_areawgt(select_idx:,:,icase)),abs(min(regimes_tagdiff_areawgt(select_idx:,:,icase))),max(scale_factor*regimes_diff2_areawgt(select_idx:,icase)),abs(min(scale_factor*regimes_diff2_areawgt(select_idx:,icase)))/))*1.1
	sres2@gsnLeftString		= "area-weighted"


	; For un-weighted
	sres@gsnXYBarChartColors = colors(0)
	sres@xyLineColors=(/colors(0),colors(0)/) 	; this is used to change the outline to be the color of the fill. I still dont know why it needs two colors.
	plots(icase,0) = gsn_csm_xy(wks,fspan(1-2*barwidth/2.,nregimes_used-2*barwidth/2,nregimes_used),regimes_tagdiff(select_idx:,0,icase),sres)

	sres@gsnLeftString		 = ""
	sres@gsnXYBarChartColors = colors(1)
	sres@xyLineColors=(/colors(1),colors(1)/)
	plots(icase,1) = gsn_csm_xy(wks,fspan(1+0*barwidth/2,nregimes_used+0*barwidth/2,nregimes_used),regimes_tagdiff(select_idx:,1,icase),sres)

	sres@gsnLeftString		 = ""
	sres@gsnXYBarChartColors = colors(2)
	sres@xyLineColors=(/colors(2),colors(2)/)
	plots(icase,2) = gsn_csm_xy(wks,fspan(1+2*barwidth/2,nregimes_used+2*barwidth/2,nregimes_used),regimes_diff2(select_idx:,icase)*scale_factor,sres)


	; For area-weighted
	sres2@gsnXYBarChartColors = colors(0)
	sres2@xyLineColors = (/colors(0),colors(0)/)
	plots_wgt(icase,0) = gsn_csm_xy(wks,fspan(1-2*barwidth/2.,nregimes_used-2*barwidth/2,nregimes_used),regimes_tagdiff_areawgt(select_idx:,0,icase),sres2)

	sres2@gsnLeftString			= ""
	sres2@gsnXYBarChartColors = colors(1)
	sres2@xyLineColors	= (/colors(1),colors(1)/)
	plots_wgt(icase,1) = gsn_csm_xy(wks,fspan(1+0*barwidth/2,nregimes_used+0*barwidth/2,nregimes_used),regimes_tagdiff_areawgt(select_idx:,1,icase),sres2)

	sres2@gsnLeftString			= ""
	sres2@gsnXYBarChartColors = colors(2)
	sres2@xyLineColors	= (/colors(2),colors(2)/)
	plots_wgt(icase,2) = gsn_csm_xy(wks,fspan(1+2*barwidth/2,nregimes_used+2*barwidth/2,nregimes_used),regimes_diff2_areawgt(select_idx:,icase)*scale_factor,sres2)

; Overlay the last three plots on the first one.
       do j=1,dimsizes(vars)-1+1
         overlay(plots(icase,0),plots(icase,j))
		 overlay(plots_wgt(icase,0),plots_wgt(icase,j))
       end do
end do

pres					= True
;pres@gsnMaximize    	= True   ; Maximize in frame
pres@gsnFrame       	= False  ; Don't advance frame.
pres@gsnPanelBottom 	= 0.13   ; Leave room for labelbar

pres@gsnPanelLeft		= 0.05
pres@gsnPanelRight		= 0.95
pres@gsnPanelXWhiteSpacePercent		= 5.0
pres@gsnPanelYWhiteSpacePercent		= 5.0

; Create some labelbars for a legend
lbres                    = True          ; labelbar only resources
lbres@vpWidthF           = 0.1           ; labelbar width
lbres@vpHeightF          = 0.1           ; labelbar height
lbres@lbBoxMajorExtentF  = 0.10          ; puts space between color boxes
lbres@lbMonoFillPattern  = True          ; Solid fill pattern
lbres@lbLabelFontHeightF = 0.010         ; font height. default is small
lbres@lbLabelJust        = "TopLeft"  ; left justify labels
lbres@lbPerimOn          = False

labels = (/"Park-RH", "Gauss-PDF","PDF-RH"/)
xpos   = (/0.15,  0.25, 0.37/)
ypos	 = 0.75

; Panel the plots.
;First plot
gsn_panel(wks,plots(:,0),(/2,2/),pres)
;drawNDCGrid(wks)
do i=0,dimsizes(case_nor)-1+1
  lbres@lbFillColors = colors(i)
  gsn_labelbar_ndc(wks,1,labels(i),xpos(i),ypos,lbres)
end do	
frame(wks)   ; Advance the frame.

; Second plot
gsn_panel(wks,plots_wgt(:,0),(/2,2/),pres)

ypos	 = 0.95
do i=0,dimsizes(case_nor)-1+1
  lbres@lbFillColors = colors(i)
  gsn_labelbar_ndc(wks,1,labels(i),xpos(i),ypos,lbres)
end do	

frame(wks)   ; Advance the frame.
exit

end if ; regime_contribution

;************************create plot***************
if(lat_lon_plot)then

do ivar=0,dimsizes(vars)-1

;wks 					= gsn_open_wks(plot_type,"figure/lat-lon-"+vars(ivar)+"-HammerProj")
;wks 					= gsn_open_wks(plot_type,"figure/lat-lon-"+vars(ivar)+"-p4K-"+vars(dimsizes(vars)-1))
wks 					= gsn_open_wks(plot_type,"figure/lat-lon-p4K-"+vars(dimsizes(vars)-1))

plot 					= new(9,graphic)
plot1 					= new(9,graphic)

;gsn_define_colormap(wks,"BlueDarkRed18")

cmap1 					= read_colormap_file("WhiteBlueGreenYellowRed")
;cmap2 					= read_colormap_file("BlueRed")
;cmap2					= read_colormap_file("BlueWhiteOrangeRed")
;cmap2				 	= read_colormap_file("NCV_roullet")
cmap2				 	= read_colormap_file("NCV_blue_red")

res                      	= True
res@gsnDraw 				= False
res@gsnFrame 				= False
;res@mpProjection			= "Hammer"
;res@mpPerimOn              = False	 

res@cnFillOn            		= True
res@cnLinesOn           		= False
res@cnLineLabelsOn      		= False
res@cnFillPalette       		= cmap1(:240,:)

res@mpLandFillColor				= "gray"

res@mpCenterLonF         		= 180.

fontheight						= 0.035
res@gsnStringFontHeightF		= fontheight
res@tmXBLabelFontHeightF		= fontheight
res@tmYLLabelFontHeightF		= fontheight

res@tmXBOn			 			= False
res@tmYLOn						= False
res@tmXTOn						= False
res@tmYROn						= False
res@tmXBTickSpacingF			= 60.
res@tmYLTickSpacingF			= 45.

res@cnLevelSelectionMode       	= "ExplicitLevels"
res@cnLevels                   	= vars_levs(ivar,:)

;res@gsnCenterString			= vars(ivar)
res@gsnRightString				= vars_units(ivar)

res@lbOrientation				= "vertical"
res@lbLabelFontHeightF			= fontheight

;;;;;; resl
resl							= res
delete(resl@cnFillPalette)
delete(resl@cnLevels)
delete(resl@mpCenterLonF)
delete(resl@mpLandFillColor)
resl@cnFillOn				 	= False
resl@cnLinesOn					= True
;resl@cnLevels					= vars_dlevs(ivar,::3)
resl@cnLevels					= vars_dlevs_sp(ivar,:)
resl@gsnLeftString				= ""
resl@gsnRightString				= ""

;;;;;;;res2
res2 							= res
delete(res2@cnFillPalette)

res2@cnFillPalette 				= cmap2(:,:)
res2@gsnSpreadColorStart		= 30
res2@gsnSpreadColorEnd			= 227

res2@cnLevelSelectionMode     	= "ExplicitLevels"
res2@cnLevels                 	= vars_dlevs_1(ivar,:)

res3							= res2
res3@cnLevels					= vars_dlevs_2(ivar,:)

var4K					= new((/9,nlat,nlon/),float)
var4K(0,:,:)				= vars_all(ivar,0,0,:,:)
var4K					= (/vars_all(0,0,0,:,:),vars_all(1,0,0,:,:),vars_all_vardiff(0,0,:,:),\ ; (var,tag,case,lat,lon)
					    vars_all(0,1,0,:,:),vars_all(1,1,0,:,:),vars_all_vardiff(1,0,:,:),\ 
					    vars_all_tagdiff(0,0,:,:),vars_all_tagdiff(1,0,:,:), vars_all_diff3(0,:,:)/) 
printVarSummary(var4K)

strings4K				= (/"RH","PDF","PDF-RH",\
					    "RH_4K","PDF_4K","PDF_4K-RH_4K",\
					    "RH_4K-RH","PDF_4K-PDF","dlt_PDF-dlt_RH"/)


plot					= new(9,graphic)
plotov					= new(9,graphic)

p = 0.5
q = -0.25


do it=0,8
if(it.eq.0.or.it.eq.1.or.it.eq.3.or.it.eq.4)then
	res@gsnLeftString			= strings4K(it)
;	plot(it)				= gsn_csm_contour_map(wks,var4K(it,:,:),res)
	plot(it)				= gsn_csm_contour_map(wks,smth9_Wrap(var4K(it,:,:),p,q,False),res)

;	plot_ov					= gsn_csm_contour(wks,var4K(it,:,:),resl)
;	overlay(plot(it),plot_ov)
else
	if(it.eq.6 .or. it.eq.7)then
		res3@gsnLeftString			= strings4K(it)
		plot(it)					= gsn_csm_contour_map(wks,smth9_Wrap(var4K(it,:,:),p,q,False),res3)
	else
		res2@gsnLeftString			= strings4K(it)
;		if(it.eq.8)then
;			res2@gsnCenterString	= sprintf("%5.2f",ra(ivar))
;		else
;			res2@gsnCenterString	= ""
;		end if
		plot(it)				= gsn_csm_contour_map(wks,smth9_Wrap(var4K(it,:,:),p,q,False),res2)
	end if
end if
end do ; do it=


resP = True
resP@gsnPanelMainString = vars(ivar)

resP@gsnPanelXWhiteSpacePercent		= 2.5
resP@gsnPanelYWhiteSpacePercent		= 5.

;if(vars(ivar).eq."SWCF".or.vars(ivar).eq."CLDLOW")then
gsn_panel(wks,plot,(/3,3/),resP)
;gsn_panel(wks,plot1,(/3,3/),resP)
;end if

end do ; do ivar=

end if ; lat_lon_plot

;=================================================
; ALPI IS DEFINE BY LTS AND PRECIPITATION.

;----- ALPI_revised
if(ALPI_revised)then
	pick_var_4d			= (/"al_st_prh","al_st_gpd"/)
	pick_var_3d			= (/"al_st_prh","al_st_gpd"/)

	use_idx_4d		= new(dimsizes(pick_var_4d),integer)
	use_idx_3d		= new(dimsizes(pick_var_3d),integer)

	do ivar=0,dimsizes(pick_var_4d)-1
		use_idx_4d(ivar)	= ind(varsA.eq.pick_var_4d(ivar))
		print(varsA(ivar))
	end do ; do ivar=
	do ivar=0,dimsizes(pick_var_3d)-1
		use_idx_3d(ivar)	= ind(varsA.eq.pick_var_3d(ivar))
		print(varsA(ivar))
	end do ; do ivar=

	Vars_3d					= varsA(use_idx_3d)
	Vars_3d_levs			= varsA_levs(use_idx_3d,:)
	Vars_3d_dlevs_1			= varsA_dlevs_1(use_idx_3d,:)
	Vars_3d_dlevs_2			= varsA_dlevs_2(use_idx_3d,:)
	Vars_3d_dlevs_sp		= varsA_dlevs_sp(use_idx_3d,:)
	Vars_3d_units			= varsA_units(use_idx_3d)
	plots_strings			= plotsA_strings

	Vars_4d					= varsA(use_idx_4d)
	Vars_4d_levs			= varsA_levs(use_idx_4d,:)
	Vars_4d_dlevs_1			= varsA_dlevs_1(use_idx_4d,:)
	Vars_4d_dlevs_2			= varsA_dlevs_2(use_idx_4d,:)
	Vars_4d_dlevs_sp		= varsA_dlevs_sp(use_idx_4d,:)
	Vars_4d_units			= varsA_units(use_idx_4d)
	plots_strings			= plotsA_strings
end if



if(ALPI_readout_data)then
; -------------------------read Vars_4D from monthly data set 

yearS				= 1980
yearE				= 1984
ntime				= (yearE-yearS+1)*12

vars_4d				= new((/dimsizes(Vars_4d),dimsizes(tag),dimsizes(case_nor),ntime,nlev,nlat,nlon/),float) ; (/time,lev,lat,lon/)
vars_4d!0			= "var"
vars_4d&var			= ispan(0,dimsizes(Vars_4d)-1,1)
vars_4d!1		 	= "tag"
vars_4d&tag			= fspan(0,dimsizes(tag)-1,dimsizes(tag))
vars_4d!2		 	= "case"
vars_4d&case		= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))

vars_3d				= new((/dimsizes(Vars_3d),dimsizes(tag),dimsizes(case_nor),ntime,nlat,nlon/),float) ; (/time,lat,lon/)
vars_3d!0			= "var"
vars_3d&var			= ispan(0,dimsizes(Vars_3d)-1,1)
vars_3d!1		 	= "tag"
vars_3d&tag			= fspan(0,dimsizes(tag)-1,dimsizes(tag))
vars_3d!2		 	= "case"
vars_3d&case		= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))


vars_mid_4d				= new((/dimsizes(Vars_mid_4d),dimsizes(tag),dimsizes(case_nor),ntime,nlev,nlat,nlon/),float) ; (/time,lev,lat,lon/)
vars_mid_4d!0			= "var"
vars_mid_4d&var			= ispan(0,dimsizes(Vars_mid_4d)-1,1)
vars_mid_4d!1		 	= "tag"
vars_mid_4d&tag			= fspan(0,dimsizes(tag)-1,dimsizes(tag))
vars_mid_4d!2		 	= "case"
vars_mid_4d&case		= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))

vars_mid_3d				= new((/dimsizes(Vars_mid_3d),dimsizes(tag),dimsizes(case_nor),ntime,nlat,nlon/),float) ; (/time,lat,lon/)
vars_mid_3d!0			= "var"
vars_mid_3d&var			= ispan(0,dimsizes(Vars_mid_3d)-1,1)
vars_mid_3d!1		 	= "tag"
vars_mid_3d&tag			= fspan(0,dimsizes(tag)-1,dimsizes(tag))
vars_mid_3d!2		 	= "case"
vars_mid_3d&case		= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))


do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1

fdir 				= "/home/lyl/WORK4/cesm1_2_1/archive/"+case(itag,icase)+"/atm/hist/"

; --- Get list of files from staryear to endyear
years                        = ispan(yearS,yearE,1)

; ------ get monthly data from specific start and end years --------- starts
do iyr=0,dimsizes(years)-1
yr_4d       				= sprinti("%04d",years(iyr))

if(iyr.eq.0)then
    tmpf         			= systemfunc("ls "+fdir+case(itag,icase)+".cam.h0.*"+yr_4d+"*.nc")
    tmpf1        			= tmpf
	delete(tmpf)
	if(iyr.eq.dimsizes(years)-1)then
		alls				= tmpf1
	end if
else
    tmpf         			= systemfunc("ls "+fdir+case(itag,icase)+".cam.h0.*"+yr_4d+"*.nc")
    tmpf2        			= array_append_record(tmpf1,tmpf,0)
	delete(tmpf)
    delete(tmpf1)
    tmpf1        			= tmpf2
    delete(tmpf2)
    if(iyr.eq.dimsizes(years)-1)then
        alls    			= tmpf1
		delete(tmpf1)
    end if
end if
end do ; do iyr=

all_files       			= alls
delete(alls)

print(all_files)
fall      					= addfiles (all_files, "r")   ; note the "s" of addfile
ListSetType (fall, "cat")             ; concatenate or "merge" (default)
printVarSummary(fall)

; get lev
lev							= fall[0]->lev
printVarSummary(lev)

tmp							= fall[0]->OCNFRAC(0,:,:)
if(horizon_interp)then
	OCNFRAC					= linint2_Wrap(tmp&lon,tmp&lat,tmp,True,klon,klat,0)
else
	OCNFRAC					= tmp
end if
delete(tmp)

; --- read vars_4d
do ivar=0,dimsizes(Vars_4d)-1
print(Vars_4d(ivar))
	if(isfilevar(fall[0],Vars_4d(ivar)))then
		tmp    						= fall[:]->$Vars_4d(ivar)$(:,:,:,:)                    ; note syntax [:]
		tmp@_FillValue  			= default_fillvalue("float")
		print("read "+Vars_4d(ivar)+" well;")
		if(horizon_interp)then
			vars_4d(ivar,itag,icase,:,:,:,:)		= linint2_Wrap(tmp&lon,tmp&lat,tmp,True,klon,klat,0)
			print("regrid "+Vars_4d(ivar)+" well;")
		else
			vars_4d(ivar,itag,icase,:,:,:,:)		= tmp
		end if
		if(Vars_4d(ivar).eq."sgm_tota" .or. Vars_4d(ivar).eq."sgm_shal" .or. Vars_4d(ivar).eq."sgm_turb")then
			vars_4d(ivar,itag,icase,:,:,:,:)	= abs(vars_4d(ivar,itag,icase,:,:,:,:))*1e6
		end if
;		if(Vars_4d(ivar).eq."CLOUD")then
;			vars_4d(ivar,itag,icase,:,:,:,:)	= vars_4d(ivar,itag,icase,:,:,:,:)*100.
;		end if
		if(Vars_4d(ivar).eq."CLDLIQ")then
			vars_4d(ivar,itag,icase,:,:,:,:)		= vars_4d(ivar,itag,icase,:,:,:,:)*1000.
		end if
		delete(tmp)
	end if
end do ; do ivar=

; --- read vars_3d
do ivar=0,dimsizes(Vars_3d)-1
	print(Vars_3d(ivar))
	if(isfilevar(fall[0],Vars_3d(ivar))) then
		tmp									= fall[:]->$Vars_3d(ivar)$
		tmp@_FillValue  					= default_fillvalue("float")
		if(dimsizes(dimsizes(tmp)).eq.4)then ; figure out variables having vertical level coordinate
			tmp1							= fall[:]->$Vars_3d(ivar)$(:,{spec_lev},:,:)
			print("read "+Vars_3d(ivar)+" well;")
			if(horizon_interp)then
				vars_3d(ivar,itag,icase,:,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
			else
				vars_3d(ivar,itag,icase,:,:,:)	= tmp1
			end if
			delete(tmp1)
		else
			tmp1						 	= fall[:]->$Vars_3d(ivar)$
			print("read "+Vars_3d(ivar)+" well;")
			if(horizon_interp)then
				vars_3d(ivar,itag,icase,:,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
			else
				vars_3d(ivar,itag,icase,:,:,:)	= tmp1
			end if
			delete(tmp1)
		end if
		delete(tmp)
		if(Vars_3d(ivar).eq."TGCLDLWP")then
			vars_3d(ivar,itag,icase,:,:,:)	= vars_3d(ivar,itag,icase,:,:,:)*1e3
		end if
		if(Vars_3d(ivar).eq."sgm_tota" .or. Vars_3d(ivar).eq."sgm_shal" .or. Vars_3d(ivar).eq."sgm_turb")then
			vars_3d(ivar,itag,icase,:,:,:)	= abs(vars_3d(ivar,itag,icase,:,:,:))*1e6
		end if
	else
			if(Vars_3d(ivar).eq."netCRF")then
				SWCF						= fall[:]->SWCF(:,:,:)
				LWCF						= fall[:]->LWCF(:,:,:)
				tmp1						= (/SWCF+LWCF/)
				copy_VarCoords(SWCF,tmp1)
				if(horizon_interp)then
					vars_3d(ivar,itag,icase,:,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
				else
					vars_3d(ivar,itag,icase,:,:,:)   = tmp1
				end if
				delete(tmp1)
				delete([/SWCF,LWCF/])
			end if
			if(Vars_3d(ivar).eq."term1_shal")then
				umf_shal_tmp			= fall[:]->umf_shal(:,{spec_lev},:,:)
				qwu_qw_tmp				= fall[:]->qwu_qw(:,{spec_lev},:,:)
				clddep2_tmp				= fall[:]->clddep2(:,:,:)
				wstarPBL_tmp			= fall[:]->wstarPBL(:,:,:)
				tmp1					= -0.5*umf_shal_tmp*qwu_qw_tmp*clddep2_tmp/where(wstarPBL_tmp.ne.0,wstarPBL_tmp,vars_4d@_FillValue)
				tmp1					= tmp1*1e6
				copy_VarCoords(umf_shal_tmp,tmp1)
				if(horizon_interp)then
					vars_3d(ivar,itag,icase,:,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
				else
					vars_3d(ivar,itag,icase,:,:,:)	= tmp1
				end if
				delete(tmp1)
				delete([/umf_shal_tmp,qwu_qw_tmp,clddep2_tmp,wstarPBL_tmp/])
			end if
			if(Vars_3d(ivar).eq."term2_shal")then
				umf_shal_tmp			= fall[:]->umf_shal(:,{spec_lev},:,:)
				thlu_thl_tmp			= fall[:]->thlu_thl(:,{spec_lev},:,:)
				clddep2_tmp				= fall[:]->clddep2(:,:,:)
				wstarPBL_tmp			= fall[:]->wstarPBL(:,:,:)
				tmp1					= -0.5*umf_shal_tmp*thlu_thl_tmp*clddep2_tmp/where(wstarPBL_tmp.ne.0,wstarPBL_tmp,vars_4d@_FillValue)
				tmp1					= tmp1*1e6
				copy_VarCoords(umf_shal_tmp,tmp1)
				if(horizon_interp)then
					vars_3d(ivar,itag,icase,:,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
				else
					vars_3d(ivar,itag,icase,:,:,:)	= tmp1
				end if 
				delete(tmp1)
				delete([/umf_shal_tmp,thlu_thl_tmp,clddep2_tmp,wstarPBL_tmp/])
			end if
			if(Vars_3d(ivar).eq."term3_shal")then
				umf_shal_tmp			= fall[:]->umf_shal(:,{spec_lev},:,:)
				qwu_thlu_tmp			= fall[:]->qwu_thlu(:,{spec_lev},:,:)
				clddep2_tmp				= fall[:]->clddep2(:,:,:)
				wstarPBL_tmp			= fall[:]->wstarPBL(:,:,:)
				tmp1					= -0.5*umf_shal_tmp*qwu_thlu_tmp*clddep2_tmp/where(wstarPBL_tmp.ne.0,wstarPBL_tmp,vars_4d@_FillValue)
				tmp1					= tmp1*1e6
				copy_VarCoords(umf_shal_tmp,tmp1)
				if(horizon_interp)then
					vars_3d(ivar,itag,icase,:,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
				else
					vars_3d(ivar,itag,icase,:,:,:)	= tmp1
				end if
				delete(tmp1)
				delete([/umf_shal_tmp,qwu_thlu_tmp,clddep2_tmp,wstarPBL_tmp/])
			end if
	end if
end do ; do ivar=

; --- read vars_mid_4d
do ivar=0,dimsizes(Vars_mid_4d)-1
print(Vars_mid_4d(ivar))
	if(isfilevar(fall[0],Vars_mid_4d(ivar)))then
		tmp    						= fall[:]->$Vars_mid_4d(ivar)$(:,:,:,:)                    ; note syntax [:]
		tmp@_FillValue  			= default_fillvalue("float")
		print("read "+Vars_mid_4d(ivar)+" well;")
		if(horizon_interp)then
			vars_mid_4d(ivar,itag,icase,:,:,:,:)		= linint2_Wrap(tmp&lon,tmp&lat,tmp,True,klon,klat,0)
			print("regrid "+Vars_mid_4d(ivar)+" well;")
		else
			vars_mid_4d(ivar,itag,icase,:,:,:,:)		= tmp
		end if
		delete(tmp)
	end if
end do ; do ivar=

; --- read vars_mid_3d
do ivar=0,dimsizes(Vars_mid_3d)-1
print(Vars_mid_3d(ivar))
	if(isfilevar(fall[0],Vars_mid_3d(ivar)))then
		tmp    						= fall[:]->$Vars_mid_3d(ivar)$(:,:,:)                    ; (time,lat,lon)
		tmp@_FillValue  			= default_fillvalue("float")
		print("read "+Vars_mid_3d(ivar)+" well;")
		if(horizon_interp)then
			vars_mid_3d(ivar,itag,icase,:,:,:)		= linint2_Wrap(tmp&lon,tmp&lat,tmp,True,klon,klat,0)
			print("regrid "+Vars_mid_3d(ivar)+" well;")
		else
			vars_mid_3d(ivar,itag,icase,:,:,:)		= tmp
		end if
		delete(tmp)
	else
		if(Vars_mid_3d(ivar).eq."PRECT")then
			precl_tmp				= fall[:]->PRECL(:,:,:)
			precc_tmp				= fall[:]->PRECC(:,:,:)
			tmp1					= precl_tmp+precc_tmp
			copy_VarCoords(precl_tmp,tmp1)
			if(horizon_interp)then
				vars_mid_3d(ivar,itag,icase,:,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
				print("regrid "+Vars_mid_3d(ivar)+" well;")
			else
				vars_mid_3d(ivar,itag,icase,:,:,:)   = tmp1
			end if
			delete(tmp1)
			vars_mid_3d(ivar,itag,icase,:,:,:)	= vars_mid_3d(ivar,itag,icase,:,:,:)*8.64e7
			delete(precl_tmp)
			delete(precc_tmp)
		end if
	end if
end do ; do ivar=


end do ; do icase=
end do ; do itag=

printVarSummary(vars_4d)
printMinMax(vars_4d,True)
printVarSummary(vars_3d)
printMinMax(vars_3d,True)
printVarSummary(vars_mid_4d)
printMinMax(vars_mid_4d,True)
printVarSummary(vars_mid_3d)
printMinMax(vars_mid_3d,True)


; ----------- 1.1 calculate LTS
LTS					= new((/dimsizes(tag),dimsizes(case_nor),ntime,nlat,nlon/),float) ; unit is degreeC
LTS!0		 		= "tag"
LTS&tag				= fspan(0,dimsizes(tag)-1,dimsizes(tag))
LTS!1		 		= "case"
LTS&case			= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))

EIS					= LTS

do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1
do itime=0,ntime-1
result							= cal_eis(vars_mid_4d(0,itag,icase,itime,:,:,:),vars_mid_4d(1,itag,icase,itime,:,:,:),lev)
LTS(itag,icase,itime,:,:)		= tofloat(result[0])
EIS(itag,icase,itime,:,:)		= tofloat(result[1])
end do ; do itime=
end do ; do icase=
end do ; do itag=

copy_VarCoords(vars_mid_4d(var|0,tag|:,case|:,time|:,lev|0,lat|:,lon|:),LTS)
copy_VarCoords(vars_mid_4d(var|0,tag|:,case|:,time|:,lev|0,lat|:,lon|:),EIS)

printVarSummary(LTS)
printMinMax(LTS,True)
printVarSummary(EIS)
printMinMax(EIS,True)

; ---------- 1.2 get precipitation
;--- separate PRECT from vars_mid_3d to avoid the plotting difficulty if I dont want to plot PRECT.
do ivar=0,dimsizes(Vars_mid_3d)-1
	if(Vars_mid_3d(ivar).eq."PRECT")then
		print(Vars_mid_3d(ivar))
		prect					= vars_mid_3d(var|ivar,tag|:,case|:,time|:,lat|:,lon|:)
		break
	end if
end do ; do ivar=
printVarSummary(prect)
printMinMax(prect,True)

; --- get OCNFRAC (time,lat,lon)
OCNFRAC_3d	= conform_dims(dimsizes(LTS(tag|0,case|0,time|:,lat|:,lon|:)),OCNFRAC,(/1,2/))
copy_VarCoords(LTS(tag|0,case|0,time|:,lat|:,lon|:),OCNFRAC_3d)
printVarSummary(OCNFRAC_3d)
printMinMax(OCNFRAC_3d,True)

if(ALPI_mask_ocn)then
	; addition: mask land region
	LTS_ocn			= LTS
	EIS_ocn			= EIS
	prect_ocn		= prect
	do itag=0,dimsizes(tag)-1
	do icase=0,dimsizes(case_nor)-1
		LTS_ocn(itag,icase,:,:,:)		= where(OCNFRAC_3d.gt.0.5, LTS(itag,icase,:,:,:),vars_4d@_FillValue)
		EIS_ocn(itag,icase,:,:,:)		= where(OCNFRAC_3d.gt.0.5, EIS(itag,icase,:,:,:),vars_4d@_FillValue)
		prect_ocn(itag,icase,:,:,:)		= where(OCNFRAC_3d.gt.0.5, prect(itag,icase,:,:,:),vars_4d@_FillValue)
	end do ; do icase=
	end do ; do itag=
	copy_VarCoords(LTS,LTS_ocn)
	copy_VarCoords(EIS,EIS_ocn)
	copy_VarCoords(prect,prect_ocn)
else
	LTS_ocn			= LTS
	EIS_ocn			= EIS 
	prect_ocn		= prect
end if ; ALPI_mask_ocn

LTS_tropics	 	= LTS_ocn(:,:,:,{latS_trop:latE_trop},:)
EIS_tropics	 	= EIS_ocn(:,:,:,{latS_trop:latE_trop},:)
prect_tropics	 = prect_ocn(:,:,:,{latS_trop:latE_trop},:)

printVarSummary(LTS_tropics)
printVarSummary(EIS_tropics)
printVarSummary(prect_tropics)

; define ALPI array
ALPI				= LTS_tropics

; --------- 1.3 get LTS_distance and Prec_distance
do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1
do itime=0,ntime-1
	;---1.3.1.1 find LTS_max and LTS_min
	LTS_max				= max(LTS_tropics(itag,icase,itime,:,:))
	LTS_min				= min(LTS_tropics(itag,icase,itime,:,:))
	;---1.3.1.2 find EIS_max and EIS_min
	EIS_max				= max(EIS_tropics(itag,icase,itime,:,:))
	EIS_min				= min(EIS_tropics(itag,icase,itime,:,:))
	;---1.3.2 find prec_max
	prect_max			= max(prect_tropics(itag,icase,itime,:,:))
	prect_min			= min(prect_tropics(itag,icase,itime,:,:))
	
	;---1.3.3 LTS_distance
	LTS_max_2D				= conform_dims(dimsizes(LTS_tropics(itag,icase,itime,:,:)),LTS_max,-1)
	copy_VarCoords(LTS_tropics(itag,icase,itime,:,:),LTS_max_2D)
	LTS_min_2D				= conform_dims(dimsizes(LTS_tropics(itag,icase,itime,:,:)),LTS_min,-1)
	copy_VarCoords(LTS_tropics(itag,icase,itime,:,:),LTS_min_2D)
	LTS_distance			= (LTS_max_2D-LTS_tropics(itag,icase,itime,:,:))/(LTS_max_2D-LTS_min_2D)+0.1
	printMinMax(LTS_distance,True)
	
	;---1.3.3.1 EIS_distance
	EIS_max_2D				= conform_dims(dimsizes(EIS_tropics(itag,icase,itime,:,:)),EIS_max,-1)
	copy_VarCoords(EIS_tropics(itag,icase,itime,:,:),EIS_max_2D)
	EIS_min_2D				= conform_dims(dimsizes(EIS_tropics(itag,icase,itime,:,:)),EIS_min,-1)
	copy_VarCoords(EIS_tropics(itag,icase,itime,:,:),EIS_min_2D)
	EIS_distance			= (EIS_max_2D-EIS_tropics(itag,icase,itime,:,:))/(EIS_max_2D-EIS_min_2D)+0.1
	printMinMax(EIS_distance,True)
	
	;---1.3.4 Prec_distance
	prect_max_2D				= conform_dims(dimsizes(prect_tropics(itag,icase,itime,:,:)),prect_max,-1)
	copy_VarCoords(prect_tropics(itag,icase,itime,:,:),prect_max_2D)
	prect_min_2D				= conform_dims(dimsizes(prect_tropics(itag,icase,itime,:,:)),prect_min,-1)
	copy_VarCoords(prect_tropics(itag,icase,itime,:,:),prect_min_2D)
	prect_distance				= (prect_max_2D-prect_tropics(itag,icase,itime,:,:))/(prect_max_2D-prect_min_2D)+0.1
	printMinMax(prect_distance,True)
	
	;------------ 1.4 get ALPI
	if(ALPIsub_revised)then
		ALPI(itag,icase,itime,:,:)		= atan2(prect_distance,EIS_distance)
	else
		ALPI(itag,icase,itime,:,:)		= atan2(prect_distance,LTS_distance)
	end if

end do ; do itime=
end do ; do icase=
end do ; do itag=

r2d							= 45.0/atan(1.0)     ; conversion factor (radians to degrees)
ALPI						= ALPI*r2d

printVarSummary(ALPI) ; (tag, case, lat, lon)
printMinMax(ALPI,True)

; =========== bin based on ALPI percentile
; ----- 1.1 sorting based on ALPI index: increasing

full_length					= dimsizes(ndtooned(ALPI(0,0,0,{latS_trop:latE_trop},:)))
nlength						= num(.not.ismissing(ndtooned(ALPI(0,0,0,{latS_trop:latE_trop},:))))
ALPI_bins					= new((/dimsizes(tag),dimsizes(case_nor),ntime,nlev,nbins/),float)
vars_4d_bin					= new((/dimsizes(Vars_4d),dimsizes(tag),dimsizes(case_nor),ntime,nlev,full_length,nbins/),float)
vars_4d_bin_mean			= new((/dimsizes(Vars_4d),dimsizes(tag),dimsizes(case_nor),ntime,nlev,nbins/),float)
vars_4d_bin_std				= new((/dimsizes(Vars_4d),dimsizes(tag),dimsizes(case_nor),ntime,nlev,nbins/),float)
vars_3d_bin					= new((/dimsizes(Vars_3d),dimsizes(tag),dimsizes(case_nor),ntime,full_length,nbins/),float)
vars_3d_bin_mean			= new((/dimsizes(Vars_3d),dimsizes(tag),dimsizes(case_nor),ntime,nbins/),float)
vars_3d_bin_std				= new((/dimsizes(Vars_3d),dimsizes(tag),dimsizes(case_nor),ntime,nbins/),float)


vars_4d_bin!0				= "var"
vars_4d_bin&var				= ispan(0,dimsizes(Vars_4d)-1,1)
vars_4d_bin!1				= "tag"
vars_4d_bin&tag				= ispan(0,dimsizes(tag)-1,1)
vars_4d_bin!2				= "case"
vars_4d_bin&case			= ispan(0,dimsizes(case_nor)-1,1)
vars_4d_bin!3				= "time"
vars_4d_bin&time			= vars_4d&time
vars_4d_bin!4				= "lev"
vars_4d_bin&lev				= vars_4d&lev
vars_4d_bin!6				= "bin"
vars_4d_bin&bin				= ispan(1,nbins,1)

copy_VarCoords(vars_4d_bin(:,:,:,:,:,0,:),vars_4d_bin_mean)
copy_VarCoords(vars_4d_bin(:,:,:,:,:,0,:),vars_4d_bin_std)

vars_3d_bin!0				= "var"
vars_3d_bin&var				= ispan(0,dimsizes(Vars_3d)-1,1)

copy_VarCoords(vars_4d_bin(0,:,:,:,0,:,:),vars_3d_bin(0,:,:,:,:,:))
copy_VarCoords(vars_3d_bin(:,:,:,:,0,:),vars_3d_bin_mean)
copy_VarCoords(vars_3d_bin(:,:,:,:,0,:),vars_3d_bin_std)


ALPI_mask					= new((/dimsizes(tag),dimsizes(case_nor),ntime,full_length/),float)

print(nlength)

do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1
do itime=0,1;ntime-1
	ALPI_1D						= ndtooned(ALPI(itag,icase,itime,{latS_trop:latE_trop},:))
	ALPI_1D_sorted				= ALPI_1D
	ip							= dim_pqsort(ALPI_1D_sorted,2) ; sorting and change the raw data
	printVarSummary(ALPI_1D_sorted)

	; ----for Vars_3d ----- 
	do ivar=0,dimsizes(Vars_3d)-1
		vars_3d_1D					= ndtooned(vars_3d(ivar,itag,icase,itime,{latS_trop:latE_trop},:))

		lat							= vars_3d&lat
		clat                		= tofloat(cos(0.01745329*lat))
		clat!0              		= "lat"
		clat&lat            		= vars_3d&lat
		clat_2D						= conform_dims(dimsizes(vars_3d(ivar,itag,icase,itime,:,:)),clat,0)
		copy_VarCoords(vars_3d(ivar,itag,icase,itime,:,:),clat_2D)
		clat_2D_1D					= ndtooned(clat_2D({latS_trop:latE_trop},:))

		do ibin=0,nbins-1
			print("itag="+itag+"	icase="+icase+"	itime="+itime+"	ivar="+ivar+"	ibin="+ibin)
			tmp1		= ALPI_1D_sorted(toint(nlength/nbins)*ibin)
			; avoid missing some points at the end
			; notion: pls dont write "tmp2 = ALPI_1D(dimsizes(ALPI_1D)-1)" when missing value exists.
			if(ibin.eq.(nbins-1))then
				tmp2		= ALPI_1D_sorted(nlength-1)
			else
				tmp2		= ALPI_1D_sorted(toint(nlength/nbins)*(ibin+1))
			end if
			vars_3d_bin(ivar,itag,icase,itime,:,ibin)	= where(ALPI_1D_sorted.gt.tmp1 .and. ALPI_1D_sorted.lt.tmp2, vars_3d_1D(ip), vars_3d@_FillValue)
			vars_3d_bin_mean(ivar,itag,icase,itime,ibin)		= dim_avg_wgt_n_Wrap(vars_3d_bin(ivar,itag,icase,itime,:,ibin),clat_2D_1D(ip),1,0) 	; the weighting effect is small. 2018.12.18
			vars_3d_bin_std(ivar,itag,icase,itime,ibin)		= dim_stddev(vars_3d_bin(ivar,itag,icase,itime,:,ibin))
	
	; ========= very brillient idea suggested by Jishi! in order to reverse to the Lat-Lon plot, we need to keep the ip index always.
			ALPI_mask(itag,icase,itime,ip(ind(ALPI_1D_sorted.gt.tmp1 .and. ALPI_1D_sorted.lt.tmp2)))	=  ALPI_mask_val(ibin) 
	; ========
		end do ; do ibin=
		delete([/tmp1,tmp2/])
	end do ; do ivar=
	
	; ----for Vars_4d ----- 
	do ivar=0,dimsizes(Vars_4d)-1
		do ilev=0,nlev-1
			vars_4d_1D					= ndtooned(vars_4d(ivar,itag,icase,itime,ilev,{latS_trop:latE_trop},:))
			printVarSummary(vars_4d_1D)
			lat							= vars_4d&lat
			clat                		= tofloat(cos(0.01745329*lat))
			clat!0              		= "lat"
			clat&lat            		= vars_4d&lat
			clat_2D						= conform_dims(dimsizes(vars_4d(ivar,itag,icase,itime,ilev,:,:)),clat,0)
			copy_VarCoords(vars_4d(ivar,itag,icase,itime,ilev,:,:),clat_2D)
			clat_2D_1D					= ndtooned(clat_2D({latS_trop:latE_trop},:))
			;	print(ALPI_1D_sorted+"		"+clat_2D_1D(ip))
	
			do ibin=0,nbins-1
				print("itag="+itag+"	icase="+icase+"	itime="+itime+"	ilev="+ilev+"	ivar="+ivar+"	ibin="+ibin)
				tmp1		= ALPI_1D_sorted(toint(nlength/nbins)*ibin)
				; avoid missing some points at the end
				; notion: pls dont write "tmp2 = ALPI_1D(dimsizes(ALPI_1D)-1)" when missing value exists.
				if(ibin.eq.(nbins-1))then 
					tmp2		= ALPI_1D_sorted(nlength-1)
				else
					tmp2		= ALPI_1D_sorted(toint(nlength/nbins)*(ibin+1))
				end if
				vars_4d_bin(ivar,itag,icase,itime,ilev,:,ibin)	= \
					where(ALPI_1D_sorted.gt.tmp1 .and. ALPI_1D_sorted.lt.tmp2, vars_4d_1D(ip), vars_4d@_FillValue)
				;print(num(.not.ismissing(vars_4d_bin(ivar,itag,icase,itime,ilev,:,ibin)))) ; 1366
	
				vars_4d_bin_mean(ivar,itag,icase,itime,ilev,ibin)	=\
					dim_avg_wgt_n_Wrap(vars_4d_bin(ivar,itag,icase,itime,ilev,:,ibin),clat_2D_1D(ip),1,0) 	; the weighting effect is small. 2018.12.18
				vars_4d_bin_std(ivar,itag,icase,itime,ilev,ibin)	=\
					dim_stddev(vars_4d_bin(ivar,itag,icase,itime,ilev,:,ibin))
	
				; ======== very brillient idea suggested by Jishi! in order to reverse to the Lat-Lon plot, we need to keep the ip index always.
				ALPI_mask(itag,icase,itime,ip(ind(ALPI_1D_sorted.gt.tmp1 .and. ALPI_1D_sorted.lt.tmp2)))	=  ALPI_mask_val(ibin) 
				; =======
			end do ; do ibin=
	
			;if(ilev.eq.21)then
			;	; -------- used to test the binned variables
			;	wks						= gsn_open_wks(plot_type,"./figure/ALPI_binned_variable_test")
			;	
			;	dims					= dimsizes(vars_4d_bin(ivar,itag,icase,itime,ilev,:,:))
			;	tmps					= transpose(vars_4d_bin(ivar,itag,icase,itime,ilev,:,:))
			;	
			;	res						= True
			;	res@vpWidthF			= 0.7
			;	res@vpHeightF			= 0.4
			;	res@xyLineColors		= ispan(2,250,250/dims(1))
			;	res@tiYAxisString		= Vars_4d(ivar)
			;	res@tiXAxisString		= "ndtooned Vars_4d [Lat-Lon]"
			;	
			;	plot					= gsn_csm_y(wks,tmps,res)
			;	delete([/plot,tmps,dims/])
			;	; ---------- 
			;end if
	
	delete([/tmp1,tmp2/])
	end do ; do ilev=
end do ; do ivar=

end do ; do itime=
end do ; do icase=
end do ; do itag=


printVarSummary(vars_4d_bin)
printMinMax(vars_4d_bin,True)
printVarSummary(vars_4d_bin_mean)
printMinMax(vars_4d_bin_mean,True)
printVarSummary(vars_3d_bin_mean)
printMinMax(vars_3d_bin_mean,True)

; ----------------- try to output the binned variables for convenience -----------------------------------

outdir						= "/home/share3/lyl/work3/qinyi/mid-data/cfmip/"
outfile						= "ALPI_binned_vars_same_input.nc"

Vars_out_4d					= Vars_4d
Vars_out_3d					= Vars_3d

system("/bin/rm -f "+outdir+outfile)
fout						= addfile(outdir+outfile,"c")

; ---- output these two variables to faciliate the ALPIsub_mask_plot
fout->ALPI					= ALPI
fout->ALPI_mask				= ALPI_mask

do itag=0,dimsizes(tag)-1
	do ivar=0,dimsizes(Vars_4d)-1
		Vars_out_4d(ivar)			= Vars_4d(ivar)+"_"+tag(itag)+"_4d_mean"
		print(Vars_out_4d(ivar))
		fout->$Vars_out_4d(ivar)$		= vars_4d_bin_mean(ivar,itag,:,:,:,:)

		Vars_out_4d(ivar)			= Vars_4d(ivar)+"_"+tag(itag)+"_4d_std"
		print(Vars_out_4d(ivar))
		fout->$Vars_out_4d(ivar)$		= vars_4d_bin_std(ivar,itag,:,:,:,:)
	end do ; do ivar=

	do ivar=0,dimsizes(Vars_3d)-1
		Vars_out_3d(ivar)			= Vars_3d(ivar)+"_"+tag(itag)+"_3d_mean"
		print(Vars_out_3d(ivar))
		fout->$Vars_out_3d(ivar)$		= vars_3d_bin_mean(ivar,itag,:,:,:)

		Vars_out_3d(ivar)			= Vars_3d(ivar)+"_"+tag(itag)+"_3d_std"
		print(Vars_out_3d(ivar))
		fout->$Vars_out_3d(ivar)$		= vars_3d_bin_std(ivar,itag,:,:,:)
	end do ; do ivar=
end do ; do itag=

date=systemfunc("date")
print(date)

exit

end if ; ALPI_readout_data



; -------------------------------------------- read data from mid-data ----------------------------------

indir						= "/home/share3/lyl/work3/qinyi/mid-data/cfmip/"
infile						= "ALPI_binned_vars_same_input.nc"

fin							= addfile(indir+infile,"r")
time						= fin->time
lev							= fin->lev
bin							= fin->bin
ntime						= dimsizes(time)
nlev						= dimsizes(lev)
nbin						= dimsizes(bin)

ALPI						= fin->ALPI
ALPI_mask					= fin->ALPI_mask
printVarSummary(ALPI)
printMinMax(ALPI,True)
printVarSummary(ALPI_mask)
printMinMax(ALPI_mask,True)

if(ALPIsub_bin_2d_plot)then
	vars_4d_bin_mean			= new((/dimsizes(Vars_4d),dimsizes(tag),dimsizes(case_nor),ntime,nlev,nbins/),float)
	vars_4d_bin_std				= new((/dimsizes(Vars_4d),dimsizes(tag),dimsizes(case_nor),ntime,nlev,nbins/),float)
	printVarSummary(vars_4d_bin_mean)

	Vars_out_4d					= Vars_4d
end if

if(ALPIsub_bin_1d_plot)then
	vars_3d_bin_mean			= new((/dimsizes(Vars_3d),dimsizes(tag),dimsizes(case_nor),ntime,nbins/),float)
	vars_3d_bin_std				= new((/dimsizes(Vars_3d),dimsizes(tag),dimsizes(case_nor),ntime,nbins/),float)

	Vars_out_3d					= Vars_3d
end if


if(ALPIsub_bin_2d_plot)then
	do itag=0,dimsizes(tag)-1
		do ivar=0,dimsizes(Vars_4d)-1
			Vars_out_4d(ivar)			= Vars_4d(ivar)+"_"+tag(itag)+"_4d_mean"
			print(Vars_out_4d(ivar))
			vars_4d_bin_mean(ivar,itag,:,:,:,:)	= fin->$Vars_out_4d(ivar)$

			Vars_out_4d(ivar)			= Vars_4d(ivar)+"_"+tag(itag)+"_4d_std"
			print(Vars_out_4d(ivar))
			vars_4d_bin_std(ivar,itag,:,:,:,:) = fin->$Vars_out_4d(ivar)$
		end do ; do ivar=
	end do ; do itag=

	printVarSummary(vars_4d_bin_mean)
	printMinMax(vars_4d_bin_mean,True)
	
	; ----- get the time average in each bin
	vars_4d_bin_mean_tm		= dim_avg_n_Wrap(vars_4d_bin_mean,3)
	vars_4d_bin_std_tm		= dim_avg_n_Wrap(vars_4d_bin_std,3)
	
	printVarSummary(vars_4d_bin_mean_tm)
	printMinMax(vars_4d_bin_mean_tm,True)
end if

if(ALPIsub_bin_1d_plot)then
	do itag=0,dimsizes(tag)-1
		do ivar=0,dimsizes(Vars_3d)-1
			Vars_out_3d(ivar)			= Vars_3d(ivar)+"_"+tag(itag)+"_3d_mean"
			print(Vars_out_3d(ivar))
		 	vars_3d_bin_mean(ivar,itag,:,:,:) = fin->$Vars_out_3d(ivar)$	

			Vars_out_3d(ivar)			= Vars_3d(ivar)+"_"+tag(itag)+"_3d_std"
			print(Vars_out_3d(ivar))
			vars_3d_bin_std(ivar,itag,:,:,:) = fin->$Vars_out_3d(ivar)$	
		end do ; do ivar=
	end do ; do itag=

; ----- get the time average in each bin
	vars_3d_bin_mean_tm		= dim_avg_n_Wrap(vars_3d_bin_mean,3)
	vars_3d_bin_std_tm		= dim_avg_n_Wrap(vars_3d_bin_std,3)
	
	printVarSummary(vars_3d_bin_mean_tm)
	printMinMax(vars_3d_bin_mean_tm,True)
end if

; ----------------------------------------------------
; -------- ALPI_binned_vars plot
; ----------------------------------------------------
if(ALPIsub_bin_2d_plot)then
	trYMinFs	= (/600,750/)
	
	do itrYMinF = 0,dimsizes(trYMinFs)-1
	
		if(ALPIsub_revised)then
			wks  	= gsn_open_wks (plot_type,"./figure/ALPI_binned_vars_nbins_"+nbins+"_revised_ALPI"+"_"+Vars_4d(dimsizes(Vars_4d)-1)+"_2d_"+trYMinFs(itrYMinF))
		else
			wks  	= gsn_open_wks (plot_type,"./figure/ALPI_binned_vars_nbins_"+nbins+"_2d_"+trYMinFs(itrYMinF))
		end if
	
	
	; -------- for lev-bin variables 
	cmap1 							= read_colormap_file("WhiteBlueGreenYellowRed")
	;cmap2 							= read_colormap_file("BlueRed")
	;cmap2							= read_colormap_file("BlueWhiteOrangeRed")
	cmap2							= read_colormap_file("NCV_blu_red")
	
	
	plot_3d							= new((/dimsizes(tag)*dimsizes(case_nor)+5/),graphic)
	
	res_con							= True
	res_con@trYReverse				= True
	res_con@cnFillOn				= True
	res_con@cnLinesOn				= False
	;res_con@cnLineLabelsOn			= False
	res_con@tiYAxisString			= "Pressure (hPa)"
	;res_con@tiXAxisString           = "ALPI percentile (%)"
	res_con@gsnFrame				= False
	res_con@gsnDraw					= False
	
	res_con@vpWidthF				= 0.7
	res_con@vpHeightF				= 0.4
	
;	res_con@trYMinF					= trYMinFs(itrYMinF)
	
	fontheight						= 0.035
	res_con@tmXBMinorOn				= False
	res_con@tmYLMinorOn				= False
	res_con@tmXTOn					= False
	res_con@tmYROn					= False
	res_con@gsnStringFontHeightF	= fontheight
	res_con@tiXAxisFontHeightF		= fontheight
	res_con@tiYAxisFontHeightF		= fontheight
	res_con@tmXBLabelFontHeightF	= fontheight*1.05
	res_con@tmYLLabelFontHeightF	= fontheight*1.05
	res_con@gsnRightStringFontHeightF	= fontheight*0.5

	;res_con@pmLabelBarOrthogonalPosF	= 0.4
	
	res_con@cnLevelSelectionMode	= "ExplicitLevels"
	
	res_con@tmXBMode				= "Explicit"
	res_con@tmXBValues				= fspan(0,nbins,nbins+1)
	res_con@tmXBLabels				= fspan(0,nbins,nbins+1)*100/nbins
	
	res_cond						= res_con
	res_cond@cnLinesOn				= False
	
	res_con@cnFillPalette       	= cmap1(:240,:)
	res_cond@cnFillPalette			= cmap2
	
	res_cond@gsnSpreadColorStart     = 40
	res_cond@gsnSpreadColorEnd       = 217
	
	; (var,tag,case,lev,bin)
	vars_4d_bin_mean_tm_casediff	= (/vars_4d_bin_mean_tm(1,:,0,:,:)-vars_4d_bin_mean_tm(0,:,0,:,:)/) ; al_st_gpd - al_st_prh
	vars_4d_bin_mean_tm_tagdiff		= (/vars_4d_bin_mean_tm(:,1,0,:,:)-vars_4d_bin_mean_tm(:,0,0,:,:)/)
	vars_4d_bin_mean_tm_diff2		= (/vars_4d_bin_mean_tm_tagdiff(1,:,:)-vars_4d_bin_mean_tm_tagdiff(0,:,:)/)
	copy_VarCoords(vars_4d_bin_mean_tm(0,:,0,:,:),vars_4d_bin_mean_tm_casediff)
	copy_VarCoords(vars_4d_bin_mean_tm(:,0,0,:,:),vars_4d_bin_mean_tm_tagdiff)
	copy_VarCoords(vars_4d_bin_mean_tm(0,0,0,:,:),vars_4d_bin_mean_tm_diff2)
	
	
;	do ivar=0,dimsizes(Vars_4d)-1
		res_con@gsnRightString = "CLOUD_"+spec_lev+" hPa"
		res_cond@gsnRightString = "CLOUD_"+spec_lev+" hPa"
	
		leftstring						= (/"RH_CTL","PDF_CTL","PDF-RH_CTL",\
											"RH_P4K","PDF_P4K","PDF-RH_P4K",\
											"RH_P4K-CTL","PDF_P4K-CTL","PDF-RH_P4K-CTL"/)
		print(dimsizes(leftstring))

		data	= (/vars_4d_bin_mean_tm(0,0,0,{trYMinFs(itrYMinF):},:),vars_4d_bin_mean_tm(1,0,0,{trYMinFs(itrYMinF):},:),vars_4d_bin_mean_tm_casediff(0,{trYMinFs(itrYMinF):},:),\
					vars_4d_bin_mean_tm(0,1,1,{trYMinFs(itrYMinF):},:),vars_4d_bin_mean_tm(1,1,1,{trYMinFs(itrYMinF):},:),vars_4d_bin_mean_tm_casediff(1,{trYMinFs(itrYMinF):},:),\
					vars_4d_bin_mean_tm_tagdiff(0,{trYMinFs(itrYMinF):},:),vars_4d_bin_mean_tm_tagdiff(1,{trYMinFs(itrYMinF):},:),vars_4d_bin_mean_tm_diff2({trYMinFs(itrYMinF):},:)\
				  /)

		copy_VarCoords(vars_4d_bin_mean_tm(0,0,0,{trYMinFs(itrYMinF):},:),data(0,:,:))
		printVarSummary(data)

		; --- attention: here the var number should be larger than 1, at least two.
		levels	= (/Vars_4d_levs(0,:),Vars_4d_levs(1,:),Vars_4d_dlevs_1(0,:),\
					Vars_4d_levs(0,:),Vars_4d_levs(1,:),Vars_4d_dlevs_1(1,:),\
					Vars_4d_dlevs_1(0,:),Vars_4d_dlevs_1(0,:),Vars_4d_dlevs_1(0,:)/)

		count = 0
		do it = 0,dimsizes(leftstring)-1
			print(it)
			if(it.eq.0 .or. it.eq.1 .or. it.eq.3 .or. it.eq.4)then
				res_con@gsnLeftString		= plots_strings(count)+" "+leftstring(it)
				res_con@gsnLeftString		= leftstring(it)

				res_con@cnLevels            = levels(it,:)
				plot_3d(it)					= gsn_csm_contour(wks,data(it,:,:),res_con)
			else
				if (it.eq.2 .or. it.eq.5)then ; case_diff
					symMinMaxPlt(data((/2,5/),:,:),15,False,res_cond)
				end if
				if (it.eq.6 .or. it.eq.7)then ; tag_diff
					symMinMaxPlt(data((/6,7/),:,:),15,False,res_cond)
				end if
				if(it.eq.8)then ; diff2
					symMinMaxPlt(data((/8/),:,:),15,False,res_cond)
				end if

				res_cond@gsnLeftString		= plots_strings(count)+" "+leftstring(it)
				res_cond@gsnLeftString		= leftstring(it)
				;res_cond@cnLevels			= levels(it,:)
				plot_3d(it)					= gsn_csm_contour(wks,data(it,:,:),res_cond)
			end if
		count = count +1
		end do ; do it=
		gsn_panel(wks,plot_3d(:),(/3,3/),False)
		delete([/data,leftstring,levels/])
;	end do ; do ivar=


;	resP								= True
;	resP@gsnPanelFigureStrings			= plots_strings(0:dimsizes(Vars_4d)*3-1)
;	resP@amJust   						= "TopLeft"
;	;resP@gsnPanelFigureStringsPerimOn	= False


	delete([/res_con, res_cond, wks, plot_3d/])
	
	end do ; do itrYMinF=
	
end if ; ALPIsub_bin_2d_plot


if(ALPIsub_bin_1d_plot)then

		; ---------------------- ploting starts here ------------------------------------

		if(ALPIsub_revised)then
			wks  	= gsn_open_wks (plot_type,"./figure/ALPI_binned_vars_nbins_"+nbins+"_revised_ALPI"+"_"+Vars_3d(dimsizes(Vars_3d)-1)+"_1d")
		else
			wks  	= gsn_open_wks (plot_type,"./figure/ALPI_binned_vars_nbins_"+nbins+"_1d")
		end if

		plot						= new(1,graphic)
		error_bar					= new((/dimsizes(tag)*dimsizes(case_nor),nbins/),graphic)
		upp_bound					= new((/dimsizes(tag)*dimsizes(case_nor),nbins/),graphic)
		dwn_bound					= new((/dimsizes(tag)*dimsizes(case_nor),nbins/),graphic)
		centers						= new((/dimsizes(tag)*dimsizes(case_nor),nbins/),graphic)
		lines						= new((/dimsizes(tag)*dimsizes(case_nor)/),graphic)
		polygons					= new((/dimsizes(tag)*dimsizes(case_nor)/),graphic)
		
		plot_diff					= new(1,graphic)
		lines_diff					= new((/dimsizes(case_nor)+1/),graphic)
		lines_ref0					= new((/dimsizes(case_nor)+1/),graphic)
		centers_diff				= new((/dimsizes(case_nor)+1,nbins/),graphic)
	
	
		res  						= True
		res@gsnDraw					= False
		res@gsnFrame				= False
		
		res@vpWidthF				= 0.7
		res@vpHeightF				= 0.4
		
		res@xyLineColor				= "white" ;(/"red","blue","lightpink","cyan"/)
		res@xyLineThicknessF		= 3.0
		res@tiXAxisString			= "ALPI percentile (%)"
		
		res@tmXBMode				= "Explicit"
		res@tmXBValues				= fspan(0,nbins,nbins+1)
		res@tmXBLabels				= fspan(0,nbins,nbins+1)*100/nbins
		res@trXMaxF					= (nbins)
		
		fontheight					= 0.035
		res@tiYAxisFontHeightF		=fontheight 
		res@tiXAxisFontHeightF		=fontheight 
		res@tmXBLabelFontHeightF	=fontheight 
		res@tmYLLabelFontHeightF	=fontheight 
		
		res@gsnStringFontHeightF	= fontheight
		
		res@tmXTOn					= False
		res@tmYROn					= False
		
		;res@tmYLFormat				= "0*+^e"
		;res@tmYLMinorOn			= False
		;res@tmXBPrecision			= 6
		
		res_diff					= res
	
				
		count_1		= 0
		count_2		= 0
	
;		do ivar=0,dimsizes(Vars_3d)-1
;			print(Vars_3d(ivar))

			tmp							= (/vars_3d_bin_mean_tm(0,0,0,:),vars_3d_bin_mean_tm(1,0,0,:),\
											vars_3d_bin_mean_tm(0,1,0,:),vars_3d_bin_mean_tm(1,1,0,:)/)
			
			tmp2						= (/vars_3d_bin_std_tm(0,0,0,:),vars_3d_bin_std_tm(1,0,0,:),\
											vars_3d_bin_std_tm(0,1,0,:),vars_3d_bin_std_tm(1,1,0,:)/)

			tmp_diff1					= (/vars_3d_bin_mean_tm(:,1,0,:)-vars_3d_bin_mean_tm(:,0,0,:)/) ; difference between P4K and CTL
			tmp_diff2					= (/tmp_diff1(1,:)-tmp_diff1(0,:)/) ; difference between PDF_P4K-CTL and RH_P4K-CTL
			tmp_diff_final				= new((/dimsizes(case_nor)+1,nbins/),float)
			tmp_diff_final(0,:)			= tmp_diff1(0,:) ; al_st_prh: P4K-CTL
			tmp_diff_final(1,:)			= tmp_diff1(1,:) ; al_st_gpd: P4K-CTL
			tmp_diff_final(2,:)			= tmp_diff2(:) ; al_st_gpd - al_st_prh: P4K-CTL

			if(change_Yaxis)then
				tmp_diff_final			= tmp_diff_final*10.
			end if
	
		
			if(.not.change_Yaxis)then
				res_diff@tiYAxisString			= "CLOUD_"+spec_lev+"hPa (10~S~-1~N~)"
				res@tiYAxisString				= "CLOUD_"+spec_lev+"hPa"
			else
				res_diff@tiYAxisString			= "CLOUD_"+spec_lev+"hPa (10~S~-1~N~)"
				res@tiYAxisString				= "CLOUD_"+spec_lev+"hPa"
			end if

				

			; -----------------------special setting for res@ -----------------------
		
			res@gsnLeftString			= "(a)"
		
			res@trYMaxF					= max((/max(tmp+tmp2)*1.5,max(tmp+tmp2)*0.5/))
			res@trYMinF					= min((/min(tmp-tmp2)*1.5,min(tmp-tmp2)*0.5/))
		
			plot	 					= gsn_csm_y(wks, tmp, res)
		
			count_1 = count_1 + 1
		
			; -----------------------special setting for res_diff@ -----------------------

			res_diff@gsnLeftString				= "(b)"
		
			res_diff@trYMaxF					= max((/max(tmp_diff_final(:1,:))*1.1,max(tmp_diff_final(:1,:))*0.9/))
			res_diff@trYMinF					= min((/min(tmp_diff_final(:1,:))*1.1,min(tmp_diff_final(:1,:))*0.9/))

			plot_diff							= gsn_csm_y(wks, tmp_diff_final(:1,:), res_diff)
		
	
			count_2	= count_2 + 1
		
			; -----------------------add error bar-----------------------
			polyres						= True
			polyres@gsMarkerIndex		= 1
			polyres@gsMarkerSizeF		= 0.02
			
		
			colors						= (/"red","blue","red","blue","black","black"/)
			;dashes						= (/0,12,0,12,0,12/)
			dashes						= (/0,0,12,12,0,12/)
			thickes						= (/1.5,1.5,1.5,1.5,3,3/) ;(/2,2,2,2,4,4/)
			
			colors_diff					= (/"red","blue","black"/)
			dashes_diff					= (/0,0,0/)
			thickes_diff				= (/1.5,1.5,3./)
			
		
			;polyres@gsLineColor		= (/"red","blue","lightpink","cyan"/)

			; ---------------------------- plot climatology value ---------------------------
			
			; loop through the points
			do idx1=0,dimsizes(tag)*dimsizes(case_nor)-1
					; ------------------- add base value line ---------------------------
					polyres@gsLineColor			= colors(idx1)
					polyres@gsLineDashPattern	= dashes(idx1)
					polyres@gsLineThicknessF	= thickes(idx1)
			
					lines(idx1)					= gsn_add_polyline(wks,plot,ispan(0,nbins-1,1)+0.5,tmp(idx1,:),polyres)
					
					; ------------------- use gsn_add_polygon to add shaded uncertainty range---------------------------
					
					gsres							= True
					;gsres@tfPolyDrawOrder			= "Predraw"
					gsres@gsFillColor				= colors(idx1)
					gsres@gsFillOpacityF			= 0.1
					
					xp								= new( (/2*nbins/),float )
					yp								= new( (/2*nbins/),float )
					
					do idx2=0,nbins-1
						yp(idx2)					= tmp(idx1,idx2)+tmp2(idx1,idx2)
						xp(idx2)					= idx2+0.5
						yp(2*nbins-1-idx2)			= tmp(idx1,idx2)-tmp2(idx1,idx2)
						xp(2*nbins-1-idx2)			= idx2+0.5
					end do ; do idx2=
			
					;polygons(idx1)			= gsn_add_polygon(wks,plot,xp,yp,gsres)

					; ------------------- add error bar, upper and lower bounds ---------------------------
					do idx2=0,nbins-1
						error_bar(idx1,idx2)	= gsn_add_polyline(wks,plot,(/idx2,idx2/)+0.5,\
													(/tmp(idx1,idx2)+tmp2(idx1,idx2),tmp(idx1,idx2)-tmp2(idx1,idx2)/),polyres)
						upp_bound(idx1,idx2)	= gsn_add_polyline(wks,plot,(/idx2-0.1,idx2+0.1/)+0.5,\
													(/tmp(idx1,idx2)+tmp2(idx1,idx2),tmp(idx1,idx2)+tmp2(idx1,idx2)/),polyres)
						dwn_bound(idx1,idx2)	= gsn_add_polyline(wks,plot,(/idx2-0.1,idx2+0.1/)+0.5,\
													(/tmp(idx1,idx2)-tmp2(idx1,idx2),tmp(idx1,idx2)-tmp2(idx1,idx2)/),polyres)
					end do ; do idx2=
					
			end do ; do idx1=
			
			; ---------------------------- plot difference value ---------------------------

			;do idx1=0,dimsizes(case_nor)
			do idx1=0,dimsizes(case_nor)-1
					polyres@gsLineColor					= colors_diff(idx1) ;colors(idx1*2)
					polyres@gsLineDashPattern 			= dashes_diff(idx1) ;dashes(idx1*2)
					polyres@gsLineThicknessF			= thickes_diff(idx1) ;thickes(idx1*2)
				
				
					; ------------------------ plot the climatological difference line --------------------
					lines_diff(idx1)	= gsn_add_polyline(wks,plot_diff,ispan(0,nbins-1,1)+0.5,\
																tmp_diff_final(idx1,:),polyres)

					; ------------------------ add line marker --------------------
					do idx2=0,nbins-1
						if(idx1.eq.2)then
						centers_diff(idx1,idx2)	= gsn_add_polymarker(wks,plot_diff,idx2+0.5,\
																			tmp_diff_final(idx1,idx2),polyres)
						end if
					end do ; do idx2=

					; ------------------------ add zero-line as reference line --------------------
					polyres@gsLineColor						= "grey"
					polyres@gsLineDashPattern				= 14
					
			
					lines_ref0(idx1)			= gsn_add_polyline(wks,plot_diff,ispan(0,nbins,1),\
																		fspan(0,0,nbins+1),polyres)
			
			end do ; do idx1=
		
;		end do ; do ivar=
	
		; ------------------------Add legends -----------------------------------
		;legends                         = vars
		legends							= (/"RH_CTL","PDF_CTL","RH_P4K","PDF_P4K"/)
		;legends_diff1					= (/"RH_P4K-CTL", "PDF_P4K-CTL","PDF-RH_P4K-CTL"/)
		legends_diff1					= (/"RH_P4K-CTL", "PDF_P4K-CTL"/)
		
	
		
		lgres							= True
		
		lgres@lgLabelFontHeightF        = 0.10
		
		lgres@vpWidthF                  = 0.25
		lgres@vpHeightF                 = 0.15
		lgres@lgPerimOn                 = False
		
		lgres1							= lgres
		lgres2							= lgres
		
		lgres1@lgLineColors				= colors
		lgres1@lgLineThicknessF			= res@xyLineThicknessF
		lgres1@lgDashIndexes			= dashes
		
		lgres2@lgLineColors				= (/"red","blue","black"/) ;colors(::2)
		lgres2@lgLineThicknessF			= res@xyLineThicknessF
		;lgres2@lgDashIndexes			= dashes(::2)
		lgres2@lgDashIndexes			= (/0,0,0/) ;dashes(::2)
		
		
		
		lbid1                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres1)
		lbid2                           = gsn_create_legend(wks,dimsizes(legends_diff1),legends_diff1,lgres2)
		
		amres                           = True
		amres@amJust					= "TopRight"
		amres@amParallelPosF            = 0.5
		amres@amOrthogonalPosF          = -0.5
		
	
		annoid1                         = gsn_add_annotation(plot,lbid1,amres)
		annoid2                         = gsn_add_annotation(plot_diff,lbid2,amres)
		
		pres							= True
		pres@gsnPanelYWhiteSpacePercent = 5.0
		pres@gsnPanelXWhiteSpacePercent = 5.0
		pres@gsnPanelLeft				= 0.05
		pres@gsnPanelRight				= 0.95
;		pres@gsnPanelFigureStrings		= plots_strings
;		pres@amJust   					= "TopLeft"

		;pres@gsnPanelRowSpec = True
		
		;gsn_panel(wks,plot_diff(1:),(/1,3/),pres)
		;gsn_panel(wks,plot(1:),(/1,3/),pres)

		gsn_panel(wks,(/plot,plot_diff/),(/1,2/),pres)
		
	
		
		delete([/count_1, count_2/])
		exit
end if ; ALPIsub_bin_1d_plot



end

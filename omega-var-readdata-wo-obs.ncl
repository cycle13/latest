load "./functions_read_file.ncl"
load "./functions_calculate.ncl"
load "./functions_plot.ncl"
;load "./functions_discretize.ncl"
load "./functions_regrid.ncl"
load "./functions_bony_method.ncl"

;;;;;processing steps
; 1. read data of model results: function read_data_monthly()
; 2. organize and output data into one file.

; 2018.09.06
; add landmask option in read_data_monthly function.
; 2018.09.09
; add ocnmask option in read_data_monthly function.
; 2018.09.24
; add sgm_tota into Vars
; add sgm_tota_uns and sgm_tota_sat into Vars
; 2018.10.02
; outSGM_2 dataset have been removed.

begin

landmask	= False
ocnmask		= False

; prescribe the start and end years --- 2019-04-26 22:57:39
yearS		= 1980
yearE		= 1988
; resolution
res			= "2deg"


fdir 		= "/home/lyl/WORK4/cesm1_2_1/archive/"
fappend 	= "/atm/hist/"

; cases with MG10 output
;case_nor	= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_mac2_amip_outSGM","FAMIPC5_f09f09_mac2_amip_outSGM_MG10"/)
;case_p4K	= (/"FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM_MG10"/)

; cases with high sensitivity
;case_nor	= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_mac2_amip","FAMIPC5_f09f09_mac2_amip_outSGM_MG10"/)
;case_p4K	= (/"FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_mac2_amip-p4K_debug","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM_MG10"/)

; cases with normal output and MG15 to compare with default simulation.
if(res.eq."1deg")then
	case_nor	= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_mac2_amip_outSGM"/)
	case_p4K	= (/"FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM"/)
	case_4xCO2	= (/"FAMIPC5_f09f09_MG15_amip-4xCO2","FAMIPC5_f09f09_mac2_amip-4xCO2"/)
else
	case_nor	= (/"FAMIPC5_f19f19_MG15_amip","FAMIPC5_f19f19_mac2_amip_outSGM"/)
	case_p4K	= (/"FAMIPC5_f19f19_MG15_amip-p4K","FAMIPC5_f19f19_mac2_amip-p4K_debug_outSGM"/)
end if


;case_all	= (/case_nor, case_p4K, case_4xCO2/)
case_all	= (/case_nor, case_p4K/)
printVarSummary(case_all) ; 3x2

Vars 		= (/"sgm_tota","sgm_tota_sat","sgm_tota_uns","N1","deltaq_sat","OMEGA","SWCF","LWCF","CLDLOW","CLDMED","CLDHGH","CLDTOT","TGCLDLWP","netCRF"/)
;Vars 		= (/"OMEGA","SWCF","LWCF","CLDLOW","CLDMED","CLDHGH","CLDTOT","TGCLDLWP","netCRF"/)

;tag 		= (/"nor","p4K","4xCO2"/)
tag			= (/"nor","p4K"/)


if(res.eq."1deg")then
	nlat					= 192
	nlon					= 288
else
	nlat					= 96
	nlon					= 144
end if 

ntime					= (yearE-yearS+1)*12
print(ntime)

;ntime					= 60	; 1980-1984 -- 12*5=60 months

vars_all				= new((/dimsizes(Vars),dimsizes(tag),dimsizes(case_nor),ntime,nlat,nlon/),float)
vars_all!0				= "var"
vars_all&var			= fspan(0,dimsizes(Vars)-1,dimsizes(Vars))
vars_all!1				= "tag"
vars_all&tag			= fspan(0,dimsizes(tag)-1,dimsizes(tag))
vars_all!2				= "case"
vars_all&case			= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))

vars_all_diff			= new((/dimsizes(Vars),dimsizes(case_nor),ntime,nlat,nlon/),float)
vars_all_diff!0			= "var"
vars_all_diff&var		= fspan(0,dimsizes(Vars)-1,dimsizes(Vars))
vars_all_diff!1			= "case"
vars_all_diff&case		= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))


do itag=0,dimsizes(tag)-1	; nor, p4K, 4xCO2
do icase=0,dimsizes(case_nor)-1 ; def, new
do ivar=0,dimsizes(Vars)-1

print(Vars(ivar))

; --- 2019-04-26 23:01:01 qinyi
; --- please add the startyear and endyear as the args to make data-reading more flexible.
;
if(Vars(ivar).eq."OMEGA")then
	OMEGA_tmp 		= read_data_monthly(fdir,fappend,case_all(itag,icase),Vars(ivar),landmask,ocnmask)
	OMEGA_tmp_500 		= OMEGA_tmp(:,{500},:,:)

	;=======change unit from Pa/s to hPa/day
	OMEGA_tmp_500 		= OMEGA_tmp_500*24*3600/100.

	vars_all(ivar,itag,icase,:,:,:)		= OMEGA_tmp_500
	delete([/OMEGA_tmp,OMEGA_tmp_500/])
else if(Vars(ivar).eq."sgm_tota")then
	if(icase.eq.1)then
		sgm_tmp			= read_data_monthly(fdir,fappend,case_all(itag,icase),Vars(ivar),landmask,ocnmask)
		sgm_tmp_900		= sgm_tmp(:,{900},:,:)
		; change unit
		sgm_tmp_900		= sqrt(abs(sgm_tmp_900))*1e3
		vars_all(ivar,itag,icase,:,:,:)		= sgm_tmp_900
		delete([/sgm_tmp,sgm_tmp_900/])
	else
		vars_all(ivar,itag,icase,:,:,:)		= vars_all@_FillValue
	end if
else if(Vars(ivar).eq."sgm_tota_sat")then
	if(icase.eq.1)then
		sgm_tmp			= read_data_monthly(fdir,fappend,case_all(itag,icase),"sgm_tota",landmask,ocnmask)
		sgm_tmp_900		= sgm_tmp(:,{900},:,:)
		deltaq_sat_tmp	= read_data_monthly(fdir,fappend,case_all(itag,icase),"deltaq_sat",landmask,ocnmask)
		deltaq_sat_tmp_900 = deltaq_sat_tmp(:,{900},:,:)*1e3
		sgm_tmp_900		= where(deltaq_sat_tmp_900.lt.1e-2,sgm_tmp@_FillValue,sgm_tmp_900)
		; change unit
		sgm_tmp_900		= sqrt(abs(sgm_tmp_900))*1e3
		vars_all(ivar,itag,icase,:,:,:)		= sgm_tmp_900
		delete([/sgm_tmp,sgm_tmp_900,deltaq_sat_tmp,deltaq_sat_tmp_900/])
	else
		vars_all(ivar,itag,icase,:,:,:)		= vars_all@_FillValue
	end if
else if(Vars(ivar).eq."sgm_tota_uns")then
	if(icase.eq.1)then
		sgm_tmp			= read_data_monthly(fdir,fappend,case_all(itag,icase),"sgm_tota",landmask,ocnmask)
		sgm_tmp_900		= sgm_tmp(:,{900},:,:)
		deltaq_sat_tmp	= read_data_monthly(fdir,fappend,case_all(itag,icase),"deltaq_sat",landmask,ocnmask)
		deltaq_sat_tmp_900 = deltaq_sat_tmp(:,{900},:,:)*1e3
		sgm_tmp_900		= where(deltaq_sat_tmp_900.gt.1e-2,sgm_tmp@_FillValue,sgm_tmp_900)
		; change unit
		sgm_tmp_900		= sqrt(abs(sgm_tmp_900))*1e3
		vars_all(ivar,itag,icase,:,:,:)		= sgm_tmp_900
		delete([/sgm_tmp,sgm_tmp_900,deltaq_sat_tmp,deltaq_sat_tmp_900/])
	else
		vars_all(ivar,itag,icase,:,:,:)		= vars_all@_FillValue
	end if
else if(Vars(ivar).eq."N1")then
	if(icase.eq.1)then
		N1_tmp			= read_data_monthly(fdir,fappend,case_all(itag,icase),Vars(ivar),landmask,ocnmask)
		N1_tmp_900		= N1_tmp(:,{900},:,:)
		vars_all(ivar,itag,icase,:,:,:)	= N1_tmp_900
		delete([/N1_tmp,N1_tmp_900/])
	else
		vars_all(ivar,itag,icase,:,:,:)	= vars_all@_FillValue
	end if
else if(Vars(ivar).eq."deltaq_sat")then
	if(icase.eq.1)then
		deltaq_sat_tmp			= read_data_monthly(fdir,fappend,case_all(itag,icase),Vars(ivar),landmask,ocnmask)
		deltaq_sat_tmp_900		= deltaq_sat_tmp(:,{900},:,:)
		; change unit
		deltaq_sat_tmp_900		= deltaq_sat_tmp_900*1e3
		vars_all(ivar,itag,icase,:,:,:)	= deltaq_sat_tmp_900
		delete([/deltaq_sat_tmp,deltaq_sat_tmp_900/])
	else
		vars_all(ivar,itag,icase,:,:,:)	= vars_all@_FillValue
	end if
else if(Vars(ivar).eq."netCRF")then
	SWCF_tmp					= read_data_monthly(fdir,fappend,case_all(itag,icase),"SWCF",landmask,ocnmask)
	LWCF_tmp					= read_data_monthly(fdir,fappend,case_all(itag,icase),"LWCF",landmask,ocnmask)
	netCRF_tmp					= SWCF_tmp + LWCF_tmp
	vars_all(ivar,itag,icase,:,:,:) = (/netCRF_tmp/)
	delete([/SWCF_tmp,LWCF_tmp,netCRF_tmp/])
else
	vars_all(ivar,itag,icase,:,:,:)		= read_data_monthly(fdir,fappend,case_all(itag,icase),Vars(ivar),landmask,ocnmask)
	if(Vars(ivar).eq."CLDLOW" .or. Vars(ivar).eq."CLDMED" .or. Vars(ivar).eq."CLDHGH" .or. Vars(ivar).eq."CLDTOT")then
		vars_all(ivar,itag,icase,:,:,:)		= vars_all(ivar,itag,icase,:,:,:)*100.
	end if
	if(Vars(ivar).eq."TGCLDLWP")then
		vars_all(ivar,itag,icase,:,:,:)		= vars_all(ivar,itag,icase,:,:,:)*1e3
	end if
end if
end if
end if
end if
end if
end if
end if

end do ; do ivar=
end do ; do icase=
end do ; do itag=

printVarSummary(vars_all)
; ======== get difference between p4K and cntl experiments
vars_all_diff				= vars_all(:,1,:,:,:,:)-vars_all(:,0,:,:,:,:)
copy_VarCoords_1(vars_all(var|:,case|:,time|:,lat|:,lon|:,tag|:),vars_all_diff)
printVarSummary(vars_all_diff)

	
;;;;; output files
setfileoption("nc","format","netcdf4") 

outdir="/home/share3/lyl/work3/qinyi/mid-data/cfmip/"
outname=outdir+"omega-var-outdata-wo-obs-"+res
print(outname)

if(landmask)then
	system("/bin/rm -f "+outname+"-ocn.nc")
	outfile = addfile(outname+"-ocn.nc","c")
else
	if(ocnmask)then
		system("/bin/rm -f "+outname+"-lnd.nc")
		outfile = addfile(outname+"-lnd.nc","c")
	else
		system("/bin/rm -f "+outname+".nc")
		outfile = addfile(outname+".nc","c")
	end if
end if

print(outfile)
	
Vars_out						= Vars
Vars_out_diff					= Vars

do ivar=0,dimsizes(Vars)-1
	; mean
	do itag=0,dimsizes(tag)-1
		Vars_out(ivar)				= Vars(ivar)+"_"+tag(itag)
		print(Vars_out(ivar))
		outfile->$Vars_out(ivar)$		= vars_all(ivar,itag,:,:,:,:)
	end do ; do itag=

	; difference
	Vars_out_diff(ivar)			= Vars(ivar)+"_diff"
	outfile->$Vars_out_diff(ivar)$		= vars_all_diff(ivar,:,:,:,:)

	end do ; do ivar=

end


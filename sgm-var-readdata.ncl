load "./functions_read_file.ncl"
load "./functions_calculate.ncl"
load "./functions_plot.ncl"
load "./functions_discretize.ncl"
load "./functions_regrid.ncl"
load "./functions_bony_method.ncl"

;;;;;processing steps
; 1. read data of model results: function read_data_monthly()
; 2. output data into ./tmp-data/sgm-var-outdata.nc

; 2018.09.06 
; deal with saturation deficit, and figure out its relation with sgm_tota
; 2018.09.07
; mask land grids by modifying functions_read_file.ncl, add option landmask


begin

landmask		= False
ocnmask			= False

fdir 					= "/home/lyl/WORK4/cesm1_2_1/archive/"
;fappend 				= "/atm/hist/"
fappend 				= "/atm/hist/nobudget/"

case					= (/"FAMIPC5_f09f09_mac2_amip_outSGM","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM"/)

printVarSummary(case)

var 					= (/"sgm_tota","SWCF","LWCF","CLDLOW","CLDMED","CLDHGH","CLDTOT","N1","deltaq_sat","TGCLDLWP","netCRF"/)

tag 					= (/"nor","p4K"/)

nlat					= 192
nlon					= 288
ntime					= 60	; 1980-1984 -- 12*5=60 months




vars_all				= new((/dimsizes(var),dimsizes(case),ntime,nlat,nlon/),float)
vars_all!0				= "var"
vars_all&var				= fspan(0,dimsizes(var)-1,dimsizes(var))
vars_all!1				= "case"
vars_all&case				= fspan(0,dimsizes(case)-1,dimsizes(case))

do icase=0,dimsizes(case)-1


do ivar=0,dimsizes(var)-1

if(var(ivar).eq."sgm_tota")then
	sgm_tota 			=  read_data_monthly(fdir,fappend,case(icase),var(ivar),landmask,ocnmask)
	sgm_tota_900			= sgm_tota(:,{900},:,:) 
	vars_all(ivar,icase,:,:,:)	= sgm_tota_900
	; change units
	vars_all(ivar,icase,:,:,:)	= sqrt(abs(vars_all(ivar,icase,:,:,:)))*1e3
else if(var(ivar).eq."N1")then
	N1				= read_data_monthly(fdir,fappend,case(icase),var(ivar),landmask,ocnmask)
	N1_900				= N1(:,{900},:,:)
	vars_all(ivar,icase,:,:,:)	= N1_900
else if(var(ivar).eq."deltaq_sat")then
	deltaq_sat			= read_data_monthly(fdir,fappend,case(icase),var(ivar),landmask,ocnmask)
	deltaq_sat_900			= deltaq_sat(:,{900},:,:)
	vars_all(ivar,icase,:,:,:)	= deltaq_sat_900
	vars_all(ivar,icase,:,:,:)	= vars_all(ivar,icase,:,:,:)*1e3	; from kg/kg to g/kg
else if(var(ivar).eq."netCRF")then
	SWCF				= read_data_monthly(fdir,fappend,case(icase),"SWCF",landmask,ocnmask)
	LWCF				= read_data_monthly(fdir,fappend,case(icase),"LWCF",landmask,ocnmask)
	netCRF				= SWCF + LWCF
	vars_all(ivar,icase,:,:,:)	= (/netCRF/)
else
	vars_all(ivar,icase,:,:,:)	= read_data_monthly(fdir,fappend,case(icase),var(ivar),landmask,ocnmask)

	if(var(ivar).eq."CLDLOW".or.var(ivar).eq."CLDMED".or.var(ivar).eq."CLDHGH".or.var(ivar).eq."CLDTOT")then
	vars_all(ivar,icase,:,:,:)= vars_all(ivar,icase,:,:,:)*100.
	end if
	if(var(ivar).eq."TGCLDLWP")then
	vars_all(ivar,icase,:,:,:) = vars_all(ivar,icase,:,:,:)*1e3	; kg/kg to g/kg
	end if
end if
end if
end if
end if

end do ; do ivar=

end do ; do icase=

printVarSummary(vars_all)
printMinMax(vars_all,False)



;;;;; output files
if(landmask)then
	system("/bin/rm -f ./tmp-data/sgm-var-outdata-with-N1-ocn.nc")
	ncdf = addfile("./tmp-data/sgm-var-outdata-with-N1-ocn.nc","c")
else
	if(ocnmask)then
		system("/bin/rm -f ./tmp-data/sgm-var-outdata-with-N1-lnd.nc")
		ncdf = addfile("./tmp-data/sgm-var-outdata-with-N1-lnd.nc","c")
	else
		system("/bin/rm -f ./tmp-data/sgm-var-outdata-with-N1.nc")
		ncdf = addfile("./tmp-data/sgm-var-outdata-with-N1.nc","c")
	end if
end if

do ivar=0,dimsizes(var)-1

ncdf->$var(ivar)$		= vars_all(ivar,:,:,:,:)

end do ; do ivar=
end


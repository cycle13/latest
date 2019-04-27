
resolution					= "2deg"

fdir						= "/home/share3/lyl/work3/qinyi/mid-data/cfmip/"
case_append					= ".ANN.climo.nc"

;case_nor					= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_mac2_amip_outSGM"/)
;case_4K						= (/"FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM"/)

;case_nor					= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_mac2_amip_outSGM_MG10"/)
;case_4K						= (/"FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM_MG10"/)

;case_nor					= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_MG10_amip","FAMIPC5_f09f09_mac2_amip","FAMIPC5_f09f09_mac2_amip_outSGM","FAMIPC5_f09f09_mac2_amip_outSGM_MG10","FAMIPC5_f09f09_mac2_amip_outSGM_MG10"/)
;case_4K						= (/"FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_MG10_amip-p4K","FAMIPC5_f09f09_mac2_amip-p4K_debugoff_outSGM_MG10","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM_MG10","FAMIPC5_f09f09_mac2_amip-p4K_debugoff_outSGM_MG10"/)

if(resolution.eq."2deg")then
	case_nor		= (/"FAMIPC5_f19f19_MG15_amip","FAMIPC5_f19f19_mac2_amip_outSGM"/)
	case_4K			= (/"FAMIPC5_f19f19_MG15_amip-p4K","FAMIPC5_f19f19_mac2_amip-p4K_debug_outSGM"/)
else
	case_nor		= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_mac2_amip_outSGM_MG10"/)
	case_4K			= (/"FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM_MG10"/)
end if


case_all					= (/case_nor,case_4K/)

tag							= (/"nor","p4K"/)
;								  0       1     2        3      4      5      6       7
Vars						= (/"FLNT","FSNT","FLNTC","FSNTC","LWCF","SWCF","netCRF","TS"/)
Vars_dlt					= "dlt_"+Vars

if(resolution.eq."2deg")then
	nlat						= 96
	nlon						= 144
else
	nlat						= 192
	nlon						= 288
end if

vars_all					= new((/dimsizes(Vars),dimsizes(tag),dimsizes(case_nor),nlat,nlon/),float)
vars_all!0					= "var"
vars_all&var				= fspan(0,dimsizes(Vars)-1,dimsizes(Vars))
vars_all!1					= "tag"
vars_all&tag				= fspan(0,dimsizes(tag)-1,dimsizes(tag))
vars_all!2					= "case"
vars_all&case				= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))

do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1
do ivar=0,dimsizes(Vars)-1
	f1						= addfile(fdir+case_all(itag,icase)+case_append,"r")
	if(isfilevar(f1,Vars(ivar)))then
		vars_all(ivar,itag,icase,:,:)			= f1->$Vars(ivar)$(0,:,:)
	else
		if(Vars(ivar).eq."netCRF")then
			SWCF_tmp							= f1->SWCF(0,:,:)
			LWCF_tmp							= f1->LWCF(0,:,:)
			vars_all(ivar,itag,icase,:,:)		= (/SWCF_tmp + LWCF_tmp/)
		end if
	end if
end do ; do ivar=
end do ; do icase=
end do ; do itag=

; =================== global mean =====================
vars_all_aave							= new((/dimsizes(Vars),dimsizes(tag),dimsizes(case_nor)/),float)
vars_all_aave!0							= "var"
vars_all_aave&var						= fspan(0,dimsizes(Vars)-1,dimsizes(Vars))
vars_all_aave!1							= "tag"
vars_all_aave&tag						= fspan(0,dimsizes(tag)-1,dimsizes(tag))
vars_all_aave!2							= "case"
vars_all_aave&case						= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))


lat_here                                = vars_all&lat
rad                                     = 4.0*atan(1.0)/180.
clat                                    = cos(lat_here*rad)
clat!0                                  = "lat"
clat&lat                                = lat_here

latS									= -30.
latE									= 30.
lonS									= 0.
lonE									= 360.

vars_all_aave                     		= wgt_areaave_Wrap(vars_all(:,:,:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1.0,1) 
printVarSummary(vars_all_aave)
printMinMax(vars_all_aave,True)

; ============ difference between P4K and Nor =============
vars_all_aave_diff						= vars_all_aave(:,1,:)-vars_all_aave(:,0,:)
copy_VarCoords_1(vars_all_aave(var|:,case|:,tag|:),vars_all_aave_diff)
printVarSummary(vars_all_aave_diff)

print("Vars_dlt		    RH_MG15			PDF_MG15")
print(Vars_dlt+"		"+sprintf("%9.4f",vars_all_aave_diff(:,0))+"	"+sprintf("%9.4f",vars_all_aave_diff(:,1)))
;

;print("Vars_dlt		    RH_MG15		RH_MG10		  PDF		PDF_MG15		PDF_MG10	PDF_MG10_debugoff")
;print(Vars_dlt+"		"+sprintf("%9.4f",vars_all_aave_diff(:,0))+"	"+sprintf("%9.4f",vars_all_aave_diff(:,1))+"	"+sprintf("%9.4f",vars_all_aave_diff(:,2))+"		"+sprintf("%9.4f",vars_all_aave_diff(:,3))+"		"+sprintf("%9.4f",vars_all_aave_diff(:,4))+"		"+sprintf("%9.4f",vars_all_aave_diff(:,5)))
;
; =========== calculate climate sensitivity (Cess, 1990) ==========
Vars_out					= (/"G","Gc","lambda","lambda_c","lambda/lambda_c","1+netCRF/G"/)
Vars_out_fix				= Vars_out
do ivar = 0,dimsizes(Vars_out)-1
	lens					= strlen(Vars_out(ivar))
	Vars_out_fix(ivar)		= str_concat((/conform_dims(20-lens, " ",-1)/))
	Vars_out_fix(ivar)		= str_concat((/Vars_out(ivar),Vars_out_fix(ivar)/))
	delete(lens)
end do ; do ivar=

; G: forcing \dlt_LW - \dlt_SW
; Gc: forcing \dlt_LWc - \dlt_SWc
G							= vars_all_aave_diff(0,:) - vars_all_aave_diff(1,:)	; FLNT - FSNT
Gc							= vars_all_aave_diff(2,:) - vars_all_aave_diff(3,:)	; FLNTC - FSNTC

; lambda: \dlt_TS / G
; lambda_c: \dlt_Ts / Gc

; 2018-12-10: change \dlt_TS from constant 4K to varied value from calculation
;dlt_Ts						= 4.
dlt_Ts						= vars_all_aave_diff(7,:)
;print(dlt_Ts)
lambda		 				= dlt_Ts/G
lambda_c					= dlt_Ts/Gc

data						= (/G, Gc, lambda, lambda_c, lambda/lambda_c, (1+vars_all_aave_diff(6,:)/G) /)
printVarSummary(data)

print(Vars_out_fix+"	"+sprintf("%9.4f",data(:,0))+"	"+sprintf("%9.4f",data(:,1)))
exit

;print(Vars_out_fix+"	"+sprintf("%9.4f",data(:,0))+"	"+sprintf("%9.4f",data(:,1))+"	"+sprintf("%9.4f",data(:,2))+"		"+sprintf("%9.4f",data(:,3))+"		"+sprintf("%9.4f",data(:,4))+"		"+sprintf("%9.4f",data(:,5)))





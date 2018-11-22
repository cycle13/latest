load "./functions_read_file.ncl"
load "./functions_calculate.ncl"
load "./functions_plot.ncl"
;load "./functions_discretize.ncl"
load "./functions_regrid.ncl"
load "./functions_bony_method.ncl"

;;;;;processing steps
; 1. read data of model results: function read_data_monthly()
; 4. do bony-type binning: function bony_method():30S~30N
; 5. organize data into specific format: function do_pdf_var()

; notes for all options
; 1. do_pdf_var = "T" and do_pdf_omega = "T", -------- form pdf_sgm_only.nc
; 2. do_pdf_var = "F" and do_pdf_omega = "F", run the script again

; 2018.09.06
; deal with N1, deltaq_sat and figure out their relations with sgm_tota
; change nbinx and nbins from 25 to 50 ----- return to 25

; 2018.09.11
; add sat and unsat dimension to further check the effect from new cloud scheme

; 2018.09.23
; get netCRF directly from reading raw data, rather than using SWCF_bin and LWCF_bin.
; this is beneficial to the standard deviation calculation.


begin

landmask				= False
ocnmask					= False

do_read_data			= False
do_pdf_var 				= True
do_pdf_omega			= True

mask_sat				= False
mask_unsat				= False

sat					= (/mask_sat,mask_unsat/)

;colors 					= (/"red","purple","grey"/) ; PDF, PDF-p4K, difference
colors					= (/"cyan4","blue","purple","red"/)

option					= (/"p4K","4xCO2"/)

Vars 		= (/"sgm_tota","SWCF","LWCF","CLDLOW","CLDMED","CLDHGH","CLDTOT","N1","deltaq_sat","TGCLDLWP","netCRF"/)
Vars_bin	= (/"SWCF","LWCF","CLDLOW","CLDMED","CLDHGH","CLDTOT","N1","deltaq_sat","TGCLDLWP","netCRF"/)
;	 	      0       1      2          3         4         5         6        7        8	      9
strings		= Vars_bin
valmin		= (/-140, -20, 0, 0, 0, 0, 0, 0, 0, -80/)
valmax		= (/20,	60, 80, 80, 80, 80, 1, 0.5, 200, 10/)
units		= (/"(W/m~S~2~N~)", "(W/m~S~2~N~)", "(fraction)", "(fraction)", "(fraction)", "(fraction)", "(fraction)", "(g/kg)", "(g/m2)", "(W/m2)"/)




nlat					= 192
nlon					= 288
ntime					= 60	; 1980-1984 -- 12*5=60 months

case					= (/"FAMIPC5_f09f09_mac2_amip_outSGM","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM"/)

nbinx 					= 25
nbiny 					= 50
nx_min 					= 0.
nx_max 					= 0.5


if(do_read_data)then
vars_all				= new((/dimsizes(Vars),dimsizes(case),ntime,nlat,nlon/),float)
vars_all!0				= "var"
vars_all&var				= fspan(0,dimsizes(Vars)-1,dimsizes(Vars))
vars_all!1				= "case"
vars_all&case				= fspan(0,dimsizes(case)-1,dimsizes(case))

; to save variable with sat and unsat array dimension
vars_all_new				= new((/dimsizes(Vars),dimsizes(sat),dimsizes(case),ntime,nlat,nlon/),float)
vars_all_new!0				= "var"
vars_all_new&var			= fspan(0,dimsizes(Vars)-1,dimsizes(Vars))
vars_all_new!1				= "sat"
vars_all_new&sat			= fspan(0,dimsizes(sat)-1,dimsizes(sat))
vars_all_new!2				= "case"
vars_all_new&case			= fspan(0,dimsizes(case)-1,dimsizes(case))


;====== read data in
;infile1 = addfile("./tmp-data/sgm-var-outdata.nc","r")
if(landmask)then
	infile1 = addfile("./tmp-data/sgm-var-outdata-with-N1-ocn.nc","r")
else
	if(ocnmask)then
		infile1 = addfile("./tmp-data/sgm-var-outdata-with-N1-lnd.nc","r")
	else
		infile1 = addfile("./tmp-data/sgm-var-outdata-with-N1.nc","r")
	end if
end if
	

do ivar=0,dimsizes(Vars)-1
	vars_all(ivar,:,:,:,:)			 = infile1->$Vars(ivar)$
	if(Vars(ivar).eq."deltaq_sat")then
		deltaq_sat_use			= infile1->$Vars(ivar)$
	end if
end do ; do ivar=
printVarSummary(vars_all)

vars_all_new(0,0,:,:,:,:)		= vars_all(0,:,:,:,:) ; copy coordinates

if(mask_sat .or. mask_unsat)then
print("qinyi")
do ivar=0,dimsizes(Vars)-1
		vars_all_new(ivar,0,:,:,:,:)	= where(deltaq_sat_use.le.1e-2,vars_all(ivar,:,:,:,:),vars_all@_FillValue)
		vars_all_new(ivar,1,:,:,:,:)	= where(deltaq_sat_use.ge.1e-2,vars_all(ivar,:,:,:,:),vars_all@_FillValue)
end do ; do ivar=
else
	vars_all_new(:,0,:,:,:,:)		= vars_all(:,:,:,:,:)
	vars_all_new(:,1,:,:,:,:)		= vars_all(:,:,:,:,:)
end if

;sgm_all					= vars_all(0,:,:,:,:)
sgm_all_new					= vars_all_new(0,:,:,:,:,:)

printVarSummary(sgm_all_new)
printVarSummary(vars_all_new)

;wks	= gsn_open_wks("x11","test")
;
;res		= True
;res@cnFillOn	= True
;res@cnLinesOn	= False
;plot	= gsn_csm_contour_map(wks,vars_all(7,1,4,:,:),res)

end if ; do_read_data


if(do_pdf_var .and. do_read_data)then

sgm_vars_1D_bin					= new((/dimsizes(Vars_bin),dimsizes(sat),dimsizes(case),nbinx/),double) ; except for sgm, but include netCRF
sgm_vars_1D_bin!0				= "var"
sgm_vars_1D_bin&var				= fspan(0,dimsizes(Vars_bin)-1,dimsizes(Vars_bin))
sgm_vars_1D_bin!1				= "sat"
sgm_vars_1D_bin&sat				= fspan(0,dimsizes(sat)-1,dimsizes(sat))

; standard deviation
sgm_vars_1D_bin_std				= sgm_vars_1D_bin


;		      0       1       2      3          4        5       6
;var 		= (/"SWCF","LWCF","CLDLOW","CLDMED","CLDHGH","CLDTOT","netCRF"/)

;outvar		= (/"sgm_SWCF","sgm_LWCF","sgm_CLDLOW","sgm_CLDMED","sgm_CLDHGH","sgm_CLDTOT","sgm_netCRF"/)
outVars		= Vars_bin
outVars_std	= Vars_bin

printMinMax(vars_all(8,:,:,:,:),True)

do ivar=0,dimsizes(Vars_bin)-1
	outVars(ivar) = "sgm_"+Vars_bin(ivar)
	outVars_std(ivar) = "sgm_"+Vars_bin(ivar)+"_std"
	print(outVars(ivar))
	result = pdf_sgm_var_wo_obs(vars_all_new(0,0,:,:,:,:),vars_all_new(ivar+1,0,:,:,:,:),nbinx,nbiny,nx_min,nx_max)	; mask_sat
	sgm_vars_1D_bin(ivar,0,:,:) = result[0]
	sgm_vars_1D_bin_std(ivar,0,:,:) = result[1]
	result = pdf_sgm_var_wo_obs(vars_all_new(0,1,:,:,:,:),vars_all_new(ivar+1,1,:,:,:,:),nbinx,nbiny,nx_min,nx_max)	; mask_unsat
	sgm_vars_1D_bin(ivar,1,:,:) = result[0]
	sgm_vars_1D_bin_std(ivar,1,:,:) = result[1]
	delete(result)
end do ; do ivar=

printVarSummary(sgm_vars_1D_bin)
printMinMax(sgm_vars_1D_bin,True)

printVarSummary(sgm_vars_1D_bin_std)
printMinMax(sgm_vars_1D_bin_std,True)

if(landmask)then
	system("/bin/rm -f ./tmp-data/pdf-sgm-var-ocn.nc")
	outfile1 = addfile("./tmp-data/pdf-sgm-var-ocn.nc","c")
else
	if(ocnmask)then
		system("/bin/rm -f ./tmp-data/pdf-sgm-var-lnd.nc")
		outfile1 = addfile("./tmp-data/pdf-sgm-var-lnd.nc","c")
	else
		system("/bin/rm -f ./tmp-data/pdf-sgm-var.nc")
		outfile1 = addfile("./tmp-data/pdf-sgm-var.nc","c")
	end if
end if
		

do ivar=0,dimsizes(outVars)-1
	outfile1->$outVars(ivar)$			= sgm_vars_1D_bin(ivar,:,:,:)
	outfile1->$outVars_std(ivar)$			= sgm_vars_1D_bin_std(ivar,:,:,:)
end do ; do ivar=

end if ; do_pdf_var

; omega pdf distribution

if(do_read_data .and. do_pdf_omega)then

sgm_PDF_bin			= new((/dimsizes(sat),dimsizes(case),nbinx/),double)
sgm_PDF_bin!0			= "sat"
sgm_PDF_bin&sat			= fspan(0,dimsizes(sat)-1,dimsizes(sat))


result = pdf_regimevar_only(sgm_all_new(0,:,:,:,:),nbinx,nx_min,nx_max)
printVarSummary(result)
sgm_PDF_bin(0,:,:)	= result
delete(result)

result = pdf_regimevar_only(sgm_all_new(1,:,:,:,:),nbinx,nx_min,nx_max)
sgm_PDF_bin(1,:,:)	= result
delete(result)
printVarSummary(sgm_PDF_bin)

if(landmask)then
	system("/bin/rm -f ./tmp-data/pdf-sgm-only-ocn.nc")
	outfile1 = addfile("./tmp-data/pdf-sgm-only-ocn.nc","c")
else
	if(ocnmask)then
		system("/bin/rm -f ./tmp-data/pdf-sgm-only-lnd.nc")
		outfile1 = addfile("./tmp-data/pdf-sgm-only-lnd.nc","c")
	else
		system("/bin/rm -f ./tmp-data/pdf-sgm-only.nc")
		outfile1 = addfile("./tmp-data/pdf-sgm-only.nc","c")
	end if
end if

outfile1->sgm_PDF_bin = sgm_PDF_bin
delete(outfile1)
exit
end if ; do_pdf_omega

;=================== Part 2: start to read data and prepare plotting ------------

;======= 2.1 SGM and Vars

if(landmask)then
	infile3 = addfile("./tmp-data/pdf-sgm-var-ocn.nc","r")
else
	if(ocnmask)then
		infile3 = addfile("./tmp-data/pdf-sgm-var-lnd.nc","r")
	else
		infile3 = addfile("./tmp-data/pdf-sgm-var.nc","r")
	end if
end if
	
Vars2				= Vars_bin
Vars2_std			= Vars_bin


varsdata			= new((/dimsizes(Vars2),dimsizes(sat),dimsizes(case),nbinx/),double)
varsdata!0			= "var"
varsdata&var			= fspan(0,dimsizes(Vars2)-1,dimsizes(Vars2))
varsdata!1			= "sat"
varsdata&sat			= fspan(0,dimsizes(sat)-1,dimsizes(sat))

varsdata_std		= varsdata

do ivar=0,dimsizes(Vars2)-1
Vars2(ivar)			= "sgm_"+Vars_bin(ivar)
varsdata(ivar,:,:,:)		= infile3->$Vars2(ivar)$
Vars2_std(ivar)		= "sgm_"+Vars_bin(ivar)+"_std"
varsdata_std(ivar,:,:,:) 	= infile3->$Vars2_std(ivar)$
end do ; do ivar=

printVarSummary(varsdata)
printVarSummary(varsdata_std)

varsdata_diff			= new((/dimsizes(Vars2),dimsizes(sat),nbinx/),double)
varsdata_diff			= varsdata(:,:,1,:)-varsdata(:,:,0,:)	; p4K-cntl

copy_VarCoords_1(varsdata(var|:,sat|:,binx_center|:,case|:),varsdata_diff)

printVarSummary(varsdata_diff)
printMinMax(varsdata_diff,True)

; ======== 2.2 SGM itself

if(landmask)then
	infile4 = addfile("./tmp-data/pdf-sgm-only-ocn.nc","r")
else
	if(ocnmask)then
		infile4 = addfile("./tmp-data/pdf-sgm-only-lnd.nc","r")
	else
		infile4 = addfile("./tmp-data/pdf-sgm-only.nc","r")
	end if
end if

sgm_PDF_bin 			= infile4->sgm_PDF_bin
printVarSummary(sgm_PDF_bin)

; get the difference to make it clear
sgm_PDF_bin_diff		= sgm_PDF_bin(:,1,:)-sgm_PDF_bin(:,0,:)
copy_VarAtts(sgm_PDF_bin,sgm_PDF_bin_diff)
copy_VarCoords_1(sgm_PDF_bin(sat|:,bin_center|:,case|:),sgm_PDF_bin_diff)
printVarSummary(sgm_PDF_bin_diff)
printMinMax(sgm_PDF_bin_diff,True)

;============= 2.3 further parts about calculation

; 2018-04-18 18:28:31
; get the weighting varsdata by multiplying the sgm_PDF_bin itself
varsdata_wgt				= varsdata
do ivar=0,dimsizes(Vars2)-1
varsdata_wgt(ivar,:,:,:)		= varsdata(ivar,:,:,:)*sgm_PDF_bin/100.
end do
printVarSummary(varsdata_wgt)
printMinMax(varsdata_wgt,False)

varsdata_wgt_diff			= new((/dimsizes(Vars2),dimsizes(sat),nbinx/),double)
varsdata_wgt_diff			= (/varsdata_wgt(:,:,1,:)-varsdata_wgt(:,:,0,:)/)
copy_VarCoords_1(varsdata_wgt(var|:,sat|:,binx_center|:,case|:),varsdata_wgt_diff)

printVarSummary(varsdata_wgt_diff)
printMinMax(varsdata_wgt_diff,True)

; 2018-04-19 11:01:20
; get Cw*\Delta_Pw and Pw*\Delta_Cw
; Cw
Cw		 			= new((/dimsizes(Vars2),dimsizes(sat),nbinx/),double)
Cw					= varsdata(:,:,0,:)
; Pw
Pw					= new((/dimsizes(sat),nbinx/),double)
Pw					= sgm_PDF_bin(:,0,:)
; \Delta_Cw
Delta_Cw				= Cw
Delta_Cw				= (/varsdata(:,:,1,:)-varsdata(:,:,0,:)/)
; \Delta_Pw
Delta_Pw				= Pw
Delta_Pw				= (/sgm_PDF_bin(:,1,:)-sgm_PDF_bin(:,0,:)/)
; PwCw, Cw*\Delta_Pw and Pw*\Delta_Cw
PwCw					= Cw
Cw_Delta_Pw				= Cw
Pw_Delta_Cw				= Cw
Delta_Pw_Delta_Cw			= Cw

do ivar=0,dimsizes(Vars2)-1
PwCw(ivar,:,:)				= Pw(:,:)*Cw(ivar,:,:)/100.
Cw_Delta_Pw(ivar,:,:)			= Cw(ivar,:,:)*Delta_Pw/100.
Pw_Delta_Cw(ivar,:,:)			= Pw(:,:)*Delta_Cw(ivar,:,:)/100.
Delta_Pw_Delta_Cw(ivar,:,:)		= Delta_Pw(:,:)*Delta_Cw(ivar,:,:)/100.
end do ; do ivar=

printVarSummary(Delta_Pw_Delta_Cw)
printMinMax(Delta_Pw_Delta_Cw,False)

CwPw_all			 	= new((/dimsizes(Vars2),3,dimsizes(sat),nbinx/),double)
CwPw_all!1				= "component"
CwPw_all&component			= (/0,1,2/)
CwPw_all(:,0,:,:)			= Cw_Delta_Pw
CwPw_all(:,1,:,:)			= Pw_Delta_Cw
CwPw_all(:,2,:,:)			= Delta_Pw_Delta_Cw
printVarSummary(CwPw_all)

; ============ averaged over all sgm
CwPw_all_name			= (/"dynamic","thermodynamic","co-variation"/)
CwPw_all_avg			= dim_sum_Wrap(CwPw_all)

Vars_out                    = Vars2
Vars_out_fix                = Vars_out

do ivar = 0,dimsizes(Vars_out)-1
    lens                    	  = strlen(Vars_out(ivar))
    Vars_out_fix(ivar)      = str_concat((/conform_dims(15-lens, " ",-1)/))
    Vars_out_fix(ivar)      = str_concat((/Vars_out(ivar),Vars_out_fix(ivar)/))
    delete(lens)
end do ; do ivar=

print(CwPw_all_name(0)+"		"+CwPw_all_name(1)+"		"+CwPw_all_name(2))
print(Vars_out_fix+"	"+sprintf("%9.5f",CwPw_all_avg(:,0,0))+"	"+sprintf("%9.5f",CwPw_all_avg(:,1,0))+"	"+sprintf("%9.5f",CwPw_all_avg(:,2,0)))


plot_type			= "pdf"

if(mask_sat.or.mask_unsat)then
wks  				= gsn_open_wks (plot_type,"./figure/sgm-var-CwPw_all-with-N1-sat-unsat-std")
else
wks  				= gsn_open_wks (plot_type,"./figure/sgm-var-CwPw_all-with-N1-std")
end if

plot 				= new(12,graphic)
plotdiff 			= new(12,graphic)

plot_wgt			= new(12,graphic)
plot_wgt_diff		= new(12,graphic)

plot_CwPw			= new(12,graphic)
plot_PwCw			= new(12,graphic)
plot_DPwDCw			= new(12,graphic)
plot_CwPw_all		= new(12,graphic)

nplots				= 7

legend_0			= new(nplots,graphic)
legend_1			= new(nplots,graphic)


res  							= True

res@gsnFrame 					= False
res@gsnDraw 					= False

res@vpWidthF 					= 0.7
res@vpHeightF 					= 0.4

res@trXMinF 					= 0.
res@trXMaxF 					= 1.

res@tiXAxisString 				= "subgrid-scale variance"

res@xyLineThicknessF 			= 4.
res@xyDashPattern 				= 0

fontheight						= 0.03
res@gsnStringFontHeightF		= fontheight
res@tmXBLabelFontHeightF		= fontheight
res@tmYLLabelFontHeightF		= fontheight
res@tiXAxisFontHeightF			= fontheight
res@tiYAxisFontHeightF			= fontheight

res@trXMaxF						= 0.5

vector							= (/0,1/)
vector_diff						= (/0,2/)

res1							= res
res1@xyLineColors 				= colors

res2							= res
res2@tiYAxisString				= ""
res2@xyLineColors				= colors(vector_diff)

if(mask_sat.or.mask_unsat)then
	; reshape: (0,0) (0,1) (1,0) (1,1)  --> (0) (1) (2) (3)
	tmp								= reshape(sgm_PDF_bin,(/dimsizes(sat)*dimsizes(case),nbinx/))
else ; dont separate sat and unsat situation
	tmp								= sgm_PDF_bin(0,:,:)
end if

plot0 							= gsn_csm_xy(wks,sgm_PDF_bin@bin_center,tmp,res1)
delete(tmp)

if(mask_sat.or.mask_unsat)then
	tmp							= sgm_PDF_bin_diff
else
	tmp							= sgm_PDF_bin_diff(0,:)
end if

plot1 							= gsn_csm_xy(wks,sgm_PDF_bin_diff@bin_center,tmp,res2)
delete(tmp)


res3							= res
;res3@xyLineColors 				= colors

res4							= res
res4@xyLineColors				= colors(vector_diff)
res4@gsnYRefLine				= 0.

res5							= res
if(mask_sat .or. mask_unsat)then
res5@xyLineColors				= (/"red","red","blue","blue","cyan4","cyan4"/)
res5@xyDashPatterns				= (/14,0,14,0,14,0/)
res5@xyLineThicknesses			= (/2,3,2,3,2,3/)
else
res5@xyLineColors				= (/"red","blue","cyan4"/)
res5@xyDashPatterns				= (/0,0,0/)
res5@xyLineThicknesses			= (/2,2,2/)
end if

res5@gsnYRefLine				= 0.


do ivar=0,dimsizes(Vars2)-1
	res3@trYMinF 					= valmin(ivar)
	res3@trYMaxF 					= valmax(ivar)
	res3@tiYAxisString 				= strings(ivar)+" "+units(ivar)
	res4@tiYAxisString 				= strings(ivar)+" "+units(ivar)
	res5@tiYAxisString 				= strings(ivar)+" "+units(ivar)
	
	if(mask_sat .or. mask_unsat)then
	tmp								= reshape(varsdata(ivar,:,:,:),(/dimsizes(sat)*dimsizes(case),nbinx/))
	tmp_std							= reshape(varsdata_std(ivar,:,:,:),(/dimsizes(sat)*dimsizes(case),nbinx/))
	else
	tmp								= reshape(varsdata(ivar,0,:,:),(/dimsizes(case),nbinx/))
	tmp_std							= reshape(varsdata_std(ivar,0,:,:),(/dimsizes(case),nbinx/))
	end if
	
	; initialize one plot as the baseline for other plots. Others plots will overlay on this one.
	res3@xyLineColors 				= colors(0)
	plot_int						= gsn_csm_xy(wks,varsdata&binx_center,tmp(0,:),res3)

	if(mask_sat .or. mask_unsat)then
	nfigs							= dimsizes(sat)*dimsizes(case)
	else
	nfigs							= dimsizes(case)
	end if
	do ifig = 0,nfigs-1
	
		res3@xyLineColors 	= colors(ifig)
		plot_tmp			= gsn_csm_xy (wks,varsdata&binx_center,tmp(ifig,:), res3)

		;************************************************
		; define a polygon centered with width of 1 sigma
		;************************************************
		xp    = new( (/2*nbinx/), float )
		yp    = new( (/2*nbinx/), float )
		
		do k=0,nbinx-1
			dx					= tmp_std(ifig,k)
			yp(k)          		= tofloat(tmp(ifig,k)+dx)
			xp(k)          		= tofloat(varsdata&binx_center(k))
			xp(2*nbinx-1-k) 	= tofloat(varsdata&binx_center(k))
			yp(2*nbinx-1-k) 	= tofloat(tmp(ifig,k)-dx)
		end do
		
		gsres                  	= True                        ; poly res
		gsres@tfPolyDrawOrder   = "Predraw"                   ; draw this first
		gsres@gsFillColor       = colors(ifig)                 ; color chosen
		gsres@gsFillOpacityF	= 0.1
		
		dummy 					= gsn_add_polygon (wks,plot_tmp,xp,yp,gsres)
		str1 					= unique_string("dummy")
		plot_int@$str1$ 		= dummy
		overlay(plot_int,plot_tmp)
	
	end do ; do ifig=
	plot(ivar) 					= plot_int
	
	delete([/tmp,tmp_std/])
	
	if(mask_sat .or. mask_unsat)then
	tmp1							= varsdata_diff(ivar,:,:)
	tmp2							= varsdata_wgt_diff(ivar,:,:)
	else
	tmp1							= varsdata_diff(ivar,0,:)
	tmp2							= varsdata_wgt_diff(ivar,0,:)
	end if

	plotdiff(ivar) 				= gsn_csm_xy (wks,varsdata_diff&binx_center,tmp1, res4)
	plot_wgt_diff(ivar) 		= gsn_csm_xy (wks,varsdata_wgt_diff&binx_center,tmp2, res4)

	if(mask_sat .or. mask_unsat)then
	tmp							= reshape(CwPw_all(ivar,:,:,:),(/dimsizes(sat)*3,nbinx/))
	else
	tmp							= reshape(CwPw_all(ivar,:,0,:),(/3,nbinx/))
	end if

	plot_CwPw_all(ivar)			= gsn_csm_xy (wks,CwPw_all&binx_center,tmp, res5)
	delete(tmp)
	
	end do ; do ivar=

;------------------------------Legend preparing---------------------------------------
;Attach a legend
lgres 							= True
lgres@lgItemType 				= "Lines"
lgres@lgLabelFontHeightF 		= .20
lgres@vpWidthF 					= 0.25
lgres@vpHeightF 				= 0.15
;lgres@lgPerimThicknessF 		= 1.0
lgres@lgMonoDashIndex 			= True
lgres@lgDashIndex 				= 0
lgres@lgLineThicknessF 			= res@xyLineThicknessF

lgres1							= lgres

lgres1@lgLineColors 			= colors
;labels_1 						= (/"unsat: PDF","sat: PDF","unsat: PDF_p4K","sat: PDF_p4K"/)
if(mask_sat .or. mask_unsat)then
labels_1 						= (/"unsat: PDF","unsat: PDF_p4K","sat: PDF","sat: PDF_p4K"/)
else
labels_1 						= (/"PDF","PDF_p4K"/)
end if


lgres2							= lgres
lgres2@lgLineColors				= colors(vector_diff)
if(mask_sat .or. mask_unsat)then
labels_2 						= (/"unsat: PDF_p4K-PDF","sat: PDF_p4K-PDF"/)
else
labels_2 						= (/"PDF_p4K-PDF"/)
end if


; for omega-PDF only
legend0							= gsn_create_legend(wks,dimsizes(labels_1),labels_1,lgres1)
legend1							= gsn_create_legend(wks,dimsizes(labels_2),labels_2,lgres2)

; for omega-var plots
do ivar=0,nplots-1
	legend_0(ivar)				= gsn_create_legend(wks,dimsizes(labels_1),labels_1,lgres1)
	legend_1(ivar)				= gsn_create_legend(wks,dimsizes(labels_2),labels_2,lgres2)
end do ; do ivar=

;---------------------------------------------
amres 							= True
amres@amJust 					= "TopRight"
amres@amParallelPosF 			= 0.5
amres@amOrthogonalPosF 			= -0.5

; for omega-PDF only
annoid0 						= gsn_add_annotation(plot0,legend0,amres)
annoid1 						= gsn_add_annotation(plot1,legend1,amres)

; for omega-var plots

; combine all components into one array "CwPw_all", so terms above aren't plotted anymore.
;tmp_plot						= (/plot(0),plotdiff(2)/)
tmp_plot						= (/plot(0),plot_wgt_diff(0)/)
legends							= (/legend_0(0),legend_1(0)/)

annoid							= new(dimsizes(tmp_plot),graphic)
do ivar=0,dimsizes(tmp_plot)-1
	annoid(ivar) 				= gsn_add_annotation(tmp_plot(ivar),legends(ivar),amres)
end do ; do ivar=

; for CwPw_all
lgres3							= lgres
lgres3@vpWidthF 				= 0.30
lgres3@vpHeightF 				= 0.20

lgres3@lgMonoDashIndex			= False
lgres3@lgMonoLineThickness		= False
lgres3@lgLineColors				= res5@xyLineColors
lgres3@lgDashIndexes			= res5@xyDashPatterns
;lgres3@lgLineThicknessF		= res5@xyLineThicknessF
lgres3@lgLineThicknesses		= res5@xyLineThicknesses

if(mask_sat .or. mask_unsat)then
labels_3 						= (/"unsat: PDF-Cw_DPw","sat: PDF-Cw_DPw",\
				    			"unsat: PDF-Pw_DCw","sat: PDF-Pw_DCw",\
				    			"unsat: PDF-DCw_DPw","sat: PDF-DCw_DPw"/)
else
labels_3 						= (/"PDF-Cw_DPw",\
				    			"PDF-Pw_DCw",\
				    			"PDF-DCw_DPw"/)
end if

legend3							= gsn_create_legend(wks,dimsizes(labels_3),labels_3,lgres3)
annoid3							= gsn_add_annotation(plot_CwPw_all(0),legend3,amres)


;------------------------Panel plots----------------------------------------------
resP 							= True
;resP@gsnPanelRowSpec 			= True
resP@gsnPanelXWhiteSpacePercent	= 3.
resP@gsnPanelXWhiteSpacePercent	= 5.
resP@gsnPanelRight				= 0.95
resP@gsnPanelLeft				= 0.05

; for omega-PDF only
gsn_panel(wks,(/plot0,plot1/),(/1,2/),resP)

; for omega-var plots
idx								= (/0,1,2,3,6,5,7,8,9/)
;
;titles							= (/"varpdf","varpdf_diff","varpdf_wgt_diff","CwPw_all"/)
;panel_plots		 				= (/plot(idx),plotdiff(idx),plot_wgt_diff(idx),plot_CwPw_all(idx)/)

titles							= (/"varpdf","varpdf_wgt_diff","CwPw_all"/)
panel_plots		 				= (/plot(idx),plot_wgt_diff(idx),plot_CwPw_all(idx)/)


do i=0,dimsizes(titles)-1
	resP@gsnPanelMainString	= titles(i)
	gsn_panel(wks,panel_plots(i,:),(/3,3/),resP)
end do

exit
end


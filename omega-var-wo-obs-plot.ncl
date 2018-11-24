load "./functions_read_file.ncl"
load "./functions_calculate.ncl"
load "./functions_plot.ncl"
;load "./functions_discretize.ncl"
load "./functions_regrid.ncl"
load "./functions_bony_method.ncl"

;;;;;processing steps
; 1. read data of model results: function read_data_monthly()
; 5. organize data into specific format: function do_pdf_var()

; 2018.09.24
; although I added the difference PDF for omega-var and omega-only, I think
; it is not reasonable to use the difference PDF for omega-only, because it
; cannot check the regime change. However, it would be reasonable to use 
; difference PDF for oemga-var.
; 2018.09.26 
; revise the method for standard deviation. weighted standard deviation.
; 2018.10.22
; the nx_min and nx_max are important to the regime-sorting result
; nx_min = -100 would make the \delta(PwCw) less than the real tropical mean
; of cloud variable. Extending it to 140 makes it closer to the real value.
; however, there is still 0.1 which cannot be matched well with real value.
; In addition, the co-variation is too large compared to the result shown 
; in Bony et al. (2005). I should double check it.


begin

do_read_data			= False
do_pdf_var 				= True
do_pdf_omega			= True

nbinx 					= 25
nbiny 					= 50
nx_min 					= -140.
nx_max 					= 140.

landmask				= False
ocnmask					= False

;tag 					= (/"nor","p4K","4xCO2"/)
tag 					= (/"nor","p4K"/)

; cases with MG10 output
;case_nor	= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_mac2_amip_outSGM","FAMIPC5_f09f09_mac2_amip_outSGM_MG10"/)
;case_p4K	= (/"FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM_MG10"/)

; cases with high sensitivity
;case_nor	= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_mac2_amip","FAMIPC5_f09f09_mac2_amip_outSGM_MG10"/)
;case_p4K	= (/"FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_mac2_amip-p4K_debug","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM_MG10"/)

; cases with normal output and MG15 to compare with default simulation.
case_nor	= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_mac2_amip_outSGM"/)
case_p4K	= (/"FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM"/)


;case_4xCO2	= (/"FAMIPC5_f09f09_MG15_amip-4xCO2","FAMIPC5_f09f09_mac2_amip-4xCO2"/)

;case_all	= (/case_nor, case_p4K, case_4xCO2/) ; 3 x 2
case_all	= (/case_nor, case_p4K/) ; 2 x 2

; ====== used for test ===============
;Vars		= (/"OMEGA","SWCF","LWCF"/)
;Vars_bin	= (/"SWCF","LWCF"/)
; ====== used for test ===============

Vars 		= (/"OMEGA","SWCF","LWCF","CLDLOW","CLDMED","CLDHGH","CLDTOT","TGCLDLWP","netCRF","sgm_tota","N1","deltaq_sat","sgm_tota_uns","sgm_tota_sat"/)
;				  0      1        2        3       4        5         6         7         8       9          10        11            12
Vars_bin	= (/"SWCF","LWCF","CLDLOW","CLDMED","CLDHGH","CLDTOT","TGCLDLWP","netCRF","sgm_tota","N1","deltaq_sat","sgm_tota_uns","sgm_tota_sat"/)

;Vars 		= (/"OMEGA","SWCF","LWCF","CLDLOW","CLDMED","CLDHGH","CLDTOT","TGCLDLWP","netCRF"/)
;;				  0      1        2        3       4        5         6         7         8       9          10        11            12
;Vars_bin	= (/"SWCF","LWCF","CLDLOW","CLDMED","CLDHGH","CLDTOT","TGCLDLWP","netCRF"/)


; ------- setting for plots
stringsA	= Vars_bin ; except for omega itself, and add Net CRF
valminA		= (/-140,	-20,	0,		  0,	   0,		0 ,		0,        -60., 	0,	       0,     0,	     0.,           0./)
valmaxA		= (/20,		100,	60,		 90,	   90,		90,		150,      0.,   	0.3,       1,     0.1,         0.3,          0.5/)
unitsA		= (/"(W/m~S~2~N~)",	"(W/m~S~2~N~)",	"(fraction)",	"(fraction)",	"(fraction)",	"(fraction)","(g/m3)",	"(W/m~S~2~N~)", " ", "(fraction)", "(g/kg)", "", ""/)

use_idx		= (/0,1,7,2,6,8/)
Vars_use	= Vars_bin(use_idx)
strings		= stringsA(use_idx)
valmin		= valminA(use_idx)
valmax		= valmaxA(use_idx)
units		= unitsA(use_idx)

colors 		= (/"blue","red","cyan","darkseagreen","deepskyblue","firebrick","cyan3","darkseagreen3"/)

; ------ setting for global grid information
nlat		= 192
nlon		= 288
ntime		= 60	; 1980-1984 -- 12*5=60 months

; ================================ All start here ========================================

if(do_read_data)then
print("read data starting here...")
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



if(landmask)then
	infile = addfile("./tmp-data/omega-var-outdata-wo-obs-ocn.nc","r")
else
	if(ocnmask)then
		infile = addfile("./tmp-data/omega-var-outdata-wo-obs-lnd.nc","r")
	else
		infile = addfile("./tmp-data/omega-var-outdata-wo-obs.nc","r")
	end if
end if
	
; variable names in input file.
Vars_out						= Vars
Vars_out_diff					= Vars

do ivar=0,dimsizes(Vars)-1
	do itag=0,dimsizes(tag)-1
			Vars_out(ivar)				= Vars(ivar)+"_"+tag(itag)
			print(Vars_out(ivar))
			vars_all(ivar,itag,:,:,:,:)		= infile->$Vars_out(ivar)$
	end do ; do itag=

	Vars_out_diff(ivar)					= Vars(ivar)+"_diff"
	vars_all_diff(ivar,:,:,:,:)			= infile->$Vars_out_diff(ivar)$
end do ; do ivar=

;printVarSummary(vars_all)
end if ; do_read_data

;=================== do_pdf_var ==========================
if(do_read_data .and. do_pdf_var)then
print("Do PDF binning of all variables....")
omega_vars_1D_bin			= new((/dimsizes(Vars_bin),dimsizes(tag),dimsizes(case_nor),nbinx/),double) ; except for omega itself
omega_vars_1D_bin!0			= "var"
omega_vars_1D_bin&var			= fspan(0,dimsizes(Vars_bin)-1,dimsizes(Vars_bin))
omega_vars_1D_bin!1			= "tag"
omega_vars_1D_bin&tag			= fspan(0,dimsizes(tag)-1,dimsizes(tag))
omega_vars_1D_bin!2			= "case"
omega_vars_1D_bin&case			= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))

omega_vars_1D_bin_std		= omega_vars_1D_bin

omega_vars_1D_bin_diff		= omega_vars_1D_bin(:,0,:,:)
omega_vars_1D_bin_diff_std	= omega_vars_1D_bin(:,0,:,:)


do ivar=1,dimsizes(Vars)-1 ; except for omega 
do itag=0,dimsizes(tag)-1
	print(Vars(ivar))
	result 					= pdf_omega_var_wo_obs(vars_all(0,itag,:,:,:,:),vars_all(ivar,itag,:,:,:,:),nbinx,nbiny,nx_min,nx_max)
	omega_vars_1D_bin(ivar-1,itag,:,:) 		= result[0]
	omega_vars_1D_bin_std(ivar-1,itag,:,:)	= result[1]
	delete(result)
end do ; do itag=
	print(Vars(ivar))
	result 					= pdf_omega_var_wo_obs(vars_all_diff(0,:,:,:,:),vars_all_diff(ivar,:,:,:,:),nbinx,nbiny,nx_min,nx_max)
	omega_vars_1D_bin_diff(ivar-1,:,:) 			= result[0]
	omega_vars_1D_bin_diff_std(ivar-1,:,:)		= result[1]
	delete(result)

end do ; do ivar=
printVarSummary(omega_vars_1D_bin_diff)

if(landmask)then
	system("/bin/rm -f ./tmp-data/pdf-omega-var-wo-obs-ocn.nc")
	outfile1 = addfile("./tmp-data/pdf-omega-var-wo-obs-ocn.nc","c")
else
	if(ocnmask)then
		system("/bin/rm -f ./tmp-data/pdf-omega-var-wo-obs-lnd.nc")
		outfile1 = addfile("./tmp-data/pdf-omega-var-wo-obs-lnd.nc","c")
	else
		system("/bin/rm -f ./tmp-data/pdf-omega-var-wo-obs.nc")
		outfile1 = addfile("./tmp-data/pdf-omega-var-wo-obs.nc","c")
	end if
end if


outVars					= Vars_bin
outVars_std				= Vars_bin
outVars_diff			= Vars_bin
outVars_diff_std		= Vars_bin

do ivar=0,dimsizes(Vars_bin)-1 ; except for omega itself
	do itag=0,dimsizes(tag)-1
		outVars(ivar)					= "omega_"+Vars_bin(ivar)+"_"+tag(itag)
		outVars_std(ivar)				= "omega_"+Vars_bin(ivar)+"_"+tag(itag)+"_std"
		outfile1->$outVars(ivar)$		= omega_vars_1D_bin(ivar,itag,:,:)
		outfile1->$outVars_std(ivar)$	= omega_vars_1D_bin_std(ivar,itag,:,:)
	end do ; do itag=
	outVars_diff(ivar)						= "omega_"+Vars_bin(ivar)+"_diff"
	outVars_diff_std(ivar)					= "omega_"+Vars_bin(ivar)+"_diff_std"
	outfile1->$outVars_diff(ivar)$			= omega_vars_1D_bin_diff(ivar,:,:)
	outfile1->$outVars_diff_std(ivar)$		= omega_vars_1D_bin_diff_std(ivar,:,:)

end do ; do ivar=
;exit
end if ; do_pdf_var

;=================== do_pdf_omega ==================================

; omega pdf distribution

if(do_read_data .and. do_pdf_omega)then
print("Do PDF for omega itself...")
omega_PDF_bin			= new((/dimsizes(tag),dimsizes(case_nor),nbinx/),double) 
omega_PDF_bin!0			= "tag"
omega_PDF_bin&tag		= fspan(0,dimsizes(tag)-1,dimsizes(tag))
omega_PDF_bin!1			= "case"
omega_PDF_bin&case		= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))

omega_PDF_bin_diff		= omega_PDF_bin(0,:,:)

if(landmask)then
	system("/bin/rm -f ./tmp-data/pdf-omega-only-wo-obs-ocn.nc")
	outfile1 = addfile("./tmp-data/pdf-omega-only-wo-obs-ocn.nc","c")
else
	if(ocnmask)then
		system("/bin/rm -f ./tmp-data/pdf-omega-only-wo-obs-lnd.nc")
		outfile1 = addfile("./tmp-data/pdf-omega-only-wo-obs-lnd.nc","c")
	else
		system("/bin/rm -f ./tmp-data/pdf-omega-only-wo-obs.nc")
		outfile1 = addfile("./tmp-data/pdf-omega-only-wo-obs.nc","c")
	end if
end if

do itag=0,dimsizes(tag)-1
	result = pdf_regimevar_only(vars_all(0,itag,:,:,:,:),nbinx,nx_min,nx_max)
	omega_PDF_bin(itag,:,:) = result
	delete(result)
	outVars		= Vars(0)+"_"+tag(itag)
	outfile1->$outVars$ = omega_PDF_bin(itag,:,:)
end do ; do itag=

	; difference
	result = pdf_regimevar_only(vars_all_diff(0,:,:,:,:),nbinx,nx_min,nx_max)
	omega_PDF_bin_diff		= result
	outVars		= Vars(0)+"_diff"
	outfile1->$outVars$		= omega_PDF_bin_diff
	delete(outVars)
exit
end if ; do_pdf_omega

; ==============================================================
;=============== Part2: read the omega-var-PDF data
if(landmask)then
	infile3 = addfile("./tmp-data/pdf-omega-var-wo-obs-ocn.nc","r")
else
	if(ocnmask)then
		infile3 = addfile("./tmp-data/pdf-omega-var-wo-obs-lnd.nc","r")
	else
		infile3 = addfile("./tmp-data/pdf-omega-var-wo-obs.nc","r")
	end if
end if
	

varsdata			= new((/dimsizes(Vars_use),dimsizes(tag),dimsizes(case_nor),nbinx/),double) ; except for omega itself
varsdata!0			= "var"
varsdata&var			= fspan(0,dimsizes(Vars_use)-1,dimsizes(Vars_use))
varsdata!1			= "tag"
varsdata&tag			= fspan(0,dimsizes(tag)-1,dimsizes(tag))
varsdata!2			= "case"
varsdata&case			= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))
; standard deviation
varsdata_std		= varsdata

varsdata_diff		= varsdata(:,0,:,:)
varsdata_diff_std	= varsdata(:,0,:,:)

outVars				= Vars_use
outVars_std			= Vars_use
outVars_diff		= Vars_use
outVars_diff_std	= Vars_use

do ivar=0,dimsizes(Vars_use)-1 ; 
	do itag=0,dimsizes(tag)-1
		outVars(ivar)						= "omega_"+Vars_use(ivar)+"_"+tag(itag)
		varsdata(ivar,itag,:,:)				= infile3->$outVars(ivar)$
		outVars_std(ivar)					= "omega_"+Vars_use(ivar)+"_"+tag(itag)+"_std"
		varsdata_std(ivar,itag,:,:) 		= infile3->$outVars_std(ivar)$
	end do ; do itag=
	outVars_diff(ivar)						= "omega_"+Vars_use(ivar)+"_diff"
	outVars_diff_std(ivar)					= "omega_"+Vars_use(ivar)+"_diff_std"
	varsdata_diff(ivar,:,:)					= infile3->$outVars_diff(ivar)$
	varsdata_diff_std(ivar,:,:)				= infile3->$outVars_diff_std(ivar)$

end do ; do ivar=

printVarSummary(varsdata_diff)
printMinMax(varsdata_diff,True)
delete(outVars)

;=============== read the omega-only-PDF data
omega_PDF_bin			= new((/dimsizes(tag),dimsizes(case_nor),nbinx/),double) 
omega_PDF_bin!0			= "tag"
omega_PDF_bin&tag		= fspan(0,dimsizes(tag)-1,dimsizes(tag))
omega_PDF_bin!1			= "case"
omega_PDF_bin&case		= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))

if(landmask)then
	infile4 = addfile("./tmp-data/pdf-omega-only-wo-obs-ocn.nc","r")
else
	if(ocnmask)then
		infile4 = addfile("./tmp-data/pdf-omega-only-wo-obs-lnd.nc","r")
	else
		infile4 = addfile("./tmp-data/pdf-omega-only-wo-obs.nc","r")
	end if
end if
	
do itag=0,dimsizes(tag)-1
	outVars		= Vars(0)+"_"+tag(itag)
	omega_PDF_bin(itag,:,:)		= infile4->$outVars$
end do ; do itag=

printVarSummary(omega_PDF_bin)

; get the difference to make it clear
; 2018.09.24 Yes, I think this is more reasonable, rather than introducing
; difference-first, and binning-second.
omega_PDF_bin_diff		= omega_PDF_bin
omega_PDF_bin_diff(0,:,:)	= (/omega_PDF_bin(0,:,:)-omega_PDF_bin(0,:,:)/) ; nor-nor
omega_PDF_bin_diff(1,:,:)	= (/omega_PDF_bin(1,:,:)-omega_PDF_bin(0,:,:)/) ; p4K-nor
;omega_PDF_bin_diff(2,:,:)	= (/omega_PDF_bin(2,:,:)-omega_PDF_bin(0,:,:)/) ; 4xCO2-nor

printVarSummary(omega_PDF_bin_diff)
printMinMax(omega_PDF_bin_diff,True)

; 2018-04-18 18:28:31
; =================== get the weighting varsdata by multiplying the OMEGA_PDF_bin itself ===================
varsdata_wgt					= varsdata
do ivar=0,dimsizes(Vars_use)-1
varsdata_wgt(ivar,:,:,:)		= varsdata(ivar,:,:,:)*omega_PDF_bin/100.
end do
printVarSummary(varsdata_wgt)
printMinMax(varsdata_wgt,False)

varsdata_wgt_diff			= varsdata_wgt
varsdata_wgt_diff(:,0,:,:)		= (/varsdata_wgt(:,0,:,:)-varsdata_wgt(:,0,:,:)/)
varsdata_wgt_diff(:,1,:,:)		= (/varsdata_wgt(:,1,:,:)-varsdata_wgt(:,0,:,:)/)
;varsdata_wgt_diff(:,2,:,:)		= (/varsdata_wgt(:,2,:,:)-varsdata_wgt(:,0,:,:)/)
printVarSummary(varsdata_wgt_diff)
printMinMax(varsdata_wgt_diff,True)

; ================= get the dynamic, thermodynamic and covariation term from Bony-2005 method ===================
; 2018-04-19 11:01:20
; get Cw*\Delta_Pw and Pw*\Delta_Cw

; Cw & Pw
Cw		 				= new((/dimsizes(Vars_use),dimsizes(tag),dimsizes(case_nor),nbinx/),double)
Pw						= new((/dimsizes(tag),dimsizes(case_nor),nbinx/),double)
do icase = 0,dimsizes(case_nor)-1
Cw(:,:,icase,:)				= varsdata(:,:,icase,:)
Pw(:,icase,:)				= omega_PDF_bin(:,icase,:)
end do
printVarSummary(Cw)
printVarSummary(Pw)

; \Delta_Cw
Delta_Cw					= Cw
Delta_Cw(:,0,:,:)			= (/varsdata(:,0,:,:)-varsdata(:,0,:,:)/)
Delta_Cw(:,1,:,:)			= (/varsdata(:,1,:,:)-varsdata(:,0,:,:)/)
;Delta_Cw(:,2,:,:)			= (/varsdata(:,2,:,:)-varsdata(:,0,:,:)/)
printVarSummary(Delta_Cw)

; \Delta_Pw
Delta_Pw					= Pw
Delta_Pw(0,:,:)				= omega_PDF_bin(0,:,:)-omega_PDF_bin(0,:,:)
Delta_Pw(1,:,:)				= omega_PDF_bin(1,:,:)-omega_PDF_bin(0,:,:)
;Delta_Pw(2,:,:)				= omega_PDF_bin(2,:,:)-omega_PDF_bin(0,:,:)
printVarSummary(Delta_Pw)

; PwCw, Cw*\Delta_Pw and Pw*\Delta_Cw
PwCw					= Cw
Cw_Delta_Pw				= Cw
Pw_Delta_Cw				= Cw
Delta_Pw_Delta_Cw		= Cw

do ivar=0,dimsizes(Vars_use)-1
PwCw(ivar,:,:,:)					= Pw(:,:,:)*Cw(ivar,:,:,:)/100.
Cw_Delta_Pw(ivar,:,:,:)				= Cw(ivar,:,:,:)*Delta_Pw(:,:,:)/100.
Pw_Delta_Cw(ivar,:,:,:)				= Pw(:,:,:)*Delta_Cw(ivar,:,:,:)/100.
Delta_Pw_Delta_Cw(ivar,:,:,:)		= Delta_Pw(:,:,:)*Delta_Cw(ivar,:,:,:)/100.
end do ; do ivar=

printVarSummary(Delta_Pw_Delta_Cw)
printMinMax(Delta_Pw_Delta_Cw,False)

; note: there is no meaning for CwPw_all when itag equals to zero.
CwPw_all			 		= new((/dimsizes(Vars_use),3,dimsizes(tag),dimsizes(case_nor),nbinx/),double)
CwPw_all!1					= "term"
CwPw_all&term				= (/0,1,2/)
CwPw_all(:,0,:,:,:)			= Cw_Delta_Pw
CwPw_all(:,1,:,:,:)			= Pw_Delta_Cw
CwPw_all(:,2,:,:,:)			= Delta_Pw_Delta_Cw
printVarSummary(CwPw_all)

; ============ averaged over all sgm
CwPw_all_name			= (/"dynamic","thermodynamic","co-variation"/)
CwPw_all_avg			= dim_sum_Wrap(CwPw_all)
;PwCw_avg				= dim_sum_Wrap(PwCw(:,1,:,:)-PwCw(:,0,:,:))
PwCw_avg				= dim_sum_Wrap(PwCw(:,1,:,:))-dim_sum_Wrap(PwCw(:,0,:,:))

Vars_out                    = Vars_use
Vars_out_fix                = Vars_out

do ivar = 0,dimsizes(Vars_out)-1
    lens                    	  = strlen(Vars_out(ivar))
    Vars_out_fix(ivar)      = str_concat((/conform_dims(15-lens, " ",-1)/))
    Vars_out_fix(ivar)      = str_concat((/Vars_out(ivar),Vars_out_fix(ivar)/))
    delete(lens)
end do ; do ivar=

print("Vars		cases	"+CwPw_all_name(0)+"		"+CwPw_all_name(1)+"		"+CwPw_all_name(2))
print(Vars_out_fix+"Park-RH 	"+sprintf("%9.5f",CwPw_all_avg(:,0,1,0))+"	"+sprintf("%9.5f",CwPw_all_avg(:,1,1,0))+"	"+sprintf("%9.5f",CwPw_all_avg(:,2,1,0))+"	"+sprintf("%9.5f",PwCw_avg(:,0)))
print(Vars_out_fix+"Gaus-PDF	"+sprintf("%9.5f",CwPw_all_avg(:,0,1,1))+"	"+sprintf("%9.5f",CwPw_all_avg(:,1,1,1))+"	"+sprintf("%9.5f",CwPw_all_avg(:,2,1,1))+"	"+sprintf("%9.5f",PwCw_avg(:,1)))

print("")

;====================== Plotting starts here .... ===============================================
plot_type			= "pdf"

;do itag=1,dimsizes(tag)-1
do itag=1,1 ; only p4K experiment

if(landmask)then
	wks  				= gsn_open_wks (plot_type,"./figure/omega-var-"+tag(itag)+"-CwPw_all-ocn")
else
	if(ocnmask)then
		wks  			= gsn_open_wks (plot_type,"./figure/omega-var-"+tag(itag)+"-CwPw_all-lnd")
	else
		wks  			= gsn_open_wks (plot_type,"./figure/omega-var-"+tag(itag)+"-CwPw_all-global")
	end if
end if


plot 				= new(dimsizes(Vars_use),graphic)
plotdiff 			= new(dimsizes(Vars_use),graphic)

plot_wgt			= new(dimsizes(Vars_use),graphic)
plot_wgt_diff		= new(dimsizes(Vars_use),graphic)

plot_CwPw			= new(dimsizes(Vars_use),graphic)
plot_PwCw			= new(dimsizes(Vars_use),graphic)
plot_DPwDCw			= new(dimsizes(Vars_use),graphic)
plot_CwPw_all		= new(dimsizes(Vars_use),graphic)

stds				= new(dimsizes(Vars_use),graphic)

nplots				= dimsizes(Vars_use)

legend_0			= new(nplots,graphic)
legend_1			= new(nplots,graphic)


res  							= True

res@gsnFrame 					= False
res@gsnDraw 					= False

res@vpWidthF 					= 0.7
res@vpHeightF 					= 0.4

res@trXMinF 					= nx_min
res@trXMaxF 					= nx_max

res@tiXAxisString 				= "vertical velocity (hPa/day)"

res@xyLineThicknessF 			= 4.0
res@xyDashPattern 				= 0

fontheight						= 0.03
res@gsnStringFontHeightF		= fontheight
res@tmXBLabelFontHeightF		= fontheight
res@tmYLLabelFontHeightF		= fontheight
res@tiXAxisFontHeightF			= fontheight
res@tiYAxisFontHeightF			= fontheight

; --------- OMEGA PDF only
res1								= res
res1@xyLineColors 					= colors

tmp					= reshape(omega_PDF_bin((/0,itag/),:,:),(/2*dimsizes(case_nor),nbinx/))
plot0 				= gsn_csm_xy(wks,omega_PDF_bin@bin_center,tmp,res1)
delete(tmp)


; -------OMEGA PDF difference only between sensitivity experiments and control exp. 
;  used to select variables we want to plot
vector								= (/3,4/)

res2								= res
res2@tiYAxisString					= ""
res2@xyLineColors					= colors;colors(vector)
res2@gsnYRefLine					= 0.

tmp					= reshape(omega_PDF_bin_diff((/itag/),:,:),(/dimsizes(case_nor),nbinx/))
plot1 				= gsn_csm_xy(wks,omega_PDF_bin_diff@bin_center,tmp,res2)
delete(tmp)

; ------- OMEGA-var plots start here ....
res3								= res
;res3@xyLineColors 					= colors

res4								= res
res4@xyLineColors					= colors ; colors((/3,4/))
res4@gsnYRefLine					= 0.


do ivar=0,dimsizes(Vars_use)-1

res3@trYMinF 						= valmin(ivar)
res3@trYMaxF 						= valmax(ivar)
res3@tiYAxisString 					= strings(ivar)+" "+units(ivar)
res4@tiYAxisString 					= strings(ivar)+" "+units(ivar)

tmp								= reshape(varsdata(ivar,(/0,itag/),:,:),(/2*dimsizes(case_nor),nbinx/))
tmp_std							= reshape(varsdata_std(ivar,(/0,itag/),:,:),(/2*dimsizes(case_nor),nbinx/))
	
; initialize one plot as the baseline for other plots. Others plots will overlay on this one.
res3@xyLineColors 				= colors(0)
plot_int						= gsn_csm_xy(wks,varsdata&binx_center,tmp(0,:),res3)

do ifig = 0,dimsizes(tag)*dimsizes(case_nor)-1

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
	
;	stds(ivar) 				= gsn_add_polygon (wks,plot_tmp,xp,yp,gsres)
	overlay(plot_int,plot_tmp)
	str1 					= unique_string("stds")
	plot_int@$str1$ 		= stds(ivar)
end do ; do ifig=
plot(ivar) 					= plot_int
delete([/tmp,tmp_std/])
delete(plot_int)


; for varsdata_diff, there is no "tag" dimension anymore. this is the revised version.
tmp									= reshape(varsdata_diff(ivar,:,:),(/dimsizes(case_nor),nbinx/))

plotdiff(ivar) 						= gsn_csm_xy (wks,varsdata_diff&binx_center,tmp, res4)
delete(tmp)

;plot_wgt(ivar)						= gsn_csm_xy (wks,varsdata_wgt&binx_center,varsdata_wgt(ivar,:,:), res3)
;plot_wgt_diff(ivar)				= gsn_csm_xy (wks,varsdata_wgt_diff&binx_center,varsdata_wgt_diff(ivar,:,:), res4)

; combine all components into one array "CwPw_all", so terms above aren't plotted anymore.
res5								= res
;res5@xyLineColors					= (/"red","red","blue","blue","cyan","cyan"/)
;res5@xyLineColors					= (/"red","red","red","red",\
;										"blue","blue","blue","blue",\
;										"cyan","cyan","cyan","cyan"/)
;res5@xyLineColors					= (/"red","red","red",\
;										"blue","blue","blue",\
;										"cyan","cyan","cyan"/)
ngroup = 2
res5@xyLineColors					= ndtooned(conform_dims((/3,ngroup/),(/"red","blue","cyan"/),0))
delete(res5@xyDashPattern)
res5@xyDashPatterns					= ndtooned(conform_dims((/3,ngroup/),(/14,0/),1))
res5@xyLineThicknesses				= ndtooned(conform_dims((/3,ngroup/),(/1,2/),1))
;res5@gsnYRefLine					= 0.
res5@tiYAxisString 					= strings(ivar)+" "+units(ivar)

;tmp									= reshape(CwPw_all(ivar,:,(/0,itag/),:,:),(/3*2*dimsizes(case_nor),nbinx/))
tmp									= reshape(CwPw_all(ivar,:,(/itag/),:,:),(/3*1*dimsizes(case_nor),nbinx/))

plot_CwPw_all(ivar)					= gsn_csm_xy (wks,CwPw_all&binx_center,tmp, res5)
;draw(plot_CwPw_all(ivar))
;frame(wks)

delete(tmp)

end do ; do ivar=

;gsn_panel(wks,plot,(/3,3/),False)

;------------------------------Legend preparing---------------------------------------
;Attach a legend
lgres 								= True
lgres@lgItemType 					= "Lines"
lgres@lgLabelFontHeightF 			= .08
lgres@vpWidthF 						= 0.2
lgres@vpHeightF 					= 0.17
lgres@lgPerimThicknessF 			= 1.0
;lgres@lgMonoDashIndex 				= True
;lgres@lgDashIndex 					= 0
lgres@lgMonoLineThickness			= True
lgres@lgMonoDashIndex				= True


lgres1								= lgres
lgres1@lgLineColors 				= colors
lgres1@lgLineThicknessF 			= res@xyLineThicknessF

;labels_1 							= (/"RH","PDF","PDF_2","PDF_MG10","RH_"+tag(itag),"PDF_"+tag(itag),"PDF_2_"+tag(itag),"PDF_MG10_"+tag(itag)/)
labels_1 							= (/"RH","PDF","RH_"+tag(itag),"PDF_"+tag(itag)/)

lgres2								= lgres
lgres2@lgLineColors					= colors ;colors(vector)
lgres2@lgLineThicknessF 			= res@xyLineThicknessF

;labels_2 							= labels_1((/2,3/)) ;labels_1(vector)
labels_2 							= labels_1((/2,3/)) ;labels_1(vector)

; for omega-PDF only
legend0					= gsn_create_legend(wks,dimsizes(labels_1),labels_1,lgres1)
legend1					= gsn_create_legend(wks,dimsizes(labels_2),labels_2,lgres2)

; for omega-var plots
do ivar=0,nplots-1
legend_0(ivar)			= gsn_create_legend(wks,dimsizes(labels_1),labels_1,lgres1)
legend_1(ivar)			= gsn_create_legend(wks,dimsizes(labels_2),labels_2,lgres2)
end do ; do ivar=

;---------------------------------------------
amres 						= True
amres@amJust 				= "TopLeft"
amres@amParallelPosF 		= -0.5
amres@amOrthogonalPosF 		= -0.5

; for omega-PDF only
annoid0 			= gsn_add_annotation(plot0,legend0,amres)
annoid1 			= gsn_add_annotation(plot1,legend1,amres)

; for omega-var plots
;tmp_plot			= (/plot(0),plotdiff(2),plot_wgt(3),plot_wgt_diff(0),plot_CwPw(3),plot_PwCw(0),plot_DPwDCw(0)/)
;legends			= (/legend_0(0),legend_1(0),legend_0(1),legend_1(1),legend_1(2),legend_1(3),legend_1(4)/)

; combine all components into one array "CwPw_all", so terms above aren't plotted anymore.
tmp_plot			= (/plot(0),plotdiff(0)/)
legends				= (/legend_0(0),legend_1(0)/)

annoid				= new(dimsizes(tmp_plot),graphic)
do ivar=0,dimsizes(tmp_plot)-1
annoid(ivar) 			= gsn_add_annotation(tmp_plot(ivar),legends(ivar),amres)
end do ; do ivar=


; for CwPw_all
lgres3								= lgres
lgres3@lgMonoLineThickness			= False
lgres3@lgMonoDashIndex				= False
lgres3@lgLineColors					= res5@xyLineColors
lgres3@lgDashIndexes				= res5@xyDashPatterns
;lgres3@lgLineThicknessF			= res5@xyLineThicknessF
lgres3@lgLineThicknesses			= res5@xyLineThicknesses

labels_3 			= (/"RH-Cw_DPw","PDF-Cw_DPw","RH-Pw_DCw","PDF-Pw_DCw","RH-DCw_DPw","PDF-DCw_DPw"/)

legend3				= gsn_create_legend(wks,dimsizes(labels_3),labels_3,lgres3)
annoid3				= gsn_add_annotation(plot_CwPw_all(0),legend3,amres)


;------------------------Panel plots----------------------------------------------
resP 				= True
;resP@gsnPanelRowSpec 		= True
resP@gsnPanelXWhiteSpacePercent	= 3.
resP@gsnPanelXWhiteSpacePercent	= 5.

; for omega-PDF only
gsn_panel(wks,(/plot0,plot1/),(/1,2/),False)

; for omega-var plots
idx					= ispan(0,dimsizes(use_idx)-1,1)
titles				= (/"varpdf","varpdf_diff","CwPw_all"/)
panel_plots		 	= (/plot(idx),plotdiff(idx),plot_CwPw_all(idx)/)

do i=0,dimsizes(titles)-1
	resP@gsnPanelMainString	= titles(i)
	gsn_panel(wks,panel_plots(i,:),(/2,3/),resP)
end do

end do ; do itag=

end

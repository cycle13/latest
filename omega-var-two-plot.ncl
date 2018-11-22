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

option					= (/"p4K","4xCO2"/)

var 					= (/"OMEGA","SWCF","LWCF","CLDLOW","CLDMED","CLDHGH","CLDTOT"/)

do_pdf_var 				= "F"
do_pdf_omega				= "F"


do ip=0,dimsizes(option)-1

;====== read data in
infile1 = addfile("./tmp-data/data-after-regrid-nor-2.5deg.nc","r")
infile2 = addfile("./tmp-data/data-after-regrid-"+option(ip)+"-2.5deg.nc","r")

OMEGA_nor = infile1->OMEGA_500_regrid
SWCF_nor = infile1->SWCF_regrid
LWCF_nor = infile1->LWCF_regrid
CLDLOW_nor = infile1->CLDLOW_regrid
CLDMED_nor = infile1->CLDMED_regrid
CLDHGH_nor = infile1->CLDHGH_regrid
CLDTOT_nor = infile1->CLDTOT_regrid

OMEGA_tmp = infile2->OMEGA_500_regrid
 SWCF_tmp = infile2->SWCF_regrid
 LWCF_tmp = infile2->LWCF_regrid
CLDLOW_tmp = infile2->CLDLOW_regrid
CLDMED_tmp = infile2->CLDMED_regrid
CLDHGH_tmp = infile2->CLDHGH_regrid
CLDTOT_tmp = infile2->CLDTOT_regrid


nbinx = 25
nbiny = 50
nx_min = -100.
nx_max = 100.

if(do_pdf_var.eq."T")then

result = pdf_omega_var(OMEGA_nor,SWCF_nor,OMEGA_tmp,SWCF_tmp,nbinx,nbiny,nx_min,nx_max)
OMEGA_SWCF_1D_bin = result
delete(result)
result = pdf_omega_var(OMEGA_nor,LWCF_nor,OMEGA_tmp,LWCF_tmp,nbinx,nbiny,nx_min,nx_max)
OMEGA_LWCF_1D_bin = result
delete(result)

result = pdf_omega_var(OMEGA_nor,CLDLOW_nor,OMEGA_tmp,CLDLOW_tmp,nbinx,nbiny,nx_min,nx_max)
OMEGA_CLDLOW_1D_bin = result
delete(result)

result = pdf_omega_var(OMEGA_nor,CLDMED_nor,OMEGA_tmp,CLDMED_tmp,nbinx,nbiny,nx_min,nx_max)
OMEGA_CLDMED_1D_bin = result
delete(result)

result = pdf_omega_var(OMEGA_nor,CLDHGH_nor,OMEGA_tmp,CLDHGH_tmp,nbinx,nbiny,nx_min,nx_max)
OMEGA_CLDHGH_1D_bin = result
delete(result)

result = pdf_omega_var(OMEGA_nor,CLDTOT_nor,OMEGA_tmp,CLDTOT_tmp,nbinx,nbiny,nx_min,nx_max)
OMEGA_CLDTOT_1D_bin = result
delete(result)


OMEGA_netCF_1D_bin = OMEGA_SWCF_1D_bin
OMEGA_netCF_1D_bin = (/OMEGA_SWCF_1D_bin + OMEGA_LWCF_1D_bin/)
printVarSummary(OMEGA_netCF_1D_bin)

system("/bin/rm -f ./tmp-data/pdf-omega-var-"+option(ip)+".nc")
outfile1 = addfile("./tmp-data/pdf-omega-var-"+option(ip)+".nc","c")

outfile1->OMEGA_SWCF_1D_bin = OMEGA_SWCF_1D_bin
outfile1->OMEGA_LWCF_1D_bin = OMEGA_LWCF_1D_bin
outfile1->OMEGA_netCF_1D_bin = OMEGA_netCF_1D_bin
outfile1->OMEGA_CLDLOW_1D_bin = OMEGA_CLDLOW_1D_bin
outfile1->OMEGA_CLDMED_1D_bin = OMEGA_CLDMED_1D_bin
outfile1->OMEGA_CLDHGH_1D_bin = OMEGA_CLDHGH_1D_bin
outfile1->OMEGA_CLDTOT_1D_bin = OMEGA_CLDTOT_1D_bin

delete(outfile1)

end if ; pdf_var

; omega pdf distribution

if(do_pdf_omega.eq."T")then
result = pdf_omega_only(OMEGA_nor,OMEGA_tmp,nbinx,nx_min,nx_max)
OMEGA_PDF_bin = result
delete(result)
printVarSummary(OMEGA_PDF_bin)

system("/bin/rm -f ./tmp-data/pdf-omega-only-"+option(ip)+".nc")
outfile1 = addfile("./tmp-data/pdf-omega-only-"+option(ip)+".nc","c")
outfile1->OMEGA_PDF_bin = OMEGA_PDF_bin

delete(outfile1)
end if

infile3 			= addfile("./tmp-data/pdf-omega-var-"+option(ip)+".nc","r")


nbins				= 25
ncases				= 5
vars				= ((/"OMEGA_SWCF_1D_bin","OMEGA_LWCF_1D_bin","OMEGA_netCF_1D_bin",\
				     "OMEGA_CLDLOW_1D_bin","OMEGA_CLDMED_1D_bin","OMEGA_CLDHGH_1D_bin",\
				     "OMEGA_CLDTOT_1D_bin"/))

colors 				= (/"grey","blue","red","deepskyblue","purple"/)

varsdata			= new((/dimsizes(vars),ncases,nbins/),double)
varsdata!0			= "vars"
varsdata&vars			= fspan(0,dimsizes(vars)-1,dimsizes(vars))

do ivar=0,dimsizes(vars)-1
varsdata(ivar,:,:)		= infile3->$vars(ivar)$
end do ; do ivar=

printVarSummary(varsdata)

varsdata_diff			= new((/dimsizes(vars),2,nbins/),double)
varsdata_diff(:,0,:)		= varsdata(:,0,:)
varsdata_diff(:,0,:)		= (/varsdata(:,3,:)-varsdata(:,1,:)/)
varsdata_diff(:,1,:)		= (/varsdata(:,4,:)-varsdata(:,2,:)/)

printVarSummary(varsdata_diff)
printMinMax(varsdata_diff,True)


infile4 			= addfile("./tmp-data/pdf-omega-only-"+option(ip)+".nc","r")
OMEGA_PDF_bin 			= infile4->OMEGA_PDF_bin

; get the difference to make it clear
OMEGA_PDF_bin_diff		= OMEGA_PDF_bin
OMEGA_PDF_bin_diff(0,:)		= (OMEGA_PDF_bin(0,:)-OMEGA_PDF_bin(0,:));/OMEGA_PDF_bin(0,:)
OMEGA_PDF_bin_diff(1,:)		= (OMEGA_PDF_bin(1,:)-OMEGA_PDF_bin(0,:));/OMEGA_PDF_bin(0,:)
OMEGA_PDF_bin_diff(2,:)		= (OMEGA_PDF_bin(2,:)-OMEGA_PDF_bin(0,:));/OMEGA_PDF_bin(0,:)
OMEGA_PDF_bin_diff(3,:)		= (OMEGA_PDF_bin(3,:)-OMEGA_PDF_bin(1,:));/OMEGA_PDF_bin(1,:)
OMEGA_PDF_bin_diff(4,:)		= (OMEGA_PDF_bin(4,:)-OMEGA_PDF_bin(2,:));/OMEGA_PDF_bin(2,:)
printVarSummary(OMEGA_PDF_bin_diff)
printMinMax(OMEGA_PDF_bin_diff,True)

; 2018-04-18 18:28:31
; get the weighting varsdata by multiplying the OMEGA_PDF_bin itself
varsdata_wgt			= varsdata
do ivar=0,dimsizes(vars)-1
varsdata_wgt(ivar,:,:)			= varsdata(ivar,:,:)*OMEGA_PDF_bin/100.
end do
printVarSummary(varsdata_wgt)
printMinMax(varsdata_wgt,False)

varsdata_wgt_diff			= new((/dimsizes(vars),2,nbins/),double)
varsdata_wgt_diff(:,0,:)		= varsdata_wgt(:,0,:)
varsdata_wgt_diff(:,0,:)		= (/varsdata_wgt(:,3,:)-varsdata_wgt(:,1,:)/)
varsdata_wgt_diff(:,1,:)		= (/varsdata_wgt(:,4,:)-varsdata_wgt(:,2,:)/)

printVarSummary(varsdata_wgt_diff)
printMinMax(varsdata_wgt_diff,True)

; 2018-04-19 11:01:20
; get Cw*\Delta_Pw and Pw*\Delta_Cw
; Cw
Cw		 			= new((/dimsizes(vars),2,nbins/),double)
Cw(:,0,:)				= varsdata(:,1,:) ; def-cntl
Cw(:,1,:)				= varsdata(:,2,:) ; new-cntl
; Pw
Pw					= new((/2,nbins/),double)
Pw(0,:)					= OMEGA_PDF_bin(1,:)
Pw(1,:)					= OMEGA_PDF_bin(2,:)
; \Delta_Cw
Delta_Cw				= Cw
Delta_Cw(:,0,:)				= (/varsdata(:,3,:)-varsdata(:,1,:)/)
Delta_Cw(:,1,:)				= (/varsdata(:,4,:)-varsdata(:,2,:)/)
; \Delta_Pw
Delta_Pw				= Pw
Delta_Pw(0,:)				= OMEGA_PDF_bin(3,:)-OMEGA_PDF_bin(1,:)
Delta_Pw(1,:)				= OMEGA_PDF_bin(4,:)-OMEGA_PDF_bin(2,:)
; PwCw, Cw*\Delta_Pw and Pw*\Delta_Cw
PwCw					= Cw
Cw_Delta_Pw				= Cw
Pw_Delta_Cw				= Cw
Delta_Pw_Delta_Cw			= Cw
do ivar=0,dimsizes(vars)-1
PwCw(ivar,:,:)				= Pw(:,:)*Cw(ivar,:,:)/100.
Cw_Delta_Pw(ivar,:,:)			= Cw(ivar,:,:)*Delta_Pw/100.
Pw_Delta_Cw(ivar,:,:)			= Pw(:,:)*Delta_Cw(ivar,:,:)/100.
Delta_Pw_Delta_Cw(ivar,:,:)		= Delta_Pw(:,:)*Delta_Cw(ivar,:,:)/100.
end do ; do ivar=

printVarSummary(Delta_Pw_Delta_Cw)
printMinMax(Delta_Pw_Delta_Cw,False)

CwPw_all			 	= new((/dimsizes(vars),6,nbins/),double)
CwPw_all(:,0:1,:)			= Cw_Delta_Pw
CwPw_all(:,2:3,:)			= Pw_Delta_Cw
CwPw_all(:,4:5,:)			= Delta_Pw_Delta_Cw
printVarSummary(CwPw_all)


plot_type			= "pdf"
wks  				= gsn_open_wks (plot_type,"./figure/omega-var-"+option(ip)+"-CwPw_all")

plot 				= new(9,graphic)
plotdiff 			= new(9,graphic)

plot_wgt			= new(9,graphic)
plot_wgt_diff			= new(9,graphic)

plot_CwPw			= new(9,graphic)
plot_PwCw			= new(9,graphic)
plot_DPwDCw			= new(9,graphic)
plot_CwPw_all			= new(9,graphic)

nplots				= 7

legend_0			= new(nplots,graphic)
legend_1			= new(nplots,graphic)


res  				= True

res@gsnFrame 			= False
res@gsnDraw 			= False

res@vpWidthF 			= 0.7
res@vpHeightF 			= 0.4

res@trXMinF 			= -100.
res@trXMaxF 			= 100.

res@tiXAxisString 		= "vertical velocity (hPa/day)"

res@xyLineThicknessF 		= 4.0
res@xyDashPattern 		= 0

fontheight			= 0.03
res@gsnStringFontHeightF	= fontheight
res@tmXBLabelFontHeightF	= fontheight
res@tmYLLabelFontHeightF	= fontheight
res@tiXAxisFontHeightF		= fontheight
res@tiYAxisFontHeightF		= fontheight


res1				= res
res1@xyLineColors 		= colors

vector				= (/3,4/)

res2				= res
res2@tiYAxisString		= ""
res2@xyLineColors		= colors(vector)

plot0 				= gsn_csm_xy(wks,OMEGA_PDF_bin@bin_center,OMEGA_PDF_bin,res1)

plot1 				= gsn_csm_xy(wks,OMEGA_PDF_bin_diff@bin_center,OMEGA_PDF_bin_diff(vector,:),res2)

;-------------------------------

strings				= (/"SWCF",		"LWCF",		"netCF",	"CLDLOW",	"CLDMED",	"CLDHGH",	"CLDTOT"/)
valmin				= (/-140,		-20,		-60.,		0,		0,		0,		0/)
valmax				= (/20,			100,		0.,		90,		90,		90,		90/)
units				= (/"(W/m~S~2~N~)",	"(W/m~S~2~N~)",	"(W/m~S~2~N~)",	"(fraction)",	"(fraction)",	"(fraction)",	"(fraction)"/)


res3				= res
res3@xyLineColors 		= colors

res4				= res
res4@xyLineColors		= colors((/3,4/))
res4@gsnYRefLine		= 0.

res5				= res
res5@xyLineColors		= toint(fspan(2,254,6))
res5@xyDashPatterns		= (/14,0,14,0,14,0/)
res5@xyLineThicknessF		= 2.0
res5@gsnYRefLine		= 0.


do ivar=0,dimsizes(vars)-1

;res3@trYMinF 			= valmin(ivar)
;res3@trYMaxF 			= valmax(ivar)
res3@tiYAxisString 		= strings(ivar)+" "+units(ivar)
res4@tiYAxisString 		= strings(ivar)+" "+units(ivar)
res5@tiYAxisString 		= strings(ivar)+" "+units(ivar)


plot(ivar) 			= gsn_csm_xy (wks,varsdata&binx_center,varsdata(ivar,:,:), res3)
plotdiff(ivar) 			= gsn_csm_xy (wks,varsdata_diff&binx_center,varsdata_diff(ivar,:,:), res4)

;plot_wgt(ivar)			= gsn_csm_xy (wks,varsdata_wgt&binx_center,varsdata_wgt(ivar,:,:), res3)
;plot_wgt_diff(ivar)		= gsn_csm_xy (wks,varsdata_wgt_diff&binx_center,varsdata_wgt_diff(ivar,:,:), res4)
;
;plot_CwPw(ivar)		= gsn_csm_xy (wks,Cw_Delta_Pw&binx_center,Cw_Delta_Pw(ivar,:,:), res4)
;plot_PwCw(ivar)		= gsn_csm_xy (wks,Pw_Delta_Cw&binx_center,Pw_Delta_Cw(ivar,:,:), res4)
;
;plot_DPwDCw(ivar)		= gsn_csm_xy (wks,Delta_Pw_Delta_Cw&binx_center,Delta_Pw_Delta_Cw(ivar,:,:), res4)

; combine all components into one array "CwPw_all", so terms above aren't plotted anymore.
plot_CwPw_all(ivar)		= gsn_csm_xy (wks,CwPw_all&binx_center,CwPw_all(ivar,:,:), res5)

end do ; do ivar=

;------------------------------Legend preparing---------------------------------------
;Attach a legend
lgres 				= True
lgres@lgItemType 		= "Lines"
lgres@lgLabelFontHeightF 	= .08
lgres@vpWidthF 			= 0.2
lgres@vpHeightF 		= 0.17
lgres@lgPerimThicknessF 	= 1.0
lgres@lgMonoDashIndex 		= True
lgres@lgDashIndex 		= 0
lgres@lgLineThicknessF 		= res@xyLineThicknessF

lgres1				= lgres
lgres1@lgLineColors 		= colors
labels_1 			= (/"OBS","RH","PDF","RH_"+option(ip),"PDF_"+option(ip)/)

lgres2				= lgres
lgres2@lgLineColors		= colors(vector)
labels_2 			= labels_1(vector)

; for omega-PDF only
legend0				= gsn_create_legend(wks,dimsizes(labels_1),labels_1,lgres1)
legend1				= gsn_create_legend(wks,dimsizes(labels_2),labels_2,lgres2)

; for omega-var plots
do ivar=0,nplots-1
legend_0(ivar)			= gsn_create_legend(wks,dimsizes(labels_1),labels_1,lgres1)
legend_1(ivar)			= gsn_create_legend(wks,dimsizes(labels_2),labels_2,lgres2)
end do ; do ivar=

;---------------------------------------------
amres 				= True
amres@amJust 			= "TopLeft"
amres@amParallelPosF 		= -0.5
amres@amOrthogonalPosF 		= -0.5

; for omega-PDF only
annoid0 			= gsn_add_annotation(plot0,legend0,amres)
annoid1 			= gsn_add_annotation(plot1,legend1,amres)

; for omega-var plots
;tmp_plot			= (/plot(0),plotdiff(2),plot_wgt(3),plot_wgt_diff(0),plot_CwPw(3),plot_PwCw(0),plot_DPwDCw(0)/)
;legends			= (/legend_0(0),legend_1(0),legend_0(1),legend_1(1),legend_1(2),legend_1(3),legend_1(4)/)

; combine all components into one array "CwPw_all", so terms above aren't plotted anymore.
tmp_plot			= (/plot(0),plotdiff(2)/)
legends				= (/legend_0(0),legend_1(0)/)

annoid				= new(dimsizes(tmp_plot),graphic)
do ivar=0,dimsizes(tmp_plot)-1
annoid(ivar) 			= gsn_add_annotation(tmp_plot(ivar),legends(ivar),amres)
end do ; do ivar=

; for CwPw_all
lgres3				= lgres
lgres3@lgLineColors		= res5@xyLineColors
lgres3@lgDashIndexes		= res5@xyDashPatterns
lgres3@lgLineThicknessF		= res5@xyLineThicknessF

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
idx				= (/0,1,2,3,5,6/)
;titles				= (/"varpdf","varpdf_diff","wgt_varPDF","wgt_varPDF_diff","Cw*Delta_Pw","Pw*Delta_Cw","Delta_Pw*Delta_Cw"/)
;panel_plots		 	= (/plot(idx),plotdiff(idx),plot_wgt(idx),plot_wgt_diff(idx),plot_CwPw(idx),plot_PwCw(idx),plot_DPwDCw(idx)/)
;
titles				= (/"varpdf","varpdf_diff","CwPw_all"/)
panel_plots		 	= (/plot(idx),plotdiff(idx),plot_CwPw_all(idx)/)


do i=0,dimsizes(titles)-1
	resP@gsnPanelMainString	= titles(i)
	gsn_panel(wks,panel_plots(i,:),(/2,3/),resP)
end do

end do ; do ip=

end


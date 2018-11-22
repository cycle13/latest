
load "./functions_regrid.ncl"
load "./functions_getvar.ncl"

begin

;***************file 1********************

;fdir 				= "/home/lyl/WORK4/cesm1_2_1/archive/"
;fdir_append			= "/atm/hist/"
plot_type			= "x11"

fdir				= "/home/share3/lyl/work3/qinyi/mid-data/cfmip/"

;case_nor 				= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_MG15_amip-p4K"/)
;case_p4K		 		= (/"FAMIPC5_f09f09_mac2_amip_outSGM","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM"/)
case_nor 				= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_mac2_amip_outSGM"/)
case_p4K		 		= (/"FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM"/)

case					= (/case_nor,case_p4K/)
tag					= (/"nor", "p4K"/)

append 				= ".ANN.climo.nc"

T_levs	= (/180,185,190,200,210,220,230,240,250,260,270,280,290,295,300/) ; T
T_dlevs	= (/-9,-7,-5,-4,-3,-2,-1,0,1,2,3,4,5,7,9/)
T_dlevs_sp	= (/-0.1,-0.05,0.05,0.1/)

Q_levs	= (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
Q_dlevs	= (/-1.8,-1.5,-1.2,-.9,-.6,-.3,-.1,0,.1,.3,.6,.9,1.2,1.5,1.8/)
Q_dlevs_sp	= (/-0.1,-0.05,0.05,0.1/)

RELHUM_levs	= (/5,10,15,20,25,30,40,50,60,70,75,80,85,90,95/)
RELHUM_dlevs = (/-50,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,50/)
RELHUM_dlevs_sp	= (/-0.1,-0.05,0.05,0.1/)

CLOUD_levs	= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
CLOUD_dlevs	= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
CLOUD_dlevs_sp	= (/-0.1,-0.05,0.05,0.1/)

;				   0            1          2           3
varsA		= (/"T",			"Q",		"RELHUM", "CLOUD"/)
varsA_levs	= (/T_levs,      Q_levs, 	RELHUM_levs,	CLOUD_levs/)
varsA_dlevs	= (/T_dlevs,     Q_dlevs,	RELHUM_dlevs,	CLOUD_dlevs/)
varsA_dlevs_sp	= (/T_dlevs_sp, Q_dlevs_sp,    RELHUM_dlevs_sp,    CLOUD_dlevs_sp/)
varsA_units	= (/"K",   "g/kg",       "%",       "fraction"/)
;varsA_amres	= (/amres_1, amres_2, amres_3/)

use_idx		= (/3,2/)
vars		= varsA(use_idx)
vars_levs	= varsA_levs(use_idx,:)
vars_dlevs	= varsA_dlevs(use_idx,:)
vars_dlevs_sp	= varsA_dlevs_sp(use_idx,:)
vars_units	= varsA_units(use_idx)
;vars_amres	= varsA_amres(use_idx)

nlat				= 192
nlon				= 288
nlev				= 30

vars_all			= new((/dimsizes(vars),dimsizes(tag),dimsizes(case_nor),nlev,nlat,nlon/),float)
vars_all!0			= "var"
vars_all&var		= fspan(0,dimsizes(vars)-1,dimsizes(vars))
vars_all!1		 	= "tag"
vars_all&tag		= fspan(0,dimsizes(tag)-1,dimsizes(tag))
vars_all!2		 	= "case"
vars_all&case		= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))


do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1
	ftmp1								= addfile(fdir+case(itag,icase)+append,"r")
	do ivar=0,dimsizes(vars)-1
		vars_all(ivar,itag,icase,:,:,:)	 	= ftmp1->$vars(ivar)$(0,:,:,:)
		if(vars(ivar).eq."Q")then
			vars_all(ivar,itag,icase,:,:,:)	= vars_all(ivar,itag,icase,:,:,:)*1e3
		end if
	end do ; do icase=
end do ; do ivar=
end do ; do itag=
printVarSummary(vars_all)

; ===================== difference between p4K and cntl ==========================
vars_all_tagdiff		= new((/dimsizes(vars),dimsizes(case_nor),nlev,nlat,nlon/),float)
vars_all_tagdiff		= (/vars_all(:,1,:,:,:,:)-vars_all(:,0,:,:,:,:)/)
copy_VarCoords_1(vars_all(var|:,case|:,lev|:,lat|:,lon|:,tag|:),vars_all_tagdiff)
printVarSummary(vars_all_tagdiff)

; ==================== difference between PDF and RH ===============================
vars_all_casediff		= new((/dimsizes(vars),dimsizes(tag),nlev,nlat,nlon/),float)
vars_all_casediff		= (/vars_all(:,:,1,:,:,:)-vars_all(:,:,0,:,:,:)/)
copy_VarCoords_1(vars_all(var|:,tag|:,lev|:,lat|:,lon|:,case|:),vars_all_casediff)
printVarSummary(vars_all_casediff)

; ==================== difference between PDF_4K-PDF and RH_4K-RH ===============================
vars_all_diff2		= new((/dimsizes(vars),nlev,nlat,nlon/),float)
vars_all_diff2		= (/vars_all_tagdiff(:,1,:,:,:)-vars_all_tagdiff(:,0,:,:,:)/)
copy_VarCoords_1(vars_all_tagdiff(var|:,lev|:,lat|:,lon|:,case|:),vars_all_diff2)
printVarSummary(vars_all_diff2)

; ==================== area mean =================
;						Global, Tropical
latS				= (/-90.,-30/)
latE				= (/90.,30/)
lonS				= (/0.,0/)
lonE				= (/360.,360/)
lat					= vars_all_diff2&lat
clat				= tofloat(cos(0.01745329*lat))
clat!0				= "lat"
clat&lat			= vars_all_diff2&lat

vars_all_aave				= new((/dimsizes(latS),dimsizes(vars),dimsizes(tag),dimsizes(case_nor),nlev/),float)
vars_all_tagdiff_aave		= new((/dimsizes(latS),dimsizes(vars),dimsizes(case_nor),nlev/),float)
vars_all_casediff_aave		= new((/dimsizes(latS),dimsizes(vars),dimsizes(tag),nlev/),float)
vars_all_diff2_aave			= new((/dimsizes(latS),dimsizes(vars),nlev/),float)

do ilat=0,dimsizes(latS)-1
	vars_all_aave(ilat,:,:,:,:)				= wgt_areaave_Wrap(vars_all(:,:,:,:,{latS(ilat):latE(ilat)},{lonS(ilat):lonE(ilat)}),clat({latS(ilat):latE(ilat)}),1.0,1)
	vars_all_tagdiff_aave(ilat,:,:,:)		= wgt_areaave_Wrap(vars_all_tagdiff(:,:,:,{latS(ilat):latE(ilat)},{lonS(ilat):lonE(ilat)}),clat({latS(ilat):latE(ilat)}),1.0,1)
	vars_all_casediff_aave(ilat,:,:,:)		= wgt_areaave_Wrap(vars_all_casediff(:,:,:,{latS(ilat):latE(ilat)},{lonS(ilat):lonE(ilat)}),clat({latS(ilat):latE(ilat)}),1.0,1)
	vars_all_diff2_aave(ilat,:,:)			= wgt_areaave_Wrap(vars_all_diff2(:,:,{latS(ilat):latE(ilat)},{lonS(ilat):lonE(ilat)}),clat({latS(ilat):latE(ilat)}),1.0,1)

end do ; do ilat=

printVarSummary(vars_all_aave)
printVarSummary(vars_all_tagdiff_aave)
printVarSummary(vars_all_casediff_aave)
printVarSummary(vars_all_diff2_aave)

; figure out the maximum value for casediff 
print("the difference between PDF and RH in vertical cloud fraction profile")
print("level			global		tropical")
print(sprintf("%9.5f",vars_all_casediff_aave&lev)+"		"+sprintf("%9.5f",vars_all_casediff_aave(0,0,0,:))+"	"+sprintf("%9.5f",vars_all_casediff_aave(1,0,0,:)))
print(max(vars_all_casediff_aave(0,0,0,:)))

;************************create plot***************

do ivar=0,dimsizes(vars)-1

;wks 					= gsn_open_wks("x11","figure/lat-hgt-"+vars(ivar)+"-HammerProj")
wks 					= gsn_open_wks(plot_type,"figure/value-hgt-"+vars(ivar)+"-p4K")

plot 					= new(9,graphic)
plot1 					= new(9,graphic)

;gsn_define_colormap(wks,"BlueDarkRed18")

cmap1 					= read_colormap_file("WhiteBlueGreenYellowRed")
cmap2 					= read_colormap_file("BlueRed")
;cmap2					= read_colormap_file("BlueWhiteOrangeRed")

res                      	= True
res@gsnDraw 				= False
res@gsnFrame 				= False

res@trYReverse				= True
res@tiYAxisString			= "Pressure (hPa)"
res@tiXAxisString			= vars(ivar)

res@gsnXRefLine				= 0.0
res@gsnXRefLineColor		= "grey"
res@gsnXRefLineThicknessF	= 3.0
res@gsnXRefLineDashPattern	= 2


fontheight						= 0.035
res@gsnStringFontHeightF		= fontheight
res@tmXBLabelFontHeightF		= fontheight
res@tmYLLabelFontHeightF		= fontheight
res@tiYAxisFontHeightF			= fontheight
res@tiXAxisFontHeightF			= fontheight

res@tmXBLabelAngleF 			= 45
res@tmXBLabelDeltaF 			= 0.4

;res@tmXBTickSpacingF			= 60.

;res@gsnCenterString			= vars(ivar)
res@gsnRightString				= vars_units(ivar)

res@xyLineColors				= (/"black","brown"/)
res@xyLineThicknessF			= 3.0


;;;;;;;res2
res2 							= res

res2@xyLineColors				= (/"blue","cyan"/)
res2@xyLineThicknessF			= 3.0


var4K					= new((/9,dimsizes(latS),nlev/),float)
var4K(0,:,:)				= vars_all_aave(:,ivar,0,0,:)
var4K					= (/vars_all_aave(:,ivar,0,0,:),vars_all_aave(:,ivar,0,1,:),vars_all_casediff_aave(:,ivar,0,:),\ ; (var,tag,case,lat,lon)
					    vars_all_aave(:,ivar,1,0,:),vars_all_aave(:,ivar,1,1,:),vars_all_casediff_aave(:,ivar,1,:),\ 
					    vars_all_tagdiff_aave(:,ivar,0,:),vars_all_tagdiff_aave(:,ivar,1,:), vars_all_diff2_aave(:,ivar,:)/) 
printVarSummary(var4K)

strings4K				= (/"RH","PDF","PDF-RH",\
					    "RH_4K","PDF_4K","PDF_4K-RH_4K",\
					    "RH_4K-RH","PDF_4K-PDF","dlt_PDF-dlt_RH"/)


plot					= new(9,graphic)
plotov					= new(9,graphic)

do it=0,8
if(it.eq.0.or.it.eq.1.or.it.eq.3.or.it.eq.4)then
	res@gsnLeftString			= strings4K(it)
	plot(it)				= gsn_csm_xy(wks,var4K(it,:,:),var4K&lev,res)
else
	res2@gsnLeftString			= strings4K(it)
	plot(it)				= gsn_csm_xy(wks,var4K(it,:,:),var4K&lev,res2)
end if
end do ; do it=

; ==================== add legends ====================================
legends                         = (/"Global","Tropical"/)

lgres                           	  = True
;lgres@lgItemType                		= res@xyMarkLineMode
lgres@lgLabelFontHeightF        = 0.20
lgres@vpWidthF                  = 0.25
lgres@vpHeightF                 = 0.15
lgres@lgPerimOn                 = False

lgres1							= lgres
lgres1@lgLineThicknessF 		= res@xyLineThicknessF
lgres1@lgLineColors			 	= res@xyLineColors

lgres2							= lgres
lgres2@lgLineThicknessF 		= res2@xyLineThicknessF
lgres2@lgLineColors			 	= res2@xyLineColors

lbid1                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres1)
lbid2                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres2)

amres_1                           = True
amres_1@amJust					= "CenterLeft"
amres_1@amParallelPosF            = -0.50
amres_1@amOrthogonalPosF          = 0.0

amres_2                           = True
amres_2@amJust					= "TopLeft"
amres_2@amParallelPosF            = -0.5
amres_2@amOrthogonalPosF          = -0.5

amres_3                           = True
amres_3@amJust					= "CenterRight"
amres_3@amParallelPosF            = 0.5
amres_3@amOrthogonalPosF          = 0.

amres							= True
if(vars(ivar).eq."T")then
amres							= amres_1
end if

if(vars(ivar).eq."Q")then
amres							= amres_2
end if

if(vars(ivar).eq."RELHUM")then
amres							= amres_3
end if


annoid1                         = gsn_add_annotation(plot(0),lbid1,amres)
annoid2                         = gsn_add_annotation(plot(2),lbid2,amres)

; ==================== add legends ====================================

resP = True
resP@gsnPanelMainString = vars(ivar)

resP@gsnPanelXWhiteSpacePercent		= 2.5
resP@gsnPanelYWhiteSpacePercent		= 5.

resP@gsnPanelRight					= 0.95
resP@gsnPanelLeft					= 0.05

;if(vars(ivar).eq."SWCF".or.vars(ivar).eq."CLDLOW")then
gsn_panel(wks,plot,(/3,3/),resP)
;gsn_panel(wks,plot1,(/3,3/),resP)
;end if

end do ; do ivar=


end

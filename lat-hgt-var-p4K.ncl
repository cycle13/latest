
load "./functions_regrid.ncl"
load "./functions_getvar.ncl"

begin

;***************file 1********************

;fdir 				= "/home/lyl/WORK4/cesm1_2_1/archive/"
;fdir_append			= "/atm/hist/"
plot_type			= "pdf"

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

use_idx		= (/3,2/)
vars		= varsA(use_idx)
vars_levs	= varsA_levs(use_idx,:)
vars_dlevs	= varsA_dlevs(use_idx,:)
vars_dlevs_sp	= varsA_dlevs_sp(use_idx,:)
vars_units	= varsA_units(use_idx)

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

; ==================== zonal mean =================
latS				= -90.
latE				= 90.
lonS				= 0.
lonE				= 360.
lat					= vars_all_diff2&lat
clat				= tofloat(cos(0.01745329*lat))
clat!0				= "lat"
clat&lat			= vars_all_diff2&lat

vars_all_zm			= dim_avg_n_Wrap(vars_all,5)
printVarSummary(vars_all_zm)

vars_all_tagdiff_zm		= dim_avg_n_Wrap(vars_all_tagdiff,4)
vars_all_casediff_zm	= dim_avg_n_Wrap(vars_all_casediff,4)
vars_all_diff2_zm		= dim_avg_n_Wrap(vars_all_diff2,3)

printVarSummary(vars_all_tagdiff_zm)
printVarSummary(vars_all_casediff_zm)
printVarSummary(vars_all_diff2_zm)


;************************create plot***************

do ivar=0,dimsizes(vars)-1

;wks 					= gsn_open_wks("x11","figure/lat-hgt-"+vars(ivar)+"-HammerProj")
wks 					= gsn_open_wks(plot_type,"figure/lat-hgt-"+vars(ivar)+"-p4K")

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

res@cnFillOn            		= True
res@cnLinesOn           		= False
res@cnLineLabelsOn      		= False
res@cnFillPalette       		= cmap1(:240,:)


fontheight						= 0.035
res@gsnStringFontHeightF		= fontheight
res@tmXBLabelFontHeightF		= fontheight
res@tmYLLabelFontHeightF		= fontheight
res@tiYAxisFontHeightF			= fontheight

;res@tmXBTickSpacingF			= 60.

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
res2@cnLevels                 	= vars_dlevs(ivar,:)

var4K					= new((/9,nlev,nlat/),float)
var4K(0,:,:)				= vars_all_zm(ivar,0,0,:,:)
var4K					= (/vars_all_zm(ivar,0,0,:,:),vars_all_zm(ivar,0,1,:,:),vars_all_casediff_zm(ivar,0,:,:),\ ; (var,tag,case,lat,lon)
					    vars_all_zm(ivar,1,0,:,:),vars_all_zm(ivar,1,1,:,:),vars_all_casediff_zm(ivar,1,:,:),\ 
					    vars_all_tagdiff_zm(ivar,0,:,:),vars_all_tagdiff_zm(ivar,1,:,:), vars_all_diff2_zm(ivar,:,:)/) 
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
	plot(it)				= gsn_csm_contour(wks,smth9_Wrap(var4K(it,:,:),p,q,False),res)
else
	if(it.eq.2.or.it.eq.5.or.it.eq.8)then
		res2@cnLevels			= vars_dlevs(ivar,:)/5.
	else
		res2@cnLevels			= vars_dlevs(ivar,:)
	end if
	res2@gsnLeftString			= strings4K(it)
;	plot(it)				= gsn_csm_contour_map(wks,var4K(it,:,:),res2)
	plot(it)				= gsn_csm_contour(wks,smth9_Wrap(var4K(it,:,:),p,q,False),res2)
end if
end do ; do it=


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

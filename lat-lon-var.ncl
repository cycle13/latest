
load "./functions_regrid.ncl"
load "./functions_getvar.ncl"

begin

;***************file 1********************

;fdir 				= "/home/lyl/WORK4/cesm1_2_1/archive/"
;fdir_append			= "/atm/hist/"

fdir				= "/home/share3/lyl/work3/qinyi/mid-data/cfmip/"

;case 				= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_MG15_amip-4xCO2",\
;		 		    "FAMIPC5_f09f09_mac2_amip","FAMIPC5_f09f09_mac2_amip-p4K_debug","FAMIPC5_f09f09_mac2_amip-4xCO2"/)
case 				= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_MG15_amip-p4K",\
		 		    "FAMIPC5_f09f09_mac2_amip","FAMIPC5_f09f09_mac2_amip-p4K_debug"/)


;append 				= ".cam.h0.ANN.nc"
append				= "-1980.1984_regrid.nc"


SWCF_levs	= (/-170,-150,-135,-120,-105,-90,-75,-60,-45,-30,-15,0,15,30,45/)
;SWCF_dlevs	= (/-80,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,80/)
SWCF_dlevs	= (/-50,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,50/)

LWCF_levs	= (/-45,-30,-20,-10,0,10,20,30,40,50,60,70,85,100,115/) 
LWCF_dlevs	= (/-35,-30,-15,-10,-6,-4,-2, 0, 2,4,6,10,15,30,35/)

PRECT_levs	= (/.2,.5,1,2,3,4,5,6,7,8,9,10,12,14,17/)
;PRECT_dlevs	= (/-6,-5,-4,-3,-2,-1,-.5,0,.5,1,2,3,4,5,6/)
PRECT_dlevs	= (/-4,-3,-2.5,-2,-1.5,-1,-.5,0,.5,1,1.5,2,2.5,3,4/)

CLDLOW_levs	= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
CLDLOW_dlevs	= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)

CLDHGH_levs	= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
CLDHGH_dlevs	= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)


varsA		= (/"PRECT",	"PRECC",	"PRECL",	"SWCF",		"LWCF",		"CLDLOW",	"CLDHGH"/)
varsA_levs	= (/PRECT_levs, PRECT_levs,     PRECT_levs,     SWCF_levs,      LWCF_levs, 	CLDLOW_levs,	CLDHGH_levs/)
varsA_dlevs	= (/PRECT_dlevs,PRECT_dlevs,    PRECT_dlevs,    SWCF_dlevs,     LWCF_dlevs,	CLDLOW_dlevs,	CLDHGH_dlevs/)
varsA_units	= (/"mm/day",   "mm/day",       "mm/day",       "W/m~S~2~N~",   "W/m~S~2~N~",	"fraction",	"fraction"/)

use_idx		= (/3,4,5/)
vars		= varsA(use_idx)
vars_levs	= varsA_levs(use_idx,:)
vars_dlevs	= varsA_dlevs(use_idx,:)
vars_units	= varsA_units(use_idx)


f1				= addfile(fdir+case(0)+append,"r")
dSizes				= getfiledimsizes(f1)
print(dSizes)
nlat				= 180
nlon				= 360
var_def				= new((/3,dimsizes(vars),nlat,nlon/),float)
var_def!0		 	= "exps"
var_def&exps			= (/0,1,2/)
var_def!1			= "vars"
var_def&vars			= fspan(0,dimsizes(vars)-1,dimsizes(vars))

var_new				= var_def

do ivar=0,dimsizes(vars)-1
do icase=0,dimsizes(case)-1

ftmp1					= addfile(fdir+case(icase)+append,"r")
if(icase.le.2)then
	var_def(icase,ivar,:,:)	 	= ftmp1->$vars(ivar)$
else
	var_new(icase-3,ivar,:,:)	= ftmp1->$vars(ivar)$
end if
printVarSummary(var_def)
printMinMax(var_def,True)

end do
end do

var_def_diff				= new((/2,dimsizes(vars),nlat,nlon/),float)
var_def_diff!0				= "exps"
var_def_diff&exps			= (/0,1/)
var_def_diff!1				= "vars"
var_def_diff&vars			= fspan(0,dimsizes(vars)-1,dimsizes(vars))

var_def_diff(0,:,:,:)			= var_def(0,:,:,:)

var_new_diff				= var_def_diff

var_def_diff(0,:,:,:)			= var_def(1,:,:,:)-var_def(0,:,:,:)
var_def_diff(1,:,:,:)			= var_def(2,:,:,:)-var_def(0,:,:,:)

var_new_diff(0,:,:,:)			= var_new(1,:,:,:)-var_new(0,:,:,:)
var_new_diff(1,:,:,:)			= var_new(2,:,:,:)-var_new(0,:,:,:)

var_new_def_diff			= var_def

var_new_def_diff(0,:,:,:)		= var_new(0,:,:,:)-var_def(0,:,:,:)
var_new_def_diff(1,:,:,:)		= var_new(1,:,:,:)-var_def(1,:,:,:)
var_new_def_diff(2,:,:,:)		= var_new(2,:,:,:)-var_def(2,:,:,:)

printVarSummary(var_new_def_diff)
printMinMax(var_new_def_diff,True)

;************************create plot***************

do ivar=0,dimsizes(vars)-1

;wks 					= gsn_open_wks("x11","figure/lat-lon-"+vars(ivar)+"-HammerProj")
wks 					= gsn_open_wks("x11","figure/lat-lon-"+vars(ivar))

plot 					= new(9,graphic)
plot1 					= new(9,graphic)

;gsn_define_colormap(wks,"BlueDarkRed18")

cmap1 					= read_colormap_file("WhiteBlueGreenYellowRed")
cmap2 					= read_colormap_file("BlueRed")

res                      		= True
res@gsnDraw 				= False
res@gsnFrame 				= False
;res@mpProjection			= "Hammer"
;res@mpPerimOn              		= False	 

res@cnFillOn            		= True
res@cnLinesOn           		= False
res@cnLineLabelsOn      		= False
res@cnFillPalette       		= cmap1(:240,:)

res@mpLandFillColor			= "gray"

res@mpCenterLonF         		= 180.

fontheight				= 0.035
res@gsnStringFontHeightF		= fontheight
res@tmXBLabelFontHeightF		= fontheight
res@tmYLLabelFontHeightF		= fontheight

res@tmXBOn			 	= False
res@tmYLOn				= False
res@tmXTOn				= False
res@tmYROn				= False
res@tmXBTickSpacingF			= 60.
res@tmYLTickSpacingF			= 45.

res@cnLevelSelectionMode        	= "ExplicitLevels"
res@cnLevels                    	= vars_levs(ivar,:)

;res@gsnCenterString			= vars(ivar)
res@gsnRightString			= vars_units(ivar)

res@lbOrientation			= "vertical"
res@lbLabelFontHeightF			= fontheight

;;;;;;;res2
res2 					= res

delete(res2@cnFillPalette)
res2@cnFillPalette 			= cmap2(:,:)

res2@cnLevelSelectionMode        	= "ExplicitLevels"
res2@cnLevels                    	= vars_dlevs(ivar,:)

var4K					= new((/8,nlat,nlon/),float)
var4K(0,:,:)				= var_def(0,ivar,:,:)
var4K					= (/var_def(0,ivar,:,:),var_new(0,ivar,:,:),var_new_def_diff(0,ivar,:,:),\
					    var_def(1,ivar,:,:),var_new(1,ivar,:,:),var_new_def_diff(1,ivar,:,:),\
					    var_def_diff(0,ivar,:,:),var_new_diff(0,ivar,:,:)/)
printVarSummary(var4K)

var4xCO2				= new((/8,nlat,nlon/),float)
var4xCO2(0,:,:)				= var_def(0,ivar,:,:)
var4xCO2				= (/var_def(0,ivar,:,:),var_new(0,ivar,:,:),var_new_def_diff(0,ivar,:,:),\
					    var_def(2,ivar,:,:),var_new(2,ivar,:,:),var_new_def_diff(2,ivar,:,:),\
					    var_def_diff(1,ivar,:,:),var_new_diff(1,ivar,:,:)/)

strings4K				= (/"RH","PDF","PDF-RH",\
					    "RH_4K","PDF_4K","PDF_4K-RH_4K",\
					    "RH_4K-RH","PDF_4K-PDF"/)

strings4xCO2				= (/"RH","PDF","PDF-RH",\
					    "RH_4xCO2","PDF_4xCO2","PDF_4xCO2-RH_4xCO2",\
					    "RH_4xCO2-RH","PDF_4xCO2-PDF"/)


plot					= new(9,graphic)
plotov					= new(9,graphic)


do it=0,7
if(it.eq.0.or.it.eq.1.or.it.eq.3.or.it.eq.4)then
	res@gsnLeftString			= strings4K(it)
	plot(it)				= gsn_csm_contour_map(wks,var4K(it,:,:),res)
	res@gsnLeftString			= strings4xCO2(it)
	plot1(it)				= gsn_csm_contour_map(wks,var4xCO2(it,:,:),res)
else
	res2@gsnLeftString			= strings4K(it)
	plot(it)				= gsn_csm_contour_map(wks,var4K(it,:,:),res2)
	res2@gsnLeftString			= strings4xCO2(it)
	plot1(it)				= gsn_csm_contour_map(wks,var4xCO2(it,:,:),res2)
end if
end do ; do it=


resP = True
resP@gsnPanelMainString = vars(ivar)

resP@gsnPanelXWhiteSpacePercent		= 2.5
resP@gsnPanelYWhiteSpacePercent		= 5.

;if(vars(ivar).eq."SWCF".or.vars(ivar).eq."CLDLOW")then
gsn_panel(wks,plot,(/3,3/),resP)
gsn_panel(wks,plot1,(/3,3/),resP)
;end if

end do ; do ivar=


end


load "./functions_regrid.ncl"
load "./functions_getvar.ncl"

begin

;***************file 1********************

;fdir 				= "/home/lyl/WORK4/cesm1_2_1/archive/"
;fdir_append			= "/atm/hist/"

fdir				= "/home/share3/lyl/work3/qinyi/mid-data/cfmip/"

case 				= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_MG15_amip-4xCO2",\
		 		    "FAMIPC5_f09f09_mac2_amip","FAMIPC5_f09f09_mac2_amip-p4K_debug","FAMIPC5_f09f09_mac2_amip-4xCO2"/)

;append 				= ".cam.h0.ANN.nc"
append				= "-1980.1984_regrid.nc"


SWCF_levs	= (/-170,-150,-135,-120,-105,-90,-75,-60,-45,-30,-15,0,15,30,45/)
SWCF_dlevs	= (/-50,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,50/)

LWCF_levs	= (/-45,-30,-20,-10,0,10,20,30,40,50,60,70,85,100,115/) 
LWCF_dlevs	= (/-35,-30,-15,-10,-6,-4,-2, 0, 2,4,6,10,15,30,35/)

PRECT_levs	= (/.2,.5,1,2,3,4,5,6,7,8,9,10,12,14,17/)
PRECT_dlevs	= (/-4,-3,-2.5,-2,-1.5,-1,-.5,0,.5,1,1.5,2,2.5,3,4/)

CLDLOW_levs	= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
CLDLOW_dlevs	= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)

CLDHGH_levs	= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
CLDHGH_dlevs	= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)

OMEGA_levs	= (/-100,-80,-60,-40,-20,-10,-5,0,5,10,20,40,60,80,100/)
OMEGA_dlevs	= (/-21.,-18,-15,-12,-9,-6,-3,0,3,6,9,12,15,18,21/)

TGCLDLWP_levs	= (/10,25,50,75,100,150,175,200,250,300,350,400,450,500,550/)
TGCLDLWP_dlevs	= (/-120,-100,-80,-60,-40,-20,-10,0,10,20,40,60,80,100,120/)

AREL_levs	= (/6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,10.5,11.,11.5,12.,12.5,13./)
AREL_dlevs	= (/-7.,-6.,-5.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,5.,6.,7./)/10.

;			0      1         2        3       4        5        6        7         8     
varsAll		= (/"OMEGA", "PRECT", "PRECC", "PRECL",	"SWCF", "LWCF", "CLDLOW", "CLDHGH","TGCLDLWP"/)
varsAll_levs	= (/OMEGA_levs, PRECT_levs, PRECT_levs, PRECT_levs, SWCF_levs, LWCF_levs, CLDLOW_levs, CLDHGH_levs,TGCLDLWP_levs/)
varsAll_dlevs	= (/OMEGA_dlevs,PRECT_dlevs,PRECT_dlevs, PRECT_dlevs, SWCF_dlevs, LWCF_dlevs, CLDLOW_dlevs, CLDHGH_dlevs,TGCLDLWP_dlevs/)/2.
varsAll_units	= (/"hPa/day", "mm/day", "mm/day", "mm/day", "W/m~S~2~N~","W/m~S~2~N~", "fraction", "fraction","g/m~S~2~N~","micron"/)

vector		= (/8,6,4,5,7/)
vars		= varsAll((/vector/))
vars_levs	= varsAll_levs((/vector/),:)
vars_dlevs	= varsAll_dlevs((/vector/),:)
vars_units	= ((/vector/))


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

OMEGA_def			= new((/3,nlat,nlon/),float)
OMEGA_def!0			= "exps"
OMEGA_def&exps			= (/0,1,2/)

OMEGA_new			= OMEGA_def

do icase=0,dimsizes(case)-1
	ftmp1					= addfile(fdir+case(icase)+append,"r")
	; read OMEGA
	if(icase.le.2)then
		OMEGA_def(icase,:,:)	 	= ftmp1->OMEGA({500},:,:)
	else
		OMEGA_new(icase-3,:,:)		= ftmp1->OMEGA({500},:,:)
	end if

	; read other variables
	do ivar=0,dimsizes(vars)-1

	if(vars(ivar).eq."OMEGA")then
		if(icase.le.2)then
		var_def(icase,ivar,:,:)	 	= ftmp1->$vars(ivar)$({500},:,:)
		var_def(icase,ivar,:,:)		= var_def(icase,ivar,:,:)*864.
		else
		var_new(icase-3,ivar,:,:)	= ftmp1->$vars(ivar)$({500},:,:)
		var_new(icase-3,ivar,:,:)	= var_new(icase-3,ivar,:,:)*864.
		end if
	else
		if(icase.le.2)then
			var_def(icase,ivar,:,:)	 	= ftmp1->$vars(ivar)$
			if(vars(ivar).eq."TGCLDLWP")then
				var_def(icase,ivar,:,:)	= var_def(icase,ivar,:,:)*1e3
			end if
		else
			var_new(icase-3,ivar,:,:)	= ftmp1->$vars(ivar)$
			if(vars(ivar).eq."TGCLDLWP")then
				var_new(icase-3,ivar,:,:)= var_new(icase-3,ivar,:,:)*1e3
			end if
		end if

	end if

;	printVarSummary(var_def)
;	printMinMax(var_def,True)
	
	end do ; do ivar=
end do ; do icase=

OMEGA_def				= OMEGA_def*864
OMEGA_new				= OMEGA_new*864.
printVarSummary(OMEGA_def)
printMinMax(OMEGA_def,True)


var_def_diff				= new((/2,dimsizes(vars),nlat,nlon/),float)
var_def_diff!0				= "exps"
var_def_diff&exps			= (/0,1/)
var_def_diff!1				= "vars"
var_def_diff&vars			= fspan(0,dimsizes(vars)-1,dimsizes(vars))
var_def_diff(0,:,:,:)			= var_def(0,:,:,:)
var_new_diff				= var_def_diff

OMEGA_def_diff				= new((/2,nlat,nlon/),float)
OMEGA_def_diff!0			= "exps"
OMEGA_def_diff&exps			= (/0,1/)
OMEGA_new_diff				= OMEGA_def_diff


var_def_diff(0,:,:,:)			= var_def(1,:,:,:)-var_def(0,:,:,:)
var_def_diff(1,:,:,:)			= var_def(2,:,:,:)-var_def(0,:,:,:)
var_new_diff(0,:,:,:)			= var_new(1,:,:,:)-var_new(0,:,:,:)
var_new_diff(1,:,:,:)			= var_new(2,:,:,:)-var_new(0,:,:,:)

OMEGA_def_diff(0,:,:)			= OMEGA_def(1,:,:)-OMEGA_def(0,:,:)
OMEGA_def_diff(1,:,:)			= OMEGA_def(2,:,:)-OMEGA_def(0,:,:)
OMEGA_new_diff(0,:,:)			= OMEGA_new(1,:,:)-OMEGA_new(0,:,:)
OMEGA_new_diff(1,:,:)			= OMEGA_new(2,:,:)-OMEGA_new(0,:,:)


var_new_def_diff			= var_def

OMEGA_new_def_diff			= OMEGA_def

var_new_def_diff(0,:,:,:)		= var_new(0,:,:,:)-var_def(0,:,:,:)
var_new_def_diff(1,:,:,:)		= var_new(1,:,:,:)-var_def(1,:,:,:)
var_new_def_diff(2,:,:,:)		= var_new(2,:,:,:)-var_def(2,:,:,:)

OMEGA_new_def_diff(0,:,:)		= OMEGA_new(0,:,:)-OMEGA_def(0,:,:)
OMEGA_new_def_diff(1,:,:)		= OMEGA_new(1,:,:)-OMEGA_def(1,:,:)
OMEGA_new_def_diff(2,:,:)		= OMEGA_new(2,:,:)-OMEGA_def(2,:,:)


printVarSummary(var_new_def_diff)
printMinMax(var_new_def_diff,True)

printVarSummary(OMEGA_new_def_diff)
printMinMax(OMEGA_new_def_diff,True)

;************************create plot***************

plot_type				= "x11"
figure_type				= "HammerProj"
do ivar=0,dimsizes(vars)-1

if(figure_type.eq."HammerProj")then
wks 					= gsn_open_wks(plot_type,"figure/lat-lon-"+vars(ivar)+"-HammerProj-land-ocean")
else
wks 					= gsn_open_wks(plot_type,"figure/lat-lon-"+vars(ivar))
end if

plot 					= new(9,graphic)
plot1 					= new(9,graphic)

plotov 					= new(9,graphic)
plot1ov 				= new(9,graphic)


;gsn_define_colormap(wks,"BlueDarkRed18")

cmap1 					= read_colormap_file("WhiteBlueGreenYellowRed")
;cmap2 					= read_colormap_file("BlueRed")
cmap2				 	= read_colormap_file("BlueWhiteOrangeRed")

res                      		= True
;res@cnFillMode 				= "RasterFill" 
res@gsnDraw 				= False
res@gsnFrame 				= False

;if(figure_type.eq."HammerProj")then
;res@mpProjection			= "Hammer"
;else
;res@mpPerimOn              		= False	 
;end if

res@cnFillOn            		= True
res@cnLinesOn           		= False
res@cnLineLabelsOn      		= False
res@cnFillPalette       		= cmap1(:240,:)

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

;res@gsnCenterString			= vars(ivar)
;res@gsnRightString			= vars_units(ivar)
res@gsnRightString			= vars(ivar)

res@lbOrientation			= "vertical"
res@lbLabelFontHeightF			= fontheight

res@cnLevelSelectionMode        	= "ExplicitLevels"
res@cnLevels                    	= vars_levs(ivar,:)

;res@cnFillDrawOrder 			= "Predraw"


res1				 	= res
if(figure_type.eq."HammerProj")then
res1@mpProjection			= "Hammer"
res1@mpPerimOn              		= False	 
end if


; please start here..
; need to mask land lines for omega

; resov
resov				 	= res
resov@cnFillOn				= False
resov@cnLinesOn				= True
resov@cnLineLabelsOn			= False
resov@gsnContourZeroLineThicknessF	= 0.

resov@cnLineDrawOrder			= "Predraw"

resov@gsnLeftString			= ""
resov@gsnRightString			= ""
resov@gsnCenterString			= ""

resov@gsnContourNegLineDashPattern	= 2
resov@gsnContourPosLineDashPattern	= 0

;resov@cnLineThicknessF			= 2.0

;resov@cnSmoothingOn			= True
;resov@cnSmoothingTensionF		= -5
;resov@cnSmoothingDistanceF		= 0.03

delete(resov@mpCenterLonF)

resov@cnLevels				= OMEGA_levs

;;;;;;;res2
res2 					= res
if(figure_type.eq."HammerProj")then
res2@mpProjection			= "Hammer"
res2@mpPerimOn              		= False	 
end if


delete(res2@cnFillPalette)
res2@cnFillPalette 			= cmap2(:,:)
;res2@gsnSpreadColorStart         	= 40
;res2@gsnSpreadColorEnd           	= 217

res2@cnLevelSelectionMode        	= "ExplicitLevels"
res2@cnLevels                    	= vars_dlevs(ivar,:)

res2ov					= resov
res2ov@cnLevels				= OMEGA_dlevs


var4K					= new((/8,nlat,nlon/),float)
var4K(0,:,:)				= var_def(0,ivar,:,:)
var4K					= (/var_def(0,ivar,:,:),var_new(0,ivar,:,:),var_new_def_diff(0,ivar,:,:),\
					    var_def(1,ivar,:,:),var_new(1,ivar,:,:),var_new_def_diff(1,ivar,:,:),\
					    var_def_diff(0,ivar,:,:),var_new_diff(0,ivar,:,:)/)

var4xCO2				= new((/8,nlat,nlon/),float)
var4xCO2(0,:,:)				= var_def(0,ivar,:,:)
var4xCO2				= (/var_def(0,ivar,:,:),var_new(0,ivar,:,:),var_new_def_diff(0,ivar,:,:),\
					    var_def(2,ivar,:,:),var_new(2,ivar,:,:),var_new_def_diff(2,ivar,:,:),\
					    var_def_diff(1,ivar,:,:),var_new_diff(1,ivar,:,:)/)

OMEGA4K					= new((/8,nlat,nlon/),float)
OMEGA4K(0,:,:)				= OMEGA_def(0,:,:)
OMEGA4K					= (/OMEGA_def(0,:,:),     OMEGA_new(0,:,:),OMEGA_new_def_diff(0,:,:),\
					    OMEGA_def(1,:,:),     OMEGA_new(1,:,:),OMEGA_new_def_diff(1,:,:),\
					    OMEGA_def_diff(0,:,:),OMEGA_new_diff(0,:,:)/)

OMEGA4xCO2				= new((/8,nlat,nlon/),float)
OMEGA4xCO2(0,:,:)			= OMEGA_def(0,:,:)
OMEGA4xCO2				= (/OMEGA_def(0,:,:),     OMEGA_new(0,:,:),OMEGA_new_def_diff(0,:,:),\
					    OMEGA_def(2,:,:),     OMEGA_new(2,:,:),OMEGA_new_def_diff(2,:,:),\
					    OMEGA_def_diff(1,:,:),OMEGA_new_diff(1,:,:)/)


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
	res1@gsnLeftString			= strings4K(it)
	plot(it)				= gsn_csm_contour_map(wks,var4K(it,:,:),res1)
	plotov(it)				= gsn_csm_contour(wks,OMEGA4K(it,:,:),resov)
;	overlay(plot(it),plotov(it))
	res1@gsnLeftString			= strings4xCO2(it)
	plot1(it)				= gsn_csm_contour_map(wks,var4xCO2(it,:,:),res1)
	plot1ov(it)				= gsn_csm_contour(wks,OMEGA4xCO2(it,:,:),resov)
;	overlay(plot1(it),plot1ov(it))

else
	res2@gsnLeftString			= strings4K(it)
	plot(it)				= gsn_csm_contour_map(wks,var4K(it,:,:),res2)
	plotov(it)				= gsn_csm_contour(wks,OMEGA4K(it,:,:),res2ov)
;	overlay(plot(it),plotov(it))
	res2@gsnLeftString			= strings4xCO2(it)
	plot1(it)				= gsn_csm_contour_map(wks,var4xCO2(it,:,:),res2)
	plot1ov(it)				= gsn_csm_contour(wks,OMEGA4xCO2(it,:,:),res2ov)
;	overlay(plot1(it),plot1ov(it))
end if
end do ; do it=

;draw(plot(0))
;frame(wks)
;exit

resP = True
resP@gsnPanelMainString = vars(ivar)

resP@gsnPanelXWhiteSpacePercent		= 2.5
resP@gsnPanelYWhiteSpacePercent		= 5.

gsn_panel(wks,plot,(/3,3/),resP)
;gsn_panel(wks,plot1,(/3,3/),resP)

;gsn_panel(wks,(/plot(6),plot(7),plot1(6),plot1(7)/),(/2,2/),resP)

end do ; do ivar=


end

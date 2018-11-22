
load "./functions_regrid.ncl"
load "./functions_getvar.ncl"
load "./functions_calculate.ncl"

begin

;***************file 1********************

;fdir 				= "/home/lyl/WORK4/cesm1_2_1/archive/"
;fdir_append			= "/atm/hist/"

pattern_correlation	= False
regime_contribution = False
lat_lon_plot		= False

fdir				= "/home/share3/lyl/work3/qinyi/mid-data/cfmip/"

;case_nor 				= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_MG15_amip-p4K"/)
;case_p4K		 		= (/"FAMIPC5_f09f09_mac2_amip_outSGM","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM"/)
case_nor 				= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_mac2_amip_outSGM"/)
case_p4K		 		= (/"FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM"/)

case					= (/case_nor,case_p4K/)
tag					= (/"nor", "p4K"/)

append 				= ".ANN.climo.nc"

SWCF_levs	= (/-170,-150,-135,-120,-105,-90,-75,-60,-45,-30,-15,0,15,30,45/)
;SWCF_dlevs	= (/-80,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,80/)
SWCF_dlevs	= (/-50,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,50/)/2.
SWCF_dlevs_sp	= (/-10,-5,5,10/)

LWCF_levs	= (/-45,-30,-20,-10,0,10,20,30,40,50,60,70,85,100,115/) 
LWCF_dlevs	= (/-35,-30,-15,-10,-6,-4,-2, 0, 2,4,6,10,15,30,35/)/2.
LWCF_dlevs_sp	= (/-10,-5,5,10/)

PRECT_levs	= (/.2,.5,1,2,3,4,5,6,7,8,9,10,12,14,17/)
;PRECT_dlevs	= (/-6,-5,-4,-3,-2,-1,-.5,0,.5,1,2,3,4,5,6/)
PRECT_dlevs	= (/-4,-3,-2.5,-2,-1.5,-1,-.5,0,.5,1,1.5,2,2.5,3,4/)
PRECT_dlevs_sp	= (/-2,-1,1,2/)

CLDLOW_levs	= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
CLDLOW_dlevs	= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
CLDLOW_dlevs_sp	= (/-0.1,-0.05,0.05,0.1/)

CLDHGH_levs	= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
CLDHGH_dlevs	= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
CLDHGH_dlevs_sp	= (/-0.1,-0.05,0.05,0.1/)

Vars_mid	= (/"T", "Q", "OMEGA"/)

;				   0            1          2           3           4            5          6
varsA		= (/"PRECT",	"PRECC",	"PRECL",	"SWCF",		"LWCF",		"CLDLOW",	"CLDHGH"/)
varsA_levs	= (/PRECT_levs, PRECT_levs,     PRECT_levs,     SWCF_levs,      LWCF_levs, 	CLDLOW_levs,	CLDHGH_levs/)
varsA_dlevs	= (/PRECT_dlevs,PRECT_dlevs,    PRECT_dlevs,    SWCF_dlevs,     LWCF_dlevs,	CLDLOW_dlevs,	CLDHGH_dlevs/)
varsA_dlevs_sp	= (/PRECT_dlevs_sp,PRECT_dlevs_sp,    PRECT_dlevs_sp,    SWCF_dlevs_sp,     LWCF_dlevs_sp,	CLDLOW_dlevs_sp,	CLDHGH_dlevs_sp/)

varsA_units	= (/"mm/day",   "mm/day",       "mm/day",       "W/m~S~2~N~",   "W/m~S~2~N~",	"fraction",	"fraction"/)

use_idx		= (/3,4,5,6/)
vars		= varsA(use_idx)
vars_levs	= varsA_levs(use_idx,:)
vars_dlevs	= varsA_dlevs(use_idx,:)
vars_dlevs_sp	= varsA_dlevs_sp(use_idx,:)
vars_units	= varsA_units(use_idx)

nlat				= 192
nlon				= 288
nlev				= 30

vars_mid			= new((/dimsizes(Vars_mid),dimsizes(tag),dimsizes(case_nor),nlev,nlat,nlon/),float)
vars_mid!0			= "var"
vars_mid&var		= fspan(0,dimsizes(Vars_mid)-1,dimsizes(Vars_mid))
vars_mid!1		 	= "tag"
vars_mid&tag		= fspan(0,dimsizes(tag)-1,dimsizes(tag))
vars_mid!2		 	= "case"
vars_mid&case		= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))


vars_all			= new((/dimsizes(vars),dimsizes(tag),dimsizes(case_nor),nlat,nlon/),float)
vars_all!0			= "var"
vars_all&var		= fspan(0,dimsizes(vars)-1,dimsizes(vars))
vars_all!1		 	= "tag"
vars_all&tag		= fspan(0,dimsizes(tag)-1,dimsizes(tag))
vars_all!2		 	= "case"
vars_all&case		= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))


do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1

	ftmp1								= addfile(fdir+case(itag,icase)+append,"r")
	lev									= ftmp1->lev
	OCNFRAC								= ftmp1->OCNFRAC(0,:,:)
	printVarSummary(OCNFRAC)
	do ivar=0,dimsizes(vars)-1
		vars_all(ivar,itag,icase,:,:)	 	= ftmp1->$vars(ivar)$(0,:,:)
	end do ; do icase=
	do ivar=0,dimsizes(Vars_mid)-1
		vars_mid(ivar,itag,icase,:,:,:)		= ftmp1->$Vars_mid(ivar)$(0,:,:,:)
	end do ; do ivar=

end do ; do icase=
end do ; do itag=
printVarSummary(vars_all)
printVarSummary(vars_mid)

; ===================== difference between p4K and cntl ==========================
vars_all_tagdiff		= new((/dimsizes(vars),dimsizes(case_nor),nlat,nlon/),float)
vars_all_tagdiff		= (/vars_all(:,1,:,:,:)-vars_all(:,0,:,:,:)/)
copy_VarCoords_1(vars_all(var|:,case|:,lat|:,lon|:,tag|:),vars_all_tagdiff)
printVarSummary(vars_all_tagdiff)

; ==================== difference between PDF and RH ===============================
vars_all_casediff		= new((/dimsizes(vars),dimsizes(tag),nlat,nlon/),float)
vars_all_casediff(:,:,:,:)		= (/vars_all(:,:,1,:,:)-vars_all(:,:,0,:,:)/)
copy_VarCoords_1(vars_all(var|:,tag|:,lat|:,lon|:,case|:),vars_all_casediff)
printVarSummary(vars_all_casediff)

; ==================== difference between PDF_4K-PDF and RH_4K-RH ===============================
vars_all_diff2		= new((/dimsizes(vars),nlat,nlon/),float)
vars_all_diff2		= (/vars_all_tagdiff(:,1,:,:)-vars_all_tagdiff(:,0,:,:)/)
copy_VarCoords_1(vars_all_tagdiff(var|:,lat|:,lon|:,case|:),vars_all_diff2)
printVarSummary(vars_all_diff2)

if(pattern_correlation)then
; ==================== pattern correlation between (dlt_PDF-dlt_RH) and (PDF-RH) =================
latS				= -90.
latE				= 90.
lonS				= 0.
lonE				= 360.
lat					= vars_all_diff2&lat
clat				= tofloat(cos(0.01745329*lat))
clat!0				= "lat"
clat&lat			= vars_all_diff2&lat

ra					= pattern_cor(vars_all_diff2(:,{latS:latE},{lonS:lonE}),vars_all_casediff(:,0,{latS:latE},{lonS:lonE}),clat({latS:latE}),1) ; uncentered 
printVarSummary(ra)
print(ra)

end if ; pattern_correlation

if (regime_contribution)then
; ================= broken down by regions and regimes =====================================
; outside of the tropics --- broken down into: 90-60S, 60-30S, 30-60N, 60-90N
; tropics --- broken down by 500hPa vertical velocity and LTS
; 1. deep convection and land (land .or. omega < 0)
; 2. trade cumulus (ocean .and. omega > 0 .and. LTS < 17K)
; 3. stratocumulus (ocean .and. omega > 0 .and. LTS > 17K)

; 1.1 calculate LTS
LTS					= new((/dimsizes(tag),dimsizes(case_nor),nlat,nlon/),float) ; unit is degreeC
LTS!0		 		= "tag"
LTS&tag				= fspan(0,dimsizes(tag)-1,dimsizes(tag))
LTS!1		 		= "case"
LTS&case			= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))

do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1
result					= cal_eis(vars_mid(0,itag,icase,:,:,:),vars_mid(1,itag,icase,:,:,:),lev)
LTS(itag,icase,:,:)		= tofloat(result[0])
end do ; do icase=
end do ; do itag=
copy_VarCoords(vars_all(0,:,:,:,:),LTS)
printVarSummary(LTS)
printMinMax(LTS,True)

; 1.2 regime separation
regimes_name		= (/"SH-high", "SH-mid", "NH-mid", "NH-high", "Trop_deep_lnd", "Trop_Shcu", "Trop_StCu"/)
nregimes			= dimsizes(regimes_name)
regimes				= new((/nregimes,dimsizes(vars),dimsizes(tag),dimsizes(case_nor)/),float) 
data2regimes		= new((/nregimes,dimsizes(vars),dimsizes(tag),dimsizes(case_nor),nlat,nlon/),float) 
data2regimes(0,:,:,:,:,:)		= vars_all

lon					= vars_all&lon
lat                 = vars_all&lat
clat                = tofloat(cos(0.01745329*lat))
clat!0              = "lat"
clat&lat            = vars_all&lat

latSs				= (/-90,-60,30,60/)
latEs				= (/-60,-30,60,90/)
lonSs				= (/0,0,0,0/)
lonEs				= (/360,360,360,360/)

do ivar=0,dimsizes(vars)-1
do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1

	omega500				= vars_mid(2,itag,icase,{500},:,:)
	lat2d					= conform_dims(dimsizes(omega500),lat,0)
	lon2d					= conform_dims(dimsizes(omega500),lon,1)

	; mid- and high- latitudes: SH-high, SH-mid, NH-mid, NH-high
	do ilat=0,dimsizes(latSs)-1
		tmp1			= where(lat2d.gt.latSs(ilat) .and. lat2d.lt.latEs(ilat) .and. lon2d.gt.lonSs(ilat) .and. lon2d.lt.lonEs(ilat),vars_all(ivar,itag,icase,:,:),vars_all@_FillValue)
		data2regimes(ilat,ivar,itag,icase,:,:)		= (/tmp1/)
;		regimes(ilat,:,:,:)		= wgt_areaave_Wrap(vars_all(:,:,:,{latSs(ilat):latEs(ilat)},{lonSs(ilat):lonEs(ilat)}),clat({latSs(ilat):latEs(ilat)}),1,0)
		regimes(ilat,ivar,itag,icase)		= wgt_areaave_Wrap(tmp1,clat,1,0)
		delete(tmp1)
	end do ; do ilat=
	
	; tropical deep convection and land regions
	crit1					= (OCNFRAC.lt.0.2).and.(lat2d.gt.-30.and.lat2d.lt.30.)
	crit2					= (omega500.lt.0).and.(lat2d.gt.-30.and.lat2d.lt.30.)
	tmp1					= where(crit1 .or. crit2, vars_all(ivar,itag,icase,:,:),vars_all@_FillValue)
	data2regimes(4,ivar,itag,icase,:,:)		= (/tmp1/)
	regimes(4,ivar,itag,icase)		= wgt_areaave_Wrap(tmp1,clat,1,0)
	delete([/crit1, crit2, tmp1/])
	
	; tropical trade cumulus
	crit1					= (OCNFRAC.gt.0.8) .and. (omega500.gt.0).and.(lat2d.gt.-30.and.lat2d.lt.30.).and.(LTS(itag,icase,:,:).lt.17)
	tmp1					= where(crit1, vars_all(ivar,itag,icase,:,:),vars_all@_FillValue)
	data2regimes(5,ivar,itag,icase,:,:)		= (/tmp1/)
	regimes(5,ivar,itag,icase)		= wgt_areaave_Wrap(tmp1,clat,1,0)
	delete([/crit1, tmp1/])
	
	; tropical stratocumulus
	crit1					= (OCNFRAC.gt.0.8) .and. (omega500.gt.0).and.(lat2d.gt.-30.and.lat2d.lt.30.).and.(LTS(itag,icase,:,:).gt.17)
	tmp1					= where(crit1, vars_all(ivar,itag,icase,:,:),vars_all@_FillValue)
	data2regimes(6,ivar,itag,icase,:,:)		= (/tmp1/)
	regimes(6,ivar,itag,icase)		= wgt_areaave_Wrap(tmp1,clat,1,0)
	delete([/crit1, tmp1/])
end do ; do itag=
end do ; do icase=
end do ; do ivar=

printVarSummary(data2regimes)

;;; ======== the following part is used to verify whether the regime-sorting is right.
; ---- yes, plots indicate that all regimes are right.
;;wks			= gsn_open_wks("x11","test")
;;
;;res			= True
;;res@cnFillOn		= True
;;res@cnLinesOn		= False
;;res@cnLineLabelsOn	= False
;;
;;res@mpCenterLonF	= 180
;;
;;do ir=0,dimsizes(regimes_name)-1
;;do ivar=0,0;dimsizes(vars)-1
;;do itag=0,dimsizes(tag)-1
;;do icase=0,0;dimsizes(case_nor)-1
;;plot		= gsn_csm_contour_map_ce(wks,data2regimes(ir,ivar,itag,icase,:,:),res)
;;end do 
;;end do
;;end do
;;end do
;;; ======================================================

printVarSummary(regimes)
;print(regimes_name+"		"+regimes(:,0,0,0))

; -================ area-weighted regime's contribution
re   = 6.37122e06
rad  = 4.0 * atan(1.0) / 180.
con  = re * rad                 
clat = tofloat(cos(lat * rad))           ; cosine of latitude

dlon = (lon(2) - lon(1))        ; assume dlon is constant
dlat = (lat(2) - lat(1))        ; assume dlat is constant

dx   = con * dlon * clat        ; dx at each latitude
dy   = con * dlat               ; dy is constant
dydx = dy * dx                  ; dydx(nlat)

wgt  = new((/dimsizes(lat),dimsizes(lon)/), float)
wgt  = tofloat(conform(wgt, dydx, 0))
opt	= 0

; change the non-missing region value into 1.0
data2regimes_1				= where(.not.ismissing(data2regimes),1.0,vars_all@_FillValue)
; change all grids into 1. to get the total earth area.
data2regimes_full			= data2regimes_1
data2regimes_full			= where(ismissing(data2regimes_full),1.0,1.0)
printVarSummary(data2regimes_full)
; get the area weigting for each regimes
data2regimes_1_percent = wgt_areasum2(data2regimes_1, wgt, opt)/wgt_areasum2(data2regimes_full,wgt,opt)
printVarSummary(data2regimes_1_percent)

; ================= difference between P4K and CNTL

regimes_tagdiff				= (/regimes(:,:,1,:)-regimes(:,:,0,:)/) ; (regime, var, case_nor)

regimes_areawgt			 	= (/data2regimes_1_percent*regimes/)
regimes_tagdiff_areawgt		= (/regimes_areawgt(:,:,1,:)-regimes_areawgt(:,:,0,:)/) ; (regime, var, case_nor)

Vars_out                    = regimes_name
Vars_out_fix                = Vars_out
do ivar = 0,dimsizes(Vars_out)-1
    lens                    = strlen(Vars_out(ivar))
    Vars_out_fix(ivar)      = str_concat((/conform_dims(20-lens, " ",-1)/))
    Vars_out_fix(ivar)      = str_concat((/Vars_out(ivar),Vars_out_fix(ivar)/))
    delete(lens)
end do ; do ivar=

print("For SWCF")
print("regimes		tagdiff_RH	tagdiff_PDF tagdiff_areawgt_RH	tagdiff_areawgt_PDF")
print(Vars_out_fix+" "+sprintf("%9.5f",regimes_tagdiff(:,0,0))+"	"+sprintf("%9.5f",regimes_tagdiff(:,0,1))+"	"+sprintf("%9.5f",regimes_tagdiff_areawgt(:,0,0))+"	"+sprintf("%9.5f",regimes_tagdiff_areawgt(:,0,1)))

; ---- check whether adding the weigting could equal to the global mean value got from direct calculation. -- the answer is yes.
print("verify the correctness of weighting:		"+dim_sum_n(regimes_tagdiff_areawgt(:,0,0),0)+"	"+dim_sum_n(regimes_tagdiff_areawgt(:,0,1),0))
; -----


wks					 	= gsn_open_wks("pdf", "figure/Global-regime-contribution-p4K")

sres 						= True
sres@gsnDraw               = False
sres@gsnFrame              = False

sres@vpWidthF              = 0.7
sres@vpHeightF             = 0.5

sres@trXMinF               = 0.
sres@trXMaxF               = nregimes+1
;sres@trYMinF               = -1.
;sres@trYMaxF               = 1.

barwidth					= 0.30
sres@gsnXYBarChart         = True
sres@gsnXYBarChartBarWidth = barwidth           ; change bar widths

sres@gsnYRefLine 			= 0

sres@tmXBMode              = "Explicit"         ; explicit labels
sres@tmXBValues            = ispan(1,nregimes,1)
sres@tmXBLabels            = regimes_name

sres@tmXBLabelAngleF		= 45.

sres@tmXBLabelFontHeightF  = 0.0205
sres@tmXTLabelFontHeightF  = 0.0205
sres@tmYLLabelFontHeightF  = 0.0225

sres@gsnRightString        = ""
;sres@tiYAxisString         = "W/m2"

sres2						= sres

plots 						= new((/dimsizes(vars),dimsizes(case_nor)/),graphic)
plots_wgt					= new((/dimsizes(vars),dimsizes(case_nor)/),graphic)

;colors 						= (/"lightpink",   "lightblue"/)
colors 						= (/"grey",   "black"/)

; The X values for each plot is different, ensuring that the bars don't overlap.
do ivar=0,dimsizes(vars)-1
	sres@tiYAxisString		= vars(ivar)+"("+vars_units(ivar)+")"
	sres@trYMinF			= min(regimes_tagdiff(:,ivar,:))*1.1
	sres@trYMaxF			= max(regimes_tagdiff(:,ivar,:))*1.1
	sres@gsnLeftString		= "un-weighted"

	sres2@tiYAxisString		= vars(ivar)+"("+vars_units(ivar)+")"
	sres2@trYMinF			= min(regimes_tagdiff_areawgt(:,ivar,:))*1.1
	sres2@trYMaxF			= max(regimes_tagdiff_areawgt(:,ivar,:))*1.1
	sres2@gsnLeftString		= "area-weighted"

	; For un-weighted
	sres@gsnXYBarChartColors = colors(0)
	sres@xyLineColors=(/colors(0),colors(0)/) 	; this is used to change the outline to be the color of the fill. I still dont know why it needs two colors.
	plots(ivar,0) = gsn_csm_xy(wks,fspan(1-barwidth/2.,nregimes-barwidth/2,nregimes),regimes_tagdiff(:,ivar,0),sres)

	sres@gsnLeftString		 = ""
	sres@gsnXYBarChartColors = colors(1)
	sres@xyLineColors=(/colors(1),colors(1)/)
	plots(ivar,1) = gsn_csm_xy(wks,fspan(1+barwidth/2,nregimes+barwidth/2,nregimes),regimes_tagdiff(:,ivar,1),sres)

	; For area-weighted
	sres2@gsnXYBarChartColors = colors(0)
	sres2@xyLineColors = (/colors(0),colors(0)/)
	plots_wgt(ivar,0) = gsn_csm_xy(wks,fspan(1-barwidth/2.,nregimes-barwidth/2,nregimes),regimes_tagdiff_areawgt(:,ivar,0),sres2)

	sres2@gsnLeftString			= ""
	sres2@gsnXYBarChartColors = colors(1)
	sres2@xyLineColors	= (/colors(1),colors(1)/)
	plots_wgt(ivar,1) = gsn_csm_xy(wks,fspan(1+barwidth/2,nregimes+barwidth/2,nregimes),regimes_tagdiff_areawgt(:,ivar,1),sres2)

; Overlay the last three plots on the first one.
       do j=1,dimsizes(case_nor)-1
         overlay(plots(ivar,0),plots(ivar,j))
		 overlay(plots_wgt(ivar,0),plots_wgt(ivar,j))
       end do
end do

pres					= True
pres@gsnMaximize    	= True   ; Maximize in frame
pres@gsnFrame       	= False  ; Don't advance frame.
pres@gsnPanelBottom 	= 0.13   ; Leave room for labelbar

pres@gsnPanelLeft		= 0.05
pres@gsnPanelRight		= 0.95
pres@gsnPanelXWhiteSpacePercent		= 5.0
pres@gsnPanelYWhiteSpacePercent		= 5.0

; Create some labelbars for a legend
lbres                    = True          ; labelbar only resources
lbres@vpWidthF           = 0.1           ; labelbar width
lbres@vpHeightF          = 0.1           ; labelbar height
lbres@lbBoxMajorExtentF  = 0.10          ; puts space between color boxes
lbres@lbMonoFillPattern  = True          ; Solid fill pattern
lbres@lbLabelFontHeightF = 0.010         ; font height. default is small
lbres@lbLabelJust        = "TopLeft"  ; left justify labels
lbres@lbPerimOn          = False

labels = (/"Park-RH", "Gauss-PDF"/)
xpos   = (/0.15,  0.25/)
ypos	 = 0.95

; Panel the plots.
;First plot
gsn_panel(wks,plots(:,0),(/2,2/),pres)
;drawNDCGrid(wks)
do i=0,dimsizes(case_nor)-1
  lbres@lbFillColors = colors(i)
  gsn_labelbar_ndc(wks,1,labels(i),xpos(i),ypos,lbres)
end do	
frame(wks)   ; Advance the frame.

; Second plot
gsn_panel(wks,plots_wgt(:,0),(/2,2/),pres)

do i=0,dimsizes(case_nor)-1
  lbres@lbFillColors = colors(i)
  gsn_labelbar_ndc(wks,1,labels(i),xpos(i),ypos,lbres)
end do	

frame(wks)   ; Advance the frame.
exit

end if ; regime_contribution

;************************create plot***************
if(lat_lon_plot)then

do ivar=0,dimsizes(vars)-1

;wks 					= gsn_open_wks("x11","figure/lat-lon-"+vars(ivar)+"-HammerProj")
wks 					= gsn_open_wks("x11","figure/lat-lon-"+vars(ivar)+"-p4K")

plot 					= new(9,graphic)
plot1 					= new(9,graphic)

;gsn_define_colormap(wks,"BlueDarkRed18")

cmap1 					= read_colormap_file("WhiteBlueGreenYellowRed")
cmap2 					= read_colormap_file("BlueRed")
;cmap2					= read_colormap_file("BlueWhiteOrangeRed")

res                      	= True
res@gsnDraw 				= False
res@gsnFrame 				= False
;res@mpProjection			= "Hammer"
;res@mpPerimOn              = False	 

res@cnFillOn            		= True
res@cnLinesOn           		= False
res@cnLineLabelsOn      		= False
res@cnFillPalette       		= cmap1(:240,:)

res@mpLandFillColor				= "gray"

res@mpCenterLonF         		= 180.

fontheight						= 0.035
res@gsnStringFontHeightF		= fontheight
res@tmXBLabelFontHeightF		= fontheight
res@tmYLLabelFontHeightF		= fontheight

res@tmXBOn			 			= False
res@tmYLOn						= False
res@tmXTOn						= False
res@tmYROn						= False
res@tmXBTickSpacingF			= 60.
res@tmYLTickSpacingF			= 45.

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
delete(resl@mpCenterLonF)
delete(resl@mpLandFillColor)
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

var4K					= new((/9,nlat,nlon/),float)
var4K(0,:,:)				= vars_all(ivar,0,0,:,:)
var4K					= (/vars_all(ivar,0,0,:,:),vars_all(ivar,0,1,:,:),vars_all_casediff(ivar,0,:,:),\ ; (var,tag,case,lat,lon)
					    vars_all(ivar,1,0,:,:),vars_all(ivar,1,1,:,:),vars_all_casediff(ivar,1,:,:),\ 
					    vars_all_tagdiff(ivar,0,:,:),vars_all_tagdiff(ivar,1,:,:), vars_all_diff2(ivar,:,:)/) 
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
;	plot(it)				= gsn_csm_contour_map(wks,var4K(it,:,:),res)
	plot(it)				= gsn_csm_contour_map(wks,smth9_Wrap(var4K(it,:,:),p,q,False),res)

;	plot_ov					= gsn_csm_contour(wks,var4K(it,:,:),resl)
;	overlay(plot(it),plot_ov)
else
	res2@gsnLeftString			= strings4K(it)
	if(it.eq.8)then
		res2@gsnCenterString	= sprintf("%5.2f",ra(ivar))
	else
		res2@gsnCenterString	= ""
	end if
;	plot(it)				= gsn_csm_contour_map(wks,var4K(it,:,:),res2)
	plot(it)				= gsn_csm_contour_map(wks,smth9_Wrap(var4K(it,:,:),p,q,False),res2)
;	plot_ov					= gsn_csm_contour(wks,var4K(it,:,:),resl)
;	plot_ov					= gsn_csm_contour(wks,smth9_Wrap(var4K(it,:,:),p,q,False),resl)
;	overlay(plot(it),plot_ov)
end if
end do ; do it=


resP = True
resP@gsnPanelMainString = vars(ivar)

resP@gsnPanelXWhiteSpacePercent		= 2.5
resP@gsnPanelYWhiteSpacePercent		= 5.

;if(vars(ivar).eq."SWCF".or.vars(ivar).eq."CLDLOW")then
gsn_panel(wks,plot,(/3,3/),resP)
;gsn_panel(wks,plot1,(/3,3/),resP)
;end if

end do ; do ivar=

end if ; lat_lon_plot
end

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

function zpsi(lev[*]:numeric,lat[*]:numeric,lon[*]:numeric,\
              U_1[*][*][*]:numeric,V_1[*][*][*]:numeric,OMEGA_1[*][*][*],\
	      gw[*]:numeric,latS:float,latE:float)

begin
;==================================================================
; set parameters...
;==================================================================
g			= 9.81		; gravity
a			= 6.37*10^6	; radius of earth
pi		 	= 3.1428	; pi

nlev			= dimsizes(lev)
nlat			= dimsizes(lat)
nlon			= dimsizes(lon)

;==================================================================
; get average over specific zonal ranges with raw wind: u and omega
;==================================================================

OMEGA_1_mm 		= dim_avg_wgt_n_Wrap(OMEGA_1(:,{latS:latE},:),gw({latS:latE}),0,1)
OMEGA_1_mm&lev@units	= "hPa"

U_1_mm 			= dim_avg_wgt_n_Wrap(U_1(:,{latS:latE},:),gw({latS:latE}),0,1)
U_1_mm&lev@units	= "hPa"

;==================================================================
; get velocity potential: u=d()/dx; v = d()/dy
; and divergent wind: u & v
;==================================================================
sfvp_1 			= new((/2,nlev,nlat,nlon/),float)

sfvp_1(0,:,:,:)		= U_1
sfvp_1 			= uv2sfvpG(U_1,V_1)

; sfvp(0,...) is stream function
; sfvp(1,...) is velocity potential
vp_1 			= sfvp_1(1,:,:,:)

DvpDx_1 		= new((/nlev,nlat,nlon/),float)
DvpDx_1			= vp_1
DvpDx_1 		= doubletofloat(center_finite_diff_n(vp_1,lon*111000,False,0,2))


DvpDx_1_mm 		= dim_avg_wgt_n_Wrap(DvpDx_1(:,{latS:latE},:),gw({latS:latE}),0,1)
DvpDx_1_mm&lev@units	= "hPa"
;printVarSummary(DvpDx_1_mm)

;==================================================================
; calculate zonal stream function: Walker Circulation
;==================================================================
dpsi 			= new(  (/dimsizes(lev),dimsizes(lon)/) , "float")
psi 			= new(  (/dimsizes(lev),dimsizes(lon)/) , "float")
dpsidphi		= new(  (/dimsizes(lev),dimsizes(lon)/) , "float")
omega			= new(  (/dimsizes(lev),dimsizes(lon)/) , "float")

; initialize these arrays at zeros
psi 			= 0.
omega 			= 0.

;since all of these values use the same level/latitude coordinates
;you can copy v to each other variable
copy_VarCoords(DvpDx_1_mm,dpsi)
copy_VarCoords(DvpDx_1_mm,psi)
copy_VarCoords(DvpDx_1_mm,dpsidphi)
copy_VarCoords(DvpDx_1_mm,omega)

;printVarSummary(dpsi)

; get dpsi
do ilev=0,dimsizes(lev)-1
	do ilon = 0,dimsizes(lon)-1
		dpsi(ilev,ilon)	= DvpDx_1_mm(ilev,ilon)*2*pi*a/g
	end do
end do

do ilev=1,dimsizes(lev)-1
	do ilon=0,dimsizes(lon)-1
		psi(ilev,ilon) = psi(ilev-1,ilon)+dpsi(ilev,ilon)*((lev(ilev)-lev(ilev-1))*100)

	end do
end do

;printVarSummary(dpsi)
;printMinMax(dpsi,True)
;printVarSummary(psi)
;printMinMax(psi,True)

; psi is calculated!
; the next step is to calculate omega which we will do a forward and backwards differencing
; for the first and last values, but centerred differencing for
; all the other calculations

; we do this forwards/bacwards/ and centered differencing on the d(psi)/d(phi) term

do ilev=0,dimsizes(lev)-1
	dpsidphi(ilev,0)		= (psi(ilev,1)-psi(ilev,0))/(lon(1)-lon(0))*(180/pi)
	dpsidphi(ilev,nlon-1)		= (psi(ilev,nlon-1)-psi(ilev,nlon-2))/(lon(1)-lon(0))*(180/pi)

	do ilon=1,nlon-2
		dpsidphi(ilev,ilon)	= (psi(ilev,ilon+1)-psi(ilev,ilon-1))/(lon(1)-lon(0))*(180/pi)
	end do
end do

;printVarSummary(dpsidphi)
;printMinMax(dpsidphi,True)

; Finally run one last set of do lops to calculate omega
do ilev=1,dimsizes(lev)-1
	do ilon=0,dimsizes(lon)-1
		omega(ilev,ilon)	= -g/(2*pi*a^2)*dpsidphi(ilev,ilon)
	end do
end do

;printVarSummary(omega)
;printMinMax(omega,True)

; This step is just to adjust the values so they can easily be plotted. Proper units are cited on the figure.
psi 					= psi*10^-10
;omega					= omega*10^3
omega					= omega*864

;printMinMax(psi,True)
;printMinMax(omega,True)

return([/psi,omega,OMEGA_1_mm,U_1_mm/])

end

;==============================================================
; Main programing starts here....
;==============================================================

begin

;***************file 1********************
fdir				= "/home/share3/lyl/work3/qinyi/mid-data/cfmip/"

case 				= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_MG15_amip-p4K",      "FAMIPC5_f09f09_MG15_amip-4xCO2",\
		 		    "FAMIPC5_f09f09_mac2_amip","FAMIPC5_f09f09_mac2_amip-p4K_debug","FAMIPC5_f09f09_mac2_amip-4xCO2"/)

;append				= "-1980.1984_regrid.nc"
append				= "-1980.1984.nc"

vars				= (/"U","V","OMEGA"/)

f1				= addfile(fdir+case(0)+append,"r")

nlat				= dimsizes(f1->lat)
nlon				= dimsizes(f1->lon)
nlev				= dimsizes(f1->lev)

gw 				= f1->gw
lev 				= tofloat(f1->lev)
lat 				= tofloat(f1->lat)
lon 				= tofloat(f1->lon)

var_def				= new((/3,dimsizes(vars),nlev,nlat,nlon/),float)
var_def!0		 	= "exps"
var_def&exps			= (/0,1,2/)
var_def!1			= "vars"
var_def&vars			= fspan(0,dimsizes(vars)-1,dimsizes(vars))

var_new				= var_def

do ivar=0,dimsizes(vars)-1
do icase=0,dimsizes(case)-1

ftmp1						= addfile(fdir+case(icase)+append,"r")
if(icase.le.2)then
	var_def(icase,ivar,:,:,:)	 	= ftmp1->$vars(ivar)$
else
	var_new(icase-3,ivar,:,:,:)		= ftmp1->$vars(ivar)$
end if
printVarSummary(var_def)

end do
end do

;1111111111111111111111111111111

latS 			= -30.
latE 			= 30.
lonS			= 140.
lonE			= 280.

g			= 9.81		; gravity
a			= 6.37*10^6	; radius of earth
pi		 	= 3.1428	; pi


;==================================================================
; get zonal stream function from subroutine
;==================================================================
psi_def			= new((/3,nlev,nlon/),float)
omega_def		= new((/3,nlev,nlon/),float)
OMEGA_def_mm		= new((/3,nlev,nlon/),float)
U_def_mm		= new((/3,nlev,nlon/),float)

psi_new			= psi_def
omega_new		= omega_def
OMEGA_new_mm		= OMEGA_def_mm
U_new_mm		= U_def_mm

do iexp=0,2

result = zpsi(lev, lat, lon, var_def(iexp,0,:,:,:), var_def(iexp,1,:,:,:), var_def(iexp,2,:,:,:), gw, latS, latE)

psi_def(iexp,:,:)			= result[0]
omega_def(iexp,:,:)			= result[1]
OMEGA_def_mm(iexp,:,:)			= result[2]
U_def_mm(iexp,:,:)			= result[3]

result = zpsi(lev, lat, lon, var_new(iexp,0,:,:,:), var_new(iexp,1,:,:,:), var_new(iexp,2,:,:,:), gw, latS, latE)
psi_new(iexp,:,:)			= result[0]
omega_new(iexp,:,:)			= result[1]
OMEGA_new_mm(iexp,:,:)			= result[2]
U_new_mm(iexp,:,:)			= result[3]

end do ; do iexp=

printVarSummary(psi_def)
printVarSummary(omega_new)

; Get difference between simulations

psi_def_diff				= new((/2,nlev,nlon/),float)
psi_def_diff!0				= "exps"
psi_def_diff&exps			= (/0,1/)
psi_def_diff(0,:,:)			= psi_def(0,:,:)

psi_new_diff				= psi_def_diff

; note: percent change could lead to very large percentage over the region where the absolute value is small.
; this is a problem to plot this kind of figure.

psi_def_diff(0,:,:)			= (psi_def(1,:,:)-psi_def(0,:,:));/mask(psi_def(0,:,:),psi_def(0,:,:).eq.0,False) ; p4K-cntl
psi_def_diff(1,:,:)			= (psi_def(2,:,:)-psi_def(0,:,:));/mask(psi_def(0,:,:),psi_def(0,:,:).eq.0,False) ; 4xCO2-cntl
psi_new_diff(0,:,:)			= (psi_new(1,:,:)-psi_new(0,:,:));/mask(psi_new(0,:,:),psi_new(0,:,:).eq.0,False)
psi_new_diff(1,:,:)			= (psi_new(2,:,:)-psi_new(0,:,:));/mask(psi_new(0,:,:),psi_new(0,:,:).eq.0,False)

psi_new_def_diff			= psi_def

psi_new_def_diff(0,:,:)			= (psi_new(0,:,:)-psi_def(0,:,:));/mask(psi_def(0,:,:),psi_def(0,:,:).eq.0,False) ; cntl
psi_new_def_diff(1,:,:)			= (psi_new(1,:,:)-psi_def(1,:,:));/mask(psi_def(1,:,:),psi_def(1,:,:).eq.0,False) ; p4K
psi_new_def_diff(2,:,:)			= (psi_new(2,:,:)-psi_def(2,:,:));/mask(psi_def(2,:,:),psi_def(2,:,:).eq.0,False) ; 4xCO2

printVarSummary(psi_new_def_diff)
printMinMax(psi_new_def_diff,True)

omega_def_diff				= new((/2,nlev,nlon/),float)
omega_def_diff!0			= "exps"
omega_def_diff&exps			= (/0,1/)
omega_def_diff(0,:,:)			= omega_def(0,:,:)

omega_new_diff				= omega_def_diff
omega_def_diff(0,:,:)			= omega_def(1,:,:)-omega_def(0,:,:) ; p4K-cntl
omega_def_diff(1,:,:)			= omega_def(2,:,:)-omega_def(0,:,:) ; 4xCO2-cntl

omega_new_diff(0,:,:)			= omega_new(1,:,:)-omega_new(0,:,:) 
omega_new_diff(1,:,:)			= omega_new(2,:,:)-omega_new(0,:,:)

omega_new_def_diff			= omega_def

omega_new_def_diff(0,:,:)		= omega_new(0,:,:)-omega_def(0,:,:) ; cntl
omega_new_def_diff(1,:,:)		= omega_new(1,:,:)-omega_def(1,:,:) ; p4K
omega_new_def_diff(2,:,:)		= omega_new(2,:,:)-omega_def(2,:,:) ; 4xCO2

printVarSummary(omega_new_def_diff)
printMinMax(omega_new_def_diff,True)


;================================================================
; Plot starts here...
;================================================================
res 					= True

res@gsnDraw 				= False
res@gsnFrame 				= False

res@trYMinF 				= 100.
res@tiYAxisString			= "Pressure (hPa)"
res@tiYAxisFontHeightF			= 0.0325

res@trXMinF				= lonS
res@trXMaxF				= lonE

res@vpWidthF				= 0.7
res@vpHeightF				= 0.4

res@gsnPresHgtHeightLabelOn		= False
res@tmYROn				= False

pmslres 				= res 
lrhes 					= res
psires 					= res
omegres 				= res

;---------------------------------------------------
;---------------------------------------------------
; This first one is for streamfunction - (psi)
;---------------------------------------------------
;---------------------------------------------------
psires 					= True

psires@gsnLeftString 			= ""
psires@gsnRightString 			= ""
psires@gsnStringFontHeightF		= 0.0325

psires@pmTickMarkDisplayMode 		= "Always"
;---------------------------------------------------
; shaded or lines
;---------------------------------------------------
psires@cnFillOn 			= False
psires@cnLinesOn 			= True                       
psires@cnInfoLabelOn         		= False          
psires@cnLineLabelsOn 			= True
;psires@cnFillDrawOrder      		= "PostDraw"

;---------------------------------------------------
; title
;---------------------------------------------------
psires@tiMainFontHeightF 		= 0.012              
;---------------------------------------------------
; label bar 
;---------------------------------------------------
psires@lbLabelAutoStride 		= True
psires@lbLabelOffsetF 			= 0.1
;psires@lbBoxLinesOn          		= False
;psires@lbTitleString 			= "[10~S~-5~N~s~S~-1~N~]"
;psires@lbTitlePosition 		= "Right"
;psires@lbTitleDirection 		= "Across"
;---------------------------------------------------
; line feature: thickness, color, et al
;---------------------------------------------------
psires@cnLineThicknessF 		= 1.5
psires@gsnContourNegLineDashPattern 	= 2
;psires@gsnContourZeroLineThicknessF 	= 6.0
; hiding the zero line
psires@gsnContourZeroLineThicknessF	= 0.0

psires@tmXBOn				= False
psires@tmYLOn				= False
psires@tmYROn				= False

;---------------------------------------------------
; line range
;---------------------------------------------------
psires@cnLineLabelPlacementMode		= "Computed"
psires@cnLineLabelFontHeightF		= 0.03

;psires@cnExplicitLineLabelsOn		= True
;psires@cnLineLabelStrings		= (/6,9,12/)

;---------------------------------------------------
;---------------------------------------------------
; This first one is for omega 
;---------------------------------------------------
;---------------------------------------------------

omegres 				= True

omegres@gsnLeftString 			= ""
omegres@gsnRightString 			= ""

omegres@pmTickMarkDisplayMode 		= "Always"

omegres@cnFillOn 			= True
omegres@cnLinesOn 			= False
omegres@cnInfoLabelOn         		= False         
omegres@cnLineLabelsOn 			= False                
omegres@tiMainFontHeightF 		= 0.012             
omegres@cnFillDrawOrder      		= "PreDraw"
;---------------------------------------------------
; label bar 
;---------------------------------------------------
omegres@lbLabelAutoStride 		= True
;omegres@lbBoxLinesOn          		= False
;omegres@lbTitleString 			= "Omega [10~S~-3~N~ hPa s~S~-1~N~]"
;omegres@lbTitlePosition 		= "Right"
;omegres@lbTitleDirection 		= "Across"

omegres@lbTitleString			= "hPa/day"
;omegres@lbTitleJust			= "BottomRight"

;omegres@lbTitlePosition			= "Right"
;omegres@lbTitleDirection		= "Across"

omegres@lbTitleFontHeightF		= 0.03
;---------------------------------------------------
; line range
;---------------------------------------------------
fontheight				= 0.0325

omegres@tmXBLabelFontHeightF		= fontheight 
omegres@tmYLLabelFontHeightF		= fontheight 
omegres@tmYRLabelFontHeightF		= fontheight 
omegres@gsnStringFontHeightF		= fontheight 
omegres@lbLabelFontHeightF		= fontheight 

omegres@tmYLMode                    	= "Explicit"
omegres@tmYLValues                  	= (/850,700,500,300,200/)
omegres@tmYLLabels                  	= (/850,700,500,300,200/)

omegres@tmXBTickSpacingF		= 60
omegres@gsnRightStringOrthogonalPosF 	= .02

;*******************((CREATE PLOTS))******************
;*****************************************************
wks_type 				= "x11"
;wks_type@wkWidth 			= 1250
;wks_type@wkHeight 			= 1250

plotName 				= "omega-u-zpsi_"+latS+"to"+latE+"_revised"
wks_name 				= "./figure/"+plotName
wks 					= gsn_open_wks(wks_type,wks_name)


gsn_define_colormap(wks,"BlueYellowRed")


psires@cnLevelSelectionMode 		= "ManualLevels"
psires@cnMinLevelValF       		= -15.
psires@cnMaxLevelValF       		= 15. 
psires@cnLevelSpacingF      		= 3. 

omegres@cnLevelSelectionMode 		= "ManualLevels"
omegres@cnMinLevelValF       		= -30.
omegres@cnMaxLevelValF       		= 30.
omegres@cnLevelSpacingF      		= 5. 


; resource change for difference plots
psires1					= psires
omegres1				= omegres

psires1@cnLevelSelectionMode 		= "ManualLevels"
psires1@cnMinLevelValF       		= -2.
psires1@cnMaxLevelValF       		= 2.
psires1@cnLevelSpacingF      		= 0.5

omegres1@cnLevelSelectionMode 		= "ManualLevels"
omegres1@cnMinLevelValF       		= -10.
omegres1@cnMaxLevelValF       		= 10.
omegres1@cnLevelSpacingF      		= 2.


psi4K					= new((/8,nlev,nlon/),float)
psi4K(0,:,:)				= psi_def(0,:,:)
psi4K					= (/psi_def(0,:,:),psi_new(0,:,:),psi_new_def_diff(0,:,:),\
					    psi_def(1,:,:),psi_new(1,:,:),psi_new_def_diff(1,:,:),\
					    psi_def_diff(0,:,:),psi_new_diff(0,:,:)/)
printVarSummary(psi4K)

omega4K					= new((/8,nlev,nlon/),float)
omega4K(0,:,:)				= omega_def(0,:,:)
omega4K					= (/omega_def(0,:,:),omega_new(0,:,:),omega_new_def_diff(0,:,:),\
					    omega_def(1,:,:),omega_new(1,:,:),omega_new_def_diff(1,:,:),\
					    omega_def_diff(0,:,:),omega_new_diff(0,:,:)/)

psi4xCO2				= new((/8,nlev,nlon/),float)
psi4xCO2(0,:,:)				= psi_def(0,:,:)
psi4xCO2				= (/psi_def(0,:,:),psi_new(0,:,:),psi_new_def_diff(0,:,:),\
					    psi_def(2,:,:),psi_new(2,:,:),psi_new_def_diff(2,:,:),\
					    psi_def_diff(1,:,:),psi_new_diff(1,:,:)/)

omega4xCO2				= new((/8,nlev,nlon/),float)
omega4xCO2(0,:,:)			= omega_def(0,:,:)
omega4xCO2				= (/omega_def(0,:,:),omega_new(0,:,:),omega_new_def_diff(0,:,:),\
					    omega_def(2,:,:),omega_new(2,:,:),omega_new_def_diff(2,:,:),\
					    omega_def_diff(1,:,:),omega_new_diff(1,:,:)/)

strings4K				= (/"RH","PDF","PDF-RH",\
					    "RH_4K","PDF_4K","PDF_4K-RH_4K",\
					    "RH_4K-RH","PDF_4K-PDF"/)

strings4xCO2				= (/"RH","PDF","PDF-RH",\
					    "RH_4xCO2","PDF_4xCO2","PDF_4xCO2-RH_4xCO2",\
					    "RH_4xCO2-RH","PDF_4xCO2-PDF"/)


plot					= new(9,graphic)
plotov					= new(9,graphic)

plot1					= new(9,graphic)
plotov1					= new(9,graphic)


do it=0,7
if(it.eq.0.or.it.eq.1.or.it.eq.3.or.it.eq.4)then
;p4K
	plot(it)				= gsn_csm_pres_hgt(wks,psi4K(it,:,:),psires)
	omegres@gsnLeftString			= strings4K(it)
	plotov(it)				= gsn_csm_pres_hgt(wks,omega4K(it,:,:),omegres)
;4xCO2
	plot1(it)				= gsn_csm_pres_hgt(wks,psi4xCO2(it,:,:),psires)
	omegres@gsnLeftString			= strings4xCO2(it)
	plotov1(it)				= gsn_csm_pres_hgt(wks,omega4xCO2(it,:,:),omegres)

else
;p4K
	plot(it)				= gsn_csm_pres_hgt(wks,psi4K(it,:,:),psires1)
	omegres1@gsnLeftString			= strings4K(it)
	plotov(it)				= gsn_csm_pres_hgt(wks,omega4K(it,:,:),omegres1)
;4xCO2
	plot1(it)				= gsn_csm_pres_hgt(wks,psi4xCO2(it,:,:),psires1)
	omegres1@gsnLeftString			= strings4xCO2(it)
	plotov1(it)				= gsn_csm_pres_hgt(wks,omega4xCO2(it,:,:),omegres1)
end if
overlay(plot(it),plotov(it))
overlay(plot1(it),plotov1(it))
end do ; do it=

resP				= True
;resP@amJust                     = "BottomLeft"
;resP@gsnPanelFigureStrings      = (/"a)","b)","c)"/)

gsn_panel(wks,plot,(/3,3/),resP)
gsn_panel(wks,plot1,(/3,3/),resP)


print( "Figures are finished")
system( "date" )

end


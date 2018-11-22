load "./functions_read_file.ncl"
load "./functions_calculate.ncl"
load "./functions_plot.ncl"
;load "./functions_discretize.ncl"
load "./functions_regrid.ncl"
load "./functions_bony_method.ncl"

;;;;;processing steps

begin

landmask			= False
ocnmask				= False

;fdir 				= "/home/lyl/WORK4/cesm1_2_1/archive/"
fdir				= "/home/lyl/WORK3/qinyi/data/cfmip/"
fappend 			= "/atm/hist/"

case_nor			= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_mac2_amip_outSGM"/)
case_p4K			= (/"FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM"/)
;case_4xCO2			= (/"FAMIPC5_f09f09_MG15_amip-4xCO2","FAMIPC5_f09f09_mac2_amip-4xCO2"/)

;case_all			= (/case_nor, case_p4K, case_4xCO2/)
case_all			= (/case_nor, case_p4K/)

printVarSummary(case_all) ; 3x2

;Vars	= (/"OMEGA", "DCCLDLIQ", "ZMDLIQ", "DPDLFLIQ", "CMFDLIQ", "SHDLFLIQ", "MACPDLIQ", "MPDLIQ", "VDCLDLIQ"/)
Vars	= (/"OMEGA",  "CMFDLIQ", "SHDLFLIQ", "MACPDLIQ", "MPDLIQ"/)

Vars_mid	= (/"T", "Q", "OMEGA"/)


tag				= (/"nor","p4K"/)

;color	= (/ "grey",	  "red",    "pink",     "blue",     "cyan",      "purple",   "green",   "brown"/)
color	= (/ "grey", "red", "blue",  "cyan" /)

thick	= (/1.5,3.0/)
dash	= (/14,0/)

colors	= conform_dims((/(dimsizes(Vars)-1),dimsizes(tag),dimsizes(case_nor)/),color,0)
thicks	= conform_dims((/(dimsizes(Vars)-1),dimsizes(tag),dimsizes(case_nor)/),thick,2)
dashes	= conform_dims((/(dimsizes(Vars)-1),dimsizes(tag),dimsizes(case_nor)/),dash,1)


nlat				= 192
nlon				= 288
nlev				= 30

; ============ read data in

vars_all			= new((/dimsizes(Vars),dimsizes(tag),dimsizes(case_nor),nlev,nlat,nlon/),float)
vars_all!0			= "var"
vars_all&var			= fspan(0,dimsizes(Vars)-1,dimsizes(Vars))
vars_all!1			= "tag"
vars_all&tag			= fspan(0,dimsizes(tag)-1,dimsizes(tag))
vars_all!2			= "case"
vars_all&case			= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))


vars_mid			= new((/dimsizes(Vars_mid),dimsizes(tag),dimsizes(case_nor),nlev,nlat,nlon/),float)
vars_mid!0			= "var"
vars_mid&var		= fspan(0,dimsizes(Vars_mid)-1,dimsizes(Vars_mid))
vars_mid!1		 	= "tag"
vars_mid&tag		= fspan(0,dimsizes(tag)-1,dimsizes(tag))
vars_mid!2		 	= "case"
vars_mid&case		= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))


do itag=0,dimsizes(tag)-1	; nor, p4K, 4xCO2
do icase=0,dimsizes(case_nor)-1 ; def, new
	ftmp					= addfile(fdir+case_all(itag,icase)+".ANN.climo.nc","r")
	lev									= ftmp->lev
	OCNFRAC								= ftmp->OCNFRAC(0,:,:)
	printVarSummary(OCNFRAC)
do ivar=0,dimsizes(Vars)-1
	vars_all(ivar,itag,icase,:,:,:)		= ftmp->$Vars(ivar)$(0,:,:,:)
end do ; do ivar=
do ivar=0,dimsizes(Vars_mid)-1
	vars_mid(ivar,itag,icase,:,:,:)		= ftmp->$Vars_mid(ivar)$(0,:,:,:)
end do ; do ivar=

end do ; do icase=
end do ; do itag=

; change unit from kg/kg/s to g/kg/day
vars_all(1:,:,:,:,:,:)				= vars_all(1:,:,:,:,:,:)*8.64e7

omega_500		= vars_all(0,:,:,{500},:,:)

printVarSummary(vars_all)
printVarSummary(omega_500)

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
copy_VarCoords(vars_all(0,:,:,0,:,:),LTS)
printVarSummary(LTS)
printMinMax(LTS,True)

; ============== situation 1: mask ascending regions
lat					= vars_all&lat
lon					= vars_all&lon

lat2d					= conform_dims(dimsizes(omega_500),lat,2)
lon2d					= conform_dims(dimsizes(omega_500),lon,3)

vars_all_omega		= vars_all
do ivar = 1,dimsizes(Vars)-1
do ilev = 1,nlev-1
;	vars_all_omega(ivar,:,:,ilev,:,:) = where(omega_500.gt.0.,vars_all(ivar,:,:,ilev,:,:),vars_all@_FillValue)
	vars_all_omega(ivar,:,:,ilev,:,:) = where(lat2d.gt.-30. .and. lat2d.lt.30. .and. omega_500.gt.0. .and. LTS.lt.17,vars_all(ivar,:,:,ilev,:,:),vars_all@_FillValue)
end do ; do ilev=
end do ; do ivar=
printVarSummary(vars_all_omega)
printMinMax(vars_all_omega,True)

;wks		= gsn_open_wks("x11","test")
;
;res		= True
;res@cnFillOn	= True
;res@cnLinesOn	= False
;
;plot		= gsn_csm_contour_map_ce(wks,vars_all_omega(2,0,0,{700},:,:),res)

; ============= global mean
rad    			= 4.0*atan(1.0)/180.0
lat			= vars_all&lat
lat!0			= "lat"
lat&lat			= vars_all&lat
clat   			= cos(lat*rad)
copy_VarCoords(lat,clat)

latS			= -30.
latE			= 0.
lonS			= 180.
lonE			= 270.

vars_all_omega_gm	= wgt_areaave_Wrap(vars_all_omega,clat,1.0,1)
printVarSummary(vars_all_omega_gm)

vars_all_omega_gm_rd	= vars_all_omega_gm
vars_all_omega_gm_rd(:,:,0,:)	= (/vars_all_omega_gm(:,:,0,:)-vars_all_omega_gm(:,:,0,:)/);/mask(vars_all_omega_gm(:,:,0,:),vars_all_omega_gm(:,:,0,:).lt.0.01,False)
vars_all_omega_gm_rd(:,:,1,:)	= (/vars_all_omega_gm(:,:,1,:)-vars_all_omega_gm(:,:,0,:)/);/mask(vars_all_omega_gm(:,:,0,:),vars_all_omega_gm(:,:,0,:).lt.0.01,False)

; ============ plotting

wks			= gsn_open_wks("x11","./figure/liquid_budget_tendency")

res			= True
res@gsnFrame		= False
res@gsnDraw		= False
res@trYReverse		= True 
res@tiYAxisString	= "Pressure (hPa)"

res@xyLineColors	= ndtooned(colors)
res@xyLineThicknesses	= ndtooned(thicks)
res@xyDashPatterns	= ndtooned(dashes)

colors	= conform_dims((/(dimsizes(Vars)-1),dimsizes(tag),dimsizes(case_nor)/),color,0) ; variables
thicks	= conform_dims((/(dimsizes(Vars)-1),dimsizes(tag),dimsizes(case_nor)/),thick,2) ; def, new
dashes	= conform_dims((/(dimsizes(Vars)-1),dimsizes(tag),dimsizes(case_nor)/),dash,1) ; cntl, p4K


tmp1			= reshape(vars_all_omega_gm(1:,:,:,:),(/(dimsizes(Vars)-1)*dimsizes(tag)*dimsizes(case_nor),nlev/))
tmp2			= reshape(vars_all_omega_gm_rd(1:,:,:,:),(/(dimsizes(Vars)-1)*dimsizes(tag)*dimsizes(case_nor),nlev/))

tmp1!0			= "test"
tmp1!1			= "lev"
tmp1&lev		= vars_all&lev
copy_VarCoords(tmp1,tmp2)

plot1	 		= gsn_csm_xy(wks,tmp1,tmp1&lev,res)
plot2	 		= gsn_csm_xy(wks,tmp2,tmp2&lev,res)

;plot			= gsn_csm_xy(wks,vars_all_omega_gm(1:,0,0,:),vars_all_omega_gm&lev,res)



lgres                          	= True
lgres@lgLabelFontHeightF        = 0.08
lgres@lgPerimOn                 = False

legends2                        = (/"CTL","P4K"/)
lgres1			 	= lgres
lgres1@vpWidthF                 = 0.150
lgres1@vpHeightF                = 0.10
lgres1@lgMonoLineColor		= True
lgres1@lgLineColor		= "black"
lgres1@lgLineThicknessF		= 3.0
lgres1@lgDashIndexes		= dash

legends3                        = (/"Park-RH","Gauss"/)
lgres2			 	= lgres
lgres2@vpWidthF                  = 0.20
lgres2@vpHeightF                = 0.10
lgres2@lgMonoLineColor		= True
lgres2@lgLineColor		= "black"
lgres2@lgMonoLineThickness 	= False
lgres2@lgLineThicknesses	= (/1.5,3.0/)
lgres2@lgMonoDashIndex		= True

legends1                        = Vars(1:)
lgres@vpWidthF                  = 0.20
lgres@vpHeightF                 = 0.20
lgres@lgLineColors		= color
lgres@lgLineThicknessF		= 3.0
lgres@lgDashIndex		= True


lbid1                           = gsn_create_legend(wks,dimsizes(legends1),legends1,lgres)
lbid2                           = gsn_create_legend(wks,dimsizes(legends2),legends2,lgres1)
lbid3                           = gsn_create_legend(wks,dimsizes(legends3),legends3,lgres2)


amres                           = True
amres@amJust			= "TopRight"
amres@amParallelPosF            = 0.5
amres@amOrthogonalPosF          = -0.50
annoid1                         = gsn_add_annotation(plot1,lbid1,amres)

amres@amParallelPosF            = 0.5
amres@amOrthogonalPosF          = -0.050
annoid3                         = gsn_add_annotation(plot1,lbid3,amres)

amres@amParallelPosF            = -0.2
amres@amOrthogonalPosF          = -0.5
annoid2                         = gsn_add_annotation(plot1,lbid2,amres)

gsn_panel(wks,(/plot1,plot2/),(/1,2/),False)
end


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

function sfvp(lev[*]:numeric,lat[*]:numeric,lon[*]:numeric,\
              U_1[*][*][*]:numeric,V_1[*][*][*]:numeric,OMEGA_1[*][*][*],\
	      gw[*]:numeric,lat1:float,lat2:float)

begin

nlev = dimsizes(lev)
nlat = dimsizes(lat)
nlon = dimsizes(lon)

sfvp_1 = new((/2,nlev,nlat,nlon/),float)
sfvp_1(0,:,:,:) = U_1
sfvp_1 = uv2sfvpG(U_1,V_1)

vp_1 = sfvp_1(0,:,:,:)

DvpDx_1 = new((/nlev,nlat,nlon/),float)
DvpDx_1 = vp_1
DvpDx_1 = doubletofloat(center_finite_diff_n(vp_1,lon*111000,False,0,2))

printVarSummary(DvpDx_1)

DvpDx_1_mm = dim_avg_wgt_n_Wrap(DvpDx_1(:,{lat1:lat2},:),gw({lat1:lat2}),0,1)
lev111 = DvpDx_1_mm&lev
lev111@units = "hPa"
DvpDx_1_mm&lev = lev111 
printVarSummary(DvpDx_1_mm)

OMEGA_1_mm = dim_avg_wgt_n_Wrap(OMEGA_1(:,{lat1:lat2},:),gw({lat1:lat2}),0,1)
lev111 = OMEGA_1_mm&lev
lev111@units = "hPa"
OMEGA_1_mm&lev = lev111 
printVarSummary(OMEGA_1_mm)

printMinMax(OMEGA_1_mm,False)

OMEGA_1_mm = (/-1.0*OMEGA_1_mm*1e2/)

return([/DvpDx_1_mm,OMEGA_1_mm/])

end


begin

;***************file 1********************
fdir = "/home/lyl/WORK3/qinyi/scripts/cfmip/data/output/"
case = (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_mac2_amip",\
         "FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_mac2_amip-p4K_debug",\
         "FAMIPC5_f09f09_MG15_amip-4xCO2","FAMIPC5_f09f09_mac2_amip-4xCO2"/)

append = ".cam.h0.ANN.nc"

f0 = addfile(fdir+case(0)+append,"r")
f1 = addfile(fdir+case(1)+append,"r")
f2 = addfile(fdir+case(2)+append,"r")
f3 = addfile(fdir+case(3)+append,"r")
f4 = addfile(fdir+case(4)+append,"r")
f5 = addfile(fdir+case(5)+append,"r")

press = f1->lev ;lev(lev)
temp = f1->T(0,:,:,:) ;T(time, lev, lat, lon))
gw = f1->gw
lev = f1->lev
lat = f1->lat
lon = f1->lon

nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = dimsizes(lev)


    U_1 = f0->U(0,:,:,:)
    V_1 = f0->V(0,:,:,:)
OMEGA_1 = f0->OMEGA(0,:,:,:)

    U_2 = f1->U(0,:,:,:)
    V_2 = f1->V(0,:,:,:)
OMEGA_2 = f1->OMEGA(0,:,:,:)

    U_3 = f2->U(0,:,:,:)
    V_3 = f2->V(0,:,:,:)
OMEGA_3 = f2->OMEGA(0,:,:,:)

    U_4 = f3->U(0,:,:,:)
    V_4 = f3->V(0,:,:,:)
OMEGA_4 = f3->OMEGA(0,:,:,:)

    U_5 = f4->U(0,:,:,:)
    V_5 = f4->V(0,:,:,:)
OMEGA_5 = f4->OMEGA(0,:,:,:)

    U_6 = f5->U(0,:,:,:)
    V_6 = f5->V(0,:,:,:)
OMEGA_6 = f5->OMEGA(0,:,:,:)



;1111111111111111111111111111111

lat1 = -10.
lat2 = 10.

result = sfvp(lev,lat,lon,U_1,V_1,OMEGA_1,gw,lat1,lat2)
DvpDx_1_mm = result[0]
OMEGA_1_mm = result[1]
delete(result)

printVarSummary(DvpDx_1_mm)

result = sfvp(lev,lat,lon,U_2,V_2,OMEGA_2,gw,lat1,lat2)
DvpDx_2_mm = result[0]
OMEGA_2_mm = result[1]
delete(result)

result = sfvp(lev,lat,lon,U_2,V_2,OMEGA_2,gw,lat1,lat2)
DvpDx_2_mm = result[0]
OMEGA_2_mm = result[1]
delete(result)

result = sfvp(lev,lat,lon,U_3,V_3,OMEGA_3,gw,lat1,lat2)
DvpDx_3_mm = result[0]
OMEGA_3_mm = result[1]
delete(result)

result = sfvp(lev,lat,lon,U_4,V_4,OMEGA_4,gw,lat1,lat2)
DvpDx_4_mm = result[0]
OMEGA_4_mm = result[1]
delete(result)

result = sfvp(lev,lat,lon,U_5,V_5,OMEGA_5,gw,lat1,lat2)
DvpDx_5_mm = result[0]
OMEGA_5_mm = result[1]
delete(result)

result = sfvp(lev,lat,lon,U_6,V_6,OMEGA_6,gw,lat1,lat2)
DvpDx_6_mm = result[0]
OMEGA_6_mm = result[1]
delete(result)

;************************create plot***************
wks = gsn_open_wks("x11","./on-height-omega-CAM5-PDF")
gsn_define_colormap(wks,"BlueRed")


res = True
res@trYReverse = True
res@gsnDraw = False
res@gsnFrame = False
;res@gsnSpreadColors = True
;
;res@tiYAxisString = "pressure(hPa)"
;
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnFillOn = True
res@cnInfoLabelOn = False
;
;res@gsnYAxisIrregular2Linear = True
;
res@vpHeightF = 0.3
res@vpWidthF = 0.7

res@trXMinF = 140.
res@trXMaxF = 280.

res@trYMinF = 100.

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -1.8
res@cnMaxLevelValF = 1.8
res@cnLevelSpacingF = 0.2


; plot = gsn_csm_contour(wks,DvpDx_mm(::3,::15),res)

vcres = True
vcres@gsnDraw = False
vcres@gsnFrame = False

vcres@trYReverse = True

vcres@trXMinF = 140.
vcres@trXMaxF = 280.
vcres@trYMinF = 100.

vcres@vpHeightF = 0.3
vcres@vpWidthF = 0.7

vcres@tiYAxisString = "Pressure (hPa)"
vcres@gsnRightString = "W and divergent U"

vcres@cnLinesOn = False
vcres@cnLineLabelsOn = False
vcres@cnFillOn = True
vcres@cnInfoLabelOn = False

vcres@cnLevelSelectionMode = "ManualLevels"
vcres@cnMinLevelValF = -1.0
vcres@cnMaxLevelValF = 1.0
vcres@cnLevelSpacingF = 0.25

vcres@vcRefMagnitudeF  = 2.0                       ; add a reference vector
vcres@vcRefLengthF     = 0.015                     ; what the ref length is
vcres@vcGlyphStyle     = "CurlyVector" 
vcres@vcLineArrowThicknessF = 3.5
vcres@vcLevelSpacingF = 8.0 
vcres@vcRefAnnoOrthogonalPosF = -1.15
vcres@vcRefAnnoArrowSpaceF = 2.0
vcres@vcRefAnnoString2On = False
; vcres@vcRefAnnoString1On = False
vcres@vcVectorDrawOrder        = "Postdraw"

ovplot = new(4, graphic)
DvpDx_21_mm = DvpDx_1_mm 
DvpDx_21_mm = (/DvpDx_2_mm-DvpDx_1_mm/)
OMEGA_21_mm = OMEGA_1_mm 
OMEGA_21_mm = (/OMEGA_2_mm-OMEGA_1_mm/)

DvpDx_43_mm = DvpDx_3_mm 
DvpDx_43_mm = (/DvpDx_4_mm-DvpDx_3_mm/)
OMEGA_43_mm = OMEGA_3_mm 
OMEGA_43_mm = (/OMEGA_4_mm-OMEGA_3_mm/)

DvpDx_65_mm = DvpDx_5_mm 
DvpDx_65_mm = (/DvpDx_6_mm-DvpDx_5_mm/)
OMEGA_65_mm = OMEGA_5_mm 
OMEGA_65_mm = (/OMEGA_6_mm-OMEGA_5_mm/)

vcres@vcRefMagnitudeF  = 0.2

;vcres@gsnLeftString = "CNTL:PDF-CAM5"
;ovplot(0) = gsn_csm_pres_hgt_vector(wks,DvpDx_21_mm, DvpDx_21_mm(::2,::4),OMEGA_21_mm(::2,::4),vcres)
;vcres@gsnLeftString = "p4K: PDF-CAM5"
;ovplot(1) = gsn_csm_pres_hgt_vector(wks,DvpDx_43_mm, DvpDx_43_mm(::2,::4),OMEGA_43_mm(::2,::4),vcres)
;vcres@gsnLeftString = "4xCO2: PDF-CAM5"
;ovplot(2) = gsn_csm_pres_hgt_vector(wks,DvpDx_65_mm, DvpDx_65_mm(::2,::4),OMEGA_65_mm(::2,::4),vcres)




;;;;;;;;;;;;;;;;;;;;;;;;;;;
DvpDx_31_mm = DvpDx_1_mm 
DvpDx_31_mm = (/DvpDx_3_mm-DvpDx_1_mm/)
OMEGA_31_mm = OMEGA_1_mm 
OMEGA_31_mm = (/OMEGA_3_mm-OMEGA_1_mm/)

DvpDx_51_mm = DvpDx_1_mm 
DvpDx_51_mm = (/DvpDx_5_mm-DvpDx_1_mm/)
OMEGA_51_mm = OMEGA_1_mm 
OMEGA_51_mm = (/OMEGA_5_mm-OMEGA_1_mm/)

DvpDx_42_mm = DvpDx_2_mm 
DvpDx_42_mm = (/DvpDx_4_mm-DvpDx_2_mm/)
OMEGA_42_mm = OMEGA_2_mm 
OMEGA_42_mm = (/OMEGA_4_mm-OMEGA_2_mm/)

DvpDx_62_mm = DvpDx_2_mm 
DvpDx_62_mm = (/DvpDx_6_mm-DvpDx_2_mm/)
OMEGA_62_mm = OMEGA_2_mm 
OMEGA_62_mm = (/OMEGA_6_mm-OMEGA_2_mm/)



vcres@vcRefMagnitudeF  = 0.4
vcres@gsnLeftString = "CAM5:p4K-cntl"
ovplot(0) = gsn_csm_pres_hgt_vector(wks,DvpDx_31_mm, DvpDx_31_mm(::2,::4),OMEGA_31_mm(::2,::4),vcres)
vcres@gsnLeftString = "CAM5:4xCO2-cntl"
ovplot(1) = gsn_csm_pres_hgt_vector(wks,DvpDx_51_mm, DvpDx_51_mm(::2,::4),OMEGA_51_mm(::2,::4),vcres)
vcres@gsnLeftString = "PDF:p4K-cntl"
ovplot(2) = gsn_csm_pres_hgt_vector(wks,DvpDx_42_mm, DvpDx_42_mm(::2,::4),OMEGA_42_mm(::2,::4),vcres)
vcres@gsnLeftString = "PDF:4xCO2-cntl"
ovplot(3) = gsn_csm_pres_hgt_vector(wks,DvpDx_62_mm, DvpDx_62_mm(::2,::4),OMEGA_62_mm(::2,::4),vcres)

; overlay(plot,ovplot)

pres = True
pres@gsnMaximize = True
gsn_panel(wks, ovplot, (/3,1/), pres)
end

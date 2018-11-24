
load "./functions_regrid.ncl"
load "./functions_getvar.ncl"
load "./functions_calculate.ncl"

begin

;***************file 1********************

;fdir 				= "/home/lyl/WORK4/cesm1_2_1/archive/"
;fdir_append			= "/atm/hist/"

pattern_correlation	= True
regime_contribution = True
regime_mask_plot	= False
lat_lon_plot		= False
Sherwood_2014		= False

sherwood_index		= (/"SS","DD","LTMI"/)

plot_type			= "pdf"

select_idx			= 0	; choose how many regimes would be plotted 
scale_factor		= 3.	; use to amplify diff to make them able to plot in non-diff plot
regime_mask_val		= (/0,1,2,3,4,5,6/)
LTS_criteria		= 17.


fdir				= "/home/share3/lyl/work3/qinyi/mid-data/cfmip/"

;case_nor 				= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_MG15_amip-p4K"/)
;case_p4K		 		= (/"FAMIPC5_f09f09_mac2_amip_outSGM","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM"/)
case_nor 				= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_mac2_amip_outSGM"/)
case_p4K		 		= (/"FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM"/)

case_name				= (/"RH","PDF"/)

case					= (/case_nor,case_p4K/)
tag					= (/"nor", "p4K"/)

append 				= ".ANN.climo.nc"

; _1 for the difference between PDF and RH; 
; _2 for the difference between p4K and cntl.
SWCF_levs			= (/-170,-150,-135,-120,-105,-90,-75,-60,-45,-30,-15,0,15,30,45/)
SWCF_dlevs_1		= (/-50,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,50/)/2.
SWCF_dlevs_2		= (/-50,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,50/)/2.
SWCF_dlevs_sp		= (/-10,-5,5,10/)

LWCF_levs			= (/-45,-30,-20,-10,0,10,20,30,40,50,60,70,85,100,115/) 
LWCF_dlevs_1			= (/-35,-30,-15,-10,-6,-4,-2, 0, 2,4,6,10,15,30,35/)/2.
LWCF_dlevs_2			= (/-35,-30,-15,-10,-6,-4,-2, 0, 2,4,6,10,15,30,35/)/2.
LWCF_dlevs_sp		= (/-10,-5,5,10/)

PRECT_levs			= (/.2,.5,1,2,3,4,5,6,7,8,9,10,12,14,17/)
PRECT_dlevs_1			= (/-4,-3,-2.5,-2,-1.5,-1,-.5,0,.5,1,1.5,2,2.5,3,4/)
PRECT_dlevs_2			= (/-4,-3,-2.5,-2,-1.5,-1,-.5,0,.5,1,1.5,2,2.5,3,4/)
PRECT_dlevs_sp		= (/-2,-1,1,2/)

CLDLOW_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
CLDLOW_dlevs_1		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
CLDLOW_dlevs_2		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
CLDLOW_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

CLDHGH_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
CLDHGH_dlevs_1		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
CLDHGH_dlevs_2		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
CLDHGH_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

LTS_levs			= (/5.,6,7,8,9,10,12,14,16,18,20,22,24,26,28/)
LTS_dlevs_1			= (/-1.2,-1.,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1,1.2/)
LTS_dlevs_2			= (/-1.2,-1.,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1,1.2/)*3.
LTS_dlevs_sp		= (/-1.0, -0.5, 0.5, 1.0/)

EIS_levs			= (/-3.,-2,-1,0,1,2,3,4,5,6,8,10,12,13,14/)
EIS_dlevs_1			= (/-1.2,-1.,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1,1.2/)
EIS_dlevs_2			= (/-1.2,-1.,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1,1.2/)*3.
EIS_dlevs_sp		= (/-1.0, -0.5, 0.5, 1.0/)

ECTEI_levs			= (/-12.,-10.,-8.0, -6.0, -4.0, -3.,-2,-1,0,1,2,4,6,8,10/)
ECTEI_dlevs_1			= (/-1.2,-1.,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1,1.2/)
ECTEI_dlevs_2			= (/-1.2,-1.,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1,1.2/)*3.
ECTEI_dlevs_sp		= (/-1.0, -0.5, 0.5, 1.0/)

TGCLDLWP_levs 		= (/10,25,50,75,100,150,175,200,250,300,350,400,450,500,550/)
TGCLDLWP_dlevs_1	 	= (/-120,-100,-80,-60,-40,-20,-10,0,10,20,40,60,80,100,120/)
TGCLDLWP_dlevs_2	 	= (/-120,-100,-80,-60,-40,-20,-10,0,10,20,40,60,80,100,120/)
TGCLDLWP_dlevs_sp	= (/-100.,-50., 50., 100. /)

sgm_tota_levs		= fspan(0.,0.84,15)
sgm_tota_dlevs_1	= fspan(-0.28,0.28,15)
sgm_tota_dlevs_2	= sgm_tota_dlevs_1
sgm_tota_dlevs_sp	= (/-1, -0.5, 0.5, 1/)
print(sgm_tota_levs+"	"+sgm_tota_dlevs_1)

N1_levs				= fspan(0,0.84,15)
N1_dlevs_1			= fspan(-0.28,0.28,15)
N1_dlevs_2			= N1_dlevs_1
N1_dlevs_sp			= (/-1, -0.5, 0., 1/)

netCRF_levs			= (/-80,-75,-70,-65,-60,-55,-50,-45,-40,-30,-20,-10,0,10,20/)
netCRF_dlevs_1		= (/-50,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,50/)/2.
netCRF_dlevs_2		= netCRF_dlevs_1
netCRF_dlevs_sp		= (/-20,-10,10,20/)
print(netCRF_levs+"		"+netCRF_dlevs_1)


Vars_mid	= (/"T", "Q", "OMEGA", "RELHUM"/)

;				   0            1          2           3           4            5          6	    7      8       9		10	      11	  12 	 13
varsA		= (/"PRECT",	"PRECC",	"PRECL",	"SWCF",		"LWCF",		"CLDLOW",	"CLDHGH", "LTS", "EIS", "ECTEI","TGCLDLWP","sgm_tota","N1", "netCRF"/)
varsA_levs	= (/PRECT_levs, PRECT_levs,     PRECT_levs,     SWCF_levs,      LWCF_levs, 	CLDLOW_levs,	CLDHGH_levs, LTS_levs, EIS_levs, ECTEI_levs,TGCLDLWP_levs,sgm_tota_levs, N1_levs, netCRF_levs/)
varsA_dlevs_1	= (/PRECT_dlevs_1,PRECT_dlevs_1,    PRECT_dlevs_1,    SWCF_dlevs_1,     LWCF_dlevs_1,	CLDLOW_dlevs_1,	CLDHGH_dlevs_1, LTS_dlevs_1, EIS_dlevs_1,	ECTEI_dlevs_1,TGCLDLWP_dlevs_1, sgm_tota_dlevs_1, N1_dlevs_1, netCRF_dlevs_1/)
varsA_dlevs_2	= (/PRECT_dlevs_2,PRECT_dlevs_2,    PRECT_dlevs_2,    SWCF_dlevs_2,     LWCF_dlevs_2,	CLDLOW_dlevs_2,	CLDHGH_dlevs_2, LTS_dlevs_2, EIS_dlevs_2,	ECTEI_dlevs_2,TGCLDLWP_dlevs_2, sgm_tota_dlevs_2, N1_dlevs_2, netCRF_dlevs_2/)
varsA_dlevs_sp	= (/PRECT_dlevs_sp,PRECT_dlevs_sp,    PRECT_dlevs_sp,    SWCF_dlevs_sp,     LWCF_dlevs_sp,	CLDLOW_dlevs_sp,	CLDHGH_dlevs_sp, LTS_dlevs_sp, EIS_dlevs_sp, ECTEI_dlevs_sp,TGCLDLWP_dlevs_sp, sgm_tota_dlevs_sp, N1_dlevs_sp, netCRF_dlevs_sp/)
varsA_units	= (/"mm/day",   "mm/day",       "mm/day",       "W/m~S~2~N~",   "W/m~S~2~N~",	"fraction",	"fraction", "K",	"K",	"K", "g/m2", "No", "No", "W/m~S~2~N~"/)

;use_idx		= (/3,4,5,10/)
;use_idx			= (/9,7,8/)
use_idx			= (/13,3,5,10,11,12/)
;use_idx			= (/13,3,5,10/)

vars			= varsA(use_idx)
vars_levs		= varsA_levs(use_idx,:)
vars_dlevs_1		= varsA_dlevs_1(use_idx,:)
vars_dlevs_2		= varsA_dlevs_2(use_idx,:)

vars_dlevs_sp	= varsA_dlevs_sp(use_idx,:)
vars_units		= varsA_units(use_idx)


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

	do ivar=0,dimsizes(Vars_mid)-1
		vars_mid(ivar,itag,icase,:,:,:)		= ftmp1->$Vars_mid(ivar)$(0,:,:,:)
	end do ; do ivar=

	do ivar=0,dimsizes(vars)-1
		if(isfilevar(ftmp1,vars(ivar))) then
			tmp									= ftmp1->$vars(ivar)$
			if(dimsizes(dimsizes(tmp)).eq.4)then ; figure out variables having vertical level coordinate
				vars_all(ivar,itag,icase,:,:)		= ftmp1->$vars(ivar)$(0,{887},:,:)
			else
				vars_all(ivar,itag,icase,:,:)	 	= ftmp1->$vars(ivar)$(0,:,:)
			end if

			if(vars(ivar).eq."TGCLDLWP")then
				vars_all(ivar,itag,icase,:,:)	= vars_all(ivar,itag,icase,:,:)*1e3
			end if
			if(vars(ivar).eq."sgm_tota")then
				vars_all(ivar,itag,icase,:,:)	= sqrt(abs(vars_all(ivar,itag,icase,:,:)))*1e3
			end if
			delete(tmp)
		else
			if(vars(ivar).eq."LTS")then
				result					= cal_ECTEI(vars_mid(0,itag,icase,:,:,:),vars_mid(1,itag,icase,:,:,:),lev)
				LTS_tmp					= result[0]
				vars_all(ivar,itag,icase,:,:)	= LTS_tmp
				delete(LTS_tmp)
			end if
			if(vars(ivar).eq."EIS")then
				result					= cal_ECTEI(vars_mid(0,itag,icase,:,:,:),vars_mid(1,itag,icase,:,:,:),lev)
				EIS_tmp					= tofloat(result[1])
				vars_all(ivar,itag,icase,:,:)	= (/EIS_tmp/)
				delete(EIS_tmp)
			end if
			if(vars(ivar).eq."ECTEI")then
				result					= cal_ECTEI(vars_mid(0,itag,icase,:,:,:),vars_mid(1,itag,icase,:,:,:),lev)
				ECTEI_tmp						= tofloat(result[2])
				vars_all(ivar,itag,icase,:,:)	= (/ECTEI_tmp/)
				delete(ECTEI_tmp)
			end if
			if(vars(ivar).eq."sgm_tota" .or. vars(ivar).eq."N1")then
				vars_all(ivar,itag,icase,:,:)	= vars_all@_FillValue
			end if
			if(vars(ivar).eq."netCRF")then
				SWCF						= ftmp1->SWCF(0,:,:)
				LWCF						= ftmp1->LWCF(0,:,:)
				netCRF						= (/SWCF+LWCF/)
				vars_all(ivar,itag,icase,:,:)	= (/netCRF/)
			end if

		end if
		printMinMax(vars_all(ivar,itag,icase,:,:),True)
	end do ; do ivar=

end do ; do icase=
end do ; do itag=
printVarSummary(vars_all)
;printMinMax(vars_all(2,:,0,:,:),True)
;printMinMax(vars_all(2,:,1,:,:),True)
printVarSummary(vars_mid)

if(Sherwood_2014)then

; define array to collect the index value
vars_shwd			= new((/dimsizes(sherwood_index),dimsizes(tag),dimsizes(case_nor)/),float)
vars_shwd!0			= "var"
vars_shwd&var		= fspan(0,dimsizes(sherwood_index)-1,dimsizes(sherwood_index))
vars_shwd!1		 	= "tag"
vars_shwd&tag		= fspan(0,dimsizes(tag)-1,dimsizes(tag))
vars_shwd!2		 	= "case"
vars_shwd&case		= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))

vars_shwd_lev		= new((/dimsizes(tag),dimsizes(case_nor),nlev/),float)
vars_shwd_lev(:,:,0)	= vars_shwd(0,:,:)

; calculate the index
do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1
	; get Sherwood index (DD, SS, LTMI)
	result					= cal_LTMI(vars_mid(2,itag,icase,:,:,:),vars_mid(0,itag,icase,:,:,:),vars_mid(1,itag,icase,:,:,:),vars_mid(3,itag,icase,:,:,:),OCNFRAC(:,:),lev)
	SS						= result[0]
	DD						= result[1]
	LTMI					= result[2]
	q_domega_dp_avg			= result[3]
	vars_shwd(0,itag,icase)	= SS
	vars_shwd(1,itag,icase)	= DD
	vars_shwd(2,itag,icase) = LTMI
	vars_shwd_lev(itag,icase,:)	= q_domega_dp_avg
end do ; do icase=
end do ; do itag=
printVarSummary(vars_shwd)
printVarSummary(vars_shwd_lev)
printMinMax(vars_shwd_lev,True)


; print index on the screen
do itag=0,dimsizes(tag)-1
print(tag(itag))
print("				SS		DD		LTMI")
print(case_nor+"		"+vars_shwd(0,itag,:)+"		"+vars_shwd(1,itag,:)+"	"+vars_shwd(2,itag,:))
end do ; do itag=

; plot them in one figure
wks						= gsn_open_wks(plot_type,"./figure/Sherwood_2014_index")
res						= True
res@gsnFrame			= False
res@gsnDraw				= False

;res@vpWidthF			= 0.7
;res@vpHeightF			= 0.4

res@trYMinF				= 0.
res@trXMinF				= -0.5
res@trXMaxF				= 1.5

res@tiYAxisString		= "Sherwood_2014_mixing_index"

; --- add grid lines
;res@tmXMajorGrid		= True
;res@tmYMajorGrid		= True
;res@tmXMajorGridLineDashPattern = 14.
;res@tmYMajorGridLineDashPattern = 14.

res@tmXBMode			= "Explicit"
res@tmXBValues			= (/0,1/)
res@tmXBLabels			= (/"Park-RH","Gauss-PDF"/)

res@xyMarkLineMode		= "Markers" ; MarkLines
res@xyMarkers		 	= (/16,16,16/)
res@xyMarkerColors		= (/"blue","red","black"/)
;res@xyLineThicknessF	= 2.0
;res@xyDashPattern		= 1.

res1					= res
res1@xyMarkerSizeF		= 0.01

res2					= res
res2@xyMarkerSizeF		= 0.015

plot_0					= gsn_csm_xy(wks,vars_shwd&case-0.1,vars_shwd(:,0,:),res1)
plot_1					= gsn_csm_xy(wks,vars_shwd&case+0.1,vars_shwd(:,1,:),res2)
overlay(plot_0,plot_1)

;----------------------- Add text
tres					= True
tres@txFontHeightF		= .015
; one learning here:
; after ndtooned function, the array "tmp" would be organized from right to left:
; that is,
;case1_nor_idx0, case2_nor_idx0, case1_4K_idx0,case1_4K_idx0
;case1_nor_idx1, case2_nor_idx1, case1_4K_idx1,case1_4K_idx1
;case1_nor_idx2, case2_nor_idx2, case1_4K_idx2,case1_4K_idx2

dum						= new(dimsizes(ndtooned(vars_shwd(:,:,:))),graphic)
tmp						= ndtooned(vars_shwd(:,:,:))
tmp_xaxis				= (/-0.2,0.8,0.2,1.2,\
							-0.2,0.8,0.2,1.2,\
							-0.2,0.8,0.2,1.2/)

do i = 0,dimsizes(ndtooned(vars_shwd(:,:,:)))-1
	dum(i)				= gsn_add_text(wks,plot_0,sprintf("%5.2f",tmp(i)),tmp_xaxis(i),tmp(i),tres)
end do ; do i=

; ------------------------Add legends
legends                         = sherwood_index

lgres							= True
lgres@lgItemType				= res@xyMarkLineMode
lgres@lgMarkerIndexes			= res@xyMarkers
lgres@lgMarkerColors            = res@xyMarkerColors

lgres@lgLabelFontHeightF        = 0.08
lgres@vpWidthF                  = 0.15
lgres@vpHeightF                 = 0.1
;lgres@lgPerimOn                 = False

lbid1                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres)

amres                           = True
amres@amJust					= "BottomRight"
amres@amParallelPosF            = 0.5
amres@amOrthogonalPosF          = 0.5

annoid1                         = gsn_add_annotation(plot_1,lbid1,amres)

draw(plot_0)
frame(wks)

delete([/res, res1, res2, tmp, legends, lgres, amres, lbid1, annoid1/])
; ------------------------------------------------------------------------
; ------------------------plot vertical moisture transport q_domega_dp_avg
res								= True
res@vpWidthF					= 0.7
res@vpHeightF					= 0.6
res@gsnFrame					= False
res@gsnDraw						= False
res@trYReverse					= True
res@xyLineColors				= (/"blue", "purple", "blue", "purple"/)
;res@xyLineColors				= (/"blue", "magenta", "blue", "magenta"/)
res@xyLineThicknessF			= 3.0
res@xyDashPatterns				= (/0, 0, 1, 1/)
res@tiYAxisString				= "Pressure (hPa)"
res@tiXAxisString				= "Large-scale source (g/kg/day)"
res@gsnXRefLine					= 0.
res@gsnXRefLineColor			= "grey"
res@gsnXRefLineDashPattern		= 14.
res@gsnXRefLineThicknessF		= 3.0
res@trYMinF						= 200.

tmp								= reshape(vars_shwd_lev,(/dimsizes(tag)*dimsizes(case_nor),nlev/))
plot							= gsn_csm_xy(wks, tmp,  vars_shwd_lev&lev, res)

; ------------------------Add legends
legends                         = (/"Park-RH", "Gauss-PDF", "Park-RH, +4K", "Gauss-PDF, +4K"/)

lgres							= True
lgres@lgLineColors				= res@xyLineColors
lgres@lgLineThicknessF			= res@xyLineThicknessF
lgres@lgDashIndexes				= res@xyDashPatterns

lgres@lgLabelFontHeightF        = 0.20
lgres@vpWidthF                  = 0.20
lgres@vpHeightF                 = 0.15
lgres@lgPerimOn                 = False

lbid1                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres)

amres                           = True
amres@amJust					= "TopLeft"
amres@amParallelPosF            = -0.3
amres@amOrthogonalPosF          = -0.3

annoid1                         = gsn_add_annotation(plot,lbid1,amres)

draw(plot)
frame(wks)

exit
end if ; Sherwood_2014

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
;regimes_name		= (/"SH-high", "SH-mid", "NH-mid", "NH-high", "deep_lnd", "Shcu", "StCu"/)

nregimes			= dimsizes(regimes_name)
regimes				= new((/nregimes,dimsizes(vars),dimsizes(tag),dimsizes(case_nor)/),float) 
data2regimes		= new((/nregimes,dimsizes(vars),dimsizes(tag),dimsizes(case_nor),nlat,nlon/),float) 
data2regimes(0,:,:,:,:,:)		= vars_all

; this array is used to store each regime together, so the dimension (nregimes) is eliminated. 
data2regimes_mask					= new((/dimsizes(vars),dimsizes(tag),dimsizes(case_nor),nlat,nlon/),float) 
data2regimes_mask(:,:,:,:,:)		= vars_all

; this array is used to store each regime in each bin, with another added dimension (nregimes,....)
data2regimes_mask_s					= new((/nregimes,dimsizes(vars),dimsizes(tag),dimsizes(case_nor),nlat,nlon/),float) 
data2regimes_mask_s(0,:,:,:,:,:)		= vars_all



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
		tmp2			= where(lat2d.gt.latSs(ilat) .and. lat2d.lt.latEs(ilat) .and. lon2d.gt.lonSs(ilat) .and. lon2d.lt.lonEs(ilat),regime_mask_val(ilat),data2regimes_mask(ivar,itag,icase,:,:))
		tmp3			= where(lat2d.gt.latSs(ilat) .and. lat2d.lt.latEs(ilat) .and. lon2d.gt.lonSs(ilat) .and. lon2d.lt.lonEs(ilat),0,1)
		data2regimes(ilat,ivar,itag,icase,:,:)		= (/tmp1/)
		data2regimes_mask(ivar,itag,icase,:,:)		= (/tmp2/)
		data2regimes_mask_s(ilat,ivar,itag,icase,:,:)		= (/tmp3/)

;		regimes(ilat,:,:,:)		= wgt_areaave_Wrap(vars_all(:,:,:,{latSs(ilat):latEs(ilat)},{lonSs(ilat):lonEs(ilat)}),clat({latSs(ilat):latEs(ilat)}),1,0)
		regimes(ilat,ivar,itag,icase)		= wgt_areaave_Wrap(tmp1,clat,1,0)
		delete([/tmp1,tmp2,tmp3/])
	end do ; do ilat=
	
	; tropical deep convection and land regions
	crit1					= (OCNFRAC.lt.0.2).and.(lat2d.gt.-30.and.lat2d.lt.30.)
	crit2					= (omega500.lt.0).and.(lat2d.gt.-30.and.lat2d.lt.30.)
	tmp1					= where(crit1 .or. crit2, vars_all(ivar,itag,icase,:,:),vars_all@_FillValue)
	tmp2					= where(crit1 .or. crit2, regime_mask_val(4),data2regimes_mask(ivar,itag,icase,:,:))
	tmp3					= where(crit1 .or. crit2, 0,1)

	data2regimes(4,ivar,itag,icase,:,:)		= (/tmp1/)
	data2regimes_mask(ivar,itag,icase,:,:)		= (/tmp2/)
	data2regimes_mask_s(4,ivar,itag,icase,:,:)		= (/tmp3/)
	regimes(4,ivar,itag,icase)		= wgt_areaave_Wrap(tmp1,clat,1,0)
	delete([/crit1, crit2, tmp1, tmp2, tmp3/])
	
	; tropical trade cumulus
	crit1					= (OCNFRAC.gt.0.8) .and. (omega500.gt.0).and.(lat2d.gt.-30.and.lat2d.lt.30.).and.(LTS(itag,icase,:,:).lt.LTS_criteria)
	tmp1					= where(crit1, vars_all(ivar,itag,icase,:,:),vars_all@_FillValue)
	tmp2					= where(crit1, regime_mask_val(5),data2regimes_mask(ivar,itag,icase,:,:))
	tmp3					= where(crit1, 0,1)

	data2regimes(5,ivar,itag,icase,:,:)		= (/tmp1/)
	data2regimes_mask(ivar,itag,icase,:,:)		= (/tmp2/)
	data2regimes_mask_s(5,ivar,itag,icase,:,:)		= (/tmp3/)

	regimes(5,ivar,itag,icase)		= wgt_areaave_Wrap(tmp1,clat,1,0)
	delete([/crit1, tmp1, tmp2, tmp3/])
	
	; tropical stratocumulus
	crit1					= (OCNFRAC.gt.0.8) .and. (omega500.gt.0).and.(lat2d.gt.-30.and.lat2d.lt.30.).and.(LTS(itag,icase,:,:).gt.LTS_criteria)
	tmp1					= where(crit1, vars_all(ivar,itag,icase,:,:),vars_all@_FillValue)
	tmp2					= where(crit1, regime_mask_val(6),data2regimes_mask(ivar,itag,icase,:,:))
	tmp3					= where(crit1, 0,1)

	data2regimes(6,ivar,itag,icase,:,:)		= (/tmp1/)
	data2regimes_mask(ivar,itag,icase,:,:)		= (/tmp2/)
	data2regimes_mask_s(6,ivar,itag,icase,:,:)		= (/tmp3/)

	regimes(6,ivar,itag,icase)		= wgt_areaave_Wrap(tmp1,clat,1,0)
	delete([/crit1, tmp1, tmp2, tmp3/])
end do ; do itag=
end do ; do icase=
end do ; do ivar=

printVarSummary(data2regimes)

if(regime_mask_plot)then

; ======== the following part is used to verify whether the regime-sorting is right.
; ======== ---- yes, plots indicate that all regimes are right.

wks			= gsn_open_wks(plot_type,"./figure/distribution-of-7-regimes-crit_LTS_"+LTS_criteria)

res							= True

res@gsnFrame				= False
res@gsnDraw					= False

res@cnFillOn				= True
res@cnLinesOn				= False
res@cnLineLabelsOn			= False
res@cnFillMode       		= "RasterFill"
res@cnLevelSelectionMode	= "ExplicitLevels"
res@cnLevels         		= regime_mask_val(1:)
res@cnFillColors			= (/"lightcyan2","lightgoldenrod2","lightyellow2","lightgrey","lightskyblue2","lightseagreen","lightpink2"/)

res@lbLabelPosition  		= "Center"           ; label position
res@lbLabelAlignment 		= "BoxCenters"       ; label orientation
;res@lbLabelStrings   		= ""+ res@cnLevels   ; trick 
res@lbLabelStrings   		= ""+ (/ 0,1,2,3,4,5,6/)

res@pmLabelBarHeightF		=  0.100
res@pmLabelBarWidthF 		=  0.80

res@mpCenterLonF			= 180

res@gsnLeftString			= "Selected regimes"

res1						= res
delete([/res1@cnLevels, res1@cnFillColors, res1@lbLabelStrings/])
;--------------------------- These settings are for regimes_mask_s -----------
; each regime is colored in single plot rather than all in one plot.
res1@cnLevels				= (/1/)
res1@cnFillColors			= (/"lightcyan2","white"/)
res1@lbLabelStrings   		= ""+ (/ 0,1/)

plot						= new((/dimsizes(tag),dimsizes(case_nor)/),graphic)
do ivar=0,0;dimsizes(vars)-1
do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1
res@gsnRightString = tag(itag)+"_"+case_name(icase)
plot(itag,icase)		= gsn_csm_contour_map_ce(wks,data2regimes_mask(ivar,itag,icase,:,:),res)

;do ir=0,dimsizes(regimes_name)-1
;plot1		= gsn_csm_contour_map_ce(wks,data2regimes_mask_s(ir,ivar,itag,icase,:,:),res1)
;end do 

end do
end do
end do

gsn_panel(wks,ndtooned(plot),(/2,2/),False)


end if ; regime_mask_plot
; ======================================================

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
; 2018.10.04 difference of tagdiff between PDF and RH
regimes_diff2				= (/regimes_tagdiff(:,:,1)-regimes_tagdiff(:,:,0)/)

regimes_areawgt			 	= (/data2regimes_1_percent*regimes/)

; use the same area for cntl and p4K simulation: both comes from cntl simulation
regimes_areawgt_1			 = (/regimes/)
regimes_areawgt_1(:,:,0,:)	= (/data2regimes_1_percent(:,:,0,:)*regimes(:,:,0,:)/)
regimes_areawgt_1(:,:,1,:)	= (/data2regimes_1_percent(:,:,0,:)*regimes(:,:,1,:)/)

regimes_tagdiff_areawgt		= (/regimes_areawgt(:,:,1,:)-regimes_areawgt(:,:,0,:)/) ; (regime, var, case_nor)
regimes_tagdiff_areawgt_1	= (/regimes_areawgt_1(:,:,1,:)-regimes_areawgt_1(:,:,0,:)/)
; 2018.10.04
regimes_diff2_areawgt		= (/regimes_tagdiff_areawgt(:,:,1)-regimes_tagdiff_areawgt(:,:,0)/)
regimes_diff2_areawgt_1		= (/regimes_tagdiff_areawgt_1(:,:,1)-regimes_tagdiff_areawgt_1(:,:,0)/)

Vars_out                    = regimes_name
Vars_out_fix                = Vars_out
do ivar = 0,dimsizes(Vars_out)-1
    lens                    = strlen(Vars_out(ivar))
    Vars_out_fix(ivar)      = str_concat((/conform_dims(20-lens, " ",-1)/))
    Vars_out_fix(ivar)      = str_concat((/Vars_out(ivar),Vars_out_fix(ivar)/))
    delete(lens)
end do ; do ivar=

do ivar = 0,dimsizes(vars)-1
print("For "+vars(ivar))
print("regimes	area-weight_RH	area-weight_PDF		tagdiff_RH	tagdiff_PDF tagdiff_PDF-RH	tagdiff_areawgt_RH	tagdiff_areawgt_PDF	tagdiff_areawgt_PDF-RH")
print(Vars_out_fix+" "+sprintf("%7.3e",data2regimes_1_percent(:,ivar,0,0))+"	"+sprintf("%7.3e",data2regimes_1_percent(:,ivar,0,1))+"	"+sprintf("%7.3e",regimes_tagdiff(:,ivar,0))+"	"+sprintf("%7.3e",regimes_tagdiff(:,ivar,1))+"	"+sprintf("%7.3e",regimes_diff2(:,ivar))+"	"+sprintf("%7.3e",regimes_tagdiff_areawgt(:,ivar,0))+"	"+sprintf("%7.3e",regimes_tagdiff_areawgt(:,ivar,1))+"	"+sprintf("%7.3e",regimes_diff2_areawgt(:,ivar)))

; ---- check whether adding the weigting could equal to the global mean value got from direct calculation. -- the answer is yes.
print("verify the correctness of weighting:		"+dim_sum_n(regimes_tagdiff_areawgt(:,ivar,0),0)+"	"+dim_sum_n(regimes_tagdiff_areawgt(:,ivar,1),0))
; -----
end do ; do ivar=

nregimes_used				= dimsizes(regimes_name(select_idx:))

wks					 	= gsn_open_wks(plot_type, "figure/Global-regime-contribution-p4K-"+vars(dimsizes(vars)-1)+"-crit_LTS_"+LTS_criteria)

sres 						= True
sres@gsnDraw               = False
sres@gsnFrame              = False

sres@vpWidthF              = 0.7
sres@vpHeightF             = 0.5

sres@trXMinF               = 0.
sres@trXMaxF               = nregimes_used+1
;sres@trYMinF               = -1.
;sres@trYMaxF               = 1.

barwidth					= 0.20
sres@gsnXYBarChart         = True
sres@gsnXYBarChartBarWidth = barwidth           ; change bar widths

sres@gsnYRefLine 			= 0

sres@tmXBMode              = "Explicit"         ; explicit labels
sres@tmXBValues            = ispan(1,nregimes_used,1)
sres@tmXBLabels            = regimes_name(select_idx:)

sres@tmXBLabelAngleF		= 45.

sres@tmXBLabelFontHeightF  = 0.0205
sres@tmXTLabelFontHeightF  = 0.0205
sres@tmYLLabelFontHeightF  = 0.0225

sres@gsnRightString        = ""
;sres@tiYAxisString         = "W/m2"

sres2						= sres
sres3						= sres

plots 						= new((/dimsizes(vars),dimsizes(case_nor)+1/),graphic)
dum0						= new((/dimsizes(vars),dimsizes(case_nor)+1/),graphic)
plots_wgt					= new((/dimsizes(vars),dimsizes(case_nor)+1/),graphic)
dum1						= new((/dimsizes(vars),dimsizes(case_nor)+1/),graphic)
plots_wgt_1					= new((/dimsizes(vars),dimsizes(case_nor)+1/),graphic)
dum2						= new((/dimsizes(vars),dimsizes(case_nor)+1/),graphic)

;colors 						= (/"lightpink",   "lightblue"/)
colors 						= (/"grey",   "black",	"blue"/)

; The X values for each plot is different, ensuring that the bars don't overlap.
do ivar=0,dimsizes(vars)-1
	sres@tiYAxisString		= vars(ivar)+"("+vars_units(ivar)+")"
	sres@trYMinF			= -1.*max((/max(regimes_tagdiff(select_idx:,ivar,:)),abs(min(regimes_tagdiff(select_idx:,ivar,:))),max(scale_factor*regimes_diff2(select_idx,ivar)),abs(min(scale_factor*regimes_diff2(select_idx,ivar)))/))*1.2
	sres@trYMaxF			= 1.*max((/max(regimes_tagdiff(select_idx:,ivar,:)),abs(min(regimes_tagdiff(select_idx:,ivar,:))),max(scale_factor*regimes_diff2(select_idx,ivar)),abs(min(scale_factor*regimes_diff2(select_idx,ivar)))/))*1.2
	sres@gsnLeftString		= "un-weighted"

	sres2@tiYAxisString		= vars(ivar)+"("+vars_units(ivar)+")"
	sres2@trYMinF			= -1.*max((/max(regimes_tagdiff_areawgt(select_idx:,ivar,:)),abs(min(regimes_tagdiff_areawgt(select_idx:,ivar,:))),max(scale_factor*regimes_diff2_areawgt(select_idx:,ivar)),abs(min(scale_factor*regimes_diff2_areawgt(select_idx:,ivar)))/))*1.2
	sres2@trYMaxF			= 1.*max((/max(regimes_tagdiff_areawgt(select_idx:,ivar,:)),abs(min(regimes_tagdiff_areawgt(select_idx:,ivar,:))),max(scale_factor*regimes_diff2_areawgt(select_idx:,ivar)),abs(min(scale_factor*regimes_diff2_areawgt(select_idx:,ivar)))/))*1.2

	sres3@tiYAxisString		= vars(ivar)+"("+vars_units(ivar)+")"
	sres3@trYMinF			= -1.*max((/max(regimes_tagdiff_areawgt_1(select_idx:,ivar,:)),abs(min(regimes_tagdiff_areawgt_1(select_idx:,ivar,:))),max(scale_factor*regimes_diff2_areawgt_1(select_idx:,ivar)),abs(min(scale_factor*regimes_diff2_areawgt_1(select_idx:,ivar)))/))*1.2
	sres3@trYMaxF			= 1.*max((/max(regimes_tagdiff_areawgt_1(select_idx:,ivar,:)),abs(min(regimes_tagdiff_areawgt_1(select_idx:,ivar,:))),max(scale_factor*regimes_diff2_areawgt_1(select_idx:,ivar)),abs(min(scale_factor*regimes_diff2_areawgt_1(select_idx:,ivar)))/))*1.2

	sres2@gsnLeftString		= "area-weighted"
	sres3@gsnLeftString		= "area-weighted-but-same"

  lnres                   = True
	lnres@gsLineThicknessF  = 3.0 

	; For un-weighted
	sres@gsnXYBarChartColors = colors(0)
	sres@xyLineColors=(/colors(0),colors(0)/) 	; this is used to change the outline to be the color of the fill. I still dont know why it needs two colors.
	plots(ivar,0) = gsn_csm_xy(wks,fspan(1-2*barwidth/2.,nregimes_used-2*barwidth/2,nregimes_used),regimes_tagdiff(select_idx:,ivar,0),sres)

	sres@gsnLeftString		 = ""
	sres@gsnXYBarChartColors = colors(1)
	sres@xyLineColors=(/colors(1),colors(1)/)
	plots(ivar,1) = gsn_csm_xy(wks,fspan(1+0*barwidth/2,nregimes_used+0*barwidth/2,nregimes_used),regimes_tagdiff(select_idx:,ivar,1),sres)

	sres@gsnLeftString		 = ""
	sres@gsnXYBarChartColors = colors(2)
	sres@xyLineColors=(/colors(2),colors(2)/)
	plots(ivar,2) = gsn_csm_xy(wks,fspan(1+2*barwidth/2,nregimes_used+2*barwidth/2,nregimes_used),regimes_diff2(select_idx:,ivar)*scale_factor,sres)

	; For area-weighted
	sres2@gsnXYBarChartColors = colors(0)
	sres2@xyLineColors = (/colors(0),colors(0)/)
	plots_wgt(ivar,0) = gsn_csm_xy(wks,fspan(1-2*barwidth/2.,nregimes_used-2*barwidth/2,nregimes_used),regimes_tagdiff_areawgt(select_idx:,ivar,0),sres2)

	sres2@gsnLeftString			= ""
	sres2@gsnXYBarChartColors = colors(1)
	sres2@xyLineColors	= (/colors(1),colors(1)/)
	plots_wgt(ivar,1) = gsn_csm_xy(wks,fspan(1+0*barwidth/2,nregimes_used+0*barwidth/2,nregimes_used),regimes_tagdiff_areawgt(select_idx:,ivar,1),sres2)

	sres2@gsnLeftString			= ""
	sres2@gsnXYBarChartColors = colors(2)
	sres2@xyLineColors	= (/colors(2),colors(2)/)
	plots_wgt(ivar,2) = gsn_csm_xy(wks,fspan(1+2*barwidth/2,nregimes_used+2*barwidth/2,nregimes_used),regimes_diff2_areawgt(select_idx:,ivar)*scale_factor,sres2)

	; For area-weighted but the same area
	sres3@gsnXYBarChartColors = colors(0)
	sres3@xyLineColors = (/colors(0),colors(0)/)
	plots_wgt_1(ivar,0) = gsn_csm_xy(wks,fspan(1-2*barwidth/2.,nregimes_used-2*barwidth/2,nregimes_used),regimes_tagdiff_areawgt_1(select_idx:,ivar,0),sres3)

	sres3@gsnLeftString			= ""
	sres3@gsnXYBarChartColors = colors(1)
	sres3@xyLineColors	= (/colors(1),colors(1)/)
	plots_wgt_1(ivar,1) = gsn_csm_xy(wks,fspan(1+0*barwidth/2,nregimes_used+0*barwidth/2,nregimes_used),regimes_tagdiff_areawgt_1(select_idx:,ivar,1),sres3)

	sres3@gsnLeftString			= ""
	sres3@gsnXYBarChartColors = colors(2)
	sres3@xyLineColors	= (/colors(2),colors(2)/)
	plots_wgt_1(ivar,2) = gsn_csm_xy(wks,fspan(1+2*barwidth/2,nregimes_used+2*barwidth/2,nregimes_used),regimes_diff2_areawgt_1(select_idx:,ivar)*scale_factor,sres3)

	lnres@gsLineColor       = colors(1)
	dum0(ivar,1)			= gsn_add_polyline(wks,plots(ivar,1),fspan(1+0*barwidth/2,nregimes_used+0*barwidth/2,nregimes_used),regimes_tagdiff(select_idx:,ivar,1),lnres)
	dum1(ivar,1)			= gsn_add_polyline(wks,plots_wgt(ivar,1),fspan(1+0*barwidth/2,nregimes_used+0*barwidth/2,nregimes_used),regimes_tagdiff_areawgt(select_idx:,ivar,1),lnres)
	dum2(ivar,1)			= gsn_add_polyline(wks,plots_wgt_1(ivar,1),fspan(1+0*barwidth/2,nregimes_used+0*barwidth/2,nregimes_used),regimes_tagdiff_areawgt_1(select_idx:,ivar,1),lnres)

	if(vars(ivar).ne."sgm_tota" .and.  vars(ivar).ne."N1")then
	lnres@gsLineColor       = colors(0)
	dum0(ivar,0)			= gsn_add_polyline(wks,plots(ivar,0),fspan(1-2*barwidth/2.,nregimes_used-2*barwidth/2,nregimes_used),regimes_tagdiff(select_idx:,ivar,0),lnres)
	dum1(ivar,0)			= gsn_add_polyline(wks,plots_wgt(ivar,0),fspan(1-2*barwidth/2.,nregimes_used-2*barwidth/2,nregimes_used),regimes_tagdiff_areawgt(select_idx:,ivar,0),lnres)
	dum2(ivar,0)			= gsn_add_polyline(wks,plots_wgt_1(ivar,0),fspan(1-2*barwidth/2.,nregimes_used-2*barwidth/2,nregimes_used),regimes_tagdiff_areawgt_1(select_idx:,ivar,0),lnres)

	lnres@gsLineColor		= colors(2)
	dum0(ivar,2) 			= gsn_add_polyline(wks,plots(ivar,2),fspan(1+2*barwidth/2,nregimes_used+2*barwidth/2,nregimes_used),regimes_diff2(select_idx:,ivar)*scale_factor,lnres)
	dum1(ivar,2) 			= gsn_add_polyline(wks,plots_wgt(ivar,2),fspan(1+2*barwidth/2,nregimes_used+2*barwidth/2,nregimes_used),regimes_diff2_areawgt(select_idx:,ivar)*scale_factor,lnres)
	dum2(ivar,2) 			= gsn_add_polyline(wks,plots_wgt_1(ivar,2),fspan(1+2*barwidth/2,nregimes_used+2*barwidth/2,nregimes_used),regimes_diff2_areawgt_1(select_idx:,ivar)*scale_factor,lnres)
	end if


; Overlay the last three plots on the first one.
       do j=1,dimsizes(case_nor)-1+1
         overlay(plots(ivar,0),plots(ivar,j))
		 overlay(plots_wgt(ivar,0),plots_wgt(ivar,j))
		 overlay(plots_wgt_1(ivar,0),plots_wgt_1(ivar,j))
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

labels = (/"Park-RH", "Gauss-PDF","PDF-RH"/)
xpos   = (/0.15,  0.25, 0.37/)

if(dimsizes(plots(:,0))/2.eq.3)then
ypos	 = 0.20
end if
if(dimsizes(plots(:,0))/2.eq.2)then
ypos	= 0.75
end if

; Panel the plots.
;First plot
;gsn_panel(wks,plots(:,0),(/2,2/),pres)
gsn_panel(wks,plots(:,0),(/dimsizes(plots(:,0))/2,2/),pres)
;drawNDCGrid(wks)

do i=0,dimsizes(case_nor)-1+1
  lbres@lbFillColors = colors(i)
  gsn_labelbar_ndc(wks,1,labels(i),xpos(i),ypos,lbres)
end do	
frame(wks)   ; Advance the frame.

; Second plot
;gsn_panel(wks,plots_wgt(:,0),(/2,2/),pres)
gsn_panel(wks,plots_wgt(:,0),(/dimsizes(plots_wgt(:,0))/2,2/),pres)

do i=0,dimsizes(case_nor)-1+1
  lbres@lbFillColors = colors(i)
  gsn_labelbar_ndc(wks,1,labels(i),xpos(i),ypos,lbres)
end do	
frame(wks)

; Third plot
gsn_panel(wks,plots_wgt_1(:,0),(/dimsizes(plots_wgt_1(:,0))/2,2/),pres)

do i=0,dimsizes(case_nor)-1+1
  lbres@lbFillColors = colors(i)
  gsn_labelbar_ndc(wks,1,labels(i),xpos(i),ypos,lbres)
end do	

frame(wks)   ; Advance the frame.
exit

end if ; regime_contribution

;************************create plot***************
if(lat_lon_plot)then

do ivar=0,dimsizes(vars)-1

;wks 					= gsn_open_wks(plot_type,"figure/lat-lon-"+vars(ivar)+"-HammerProj")
wks 					= gsn_open_wks(plot_type,"figure/lat-lon-"+vars(ivar)+"-p4K-"+vars(dimsizes(vars)-1))

plot 					= new(9,graphic)
plot1 					= new(9,graphic)

;gsn_define_colormap(wks,"BlueDarkRed18")

cmap1 					= read_colormap_file("WhiteBlueGreenYellowRed")
;cmap2 					= read_colormap_file("BlueRed")
;cmap2					= read_colormap_file("BlueWhiteOrangeRed")
;cmap2				 	= read_colormap_file("NCV_roullet")
cmap2				 	= read_colormap_file("NCV_blue_red")

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
res2@cnLevels                 	= vars_dlevs_1(ivar,:)

res3							= res2
res3@cnLevels					= vars_dlevs_2(ivar,:)

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
	if(it.eq.6 .or. it.eq.7)then
		res3@gsnLeftString			= strings4K(it)
		plot(it)					= gsn_csm_contour_map(wks,smth9_Wrap(var4K(it,:,:),p,q,False),res3)
	else
		res2@gsnLeftString			= strings4K(it)
;		if(it.eq.8)then
;			res2@gsnCenterString	= sprintf("%5.2f",ra(ivar))
;		else
;			res2@gsnCenterString	= ""
;		end if
		plot(it)				= gsn_csm_contour_map(wks,smth9_Wrap(var4K(it,:,:),p,q,False),res2)
	end if
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

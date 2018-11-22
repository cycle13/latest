
begin 

fdir				= "/home/share3/lyl/work3/qinyi/mid-data/cfmip/"

case 				= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_MG15_amip-4xCO2",\
		 		    "FAMIPC5_f09f09_mac2_amip","FAMIPC5_f09f09_mac2_amip-p4K_debug","FAMIPC5_f09f09_mac2_amip-4xCO2"/)

append				= "-1980.1984_regrid.nc"

vars				= (/"T","Q"/)

nlat				= 180
nlon				= 360
nlev			 	= 30

calculation			= False 

if(calculation)then

f1				= addfile(fdir+case(0)+append,"r")

lev				= f1->lev

var_def				= new((/3,dimsizes(vars),nlev,nlat,nlon/),float)
var_def!0		 	= "exps"
var_def&exps			= (/0,1,2/)
var_def!1			= "vars"
var_def&vars			= fspan(0,dimsizes(vars)-1,dimsizes(vars))

var_new				= var_def

do icase=0,dimsizes(case)-1
	ftmp1						= addfile(fdir+case(icase)+append,"r")
	do ivar=0,dimsizes(vars)-1
		if(icase.le.2)then
			var_def(icase,ivar,:,:,:)	 = ftmp1->$vars(ivar)$
		else
			var_new(icase-3,ivar,:,:,:)	= ftmp1->$vars(ivar)$
		end if
	end do ; do ivar=
end do ; do icase=

printVarSummary(var_def)

; Start to calculate EIS...
; 1. define constant parameters
RV 			= 461.5 ; J/kg/K
RD 			= 287.04 ; J/kg/K
EPS 			= RD/RV
CPD 			= 1005.7 ; J/kg/K
CPV 			= 1870.0 ; J/kg/K
P0 			= 1000.0 ; hPa
LV 			= 2501000 ; J/kg
gravit 			= 9.8 ; m/s2

temp_def		= var_def(:,0,:,:,:)
;temp_def		= where(var_def(:,0,:,:,:).gt.270,var_def(:,0,:,:,:),var_def@_FillValue) ; K
temp_new		= var_new(:,0,:,:,:)
;temp_new		= where(var_new(:,0,:,:,:).gt.270,var_new(:,0,:,:,:),var_new@_FillValue)
qv_def			= var_def(:,1,:,:,:) ; kg/kg
qv_new			= var_new(:,1,:,:,:)

; 1. calculate potential temperature ------------------------------------
P 			= conform(var_def(:,0,:,:,:),lev,1)	; hPa
theta_def 		= temp_def*(P0/P)^0.286 ; K
theta_new 		= temp_new*(P0/P)^0.286 
copy_VarCoords(var_def(:,0,:,:,:),theta_def)
copy_VarCoords(var_new(:,0,:,:,:),theta_new)
printVarSummary(theta_def)
printMinMax(theta_def,True)

; 2. calculate moist-adiabatic potential temp. gradient Gamma_m ---------------------
; 2.1 convert spc humidty to mixing ratio first
mr_def			= qv_def/(1.-qv_def)
mr_new			= qv_new/(1.-qv_new)
copy_VarCoords(qv_def,mr_def)
copy_VarCoords(qv_new,mr_new)
printVarSummary(mr_def)
printMinMax(mr_def,True)

; 2.2 get saturation mixing ratio
; water vapor pressure
ev_def			 = mr_def*P/(EPS+mr_def) 
ev_new			 = mr_new*P/(EPS+mr_new) 
copy_VarCoords(mr_def,ev_def)
copy_VarCoords(mr_new,ev_new)
printVarSummary(ev_def)
; saturation water vapor pressure (hPa)
es_def			=6.112*exp(17.67*(temp_def-273.15)/(243.5+(temp_def-273.15)))
es_new			=6.112*exp(17.67*(temp_new-273.15)/(243.5+(temp_new-273.15)))
copy_VarCoords(temp_def,es_def)
copy_VarCoords(temp_new,es_new)
printVarSummary(es_def)
printMinMax(es_def,True)
; saturation mixing ratio
mrs_def			= EPS*es_def/(P-es_def)
mrs_new			= EPS*es_new/(P-es_new)
copy_VarCoords(es_def,mrs_def)
copy_VarCoords(es_new,mrs_new)
printVarSummary(mrs_def)
; saturation spc humidity
qs_def			= mrs_def/(1.+mrs_def)
qs_new			= mrs_new/(1.+mrs_new)
copy_VarCoords(mrs_def,qs_def)
copy_VarCoords(mrs_new,qs_new)
printVarSummary(qs_new)
printMinMax(qs_new,False)

; 2.3 get moist-adiabatic potential temperature gradient 
Gamma_def		= gravit/CPD*(1.0-(1.0+LV*qs_def/RD/temp_def)/(1.0+LV^2*qs_def/CPD/RV/temp_def^2))
Gamma_new		= gravit/CPD*(1.0-(1.0+LV*qs_new/RD/temp_new)/(1.0+LV^2*qs_new/CPD/RV/temp_new^2))
copy_VarCoords(qs_def,Gamma_def)
copy_VarCoords(qs_new,Gamma_new)
printVarSummary(Gamma_def)
printMinMax(Gamma_def,True)

; 3. calculate lifting condensation level: Bolton (1980) Eq.(21) ----------------------------
TC_def			= 2840.0/(3.5*log(temp_def(:,29,:,:))-log(ev_def(:,29,:,:))-4.805)+55.0
TC_new			= 2840.0/(3.5*log(temp_new(:,29,:,:))-log(ev_new(:,29,:,:))-4.805)+55.0
copy_VarCoords(temp_def(:,0,:,:),TC_def)
copy_VarCoords(temp_new(:,0,:,:),TC_new)
plcl_def		= P(:,29,:,:)*((TC_def/temp_def(:,29,:,:))^3.48)
plcl_new		= P(:,29,:,:)*((TC_new/temp_new(:,29,:,:))^3.48)
copy_VarCoords(temp_def(:,0,:,:),plcl_def)
copy_VarCoords(temp_new(:,0,:,:),plcl_new)

printMinMax(TC_def,True)
printMinMax(TC_new,True)
printMinMax(plcl_def,False)

; 4. get z700 and zlcl ---------------------------------------------------------------------
z700_def		= (RD*temp_def(:,29,:,:)/gravit)*log(P0/700.)
z700_new		= (RD*temp_new(:,29,:,:)/gravit)*log(P0/700.)
zlcl_def 		= (RD*temp_def(:,29,:,:)/gravit)*log(P0/plcl_def)
zlcl_new 		= (RD*temp_new(:,29,:,:)/gravit)*log(P0/plcl_new)

copy_VarCoords(temp_def(:,29,:,:),z700_def)
copy_VarCoords(temp_def(:,29,:,:),z700_new)
copy_VarCoords(temp_def(:,29,:,:),zlcl_def)
copy_VarCoords(temp_def(:,29,:,:),zlcl_new)

printMinMax(z700_def,True)
printMinMax(zlcl_new,True)

; 5. calculate estimated inversion strength (EIS) ------------------------------------------
LTS_def			= theta_def(:,{700},:,:)-theta_def(:,29,:,:)
Gamma_def_850		= Gamma_def(:,{850},:,:)
dz_def			= z700_def-zlcl_def
EIS_def			= where(dz_def.lt.0,LTS_def,LTS_def-Gamma_def(:,{850},:,:)*dz_def)
copy_VarCoords(temp_def(:,0,:,:),LTS_def)
copy_VarCoords(temp_def(:,0,:,:),EIS_def)

LTS_new			= theta_new(:,{700},:,:)-theta_new(:,29,:,:)
Gamma_new_850		= Gamma_new(:,{850},:,:)
dz_new			= z700_new-zlcl_new
EIS_new			= where(dz_new.lt.0,LTS_new,LTS_new-Gamma_new(:,{850},:,:)*dz_new)
copy_VarCoords(temp_new(:,0,:,:),LTS_new)
copy_VarCoords(temp_new(:,0,:,:),EIS_new)

printMinMax(EIS_def,False)
printMinMax(EIS_new,False)


system("/bin/rm -f ./tmp-data/EIS.nc")
ncdf 			= addfile("./tmp-data/EIS.nc","c")

ncdf->EIS_def		= EIS_def
ncdf->EIS_new		= EIS_new
ncdf->LTS_def		= LTS_def
ncdf->LTS_new		= LTS_new

end if ; calculation

infile			= addfile("./tmp-data/EIS.nc","r")

vardef			= new((/3,2,nlat,nlon/),double)
vardef!1		= "vars"
vardef&vars		= (/0,1/)
vardef(:,0,:,:)		= infile->EIS_def
vardef(:,1,:,:)		= infile->LTS_def

varnew			= vardef
varnew(:,0,:,:)		= infile->EIS_new
varnew(:,1,:,:)		= infile->LTS_new

printVarSummary(vardef)

; get the difference
vardef_diff			= new((/2,2,nlat,nlon/),double)
vardef_diff(0,:,:,:)		= vardef(0,:,:,:)
vardef_diff(0,:,:,:)		= vardef(1,:,:,:)-vardef(0,:,:,:) ; RH_4K-RH_def
vardef_diff(1,:,:,:)		= vardef(2,:,:,:)-vardef(0,:,:,:) ; RH_4xCO2-RH_def

varnew_diff			= vardef_diff
varnew_diff(0,:,:,:)		= varnew(1,:,:,:)-varnew(0,:,:,:) ; PDF_4K-PDF_def
varnew_diff(1,:,:,:)		= varnew(2,:,:,:)-varnew(0,:,:,:) ; PDF_4xCO2-PDF_def

var_newdef_diff			= new((/5,2,nlat,nlon/),double)
var_newdef_diff(0,:,:,:)	= varnew(0,:,:,:)
var_newdef_diff(0,:,:,:)	= varnew(0,:,:,:)-vardef(0,:,:,:) ; PDF-RH
var_newdef_diff(1,:,:,:)	= varnew(1,:,:,:)-vardef(1,:,:,:) ; PDF_4K-RH_4K
var_newdef_diff(2,:,:,:)	= varnew(2,:,:,:)-vardef(2,:,:,:) ; PDF_4xCO2-RH_4xCO2

var_newdef_diff(3,:,:,:) 	= var_newdef_diff(1,:,:,:)-var_newdef_diff(0,:,:,:)
var_newdef_diff(4,:,:,:) 	= var_newdef_diff(2,:,:,:)-var_newdef_diff(0,:,:,:)


printVarSummary(vardef_diff)
printVarSummary(var_newdef_diff)

; -----------------------------------------------------------------------
;create plot
; -----------------------------------------------------------------------
varout				= (/"EIS",	"LTS"/)
vars_units			= (/"K",	"K"/)
EIS_levs			= (/-3.,-2,-1,0,1,2,3,4,5,6,8,10,12/)
LTS_levs			= (/5.,6,7,8,9,10,12,14,16,18,20,22,24/)

EIS_dlevs_1			= (/-1.,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1/)
EIS_dlevs_2			= (/-1.,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1/)

LTS_dlevs_1			= (/-1.,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1/)
LTS_dlevs_2			= (/-3,-2.6,-2.2,-1.8,-1.4,-1.0,0,1.0,1.4,1.8,2.2,2.6,3.0/)


vars_levs			= (/EIS_levs,		LTS_levs/)
vars_dlevs_1			= (/EIS_dlevs_1,	LTS_dlevs_1/)
vars_dlevs_2			= (/EIS_dlevs_2,	LTS_dlevs_2/)


;do ivar=0,dimsizes(varout)-1
do ivar=0,1

wks 				= gsn_open_wks("pdf","./figure/"+varout(ivar))

cmap1 					= read_colormap_file("WhiteBlueGreenYellowRed")
;cmap2 					= read_colormap_file("BlueRed")
cmap2					= read_colormap_file("BlueDarkRed18")
;cmap2					= read_colormap_file("nrl_sirkes")
;cmap2					= read_colormap_file("hotcold_18lev")


res                      		= True
res@gsnDraw 				= False
res@gsnFrame 				= False
;res@mpProjection			= "Hammer"
;res@mpPerimOn              		= False	 

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

;res@gsnCenterString			= varout(ivar)
res@gsnRightString			= vars_units(ivar)

res@lbOrientation			= "vertical"
res@lbLabelFontHeightF			= fontheight

res@cnLevelSelectionMode        	= "ExplicitLevels"
res@cnLevels                    	= vars_levs(ivar,:)


; please start here..
; need to mask land lines for omega

;;;;;;;res2
res2 					= res

delete(res2@cnFillPalette)
res2@cnFillPalette 			= cmap2(:,:)

res2@cnLevelSelectionMode        	= "ExplicitLevels"
res2@cnLevels                    	= vars_dlevs_1(ivar,:)

res3					= res2
res3@cnLevels				= vars_dlevs_2(ivar,:)


var4K					= new((/9,nlat,nlon/),double)
var4K(0,:,:)				= vardef(0,ivar,:,:)
var4K					= (/vardef(0,ivar,:,:),		varnew(0,ivar,:,:),		var_newdef_diff(0,ivar,:,:),\
					    vardef(1,ivar,:,:),		varnew(1,ivar,:,:),		var_newdef_diff(1,ivar,:,:),\
					    vardef_diff(0,ivar,:,:),	varnew_diff(0,ivar,:,:),	var_newdef_diff(3,ivar,:,:)/)

var4xCO2				= new((/9,nlat,nlon/),double)
var4xCO2(0,:,:)				= vardef(0,ivar,:,:)
var4xCO2				= (/vardef(0,ivar,:,:),		varnew(0,ivar,:,:),		var_newdef_diff(0,ivar,:,:),\
					    vardef(2,ivar,:,:),		varnew(2,ivar,:,:),		var_newdef_diff(2,ivar,:,:),\
					    vardef_diff(1,ivar,:,:),	varnew_diff(1,ivar,:,:),	var_newdef_diff(4,ivar,:,:)/)


strings4K				= (/"RH",		"PDF",			"PDF-RH",\
					    "RH_4K",		"PDF_4K",		"PDF_4K-RH_4K",\
					    "RH_4K-RH",		"PDF_4K-PDF",		"PDF_RH_4K-PDF_RH_def"/)

strings4xCO2				= (/"RH",		"PDF",			"PDF-RH",\
					    "RH_4xCO2",		"PDF_4xCO2",		"PDF_4xCO2-RH_4xCO2",\
					    "RH_4xCO2-RH",	"PDF_4xCO2-PDF",	"PDF_RH_4xCO2-PDF_RH_def"/)


plot					= new(9,graphic)
plot1					= new(9,graphic)


do it=0,8
if(it.eq.0.or.it.eq.1.or.it.eq.3.or.it.eq.4)then
	res@gsnLeftString			= strings4K(it)
	plot(it)				= gsn_csm_contour_map(wks,var4K(it,:,:),res)
	res@gsnLeftString			= strings4xCO2(it)
	plot1(it)				= gsn_csm_contour_map(wks,var4xCO2(it,:,:),res)

else
	if((it.eq.6.or.it.eq.7).and.varout(ivar).eq."LTS")then
		res3@gsnLeftString			= strings4K(it)
		plot(it)				= gsn_csm_contour_map(wks,var4K(it,:,:),res3)
		res2@gsnLeftString			= strings4xCO2(it)
		plot1(it)				= gsn_csm_contour_map(wks,var4xCO2(it,:,:),res2)
	else
		res2@gsnLeftString			= strings4K(it)
		plot(it)				= gsn_csm_contour_map(wks,var4K(it,:,:),res2)
		res2@gsnLeftString			= strings4xCO2(it)
		plot1(it)				= gsn_csm_contour_map(wks,var4xCO2(it,:,:),res2)
	end if


end if
end do ; do it=

resP = True
resP@gsnPanelMainString = varout(ivar)

resP@gsnPanelXWhiteSpacePercent		= 2.5
resP@gsnPanelYWhiteSpacePercent		= 5.

;if(vars(ivar).eq."SWCF".or.vars(ivar).eq."CLDLOW")then
gsn_panel(wks,plot,(/3,3/),resP)
gsn_panel(wks,plot1,(/3,3/),resP)
;gsn_panel(wks,plot((/2,5,6,7/)),(/2,2/),resP)
gsn_panel(wks,(/plot(6),plot(7),plot1(6),plot1(7)/),(/2,2/),resP)

end do ; doivar=

end

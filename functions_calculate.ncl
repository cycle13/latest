
function ave_std(data[*][*][*]:numeric,nVar:integer,total_num:integer,nBin:integer)

begin

; data(nVar,total_num,nBin)
;print(data(0,200000:200100,49))
;exit


yy_data = new((/2,nVar,nBin/),float)

; ==========0: average mean; 1: standard deviation
yy_data(0,:,:) = dim_avg_n_Wrap(data,1)
yy_data(1,:,:) = dim_stddev_n_Wrap(data,1)

return(yy_data)

end

;=======================================================
function cal_theta(T[*][*][*][*]:numeric,P[*][*][*][*]:numeric,P0:numeric)

begin

theta = T*(P0/P)^0.286
copy_VarCoords(T,theta)

return(theta)

end

; 2018-10-20
function cal_eis(temp[*][*][*]:numeric,qv[*][*][*]:numeric,lev[*]:numeric)
begin

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

; 1. calculate potential temperature ------------------------------------
P 			= conform(temp,lev,0)	; hPa
theta			= temp*(P0/P)^0.286 ; K
copy_VarCoords(temp,theta)
;printVarSummary(theta)

; 2. calculate moist-adiabatic potential temp. gradient Gamma_m ---------------------
; 2.1 convert spc humidty to mixing ratio first
mr			= qv/(1.-qv)
copy_VarCoords(qv,mr)
;printVarSummary(mr)

; 2.2 get saturation mixing ratio
; water vapor pressure
ev			 = mr*P/(EPS+mr) 
copy_VarCoords(mr,ev)
;printVarSummary(ev)
; saturation water vapor pressure (hPa)
es			=6.112*exp(17.67*(temp-273.15)/(243.5+(temp-273.15)))
copy_VarCoords(temp,es)
;printVarSummary(es)
; saturation mixing ratio
mrs			= EPS*es/(P-es)
copy_VarCoords(es,mrs)
;printVarSummary(mrs)
; saturation spc humidity
qs			= mrs/(1.+mrs)
copy_VarCoords(mrs,qs)
;printVarSummary(qs)

; 2.3 get moist-adiabatic potential temperature gradient 
Gamma		= gravit/CPD*(1.0-(1.0+LV*qs/RD/temp)/(1.0+LV^2*qs/CPD/RV/temp^2))
copy_VarCoords(qs,Gamma)
;printVarSummary(Gamma)

; 3. calculate lifting condensation level: Bolton (1980) Eq.(21) ----------------------------
TC			= 2840.0/(3.5*log(temp(29,:,:))-log(ev(29,:,:))-4.805)+55.0
copy_VarCoords(temp(29,:,:),TC)

plcl		= P(29,:,:)*((TC/temp(29,:,:))^3.48)
copy_VarCoords(temp(29,:,:),plcl)

;printMinMax(TC,True)
;printMinMax(plcl,False)

; 4. get z700 and zlcl ---------------------------------------------------------------------
z700		= (RD*temp(29,:,:)/gravit)*log(P0/700.)
zlcl 		= (RD*temp(29,:,:)/gravit)*log(P0/plcl)

copy_VarCoords(temp(29,:,:),z700)
copy_VarCoords(temp(29,:,:),zlcl)

; 5. calculate estimated inversion strength (EIS) ------------------------------------------
LTS			= theta({700},:,:)-theta(29,:,:)
Gamma_850		= Gamma({850},:,:)
dz			= z700-zlcl
EIS			= where(dz.lt.0,LTS,LTS-Gamma({850},:,:)*dz)
copy_VarCoords(temp(0,:,:),LTS)
copy_VarCoords(temp(0,:,:),EIS)

return ([/LTS,EIS/])

end


; 2018-11-03
; refer to Kawai et al (2017): ECTEI = EIS - \beta*(L/cp)*(q_surf-q_700)
; \beta = (1-k)*C_qgap ; k = 0.70, C_qgap = 0.76, so beta = 0.23
function cal_ECTEI(temp[*][*][*]:numeric,qv[*][*][*]:numeric,lev[*]:numeric)
begin

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

; 1. calculate potential temperature ------------------------------------
P 			= conform(temp,lev,0)	; hPa
theta			= temp*(P0/P)^0.286 ; K
copy_VarCoords(temp,theta)
;printVarSummary(theta)

; 2. calculate moist-adiabatic potential temp. gradient Gamma_m ---------------------
; 2.1 convert spc humidty to mixing ratio first
mr			= qv/(1.-qv)
copy_VarCoords(qv,mr)
;printVarSummary(mr)

; 2.2 get saturation mixing ratio
; water vapor pressure
ev			 = mr*P/(EPS+mr) 
copy_VarCoords(mr,ev)
;printVarSummary(ev)
; saturation water vapor pressure (hPa)
es			=6.112*exp(17.67*(temp-273.15)/(243.5+(temp-273.15)))
copy_VarCoords(temp,es)
;printVarSummary(es)
; saturation mixing ratio
mrs			= EPS*es/(P-es)
copy_VarCoords(es,mrs)
;printVarSummary(mrs)
; saturation spc humidity
qs			= mrs/(1.+mrs)
copy_VarCoords(mrs,qs)
;printVarSummary(qs)

; 2.3 get moist-adiabatic potential temperature gradient 
Gamma		= gravit/CPD*(1.0-(1.0+LV*qs/RD/temp)/(1.0+LV^2*qs/CPD/RV/temp^2))
copy_VarCoords(qs,Gamma)
;printVarSummary(Gamma)

; 3. calculate lifting condensation level: Bolton (1980) Eq.(21) ----------------------------
TC			= 2840.0/(3.5*log(temp(29,:,:))-log(ev(29,:,:))-4.805)+55.0
copy_VarCoords(temp(29,:,:),TC)

plcl		= P(29,:,:)*((TC/temp(29,:,:))^3.48)
copy_VarCoords(temp(29,:,:),plcl)

;printMinMax(TC,True)
;printMinMax(plcl,False)

; 4. get z700 and zlcl ---------------------------------------------------------------------
z700		= (RD*temp(29,:,:)/gravit)*log(P0/700.)
zlcl 		= (RD*temp(29,:,:)/gravit)*log(P0/plcl)

copy_VarCoords(temp(29,:,:),z700)
copy_VarCoords(temp(29,:,:),zlcl)

; 5. calculate estimated inversion strength (EIS) ------------------------------------------
LTS			= tofloat(theta({700},:,:)-theta(29,:,:))
Gamma_850		= Gamma({850},:,:)
dz			= z700-zlcl
EIS			= tofloat(where(dz.lt.0,LTS,LTS-Gamma({850},:,:)*dz))
copy_VarCoords(temp(0,:,:),LTS)
copy_VarCoords(temp(0,:,:),EIS)

; 6. get q700 and q_surface and their difference (q_surface - q_700)
q700		=  qv({700},:,:)
qsfc		= qv(29,:,:)

qv_diff		= qsfc - q700
copy_VarCoords(q700, qv_diff)

; 7. define beta value
kk			= 0.70
C_qgap		= 0.76
beta		= (1.-kk)*C_qgap

; 8. get ECTEI
ECTEI		= EIS - beta*(LV/CPD)*qv_diff
copy_VarCoords(EIS, ECTEI)
;printVarSummary(LTS)
;printVarSummary(EIS)
;printVarSummary(ECTEI)

return ([/LTS,EIS,ECTEI/])

end

; 2018-11-06
; refer to Sherwood et al. (2014) to get low-tropospheric mixing indexes (LTMI)
function cal_LTMI(omega[*][*][*]:numeric,temp[*][*][*]:numeric,qv[*][*][*]:numeric,rh[*][*][*]:numeric,OCNFRAC[*][*]:numeric,lev[*]:numeric)
begin

; ==============================1. calculate SS ===============================
; S within a region where convective effects are a leading term in thermodynamic budgets, 
; defined by the upper quartile of the annual-mean midtropospheric ascent rate where it is upward, -omega_500 (omega the pressure velocity).
; ==============================================================================
omega_500			= omega({500},:,:)
omega_500_upward	= where((-1.*omega_500).gt.0,-1.*omega_500,omega_500@_FillValue)
;printMinMax(omega_500_upward,True)

; get the upper quartile of the annual-mean midtropospheric ascent rate
opt 				= True
opt@PrintStat 		= False
statb 				= stat_dispersion(ndtooned(omega_500_upward), opt )
print(statb(10))

; relative humidity gradient between 700hPa and 850hPa
rh_diff				= rh({700},:,:)-rh({850},:,:)
copy_VarCoords(rh(0,:,:),rh_diff)
;printVarSummary(rh_diff)
;printMinMax(rh_diff,True)

; temperature gradient between 700hPa and 850hPa
temp_diff			= temp({700},:,:)-temp({850},:,:)
copy_VarCoords(temp(0,:,:),temp_diff)
;printVarSummary(temp_diff)
;printMinMax(temp_diff,True)

; small-scale convective mixing strength --- S
SS_tmp					= (rh_diff/100.0 - temp_diff/9)/2.
copy_VarCoords(rh_diff,SS_tmp)
;printVarSummary(SS_tmp)
;printMinMax(SS_tmp,True)

; mask over the region where convective effects are a leading term in thermodynamic budgets
SS						= SS_tmp
SS						= where(omega_500_upward.gt.statb(10),SS,SS@_FillValue)
;printVarSummary(SS)
;printMinMax(SS,True)

; test the rightness of convective region for average
;wks						= gsn_open_wks("x11","test")
;
;res						= True
;res@cnFillOn			= True
;res@cnLinesOn			= False
;
;plot					= gsn_csm_contour(wks,SS,res)

; ==============================2. calculate DD ===============================

; get omega_1 (average of omega at 850hPa and 700 hPa)
omega_1					= (omega({700},:,:)+omega({850},:,:))/2.
copy_VarCoords(omega(0,:,:),omega_1)

; get omega_2 (average of omega at 600, 500 an 400 hPa)
omega_2					= (omega({600},:,:)+omega({500},:,:)+omega({400},:,:))/3.
copy_VarCoords(omega(0,:,:),omega_2)

; Delta measures the local horizontal outflow in the lower troposphere above the boundary layer.
; Moisture is transported upward and outward wherever Delta > 0 and omega_1 < 0.
; We restrict measurement to tropical ocean regions from 160 W to 30 E
; Delta = omega_2 - omega_1
Delta					= (omega_2 - omega_1)
copy_VarCoords(omega(0,:,:),Delta)
;printVarSummary(Delta)
;printMinMax(Delta,True)

;  the moisture supplied to the environment ===================
; M_(LT,Large) = -< qd_omega/d_p H(Delta) H (-omega_1) >
levels					= tofloat(omega&lev)
q_domega_dp				= 24.*3600.*qv*1e3*center_finite_diff_n(omega,100.*levels,False,0,0) ; g/kg/day

Delta_3D				= conform(qv,Delta,(/1,2/))
omega_1_3D				= conform(qv,omega_1,(/1,2/))
OCNFRAC_3D				= conform(qv,OCNFRAC,(/1,2/))
lon						= qv&lon
lat						= qv&lat
lon_3D					= conform(qv,lon,(/2/))
lat_3D					= conform(qv,lat,(/1/))
lon_2D					= conform(Delta,lon,(/1/))
lat_2D					= conform(Delta,lat,(/0/))

clat					= tofloat(cos(0.01745329*lat))
clat!0					= "lat"
clat&lat				= qv&lat

q_domega_dp				= where((Delta_3D.gt.0.and.omega_1_3D.lt.0.and.\
								 OCNFRAC_3D.gt.0.8.and.\
								 ((lon_3D.gt.200.and.lon_3D.lt.360).or.(lon_3D.gt.0.and.lon_3D.lt.30)).and.\
								 (lat_3D.gt.-30.and.lat_3D.lt.30)),q_domega_dp,qv@_FillValue)
copy_VarCoords(qv,q_domega_dp)
;printVarSummary(q_domega_dp)
;printMinMax(q_domega_dp,True)

q_domega_dp_avg			= wgt_areaave_Wrap(q_domega_dp,clat,1,0)
q_domega_dp_avg			= (/-1.*q_domega_dp_avg/)
printVarSummary(q_domega_dp_avg)
;print(q_domega_dp_avg&lev+"		"+q_domega_dp_avg)

; test the rightness of moisture supply
;wks						= gsn_open_wks("x11","test")
;
;res						= True
;res@cnFillOn			= True
;res@cnLinesOn			= False
;
;plot					= gsn_csm_contour(wks,q_domega_dp({850},:,:),res)



; the DD
D_upp				= where(Delta.gt.0.and.omega_1.lt.0.and.\
							OCNFRAC.gt.0.8.and.\
							((lon_2D.gt.200.and.lon_2D.lt.360.).or.(lon_2D.gt.0.and.lon_2D.lt.30)).and.\
							(lat_2D.gt.-30.and.lat_2D.lt.30),Delta,qv@_FillValue)
D_low				= where(omega_2.lt.0.and.\
							OCNFRAC.gt.0.8.and.\
							((lon_2D.gt.200.and.lon_2D.lt.360.).or.(lon_2D.gt.0.and.lon_2D.lt.30)).and.\
							(lat_2D.gt.-30.and.lat_2D.lt.30),-1.*omega_2,qv@_FillValue)
copy_VarCoords(Delta,D_upp)
copy_VarCoords(Delta,D_low)
;printMinMax(D_upp,True)
;printMinMax(D_low,True)

;wks						= gsn_open_wks("x11","test")
;
;res						= True
;res@cnFillOn			= True
;res@cnLinesOn			= False
;
;plot					= gsn_csm_contour(wks,D_upp,res)
;plot					= gsn_csm_contour(wks,D_low,res)
;

D_upp_avg			= wgt_areaave_Wrap(D_upp,clat,1,0)
D_low_avg			= wgt_areaave_Wrap(D_low,clat,1,0)

DD					= D_upp_avg/D_low_avg
print("D_upp_avg="+D_upp_avg+"		D_low_avg="+D_low_avg+"		DD="+DD)

; ==============================3. calculate LTMI ===============================
; get region average of SS and DD, and then get the sum of them.

SS_avg					= wgt_areaave_Wrap(SS,clat,1,0)
DD_avg					= DD

LTMI				= SS_avg + DD_avg
print("SS="+SS_avg+"		DD="+DD_avg+"		LTMI="+LTMI)


return ([/SS_avg,DD_avg,LTMI,q_domega_dp_avg/])

end



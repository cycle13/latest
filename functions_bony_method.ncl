
; 2018-09-23 
; add standard deviation1

function bony_method_std(omega[*]:numeric,var[*]:numeric,nbinsx_in:integer,nbinsy_in:integer,binx_min_in:float,binx_max_in:float)

begin

opt = True
opt@binx_min = binx_min_in
opt@binx_max = binx_max_in

;===== pdfxy
ap = pdfxy(omega, var, nbinsx_in, nbinsy_in, opt)

nbinsx = ap@nbinsx
binx_center = ap@binx_center

;===== do column integral of other variables based on vertical velocity bin
;===== use the biny center value to represent the mean value of each y-bin
biny_center = conform_dims(dimsizes(ap),ap@biny_center,0)

copy_VarCoords(ap,biny_center)

;printVarSummary(biny_center)

ap_x = new((/nbinsx/),double)
ap_x_std = new((/nbinsx/),double)

do ibinx = 0,dimsizes(binx_center)-1

wgt = ap(:,ibinx)

ap_x(ibinx) = dim_avg_wgt_n_Wrap(biny_center(:,ibinx),wgt,0,0)
; calculate standard deviation
; 1.1 get sum of weights in ibinx
sumwgt	 	= dim_sum_Wrap(wgt)
; 1.2 get the weight for each biny_center value
newwgt		= wgt/sumwgt
; 1.3 get standard deviation
tmp1			= ap_x
tmp2			= conform_dims(dimsizes(ap),tmp1,1)
tmp3		= sum(newwgt*(biny_center(:,ibinx)-tmp2(:,ibinx))^2)
ap_x_std(ibinx)	= (dimsizes(biny_center(:,0)-1))/dimsizes(biny_center(:,0))*sqrt(tmp3) ; (N-1)/N

;ap_x_std(ibinx) = dim_stddev_n_Wrap(biny_center(:,ibinx)*wgt,0)
end do

ap_x!0 = "binx_center"
ap_x&binx_center = binx_center
ap_x_std!0 = "binx_center"
ap_x_std&binx_center = binx_center

return([/ap,ap_x,ap_x_std/])

delete(ap)
delete(ap_x)
delete(ap_x_std)
delete(nbinsx)
delete(binx_center)
delete(biny_center)
delete(wgt)

end


function bony_method(omega[*]:numeric,var[*]:numeric,nbinsx_in:integer,nbinsy_in:integer,binx_min_in:float,binx_max_in:float)

begin

opt = True
opt@binx_min = binx_min_in
opt@binx_max = binx_max_in

;===== pdfxy
ap = pdfxy(omega, var, nbinsx_in, nbinsy_in, opt)

nbinsx = ap@nbinsx
binx_center = ap@binx_center

;===== do column integral of other variables based on vertical velocity bin
;===== use the biny center value to represent the mean value of each y-bin
biny_center = conform_dims(dimsizes(ap),ap@biny_center,0)

copy_VarCoords(ap,biny_center)

;printVarSummary(biny_center)

ap_x = new((/nbinsx/),double)

do ibinx = 0,dimsizes(binx_center)-1

wgt = ap(:,ibinx)

ap_x(ibinx) = dim_avg_wgt_n_Wrap(biny_center(:,ibinx),wgt,0,0)

end do

ap_x!0 = "binx_center"
ap_x&binx_center = binx_center

return([/ap,ap_x/])

delete(ap)
delete(nbinsx)
delete(binx_center)
delete(biny_center)
delete(ap_x)
delete(wgt)

end

function omega_pdf(omega[*]:numeric,nbins_in:integer,binx_min_in:float,binx_max_in:float)

begin

opt = True
opt@bin_min = binx_min_in
opt@bin_max = binx_max_in

;===== pdfxy
ap = pdfx(omega, nbins_in, opt)
ap!0	= "bin_center"
ap&bin_center	= ap@bin_center

return(ap)

delete(ap)

end

; -------------------------------------------------
; 2018.09.16
function pdf_omega_var_wo_obs(omega_nor[*][*][*][*]:numeric,var_nor[*][*][*][*]:numeric,\
		       	      nbinx,nbiny,nx_min,nx_max)

begin

tmp_omega_tropic = omega_nor(:,:,{-30:30},:)
tmp_var_tropic = var_nor(:,:,{-30:30},:)

;omega_nor_tropic = omega_nor(:,:,{-30:30},:)
;var_nor_tropic = var_nor(:,:,{-30:30},:)

dims	= dimsizes(tmp_omega_tropic)

omega_var_2D_all = new((/dims(0),nbiny,nbinx/),double)
omega_var_2D_all!0 = "case"
omega_var_2D_all&case = fspan(0,dims(0)-1,dims(0))

omega_var_1D_all = new((/dims(0),nbinx/),double)
omega_var_1D_all!0 = "case"
omega_var_1D_all&case = fspan(0,dims(0)-1,dims(0))

omega_var_1D_std 	= new((/dims(0),nbinx/),double)
omega_var_1D_std!0 	= "case"
omega_var_1D_std&case = fspan(0,dims(0)-1,dims(0))


do icase = 0,dims(0)-1
;=======change from multi-dimensions to one-dimension
tmp_omega_oned		= ndtooned(tmp_omega_tropic(icase,:,:,:))
tmp_var_oned		= ndtooned(tmp_var_tropic(icase,:,:,:))

;============new method try

;;;;;;;bony method starts here
; control result
if(.not.all(ismissing(tmp_var_oned)))then
	result = bony_method_std(tmp_omega_oned, tmp_var_oned,nbinx,nbiny,nx_min,nx_max)
	omega_var_2D_all(icase,:,:) = result[0]
	omega_var_1D_all(icase,:) = result[1]
	omega_var_1D_std(icase,:) = result[2] ; 2018.09.23
	delete(result)
else
	omega_var_2D_all(icase,:,:)		= tmp_omega_tropic@_FillValue
	omega_var_1D_all(icase,:)		= tmp_omega_tropic@_FillValue
	omega_var_1D_std(icase,:)		= tmp_omega_tropic@_FillValue
end if

delete([/tmp_omega_oned, tmp_var_oned/])
end do ; do icase =

return([/omega_var_1D_all,omega_var_1D_std/])

delete(omega_var_1D_all)
delete(oemga_var_1D_std)
delete(omega_var_2D_all)

end


function pdf_sgm_var_wo_obs(sgm_nor[*][*][*][*]:numeric,var_nor[*][*][*][*]:numeric,\
		       	      nbinx,nbiny,nx_min,nx_max)

begin

sgm_nor_tropic = sgm_nor(:,:,{-30:30},:)
var_nor_tropic = var_nor(:,:,{-30:30},:)

;=======change from multi-dimensions to one-dimension
sgm_nor_def_oned = ndtooned(sgm_nor_tropic(0,:,:,:))
sgm_nor_new_oned = ndtooned(sgm_nor_tropic(1,:,:,:))

var_nor_def_oned = ndtooned(var_nor_tropic(0,:,:,:))
var_nor_new_oned = ndtooned(var_nor_tropic(1,:,:,:))

;============
;============new method try

;;;;;;;bony method starts here
; control result
result = bony_method_std(sgm_nor_def_oned, var_nor_def_oned,nbinx,nbiny,nx_min,nx_max)
sgm_var_2D_nor_def = result[0]
sgm_var_1D_nor_def = result[1]
sgm_var_1D_nor_def_std = result[2]
delete(result)
result = bony_method_std(sgm_nor_new_oned, var_nor_new_oned,nbinx,nbiny,nx_min,nx_max)
sgm_var_2D_nor_new = result[0]
sgm_var_1D_nor_new = result[1]
sgm_var_1D_nor_new_std = result[2]
delete(result)

; merge them into one
sgm_var_1D_all = new((/2,dimsizes(sgm_var_1D_nor_new)/),double)

sgm_var_1D_all!0 = "case"
sgm_var_1D_all&case = ((/0,1/))

sgm_var_1D_all(0,:) = sgm_var_1D_nor_def
sgm_var_1D_all(1,:) = sgm_var_1D_nor_new

sgm_var_1D_std = sgm_var_1D_all
sgm_var_1D_std(0,:) = sgm_var_1D_nor_def_std
sgm_var_1D_std(1,:) = sgm_var_1D_nor_new_std

return([/sgm_var_1D_all,sgm_var_1D_std/])

delete(sgm_var_1D_all)

end


function pdf_omega_var(OMEGA_nor[*][*][*][*]:numeric,var_nor[*][*][*][*]:numeric,\
                       OMEGA_4K[*][*][*][*]:numeric,var_4K[*][*][*][*]:numeric,\
		       nbinx,nbiny,nx_min,nx_max)

begin

OMEGA_nor_tropic = OMEGA_nor(:,:,{-30:30},:)
 var_nor_tropic =  var_nor(:,:,{-30:30},:)

OMEGA_4K_tropic = OMEGA_4K(:,:,{-30:30},:)
 var_4K_tropic =  var_4K(:,:,{-30:30},:)

;=======change from multi-dimensions to one-dimension
OMEGA_nor_obs_oned = ndtooned(OMEGA_nor_tropic(0,:,:,:))
OMEGA_nor_def_oned = ndtooned(OMEGA_nor_tropic(1,:,:,:))
OMEGA_nor_new_oned = ndtooned(OMEGA_nor_tropic(2,:,:,:))

OMEGA_4K_obs_oned = ndtooned(OMEGA_4K_tropic(0,:,:,:))
OMEGA_4K_def_oned = ndtooned(OMEGA_4K_tropic(1,:,:,:))
OMEGA_4K_new_oned = ndtooned(OMEGA_4K_tropic(2,:,:,:))

var_nor_obs_oned = ndtooned(var_nor_tropic(0,:,:,:))
var_nor_def_oned = ndtooned(var_nor_tropic(1,:,:,:))
var_nor_new_oned = ndtooned(var_nor_tropic(2,:,:,:))

var_4K_obs_oned = ndtooned(var_4K_tropic(0,:,:,:))
var_4K_def_oned = ndtooned(var_4K_tropic(1,:,:,:))
var_4K_new_oned = ndtooned(var_4K_tropic(2,:,:,:))

;============
;============new method try

;;;;;;;bony method starts here
; control result
result = bony_method(OMEGA_nor_obs_oned, var_nor_obs_oned,nbinx,nbiny,nx_min,nx_max)
OMEGA_var_2D_nor_obs = result[0]
OMEGA_var_1D_nor_obs = result[1]
delete(result)
result = bony_method(OMEGA_nor_def_oned, var_nor_def_oned,nbinx,nbiny,nx_min,nx_max)
OMEGA_var_2D_nor_def = result[0]
OMEGA_var_1D_nor_def = result[1]
delete(result)
result = bony_method(OMEGA_nor_new_oned, var_nor_new_oned,nbinx,nbiny,nx_min,nx_max)
OMEGA_var_2D_nor_new = result[0]
OMEGA_var_1D_nor_new = result[1]
delete(result)

; 4K result
result = bony_method(OMEGA_4K_obs_oned, var_4K_obs_oned,nbinx,nbiny,nx_min,nx_max)
OMEGA_var_2D_4K_obs = result[0]
OMEGA_var_1D_4K_obs = result[1]
delete(result)
result = bony_method(OMEGA_4K_def_oned, var_4K_def_oned,nbinx,nbiny,nx_min,nx_max)
OMEGA_var_2D_4K_def = result[0]
OMEGA_var_1D_4K_def = result[1]
delete(result)
result = bony_method(OMEGA_4K_new_oned, var_4K_new_oned,nbinx,nbiny,nx_min,nx_max)
OMEGA_var_2D_4K_new = result[0]
OMEGA_var_1D_4K_new = result[1]
delete(result)

; merge them into one
OMEGA_var_1D_all = new((/5,dimsizes(OMEGA_var_1D_nor_new)/),double)

OMEGA_var_1D_all!0 = "ncases"
OMEGA_var_1D_all&ncases = ((/0,1,2,3,4/))

OMEGA_var_1D_all(0,:) = OMEGA_var_1D_nor_obs
OMEGA_var_1D_all(1,:) = OMEGA_var_1D_nor_def
OMEGA_var_1D_all(2,:) = OMEGA_var_1D_nor_new
OMEGA_var_1D_all(3,:) = OMEGA_var_1D_4K_def
OMEGA_var_1D_all(4,:) = OMEGA_var_1D_4K_new

return(OMEGA_var_1D_all)

delete(OMEGA_var_1D_all)

end

; 2018-07-13
; 2018-10-02 change into loop format 
function pdf_regimevar_only(regimevar_nor[*][*][*][*]:numeric,\
		       nbinx,nx_min,nx_max)

begin

regimevar_tropic = regimevar_nor(:,:,{-30:30},:)
dims	= dimsizes(regimevar_tropic)

regimevar_PDF_all = new((/dims(0),nbinx/),double)
regimevar_PDF_all!0	= "case"
regimevar_PDF_all&case	= fspan(0,dims(0)-1,dims(0))

do icase=0,dims(0)-1
	;=======change from multi-dimensions to one-dimension
	regimevar_oned = ndtooned(regimevar_tropic(icase,:,:,:))
	
	;;; omega pdf
	result = omega_pdf(regimevar_oned,nbinx,nx_min,nx_max)
	regimevar_PDF_all(icase,:) = result
	
	delete([/result,regimevar_oned/])

end do ; do icase=

return(regimevar_PDF_all)

delete(regimevar_PDF_all)

end




function pdf_omega_only(OMEGA_nor[*][*][*][*]:numeric,\
                       OMEGA_4K[*][*][*][*]:numeric,\
		       nbinx,nx_min,nx_max)

begin

OMEGA_nor_tropic = OMEGA_nor(:,:,{-30:30},:)

OMEGA_4K_tropic = OMEGA_4K(:,:,{-30:30},:)

;=======change from multi-dimensions to one-dimension
OMEGA_nor_obs_oned = ndtooned(OMEGA_nor_tropic(0,:,:,:))
OMEGA_nor_def_oned = ndtooned(OMEGA_nor_tropic(1,:,:,:))
OMEGA_nor_new_oned = ndtooned(OMEGA_nor_tropic(2,:,:,:))

OMEGA_4K_obs_oned = ndtooned(OMEGA_4K_tropic(0,:,:,:))
OMEGA_4K_def_oned = ndtooned(OMEGA_4K_tropic(1,:,:,:))
OMEGA_4K_new_oned = ndtooned(OMEGA_4K_tropic(2,:,:,:))


;;; omega pdf
result = omega_pdf(OMEGA_nor_obs_oned,nbinx,nx_min,nx_max)
OMEGA_PDF_nor_obs = result
delete(result)
printVarSummary(OMEGA_PDF_nor_obs)

result = omega_pdf(OMEGA_nor_def_oned,nbinx,nx_min,nx_max)
OMEGA_PDF_nor_def = result
delete(result)

result = omega_pdf(OMEGA_nor_new_oned,nbinx,nx_min,nx_max)
OMEGA_PDF_nor_new = result
delete(result)

;;;;
result = omega_pdf(OMEGA_4K_obs_oned,nbinx,nx_min,nx_max)
OMEGA_PDF_4K_obs = result
delete(result)

result = omega_pdf(OMEGA_4K_def_oned,nbinx,nx_min,nx_max)
OMEGA_PDF_4K_def = result
delete(result)

result = omega_pdf(OMEGA_4K_new_oned,nbinx,nx_min,nx_max)
OMEGA_PDF_4K_new = result
delete(result)


OMEGA_PDF_all = new((/5,dimsizes(OMEGA_PDF_nor_obs)/),double)
OMEGA_PDF_all(0,:) = OMEGA_PDF_nor_obs
OMEGA_PDF_all(1,:) = OMEGA_PDF_nor_def
OMEGA_PDF_all(2,:) = OMEGA_PDF_nor_new
OMEGA_PDF_all(3,:) = OMEGA_PDF_4K_def
OMEGA_PDF_all(4,:) = OMEGA_PDF_4K_new

return(OMEGA_PDF_all)

delete(OMEGA_PDF_all)

end






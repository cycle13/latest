
load "./functions_regrid.ncl"
load "./functions_getvar.ncl"
load "./functions_calculate.ncl"

load "./functions_matrix.ncl"

begin

;***************file 1********************

horizon_interp			= False ; used for different input data resolution
resolution				= "2deg"

ALPI_revised			= True ; overall control
ALPI_mask_ocn			= True
ALPIsub_revised			= True ; from LTS to EIS

ALPIsub_mask_plot		= False
ALPIsub_bin_2d_plot		= False
ALPIsub_bin_1d_plot		= False 
sgm_diff_plot			= False	; plot sgm_diff for sgm_tota, sgm_shal, sgm_turb in one figure by overlay
sgm_shal_4terms			= True  ; term1, term2 and term3 for variance of shallow convection in one figure
								; sgm_shal_4terms should work with sgm_diff_plot = True.
plot_type			= "pdf"

; add date to each plot 
date=systemfunc("date")
print(date)

latS_trop			= -30.
latE_trop			= 30.

spec_lev			= 887


nbins				= 10
ALPI_mask_val		= fspan(0,nbins-1,nbins)+0.5 ;-- the addition of 0.5 is to make sure that 
												 ;-- the value is in the center of cnLevels rather 
												 ;-- than the edge. After this change, RasterFill
												 ;-- is not needed and masking plot would be right as well.


if(resolution.eq."2deg")then
	case_nor 				= (/"FAMIPC5_f19f19_MG15_amip","FAMIPC5_f19f19_mac2_amip_outSGM"/)
	case_p4K		 		= (/"FAMIPC5_f19f19_MG15_amip-p4K","FAMIPC5_f19f19_mac2_amip-p4K_debug_outSGM"/)
else
	case_nor 				= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_mac2_amip_outSGM"/)
	case_p4K		 		= (/"FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM"/)
end if


case_name				= (/"RH","PDF"/)

case					= (/case_nor,case_p4K/)
tag					= (/"nor", "p4K"/)

append 				= ".ANN.climo.nc"

; _1 for the difference between PDF and RH; 
; _2 for the difference between p4K and CTL.
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

CLOUD_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)*0.5
CLOUD_dlevs_1		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)*0.2
CLOUD_dlevs_2		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)*0.1
CLOUD_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

CLDLIQ_levs			= (/0.01,0.02,0.03,0.04,0.05,0.08,0.1,0.12,0.14,0.16,0.18,0.20,0.22,0.24,0.26/)
CLDLIQ_dlevs_1		= (/-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,-0.01,0,0.01,0.02,0.04,0.06,0.08,0.1,0.12/)
CLDLIQ_dlevs_2		= (/-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,-0.01,0,0.01,0.02,0.04,0.06,0.08,0.1,0.12/)
CLDLIQ_dlevs_sp		= (/-0.12, -0.02, 0.02, 0.12/)

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

sgm_tota_levs		= fspan(0.,0.3,15)*0.2
sgm_tota_dlevs_1	= fspan(-0.1,0.1,15)*0.05
sgm_tota_dlevs_2	= sgm_tota_dlevs_1
sgm_tota_dlevs_sp	= (/-1, -0.5, 0.5, 1/)

sgm_shal_levs		= fspan(0.,0.3,15)
sgm_shal_dlevs_1	= fspan(-0.1,0.1,15)
sgm_shal_dlevs_2	= sgm_tota_dlevs_1
sgm_shal_dlevs_sp	= (/-1, -0.5, 0.5, 1/)

sgm_turb_levs		= fspan(0.,0.84,15)
sgm_turb_dlevs_1	= fspan(-0.28,0.28,15)
sgm_turb_dlevs_2	= sgm_tota_dlevs_1
sgm_turb_dlevs_sp	= (/-1, -0.5, 0.5, 1/)

N1_levs				= fspan(0,0.84,15)*0.25
N1_dlevs_1			= fspan(-0.28,0.28,15)*0.05
N1_dlevs_2			= N1_dlevs_1
N1_dlevs_sp			= (/-1, -0.5, 0., 1/)

netCRF_levs			= (/-80,-75,-70,-65,-60,-55,-50,-45,-40,-30,-20,-10,0,10,20/)
netCRF_dlevs_1		= (/-50,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,50/)/2.
netCRF_dlevs_2		= netCRF_dlevs_1
netCRF_dlevs_sp		= (/-20,-10,10,20/)

RELHUM_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)*100.
RELHUM_dlevs_1		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)*100.*0.5
RELHUM_dlevs_2		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)*100.*0.5
RELHUM_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)*100.

FREQSH_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
FREQSH_dlevs_1		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
FREQSH_dlevs_2		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
FREQSH_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

TREFHT_levs			= (/210,220,230,240,250,260,270,275,280,285,290,295,300,305,310/)
TREFHT_dlevs_1		= (/-5,-4,-3,-2,-1,-.5,-.2,0,.2,.5,1,2,3,4,5/)
TREFHT_dlevs_2		= (/-5,-4,-3,-2,-1,-.5,-.2,0,.2,.5,1,2,3,4,5/)
TREFHT_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

CMFMC_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
CMFMC_dlevs_1		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
CMFMC_dlevs_2		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
CMFMC_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

umf_shal_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
umf_shal_dlevs_1		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
umf_shal_dlevs_2		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
umf_shal_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

clddep2_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
clddep2_dlevs_1		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
clddep2_dlevs_2		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
clddep2_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

wstarPBL_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
wstarPBL_dlevs_1		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
wstarPBL_dlevs_2		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
wstarPBL_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

lengi_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
lengi_dlevs_1		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
lengi_dlevs_2		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
lengi_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

shi_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
shi_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
shi_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
shi_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

aa_gp_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
aa_gp_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
aa_gp_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
aa_gp_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

bb_gp_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
bb_gp_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
bb_gp_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
bb_gp_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

beta_gp_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
beta_gp_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
beta_gp_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
beta_gp_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

dqwdz_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
dqwdz_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
dqwdz_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
dqwdz_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

dthldz_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
dthldz_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
dthldz_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
dthldz_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

qw2_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
qw2_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
qw2_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
qw2_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

thl2_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
thl2_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
thl2_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
thl2_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

qwthl_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
qwthl_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
qwthl_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
qwthl_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

qwu_qw_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
qwu_qw_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
qwu_qw_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
qwu_qw_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

thlu_thl_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
thlu_thl_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
thlu_thl_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
thlu_thl_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

qwu_thlu_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
qwu_thlu_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
qwu_thlu_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
qwu_thlu_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

qtu_shal_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
qtu_shal_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
qtu_shal_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
qtu_shal_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

vert_qwflx_shal_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
vert_qwflx_shal_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
vert_qwflx_shal_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
vert_qwflx_shal_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

vert_thlflx_shal_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
vert_thlflx_shal_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
vert_thlflx_shal_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
vert_thlflx_shal_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

vert_qwthlflx_shal_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
vert_qwthlflx_shal_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
vert_qwthlflx_shal_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
vert_qwthlflx_shal_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

tau_shal_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
tau_shal_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
tau_shal_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
tau_shal_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

vert_qwflx_turb_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
vert_qwflx_turb_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
vert_qwflx_turb_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
vert_qwflx_turb_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

vert_thlflx_turb_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
vert_thlflx_turb_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
vert_thlflx_turb_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
vert_thlflx_turb_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

vert_qwthlflx_turb_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
vert_qwthlflx_turb_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
vert_qwthlflx_turb_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
vert_qwthlflx_turb_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

lengi2_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
lengi2_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
lengi2_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
lengi2_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

dqwdz_2_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
dqwdz_2_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
dqwdz_2_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
dqwdz_2_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

dthldz_2_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
dthldz_2_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
dthldz_2_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
dthldz_2_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

dqwdz_dthldz_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
dqwdz_dthldz_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
dqwdz_dthldz_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
dqwdz_dthldz_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

aa2_gp_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
aa2_gp_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
aa2_gp_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
aa2_gp_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

bb2_gp_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
bb2_gp_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
bb2_gp_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
bb2_gp_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

aa2_dqwdz2_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
aa2_dqwdz2_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
aa2_dqwdz2_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
aa2_dqwdz2_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

bb2_dthldz2_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
bb2_dthldz2_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
bb2_dthldz2_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
bb2_dthldz2_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

ab_dqwdz_dthldz_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
ab_dqwdz_dthldz_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
ab_dqwdz_dthldz_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
ab_dqwdz_dthldz_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

shi_lengi2_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
shi_lengi2_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
shi_lengi2_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
shi_lengi2_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

term1_turb_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
term1_turb_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
term1_turb_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
term1_turb_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

term2_turb_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
term2_turb_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
term2_turb_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
term2_turb_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

term3_turb_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
term3_turb_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
term3_turb_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
term3_turb_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

term1_shal_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
term1_shal_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
term1_shal_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
term1_shal_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

term2_shal_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
term2_shal_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
term2_shal_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
term2_shal_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

term3_shal_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
term3_shal_dlevs_1			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
term3_shal_dlevs_2			= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
term3_shal_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)


Vars_mid	= (/"T", "Q", "OMEGA", "RELHUM"/)

Vars_pr		= (/"PRECT", "PRECL", "PRECC"/)

;float FISCCP1_COSP(time, cosp_prs, cosp_tau, lat, lon)
Vars_cosp	= (/"FISCCP1_COSP"/)

; --- 4D variables for ALPI_revised (time,lev,lat,lon)
;Vars_4d		= (/"CLOUD","CLDLIQ","sgm_tota","sgm_shal","sgm_turb"/)
;Vars_4d			= (/"CLOUD"/)
;Vars_3d			= (/"netCRF","SWCF","CLDLOW","TGCLDLWP"/)

Vars_mid_4d		= (/"T", "Q", "OMEGA", "RELHUM"/)
Vars_mid_3d		= (/"PRECT", "PRECL", "PRECC"/)


; --- mention:
; vert_qwflx_xxx and vert_qwflx_xxx is not the same in shal and turb.
; vert_qwflx_shal includes the effect of aa_gp, while vert_qwflx_turb not.
; so pls dont use vert_qwflx_turb later if possible.

;				   0            1          2           3           4            5          6	    7      8       9		10	      11	  12 	 13        14         15           16        17         18      19			20		   21         22      23     24       25         26       27      28        29      30     31       32         33          34			35               36				    37				   38			   39			   40			    41					42	            43       44        45           46           47        48          49         50              51               52            53            54           55               56           57           58          59
varsA		= (/"PRECT",	"PRECC",	"PRECL",	"SWCF",		"LWCF",		"CLDLOW",	"CLDHGH", "LTS", "EIS", "ECTEI","TGCLDLWP","sgm_tota","N1", "netCRF","sgm_shal","sgm_turb", "RELHUM", "FREQSH", "CMFMC","umf_shal", "clddep2", "wstarPBL", "lengi", "shi", "aa_gp",	"bb_gp", "beta_gp", "dqwdz", "dthldz", "qw2", "thl2", "qwthl", "qwu_qw", "thlu_thl", "qwu_thlu", "qtu_shal", "vert_qwflx_shal", "vert_thlflx_shal","vert_qwthlflx_shal","tau_shal","vert_qwflx_turb","vert_thlflx_turb","vert_qwthlflx_turb","lengi2","dqwdz_2","dthldz_2","dqwdz_dthldz","aa2_gp", "bb2_gp", "aa2_dqwdz2","bb2_dthldz2","ab_dqwdz_dthldz","shi_lengi2", "term1_turb", "term2_turb", "term3_turb", "term1_shal", "term2_shal", "term3_shal","TREFHT","CLOUD", "CLDLIQ"/)
varsA_levs	= (/PRECT_levs, PRECT_levs,     PRECT_levs,     SWCF_levs,      LWCF_levs, 	CLDLOW_levs,	CLDHGH_levs, LTS_levs, EIS_levs, ECTEI_levs,TGCLDLWP_levs,sgm_tota_levs, N1_levs, netCRF_levs,sgm_shal_levs, sgm_turb_levs, RELHUM_levs, FREQSH_levs, CMFMC_levs, umf_shal_levs, clddep2_levs, wstarPBL_levs, lengi_levs, shi_levs, aa_gp_levs, bb_gp_levs, beta_gp_levs, dqwdz_levs, dthldz_levs, qw2_levs, thl2_levs, qwthl_levs, qwu_qw_levs, thlu_thl_levs, qwu_thlu_levs, qtu_shal_levs,vert_qwflx_shal_levs, vert_thlflx_shal_levs, vert_qwthlflx_shal_levs, tau_shal_levs, vert_qwflx_turb_levs, vert_thlflx_turb_levs, vert_qwthlflx_turb_levs, lengi2_levs,dqwdz_2_levs, dthldz_2_levs, dqwdz_dthldz_levs, aa2_gp_levs, bb2_gp_levs, aa2_dqwdz2_levs, bb2_dthldz2_levs, ab_dqwdz_dthldz_levs, shi_lengi2_levs, term1_turb_levs, term2_turb_levs, term3_turb_levs, term1_shal_levs, term2_shal_levs, term3_shal_levs, TREFHT_levs, CLOUD_levs, CLDLIQ_levs/)
varsA_dlevs_1	= (/PRECT_dlevs_1,PRECT_dlevs_1,    PRECT_dlevs_1,    SWCF_dlevs_1,     LWCF_dlevs_1,	CLDLOW_dlevs_1,	CLDHGH_dlevs_1, LTS_dlevs_1, EIS_dlevs_1,	ECTEI_dlevs_1,TGCLDLWP_dlevs_1, sgm_tota_dlevs_1, N1_dlevs_1, netCRF_dlevs_1, sgm_shal_dlevs_1, sgm_turb_dlevs_1, RELHUM_dlevs_1, FREQSH_dlevs_1, CMFMC_dlevs_1, umf_shal_dlevs_1, clddep2_dlevs_1, wstarPBL_dlevs_1, lengi_dlevs_1, shi_dlevs_1, aa_gp_dlevs_1, bb_gp_dlevs_1, beta_gp_dlevs_1, dqwdz_dlevs_1, dthldz_dlevs_1, qw2_dlevs_1, thl2_dlevs_1, qwthl_dlevs_1, qwu_qw_dlevs_1, thlu_thl_dlevs_1, qwu_thlu_dlevs_1, qtu_shal_dlevs_1,vert_qwflx_shal_dlevs_1, vert_thlflx_shal_dlevs_1, vert_qwthlflx_shal_dlevs_1, tau_shal_dlevs_1, vert_qwflx_turb_dlevs_1, vert_thlflx_turb_dlevs_1, vert_qwthlflx_turb_dlevs_1, lengi2_dlevs_1, dqwdz_2_dlevs_1, dthldz_2_dlevs_1, dqwdz_dthldz_dlevs_1, aa2_gp_dlevs_1, bb2_gp_dlevs_1,aa2_dqwdz2_dlevs_1, bb2_dthldz2_dlevs_1, ab_dqwdz_dthldz_dlevs_1, shi_lengi2_dlevs_1, term1_turb_dlevs_1, term2_turb_dlevs_1, term3_turb_dlevs_1, term1_shal_dlevs_1, term2_shal_dlevs_1, term3_shal_dlevs_1, TREFHT_dlevs_1, CLOUD_dlevs_1, CLDLIQ_dlevs_1/)
varsA_dlevs_2	= (/PRECT_dlevs_2,PRECT_dlevs_2,    PRECT_dlevs_2,    SWCF_dlevs_2,     LWCF_dlevs_2,	CLDLOW_dlevs_2,	CLDHGH_dlevs_2, LTS_dlevs_2, EIS_dlevs_2,	ECTEI_dlevs_2,TGCLDLWP_dlevs_2, sgm_tota_dlevs_2, N1_dlevs_2, netCRF_dlevs_2, sgm_shal_dlevs_2, sgm_turb_dlevs_2, RELHUM_dlevs_2, FREQSH_dlevs_2, CMFMC_dlevs_2, umf_shal_dlevs_2, clddep2_dlevs_2, wstarPBL_dlevs_2, lengi_dlevs_2, shi_dlevs_2, aa_gp_dlevs_2, bb_gp_dlevs_2, beta_gp_dlevs_2, dqwdz_dlevs_2, dthldz_dlevs_2, qw2_dlevs_2, thl2_dlevs_2, qwthl_dlevs_2, qwu_qw_dlevs_2, thlu_thl_dlevs_2, qwu_thlu_dlevs_2, qtu_shal_dlevs_2,vert_qwflx_shal_dlevs_2, vert_thlflx_shal_dlevs_2, vert_qwthlflx_shal_dlevs_2, tau_shal_dlevs_2, vert_qwflx_turb_dlevs_2, vert_thlflx_turb_dlevs_2, vert_qwthlflx_turb_dlevs_2,lengi2_dlevs_2, dqwdz_2_dlevs_2, dthldz_2_dlevs_2, dqwdz_dthldz_dlevs_2, aa2_gp_dlevs_2, bb2_gp_dlevs_2,aa2_dqwdz2_dlevs_2, bb2_dthldz2_dlevs_2, ab_dqwdz_dthldz_dlevs_2, shi_lengi2_dlevs_2, term1_turb_dlevs_2, term2_turb_dlevs_2, term3_turb_dlevs_2, term1_shal_dlevs_2, term2_shal_dlevs_2, term3_shal_dlevs_2, TREFHT_dlevs_2, CLOUD_dlevs_2, CLDLIQ_dlevs_2/)
varsA_dlevs_sp	= (/PRECT_dlevs_sp,PRECT_dlevs_sp,    PRECT_dlevs_sp,    SWCF_dlevs_sp,     LWCF_dlevs_sp,	CLDLOW_dlevs_sp,	CLDHGH_dlevs_sp, LTS_dlevs_sp, EIS_dlevs_sp, ECTEI_dlevs_sp,TGCLDLWP_dlevs_sp, sgm_tota_dlevs_sp, N1_dlevs_sp, netCRF_dlevs_sp, sgm_shal_dlevs_sp, sgm_turb_dlevs_sp, RELHUM_dlevs_sp, FREQSH_dlevs_sp, CMFMC_dlevs_sp, umf_shal_dlevs_sp, clddep2_dlevs_sp, wstarPBL_dlevs_sp, lengi_dlevs_sp, shi_dlevs_sp, aa_gp_dlevs_sp, bb_gp_dlevs_sp, beta_gp_dlevs_sp, dqwdz_dlevs_sp, dthldz_dlevs_sp, qw2_dlevs_sp, thl2_dlevs_sp, qwthl_dlevs_sp, qwu_qw_dlevs_sp, thlu_thl_dlevs_sp, qwu_thlu_dlevs_sp, qtu_shal_dlevs_sp, vert_qwflx_shal_dlevs_sp, vert_thlflx_shal_dlevs_sp, vert_qwthlflx_shal_dlevs_sp, tau_shal_dlevs_sp, vert_qwflx_turb_dlevs_sp, vert_thlflx_turb_dlevs_sp, vert_qwthlflx_turb_dlevs_sp, lengi2_dlevs_sp, dqwdz_2_dlevs_sp, dthldz_2_dlevs_sp, dqwdz_dthldz_dlevs_sp, aa2_gp_dlevs_sp, bb2_gp_dlevs_sp, aa2_dqwdz2_dlevs_sp, bb2_dthldz2_dlevs_sp, ab_dqwdz_dthldz_dlevs_sp, shi_lengi2_dlevs_sp, term1_turb_dlevs_sp, term2_turb_dlevs_sp, term3_turb_dlevs_sp, term1_shal_dlevs_sp, term2_shal_dlevs_sp, term3_shal_dlevs_sp, TREFHT_dlevs_sp, CLOUD_dlevs_sp, CLDLIQ_dlevs_sp/)

varsA_units	= (/"mm/day",   "mm/day",       "mm/day",       "W/m~S~2~N~",   "W/m~S~2~N~",	"fraction",	"fraction", "K",	"K",	"K", "g/m~S~2~N~", "-", "fraction", "W/m~S~2~N~", "-", "-", "%", "fraction","kg/m~S~2~N~/s", "kg/m~S~2~N~/s","m", "m/s", "m", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","K", "fraction", "g/kg"/)
plotsA_strings	= (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)","(m)","(n)","(o)","(p)","(q)","(r)","(s)","(t)","(u)","(v)","(w)","(x)","(y)","(z)"/)


; --- ALPI for plot
if(sgm_diff_plot .and. .not.sgm_shal_4terms)then
	ALPIsub_bin_1d_plot	= True
	pick_var_3d			= (/ "sgm_tota","sgm_shal","sgm_turb"/)	; sgm_diff_plot
end if

if(sgm_shal_4terms)then
	ALPIsub_bin_1d_plot	= True
	sgm_diff_plot		= True
	pick_var_3d			= (/ "sgm_tota","sgm_shal","sgm_turb","term1_shal","term2_shal","term3_shal"/)	; sgm_shal_4terms
;	pick_var_3d			= (/"sgm_shal","term1_shal","term2_shal","term3_shal"/) ; sgm_shal_4terms
end if

if(ALPIsub_bin_2d_plot)then
	pick_var_4d			= (/"CLOUD","CLDLIQ","RELHUM","sgm_tota","N1"/)
end if

if(ALPIsub_bin_1d_plot .and. .not. sgm_shal_4terms .and. .not.sgm_diff_plot)then
	pick_var_3d			= (/"netCRF","SWCF","CLDLOW","TGCLDLWP","sgm_tota","N1"/) ; normal
end if


if(ALPIsub_bin_2d_plot)then
	use_idx_4d		= new(dimsizes(pick_var_4d),integer)
	do ivar=0,dimsizes(pick_var_4d)-1
		use_idx_4d(ivar)	= ind(varsA.eq.pick_var_4d(ivar))
		print(varsA(ivar))
	end do ; do ivar=

	Vars_4d					= varsA(use_idx_4d)
	Vars_4d_levs			= varsA_levs(use_idx_4d,:)
	Vars_4d_dlevs_1			= varsA_dlevs_1(use_idx_4d,:)
	Vars_4d_dlevs_2			= varsA_dlevs_2(use_idx_4d,:)
	Vars_4d_dlevs_sp		= varsA_dlevs_sp(use_idx_4d,:)
	Vars_4d_units			= varsA_units(use_idx_4d)
	plots_strings			= plotsA_strings
end if

if(ALPIsub_bin_1d_plot)then
	use_idx_3d		= new(dimsizes(pick_var_3d),integer)
	do ivar=0,dimsizes(pick_var_3d)-1
		use_idx_3d(ivar)	= ind(varsA.eq.pick_var_3d(ivar))
		print(varsA(ivar))
	end do ; do ivar=

	Vars_3d					= varsA(use_idx_3d)
	Vars_3d_levs			= varsA_levs(use_idx_3d,:)
	Vars_3d_dlevs_1			= varsA_dlevs_1(use_idx_3d,:)
	Vars_3d_dlevs_2			= varsA_dlevs_2(use_idx_3d,:)
	Vars_3d_dlevs_sp		= varsA_dlevs_sp(use_idx_3d,:)
	Vars_3d_units			= varsA_units(use_idx_3d)
	plots_strings			= plotsA_strings
end if

; --- for variables only in Gauss-PDF
use_idx_pdf			= (/11,12,14,15,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58/)

;---- variables with the same Y axis
use_idx_yaxis_sgm		= (/11,14,15,56,57,58/) ; sgm_tota, sgm_turb, sgm_shal
use_idx_yaxis_crf		= (/13,3,4/)

vars_pdf			= varsA(use_idx_pdf)

vars_yaxis = NewList("lifo")

ListAppend(vars_yaxis,(/"sgm_tota","sgm_shal","sgm_turb",\
				    	"term1_shal","term2_shal","term3_shal"/))
ListAppend(vars_yaxis,(/"netCRF","SWCF","LWCF"/))
print(vars_yaxis)

nlev						= 30

if(.not.horizon_interp)then
	if(resolution.eq."2deg")then
		nlat				= 96
		nlon				= 144
	else
		nlat				= 192
		nlon				= 288
	end if
else
	nlat				= 73		; number of output latitudes
	nlon				= 144		; number of output longitudes
	klon = ispan(0,nlon-1,1)*2.5	; to 2.5 x 2.5
	klat = ispan(0,nlat-1,1)*2.5-90.
	
	;nlat				= 181		; number of output latitudes
	;nlon				= 360		; number of output longitudes
	;klon = ispan(0,nlon-1,1)*1.0	; to 1.0 x 1.0
	;klat = ispan(0,nlat-1,1)*1.0-90.
end if



;--------------------------------------------------------ALPIS-------------------------------------
if (ALPI_revised)then ; revised -- the overall switch for ALPI kind of plots and calculation

; -------------------------------------------- read data from mid-data ----------------------------------

indir						= "/home/share3/lyl/work3/qinyi/mid-data/cfmip/"
infile						= "ALPI_binned_vars.nc"

fin							= addfile(indir+infile,"r")
time						= fin->time
lev							= fin->lev
bin							= fin->bin
ntime						= dimsizes(time)
nlev						= dimsizes(lev)
nbin						= dimsizes(bin)

ALPI						= fin->ALPI
ALPI_mask					= fin->ALPI_mask
printVarSummary(ALPI)
printMinMax(ALPI,True)
printVarSummary(ALPI_mask)
printMinMax(ALPI_mask,True)

if(ALPIsub_bin_2d_plot)then
	vars_4d_bin_mean			= new((/dimsizes(Vars_4d),dimsizes(tag),dimsizes(case_nor),ntime,nlev,nbins/),float)
	vars_4d_bin_std				= new((/dimsizes(Vars_4d),dimsizes(tag),dimsizes(case_nor),ntime,nlev,nbins/),float)
	printVarSummary(vars_4d_bin_mean)

	Vars_out_4d					= Vars_4d
end if

if(ALPIsub_bin_1d_plot)then
	vars_3d_bin_mean			= new((/dimsizes(Vars_3d),dimsizes(tag),dimsizes(case_nor),ntime,nbins/),float)
	vars_3d_bin_std				= new((/dimsizes(Vars_3d),dimsizes(tag),dimsizes(case_nor),ntime,nbins/),float)

	Vars_out_3d					= Vars_3d
end if




if(ALPIsub_bin_2d_plot)then
	do itag=0,dimsizes(tag)-1
		do ivar=0,dimsizes(Vars_4d)-1
			Vars_out_4d(ivar)			= Vars_4d(ivar)+"_"+tag(itag)+"_4d_mean"
			print(Vars_out_4d(ivar))
			vars_4d_bin_mean(ivar,itag,:,:,:,:)	= fin->$Vars_out_4d(ivar)$

			Vars_out_4d(ivar)			= Vars_4d(ivar)+"_"+tag(itag)+"_4d_std"
			print(Vars_out_4d(ivar))
			vars_4d_bin_std(ivar,itag,:,:,:,:) = fin->$Vars_out_4d(ivar)$
		end do ; do ivar=
	end do ; do itag=

	printVarSummary(vars_4d_bin_mean)
	printMinMax(vars_4d_bin_mean,True)
	
	; ----- get the time average in each bin
	vars_4d_bin_mean_tm		= dim_avg_n_Wrap(vars_4d_bin_mean,3)
	vars_4d_bin_std_tm		= dim_avg_n_Wrap(vars_4d_bin_std,3)
	
	printVarSummary(vars_4d_bin_mean_tm)
	printMinMax(vars_4d_bin_mean_tm,True)
end if

if(ALPIsub_bin_1d_plot)then
	do itag=0,dimsizes(tag)-1
		do ivar=0,dimsizes(Vars_3d)-1
			Vars_out_3d(ivar)			= Vars_3d(ivar)+"_"+tag(itag)+"_3d_mean"
			print(Vars_out_3d(ivar))
		 	vars_3d_bin_mean(ivar,itag,:,:,:) = fin->$Vars_out_3d(ivar)$	

			Vars_out_3d(ivar)			= Vars_3d(ivar)+"_"+tag(itag)+"_3d_std"
			print(Vars_out_3d(ivar))
			vars_3d_bin_std(ivar,itag,:,:,:) = fin->$Vars_out_3d(ivar)$	
		end do ; do ivar=
	end do ; do itag=

; ----- get the time average in each bin
	vars_3d_bin_mean_tm		= dim_avg_n_Wrap(vars_3d_bin_mean,3)
	vars_3d_bin_std_tm		= dim_avg_n_Wrap(vars_3d_bin_std,3)
	
	printVarSummary(vars_3d_bin_mean_tm)
	printMinMax(vars_3d_bin_mean_tm,True)
end if


; ---------------------------------------ALPIsub_mask_plot ------------------------------

if(ALPIsub_mask_plot)then

	; get mask variable to verify the correctness of ALPI bin
	dims			= dimsizes(ALPI(0,0,:,:,:))
	; reshape into Lat-Lon dimension variable
	ALPI_mask_latlon								= reshape(ALPI_mask,(/dimsizes(tag),dimsizes(case_nor),dims(0),dims(1),dims(2)/))
	copy_VarCoords(ALPI,ALPI_mask_latlon(:,:,:,:,:))
	printVarSummary(ALPI_mask_latlon)
	; get the difference to identify the distinction of Gauss-PDF and Park-RH
	ALPI_mask_latlon_diff							= (/ALPI_mask_latlon(:,1,:,:,:)-ALPI_mask_latlon(:,0,:,:,:)/)
	copy_VarCoords(ALPI_mask_latlon(:,0,:,:,:),ALPI_mask_latlon_diff)
	ALPI_mask_latlon_diff2							= (/ALPI_mask_latlon(1,:,:,:,:)-ALPI_mask_latlon(0,:,:,:,:)/)
	copy_VarCoords(ALPI_mask_latlon(0,:,:,:,:),ALPI_mask_latlon_diff2)
	
	
	; ======== the following part is used to verify whether the ALPI index sorting is right. (only plot the 1st time in the array)
	; ======== ---- yes, plots indicate that all regimes are right.
	
	;wks			= gsn_open_wks(plot_type,"./figure/ALPI_mask_plot_nbins_"+nbins+"_Raster")
	wks			= gsn_open_wks(plot_type,"./figure/ALPI_mask_plot_nbins_"+nbins+"_"+date)
	
	;gsn_define_colormap(wks,"NCV_manga")
	;gsn_define_colormap(wks,"NCV_blue_red")
	;gsn_define_colormap(wks,"radar_1")
	gsn_define_colormap(wks,"perc2_9lev")
	
	
	res							= True
	
	;res@gsnSpreadColors		= True
	res@gsnFrame				= False
	res@gsnDraw					= False
	
	res@cnFillOn				= True
	res@cnLinesOn				= False
	res@cnLineLabelsOn			= False
	;res@cnFillMode       		= "RasterFill"
	res@cnLevelSelectionMode	= "ExplicitLevels"
	res@cnLevels         		= ispan(1,nbins-1,1);ALPI_mask_val(1:)
	;res@cnFillColors			= ispan(2,250,250/dimsizes(ALPI_mask_val)) ;(/"lightcyan2","lightgoldenrod2","lightyellow2","lightgrey","lightskyblue2","lightseagreen","lightpink2"/)
	
	res@lbLabelBarOn			= False
	res@lbLabelPosition  		= "Center"           ; label position
	res@lbLabelAlignment 		= "BoxCenters"       ; label orientation
	;res@lbLabelStrings   		= ""+ res@cnLevels   ; trick 
	res@lbLabelStrings   		= ""+ (ALPI_mask_val(:)-0.5) ; --- slightly tricky here: to make the labeled value look clear.
	res@gsnLeftString			= "ALPI regimes"
	
	res@mpCenterLonF			= 180
	res@mpMinLatF				= -30.
	res@mpMaxLatF				= 30.
	
	fontheight					= 0.014
	res@tmXBLabelStride			= 2.
	res@gsnStringFontHeightF	= fontheight
	res@tmXBLabelFontHeightF	= fontheight
	res@tmYLLabelFontHeightF	= fontheight
	
	
	plot_tmp					= gsn_csm_contour_map_ce(wks, ALPI_mask_latlon(0,0,0,:,:),res)
	
	fontheight					= 0.025
	res@lbLabelFontHeightF		= fontheight
	
	res@pmLabelBarHeightF		=  0.100
	res@pmLabelBarWidthF 		=  0.80
	
	res@tmXBLabelStride			= 2.
	res@gsnStringFontHeightF	= fontheight
	res@tmXBLabelFontHeightF	= fontheight
	res@tmYLLabelFontHeightF	= fontheight
	
	; -------- resource for difference plot
	resd						= True
	resd@gsnFrame				= False
	resd@gsnDraw				= False
	
	;resd@cnFillMode       		= "RasterFill"
	
	resd@cnFillOn				= True
	resd@cnLinesOn				= False
	resd@cnLineLabelsOn			= False
	resd@mpCenterLonF			= 180
	resd@mpMinLatF				= -30.
	resd@mpMaxLatF				= 30.
	resd@tmXBLabelFontHeightF	= fontheight/2.
	resd@tmYLLabelFontHeightF	= fontheight/2.
	;resd@pmLabelBarParallelPosF	= 
	resd@pmLabelBarOrthogonalPosF	= 0.25
	;resd@gsnLeftString			= "ALPI regimes"
	resd@cnLevelSelectionMode	= "ManualLevels"
	resd@cnMinLevelValF			= -2.
	resd@cnMaxLevelValF			= 2.
	resd@cnLevelSpacingF		= 0.5
	
	resd@lbLabelBarOn			= False
	
	
	plot						= new((/dimsizes(tag),dimsizes(case_nor)/),graphic)
	plot_diff					= new((/dimsizes(tag)/),graphic)
	plot_diff2					= new((/dimsizes(tag)/),graphic)
	
	count_1	= 0
	count_2	= 0
	count_3	= 0
	do itag=0,dimsizes(tag)-1
	resd@gsnLeftString			= plots_strings(count_1)
	resd@gsnRightString			= tag(itag)+"_"+"PDF-RH"
	plot_diff(itag)				= gsn_csm_contour_map_ce(wks, ALPI_mask_latlon_diff(itag,0,:,:),resd)
	count_1 = count_1 + 1
	
	do icase=0,dimsizes(case_nor)-1
	res@gsnLeftString	= plots_strings(count_2)
	res@gsnRightString = tag(itag)+"_"+case_name(icase)
	plot(itag,icase)		= gsn_csm_contour_map_ce(wks,ALPI_mask_latlon(itag,icase,0,:,:),res)
	count_2	= count_2 + 1
	end do ; do itag=
	end do ; do icase=
	
	do icase=0,dimsizes(case_nor)-1
	resd@gsnLeftString			= plots_strings(count_3)
	resd@gsnRightString			= "P4K-CTL_"+case_name(icase)
	plot_diff2(icase)			= gsn_csm_contour_map_ce(wks, ALPI_mask_latlon_diff2(icase,0,:,:),resd)
	count_3	= count_3 + 1
	end do ; do icase=
	
	pres						= True
	pres@gsnPanelLabelBar		= True
	pres@pmLabelBarWidthF		= 0.8
	pres@pmLabelBarHeightF		= 0.05
	pres@lbLabelFontHeightF		= 0.015
	pres@lbLabelPosition  		= "Center"           ; label position
	pres@lbLabelAlignment 		= "BoxCenters"       ; label orientation
	pres@lbLabelStrings   		= ""+ (ALPI_mask_val(:)-0.5) ; --- slightly tricky here: to make the labeled value look clear.
	
	presd						= True
	presd@gsnPanelLabelBar		= True
	presd@pmLabelBarWidthF		= 0.8
	presd@lbLabelFontHeightF	= 0.015
	
	gsn_panel(wks,plot_tmp,(/1,1/),pres)
	gsn_panel(wks,ndtooned(plot),(/2,2/),pres)
	gsn_panel(wks,plot_diff,(/2,1/),presd)
	gsn_panel(wks,plot_diff2,(/2,1/),presd)
	
	delete([/count_1, count_2, count_3 /])
	exit
end if ; ALPIsub_mask_plot


; ----------------------------------------------------
; -------- ALPI_binned_vars plot
; ----------------------------------------------------
if(ALPIsub_bin_2d_plot)then
	trYMinFs	= (/600,750/)
	
	do itrYMinF = 0,dimsizes(trYMinFs)-1
	
		if(ALPIsub_revised)then
			wks  	= gsn_open_wks (plot_type,"./figure/ALPI_binned_vars_nbins_"+nbins+"_revised_ALPI"+"_"+Vars_4d(dimsizes(Vars_4d)-1)+"_2d_"+trYMinFs(itrYMinF))
		else
			wks  	= gsn_open_wks (plot_type,"./figure/ALPI_binned_vars_nbins_"+nbins+"_2d_"+trYMinFs(itrYMinF))
		end if
	
	
	; -------- for lev-bin variables 
	cmap1 							= read_colormap_file("WhiteBlueGreenYellowRed")
	;cmap2 							= read_colormap_file("BlueRed")
	;cmap2							= read_colormap_file("BlueWhiteOrangeRed")
	cmap2							= read_colormap_file("NCV_blu_red")
	
	
	plot_3d							= new((/dimsizes(Vars_4d),dimsizes(tag)*dimsizes(case_nor)+5/),graphic)
	
	res_con							= True
	res_con@trYReverse				= True
	res_con@cnFillOn				= True
	res_con@cnLinesOn				= False
	;res_con@cnLineLabelsOn			= False
	res_con@tiYAxisString			= "Pressure (hPa)"
	;res_con@tiXAxisString           = "ALPI percentile (%)"
	res_con@gsnFrame				= False
	res_con@gsnDraw					= False
	
	res_con@vpWidthF				= 0.7
	res_con@vpHeightF				= 0.4
	
	res_con@trYMinF					= trYMinFs(itrYMinF)
	
	fontheight						= 0.035
	res_con@tmXBMinorOn				= False
	res_con@tmYLMinorOn				= False
	res_con@tmXTOn					= False
	res_con@tmYROn					= False
	res_con@gsnStringFontHeightF	= fontheight
	res_con@tiXAxisFontHeightF		= fontheight
	res_con@tiYAxisFontHeightF		= fontheight
	res_con@tmXBLabelFontHeightF	= fontheight*1.05
	res_con@tmYLLabelFontHeightF	= fontheight*1.05
	res_con@gsnRightStringFontHeightF	= fontheight*0.5

	;res_con@pmLabelBarOrthogonalPosF	= 0.4
	
	res_con@cnLevelSelectionMode	= "ExplicitLevels"
	
	res_con@tmXBMode				= "Explicit"
	res_con@tmXBValues				= fspan(0,nbins,nbins+1)
	res_con@tmXBLabels				= fspan(0,nbins,nbins+1)*100/nbins
	
	res_cond						= res_con
	res_cond@cnLinesOn				= False
	
	res_con@cnFillPalette       	= cmap1(:240,:)
	res_cond@cnFillPalette			= cmap2
	
	res_cond@gsnSpreadColorStart     = 40
	res_cond@gsnSpreadColorEnd       = 217
	
	vars_4d_bin_mean_tm_casediff	= (/vars_4d_bin_mean_tm(:,:,1,:,:)-vars_4d_bin_mean_tm(:,:,0,:,:)/)
	vars_4d_bin_mean_tm_tagdiff		= (/vars_4d_bin_mean_tm(:,1,:,:,:)-vars_4d_bin_mean_tm(:,0,:,:,:)/)
	vars_4d_bin_mean_tm_diff2		= (/vars_4d_bin_mean_tm_tagdiff(:,1,:,:)-vars_4d_bin_mean_tm_tagdiff(:,0,:,:)/)
	copy_VarCoords(vars_4d_bin_mean_tm(:,:,0,:,:),vars_4d_bin_mean_tm_casediff)
	copy_VarCoords(vars_4d_bin_mean_tm(:,0,:,:,:),vars_4d_bin_mean_tm_tagdiff)
	copy_VarCoords(vars_4d_bin_mean_tm(:,0,0,:,:),vars_4d_bin_mean_tm_diff2)
	
	
	do ivar=0,dimsizes(Vars_4d)-1
		res_con@gsnRightString = Vars_4d(ivar)+" ("+Vars_4d_units(ivar)+")"
		res_cond@gsnRightString = Vars_4d(ivar)+" ("+Vars_4d_units(ivar)+")"
	
	  if(.not.ismissing(ind(vars_pdf.eq.Vars_4d(ivar))))then ; for vars_pdf 
			leftstring						= (/"PDF_CTL",\
												"PDF_P4K",\
												"PDF_P4K-CTL"/)
			print(dimsizes(leftstring))
	
			data	= (/vars_4d_bin_mean_tm(ivar,0,1,{trYMinFs(itrYMinF):},:),\
						vars_4d_bin_mean_tm(ivar,1,1,{trYMinFs(itrYMinF):},:),\
						vars_4d_bin_mean_tm_tagdiff(ivar,1,{trYMinFs(itrYMinF):},:)/)
			copy_VarCoords(vars_4d_bin_mean_tm(ivar,0,0,{trYMinFs(itrYMinF):},:),data(0,:,:))
			printVarSummary(data)
			; --- attention: here the var number should be larger than 1, at least two.
			levels	= (/Vars_4d_levs(ivar,:),\
						Vars_4d_levs(ivar,:),\
						Vars_4d_dlevs_1(ivar,:)/)
	
			count = 0
			do it = 0,dimsizes(leftstring)-1
				print(it)
				if(it.eq.0 .or. it.eq.1)then
;					res_con@gsnLeftString		= plots_strings(count)+" "+leftstring(it)
					res_con@gsnLeftString		= leftstring(it)

					 res_con@cnLevels            = levels(it,:)
					plot_3d(ivar,it)			= gsn_csm_contour(wks,data(it,:,:),res_con)
				else
					if(it.eq.2)then
					symMinMaxPlt(data(it,:,:),15,False,res_cond)
					end if
;					res_cond@gsnLeftString		= plots_strings(count)+" "+leftstring(it)
					res_cond@gsnLeftString		= leftstring(it)

					; res_cond@cnLevels			= levels(it,:)
					plot_3d(ivar,it)			= gsn_csm_contour(wks,data(it,:,:),res_cond)
				end if
			count = count +1
			end do ; do it=
			gsn_panel(wks,plot_3d(ivar,:dimsizes(leftstring)-1),(/1,3/),False)
			delete([/data,leftstring,levels/])
	
		else ; not for vars_pdf
			leftstring						= (/"RH_CTL","PDF_CTL","PDF-RH_CTL",\
												"RH_P4K","PDF_P4K","PDF-RH_P4K",\
												"RH_P4K-CTL","PDF_P4K-CTL","PDF-RH_P4K-CTL"/)
			print(dimsizes(leftstring))
	
			data	= (/vars_4d_bin_mean_tm(ivar,0,0,{trYMinFs(itrYMinF):},:),vars_4d_bin_mean_tm(ivar,0,1,{trYMinFs(itrYMinF):},:),vars_4d_bin_mean_tm_casediff(ivar,0,{trYMinFs(itrYMinF):},:),\
						vars_4d_bin_mean_tm(ivar,1,0,{trYMinFs(itrYMinF):},:),vars_4d_bin_mean_tm(ivar,1,1,{trYMinFs(itrYMinF):},:),vars_4d_bin_mean_tm_casediff(ivar,1,{trYMinFs(itrYMinF):},:),\
						vars_4d_bin_mean_tm_tagdiff(ivar,0,{trYMinFs(itrYMinF):},:),vars_4d_bin_mean_tm_tagdiff(ivar,1,{trYMinFs(itrYMinF):},:),vars_4d_bin_mean_tm_diff2(ivar,{trYMinFs(itrYMinF):},:)\
					  /)
	
			copy_VarCoords(vars_4d_bin_mean_tm(ivar,0,0,{trYMinFs(itrYMinF):},:),data(0,:,:))
			printVarSummary(data)
	
			; --- attention: here the var number should be larger than 1, at least two.
			levels	= (/Vars_4d_levs(ivar,:),Vars_4d_levs(ivar,:),Vars_4d_dlevs_1(ivar,:),\
						Vars_4d_levs(ivar,:),Vars_4d_levs(ivar,:),Vars_4d_dlevs_1(ivar,:),\
						Vars_4d_dlevs_1(ivar,:),Vars_4d_dlevs_1(ivar,:),Vars_4d_dlevs_1(ivar,:)/)
	
			count = 0
			do it = 0,dimsizes(leftstring)-1
				print(it)
				if(it.eq.0 .or. it.eq.1 .or. it.eq.3 .or. it.eq.4)then
;					res_con@gsnLeftString		= plots_strings(count)+" "+leftstring(it)
					res_con@gsnLeftString		= leftstring(it)

					res_con@cnLevels            = levels(it,:)
					plot_3d(ivar,it)			= gsn_csm_contour(wks,data(it,:,:),res_con)
				else
					if (it.eq.2 .or. it.eq.5)then ; case_diff
						symMinMaxPlt(data((/2,5/),:,:),15,False,res_cond)
					end if
					if (it.eq.6 .or. it.eq.7)then ; tag_diff
						symMinMaxPlt(data((/6,7/),:,:),15,False,res_cond)
					end if
					if(it.eq.8)then ; diff2
						symMinMaxPlt(data((/8/),:,:),15,False,res_cond)
					end if
	
;					res_cond@gsnLeftString		= plots_strings(count)+" "+leftstring(it)
					res_cond@gsnLeftString		= leftstring(it)
					;res_cond@cnLevels			= levels(it,:)
					plot_3d(ivar,it)			= gsn_csm_contour(wks,data(it,:,:),res_cond)
				end if
			count = count +1
			end do ; do it=
			gsn_panel(wks,plot_3d(ivar,:),(/3,3/),False)
			delete([/data,leftstring,levels/])
		end if
	end do ; do ivar=

	resP								= True
	resP@gsnPanelFigureStrings			= plots_strings(0:dimsizes(Vars_4d)*3-1)
	resP@amJust   						= "TopLeft"
	;resP@gsnPanelFigureStringsPerimOn	= False

	plots_1 = (/plot_3d(0,6),plot_3d(0,7),plot_3d(0,8),\  ; CLOUD
			    plot_3d(1,6),plot_3d(1,7),plot_3d(1,8),\  ; CLDLIQ
			    plot_3d(2,6),plot_3d(2,7),plot_3d(2,8)/)  ; RELHUM

	gsn_panel(wks,plots_1,(/3,3/),resP)

	resP@gsnPanelFigureStrings			= plots_strings(0:dimsizes(Vars_4d)*3-1)

	plots_2 = (/plot_3d(0,7),\  ; CLOUD
			  	plot_3d(1,7),\  ; CLDLIQ
			  	plot_3d(3,2),\  ; sgm_tota
				plot_3d(4,2)/)  ; N1
	
	gsn_panel(wks,plots_2,(/2,2/),resP)

	plots_3 = (/plot_3d(3,0), plot_3d(3,1), plot_3d(3,2),\ ; sgm_tota
				plot_3d(4,0), plot_3d(4,1), plot_3d(4,2)/) ; N1
	
	gsn_panel(wks,plots_3,(/2,3/),resP)

	exit
	
	delete([/res_con, res_cond, wks, plot_3d/])
	
	end do ; do itrYMinF=
	
	exit
end if ; ALPIsub_bin_2d_plot


if(ALPIsub_bin_1d_plot)then

		;------------------------------------------------------------------------
		; ---- zhangjs---
		; ---- get the common min and max value for some sgm related variables to be used as Y axis --- starts 

		; yspan: include each term in each ilist (sgm, crf) of list[vars_yaxis]
		; 		 /ilist,iterm,2+2 (max, min, diff-max, diff-min)/
		; all_yspan: include the max/min of all terms in each ilist of list[vars_yaxis]
		; 		 /ilist,2+2 (max, min, diff-max, diff-min)/
		yspan 					= NewList("lifo")
		all_yspan 				= NewList("lifo")
		do ilist=0,ListCount(vars_yaxis)-1
			ListAppend(yspan,new((/dimsizes(vars_yaxis[ilist]),2+2/), float))
			ListAppend(all_yspan,new(2+2, float))
		end do
		;print(all_yspan)
			
		do ilist=0,ListCount(vars_yaxis)-1
			do iy=0,dimsizes(vars_yaxis[ilist])-1
				do ivar=0,dimsizes(Vars_3d)-1
					if(Vars_3d(ivar).eq.vars_yaxis[ilist](iy))then
						yspan[ilist](iy,0) 		= max(vars_3d_bin_mean_tm(ivar,:,:,:)+vars_3d_bin_std_tm(ivar,:,:,:)) ; max var
						yspan[ilist](iy,1) 		= min(vars_3d_bin_mean_tm(ivar,:,:,:)-vars_3d_bin_std_tm(ivar,:,:,:)) ; min var
						yspan[ilist](iy,2) 		= max(vars_3d_bin_mean_tm(ivar,1,:,:)-vars_3d_bin_mean_tm(ivar,0,:,:)) ; max casediff
						yspan[ilist](iy,3) 		= min(vars_3d_bin_mean_tm(ivar,1,:,:)-vars_3d_bin_mean_tm(ivar,0,:,:)) ; min casediff
					end if
				end do ; do ivar=
			end do ; do iy=
			all_yspan[ilist](0) 	= max(yspan[ilist](:,0))	
			all_yspan[ilist](1) 	= min(yspan[ilist](:,1))	
			all_yspan[ilist](2) 	= max(yspan[ilist](:,2))	
			all_yspan[ilist](3) 	= min(yspan[ilist](:,3))	
			
			print("ilist="+ilist+"	"+yspan[ilist](:,2))
			print("ilist="+ilist+"	"+yspan[ilist](:,3))
			print(all_yspan[ilist](0)+"		"+all_yspan[ilist](1))
		end do ;ilist
		
		; ---- get the common min and max value for some sgm related variables to be used as Y axis --- ends
		;------------------------------------------------------------------------

		; --------- 2019-07-24 22:13:51 qinyi----- 
		; vars_change_Yaxis: save all variables that need to change Y axis value to make it easier to read.
		; vars_change_Yaxis_units: save all variables' changed units
		; vars_change_Yaxis_scale: save all variables' needed scale to make it easier to read.
		
		vars_change_Yaxis				= (/"sgm_tota","N1","sgm_shal","sgm_turb","term1_shal","term2_shal","term3_shal","CLDLOW"/)
		vars_change_Yaxis_units			= (/"10~S~-2~N~","10~S~-1~N~","10~S~-2~N~","10~S~-2~N~","10~S~-2~N~","10~S~-2~N~","10~S~-2~N~","fraction"/)
		vars_change_Yaxis_scale			= (/1e2, 1e1, 1e2, 1e2, 1e2, 1e2, 1e2, 1/)
		
		vars_change_Yaxis_unitsd		= (/"10~S~-3~N~","10~S~-2~N~","10~S~-3~N~","10~S~-3~N~","10~S~-3~N~","10~S~-3~N~","10~S~-3~N~","10~S~-2~N~"/)
		vars_change_Yaxis_scaled		= (/1e3, 1e2, 1e3, 1e3, 1e3, 1e3, 1e3, 1e2/)
		

		; ---------------------- ploting starts here ------------------------------------

		if(ALPIsub_revised)then
			wks  	= gsn_open_wks (plot_type,"./figure/ALPI_binned_vars_nbins_"+nbins+"_revised_ALPI"+"_"+Vars_3d(dimsizes(Vars_3d)-1)+"_1d")
		else
			wks  	= gsn_open_wks (plot_type,"./figure/ALPI_binned_vars_nbins_"+nbins+"_1d")
		end if

		plot						= new(dimsizes(Vars_3d),graphic)
		error_bar					= new((/dimsizes(Vars_3d),dimsizes(tag)*dimsizes(case_nor),nbins/),graphic)
		upp_bound					= new((/dimsizes(Vars_3d),dimsizes(tag)*dimsizes(case_nor),nbins/),graphic)
		dwn_bound					= new((/dimsizes(Vars_3d),dimsizes(tag)*dimsizes(case_nor),nbins/),graphic)
		centers						= new((/dimsizes(Vars_3d),dimsizes(tag)*dimsizes(case_nor),nbins/),graphic)
		lines						= new((/dimsizes(Vars_3d),dimsizes(tag)*dimsizes(case_nor)/),graphic)
		polygons					= new((/dimsizes(Vars_3d),dimsizes(tag)*dimsizes(case_nor)/),graphic)
		
		plot_diff					= new(dimsizes(Vars_3d),graphic)
		lines_diff					= new((/dimsizes(Vars_3d),dimsizes(case_nor)+1/),graphic)
		lines_ref0					= new((/dimsizes(Vars_3d),dimsizes(case_nor)+1/),graphic)
		centers_diff				= new((/dimsizes(Vars_3d),dimsizes(case_nor)+1,nbins/),graphic)
		
		plot_diff_sgm				= new(dimsizes(Vars_3d),graphic)
		lines_diff_sgm				= new((/dimsizes(Vars_3d),dimsizes(case_nor)+1/),graphic)
		lines_ref0_sgm				= new((/dimsizes(Vars_3d),dimsizes(case_nor)+1/),graphic)
		centers_diff_sgm			= new((/dimsizes(Vars_3d),dimsizes(case_nor)+1,nbins/),graphic)
		
	
	
		res  						= True
		res@gsnDraw					= False
		res@gsnFrame				= False
		
		res@vpWidthF				= 0.7
		res@vpHeightF				= 0.4
		
		res@xyLineColor				= "white" ;(/"red","blue","lightpink","cyan"/)
		res@xyLineThicknessF		= 3.0
		res@tiXAxisString			= "ALPI percentile (%)"
		
		res@tmXBMode				= "Explicit"
		res@tmXBValues				= fspan(0,nbins,nbins+1)
		res@tmXBLabels				= fspan(0,nbins,nbins+1)*100/nbins
		res@trXMaxF					= (nbins)
		
		fontheight						= 0.035
		res@tiYAxisFontHeightF			=fontheight 
		res@tiXAxisFontHeightF			=fontheight 
		res@tmXBLabelFontHeightF		=fontheight 
		res@tmYLLabelFontHeightF		=fontheight 
		
		res@gsnStringFontHeightF	= 0.03
		
		res@tmXTOn					= False
		res@tmYROn					= False
		
		;res@tmYLFormat				= "0*+^e"
		;res@tmYLMinorOn			= False
		;res@tmXBPrecision			= 6
		
		res_diff					= res
	
				
		count_1		= 0
		count_2		= 0
	
		do ivar=0,dimsizes(Vars_3d)-1
			print(Vars_3d(ivar))

			tmp							= reshape(vars_3d_bin_mean_tm(ivar,:,:,:),(/dimsizes(tag)*dimsizes(case_nor),nbins/))
			tmp2						= reshape(vars_3d_bin_std_tm(ivar,:,:,:),(/dimsizes(tag)*dimsizes(case_nor),nbins/))
			
			tmp_diff1					= (/vars_3d_bin_mean_tm(ivar,1,:,:)-vars_3d_bin_mean_tm(ivar,0,:,:)/) ; difference between P4K and CTL
			tmp_diff2					= (/tmp_diff1(1,:)-tmp_diff1(0,:)/) ; difference between PDF_P4K-CTL and RH_P4K-CTL
			tmp_diff_final				= new((/dimsizes(case_nor)+1,nbins/),float)
			tmp_diff_final(0,:)			= tmp_diff1(0,:)
			tmp_diff_final(1,:)			= tmp_diff1(1,:)
			tmp_diff_final(2,:)			= tmp_diff2(:)
	
		
			; --------- 2019-07-24 22:13:57
			iyaxis							= str_match_ind(vars_change_Yaxis,Vars_3d(ivar))

			if(.not.ismissing(iyaxis))then
				tmp						= tmp*vars_change_Yaxis_scale(iyaxis)
				tmp2					= tmp2*vars_change_Yaxis_scale(iyaxis)
				tmp_diff_final			= tmp_diff_final*vars_change_Yaxis_scaled(iyaxis)
			end if
	
			if(.not.ismissing(iyaxis))then
				print(vars_change_Yaxis(iyaxis)+"	"+vars_change_Yaxis_units(iyaxis))
				res_diff@tiYAxisString			= Vars_3d(ivar)+" ("+vars_change_Yaxis_unitsd(iyaxis)+")"
				res@tiYAxisString				= Vars_3d(ivar)+" ("+vars_change_Yaxis_units(iyaxis)+")"
			else
				res_diff@tiYAxisString			= Vars_3d(ivar)+" ("+Vars_3d_units(ivar)+")"
				res@tiYAxisString				= Vars_3d(ivar)+" ("+Vars_3d_units(ivar)+")"
			end if

			; -----------------------special setting for res@ -----------------------
		
			res@gsnLeftString			= plots_strings(count_1)
		
			do ilist=0,0 ;ListCount(vars_yaxis)-1
				if((.not.ismissing(ind(vars_yaxis[ilist].eq.Vars_3d(ivar)))))then
					if(.not.ismissing(iyaxis))then
						res@trYMaxF					= max((/all_yspan[ilist](0)*1.1,all_yspan[ilist](0)*0.9/))*vars_change_Yaxis_scale(iyaxis)
						res@trYMinF					= min((/all_yspan[ilist](1)*1.1,all_yspan[ilist](1)*0.9/))*vars_change_Yaxis_scale(iyaxis)
					else
						res@trYMaxF					= max((/all_yspan[ilist](0)*1.1,all_yspan[ilist](0)*0.9/))
						res@trYMinF					= min((/all_yspan[ilist](1)*1.1,all_yspan[ilist](1)*0.9/))
					end if
					break	; once find it, exit the loop
				else
					res@trYMaxF					= max((/max(tmp+tmp2)*1.2,max(tmp+tmp2)*0.8/))
					res@trYMinF					= min((/min(tmp-tmp2)*1.2,min(tmp-tmp2)*0.8/))
					break 	; once find it, exit the loop
				end if
			end do  ; do ilist=
		
			plot(ivar) 					= gsn_csm_y(wks, tmp, res)
		
			count_1 = count_1 + 1
		
			; -----------------------special setting for res_diff@ -----------------------

			res_diff@gsnLeftString			= plots_strings(count_2)
		
			; --- make the sgm related variables share the same Y axis range
			do ilist=0, ListCount(vars_yaxis)-1
				;print("ilist="+ilist+"	"+all_yspan[ilist](2)+"		"+all_yspan[ilist](3))
				
				if((.not.ismissing(ind(vars_yaxis[ilist].eq.Vars_3d(ivar)))))then
					if(.not.ismissing(iyaxis))then
						res_diff@trYMaxF					= all_yspan[ilist](2)*1.1*vars_change_Yaxis_scaled(iyaxis)
						res_diff@trYMinF					= all_yspan[ilist](3)*1.1*vars_change_Yaxis_scaled(iyaxis)
					else
						res_diff@trYMaxF					= all_yspan[ilist](2)*1.1
						res_diff@trYMinF					= all_yspan[ilist](3)*1.1
					end if
					break ; once find it, exit the loop
				else
					res_diff@trYMaxF					= max(tmp_diff_final(:1,:))*1.1
					res_diff@trYMinF					= min(tmp_diff_final(:1,:))*1.1
					break ; once find it, exit the loop
				end if
			end do ; do ilist=

			plot_diff(ivar)						= gsn_csm_y(wks, tmp_diff_final(:1,:), res_diff)
		
			; -----------------------special setting for res_diff_sgm@ -----------------------

			; resource for sgm variables, because I will overlay them on one plot
			; one switch for sgm difference from sgm_tota, sgm_shal and sgm_turb
			; sgm_diff_plot
			
			res_diff_sgm								= res_diff

			if(sgm_diff_plot)then
				fontheight								= 0.020
				res_diff_sgm@tiYAxisFontHeightF			=fontheight 
				res_diff_sgm@tiXAxisFontHeightF			=fontheight 
				res_diff_sgm@tmXBLabelFontHeightF		=fontheight 
				res_diff_sgm@tmYLLabelFontHeightF		=fontheight 
				res_diff_sgm@gsnStringFontHeightF		=fontheight
				res_diff_sgm@gsnLeftString				= ""
			
				if(ivar.eq.0)then
					res_diff_sgm@gsnLeftString				= spec_lev+" hPa"
					res_diff_sgm@tiYAxisString			= "SGM ("+vars_change_Yaxis_unitsd(iyaxis)+")"
				end if
			else
				res_diff_sgm@gsnLeftString				= ""
			end if
		
		
			if(sgm_diff_plot)then
				plot_diff_sgm(ivar)				= gsn_csm_y(wks, tmp_diff_final(:1,:), res_diff_sgm)
			end if
		
			count_2	= count_2 + 1
		
			; -----------------------add error bar-----------------------
			polyres						= True
			polyres@gsMarkerIndex		= 1
			polyres@gsMarkerSizeF		= 0.02
			
			polyres_sgm					= polyres
		
			colors						= (/"red","blue","red","blue","black","black"/)
			;dashes						= (/0,12,0,12,0,12/)
			dashes						= (/0,0,12,12,0,12/)
			thickes						= (/1.5,1.5,1.5,1.5,3,3/) ;(/2,2,2,2,4,4/)
			
			colors_diff					= (/"red","blue","black"/)
			dashes_diff					= (/0,0,0/)
			thickes_diff				= (/1.5,1.5,3./)
			
										  ;sgm_tota, sgm_shal, sgm_turb
										  ; or
										  ; sgm_shal, term1, term2, term3
;			colors_sgm					= (/"black","red","blue","cyan"/)
			colors_sgm					= (/"black","red","blue",\
											"deeppink","orange","darkorange3"/)
			dashes_sgm					= (/0,12,0,12,0,12,\
											16,12,16,12,16,12/)
			thickes_sgm					= (/2,2,2,2,2,2,\
											2,2,2,2,2,2/)
			
			;polyres@gsLineColor		= (/"red","blue","lightpink","cyan"/)

			; ---------------------------- plot climatology value ---------------------------
			
			; loop through the points
			do idx1=0,dimsizes(tag)*dimsizes(case_nor)-1
			; -- avoid missing value for sgm_tota and N1 to be plotted
				if((.not.ismissing(ind(vars_pdf.eq.Vars_3d(ivar)))) .and. (idx1.eq.0 .or. idx1.eq.2))then
					continue
				else
					; ------------------- add base value line ---------------------------
					polyres@gsLineColor			= colors(idx1)
					polyres@gsLineDashPattern	= dashes(idx1)
					polyres@gsLineThicknessF	= thickes(idx1)
			
					lines(ivar,idx1)				= gsn_add_polyline(wks,plot(ivar),ispan(0,nbins-1,1)+0.5,tmp(idx1,:),polyres)
					
					; ------------------- use gsn_add_polygon to add shaded uncertainty range---------------------------
					
					gsres							= True
					;gsres@tfPolyDrawOrder			= "Predraw"
					gsres@gsFillColor				= colors(idx1)
					gsres@gsFillOpacityF			= 0.1
					
					xp								= new( (/2*nbins/),float )
					yp								= new( (/2*nbins/),float )
					
					do idx2=0,nbins-1
						yp(idx2)					= tmp(idx1,idx2)+tmp2(idx1,idx2)
						xp(idx2)					= idx2+0.5
						yp(2*nbins-1-idx2)			= tmp(idx1,idx2)-tmp2(idx1,idx2)
						xp(2*nbins-1-idx2)			= idx2+0.5
					end do ; do idx2=
			
					;polygons(ivar,idx1)			= gsn_add_polygon(wks,plot(ivar),xp,yp,gsres)

					; ------------------- add error bar, upper and lower bounds ---------------------------
					do idx2=0,nbins-1
						error_bar(ivar,idx1,idx2)	= gsn_add_polyline(wks,plot(ivar),(/idx2,idx2/)+0.5,\
													(/tmp(idx1,idx2)+tmp2(idx1,idx2),tmp(idx1,idx2)-tmp2(idx1,idx2)/),polyres)
						upp_bound(ivar,idx1,idx2)	= gsn_add_polyline(wks,plot(ivar),(/idx2-0.1,idx2+0.1/)+0.5,\
													(/tmp(idx1,idx2)+tmp2(idx1,idx2),tmp(idx1,idx2)+tmp2(idx1,idx2)/),polyres)
						dwn_bound(ivar,idx1,idx2)	= gsn_add_polyline(wks,plot(ivar),(/idx2-0.1,idx2+0.1/)+0.5,\
													(/tmp(idx1,idx2)-tmp2(idx1,idx2),tmp(idx1,idx2)-tmp2(idx1,idx2)/),polyres)
					end do ; do idx2=
					
				end if ; Vars_3d(ivar).eq."sgm_tota"...
			end do ; do idx1=
			
			; ---------------------------- plot difference value ---------------------------

			;do idx1=0,dimsizes(case_nor)
			do idx1=0,dimsizes(case_nor)-1
				if((.not.ismissing(ind(vars_pdf.eq.Vars_3d(ivar)))).and.(idx1.eq.0 .or. idx1.eq.2))then
					continue
				else
					polyres@gsLineColor					= colors_diff(idx1) ;colors(idx1*2)
					polyres@gsLineDashPattern 			= dashes_diff(idx1) ;dashes(idx1*2)
					polyres@gsLineThicknessF			= thickes_diff(idx1) ;thickes(idx1*2)
				
					if(sgm_diff_plot)then
						polyres_sgm@gsLineColor					= colors_sgm(ivar)
						polyres_sgm@gsLineDashPattern 			= dashes_sgm(ivar*idx1*2)
						polyres_sgm@gsLineThicknessF			= thickes_sgm(ivar*idx1*2)
					end if
					
					; ------------------------ plot the climatological difference line --------------------
					lines_diff(ivar,idx1)	= gsn_add_polyline(wks,plot_diff(ivar),ispan(0,nbins-1,1)+0.5,\
																tmp_diff_final(idx1,:),polyres)

					if(sgm_diff_plot)then
						lines_diff_sgm(ivar,idx1)	= gsn_add_polyline(wks,plot_diff_sgm(ivar),ispan(0,nbins-1,1)+0.5,\
																		tmp_diff_final(idx1,:),polyres_sgm)
					end if

					; ------------------------ add line marker --------------------
					do idx2=0,nbins-1
						if(idx1.eq.2)then
						centers_diff(ivar,idx1,idx2)	= gsn_add_polymarker(wks,plot_diff(ivar),idx2+0.5,\
																			tmp_diff_final(idx1,idx2),polyres)
						end if
						if(sgm_diff_plot)then
							centers_diff_sgm(ivar,idx1,idx2)	= gsn_add_polymarker(wks,plot_diff_sgm(ivar),idx2+0.5,\
																			tmp_diff_final(idx1,idx2),polyres_sgm)
						end if
					end do ; do idx2=

					; ------------------------ add zero-line as reference line --------------------
					polyres@gsLineColor						= "grey"
					polyres@gsLineDashPattern				= 14
					
					polyres_sgm@gsLineColor					= "grey"
					polyres_sgm@gsLineDashPattern			= 14
				
					lines_ref0(ivar,idx1)			= gsn_add_polyline(wks,plot_diff(ivar),ispan(0,nbins,1),\
																		fspan(0,0,nbins+1),polyres)
					if(sgm_diff_plot)then
						lines_ref0_sgm(ivar,idx1)	= gsn_add_polyline(wks,plot_diff_sgm(ivar),ispan(0,nbins,1),\
																		fspan(0,0,nbins+1),polyres_sgm)
					end if
				
				end if ; Vars_3d(ivar).eq."sgm_tota"...
			end do ; do idx1=
		
		end do ; do ivar=
		
		if(sgm_diff_plot)then
			if(.not.sgm_shal_4terms)then
				overlay(plot_diff_sgm(0),plot_diff_sgm(1))
				overlay(plot_diff_sgm(0),plot_diff_sgm(2))
			else
				overlay(plot_diff_sgm(0),plot_diff_sgm(1))
				overlay(plot_diff_sgm(0),plot_diff_sgm(2))
				overlay(plot_diff_sgm(0),plot_diff_sgm(3))
				overlay(plot_diff_sgm(0),plot_diff_sgm(4))
				overlay(plot_diff_sgm(0),plot_diff_sgm(5))
			end if
		end if
		
		; ------------------------Add legends -----------------------------------
		;legends                         = vars
		legends							= (/"RH_CTL","PDF_CTL","RH_P4K","PDF_P4K"/)
		;legends_diff1					= (/"RH_P4K-CTL", "PDF_P4K-CTL","PDF-RH_P4K-CTL"/)
		legends_diff1					= (/"RH_P4K-CTL", "PDF_P4K-CTL"/)
		
		if(.not.sgm_shal_4terms)then
			legends_sgm						= (/"SGM_tota", "SGM_shal", "SGM_turb"/)
		else
			;legends_sgm						= (/"SGM_shal","term1_shal","term2_shal","term3_shal"/)
			legends_sgm						= (/"SGM_tota", "SGM_shal","SGM_turb",\
												"term1_shal","term2_shal","term3_shal"/)
		end if
		
		
		lgres							= True
		
		if(.not.sgm_diff_plot)then
			lgres@lgLabelFontHeightF        = 0.10
		else
			lgres@lgLabelFontHeightF        = 0.05
		end if
		
		lgres@vpWidthF                  = 0.25
		lgres@vpHeightF                 = 0.15
		lgres@lgPerimOn                 = False
		
		lgres1							= lgres
		lgres2							= lgres
		lgres_sgm						= lgres
		
		lgres1@lgLineColors				= colors
		lgres1@lgLineThicknessF			= res@xyLineThicknessF
		lgres1@lgDashIndexes			= dashes
		
		lgres2@lgLineColors				= (/"red","blue","black"/) ;colors(::2)
		lgres2@lgLineThicknessF			= res@xyLineThicknessF
		;lgres2@lgDashIndexes			= dashes(::2)
		lgres2@lgDashIndexes			= (/0,0,0/) ;dashes(::2)
		
		
		lgres_sgm@lgLineColors				= colors_sgm;((/1,3,4/))
		lgres_sgm@lgLineThicknessF			= res@xyLineThicknessF
		lgres_sgm@lgDashIndexes				= dashes_sgm(::2)
		
		lbid1                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres1)
		lbid2                           = gsn_create_legend(wks,dimsizes(legends_diff1),legends_diff1,lgres2)
		if(sgm_diff_plot)then
			lbid3                           = gsn_create_legend(wks,dimsizes(legends_sgm),legends_sgm,lgres_sgm)
		end if
		
		amres                           = True
		amres@amJust					= "BottomRight"
		amres@amParallelPosF            = 0.5
		amres@amOrthogonalPosF          = 0.5
		
		amres_sgm						= amres
		amres_sgm@amJust				= "BottomRight"
		amres_sgm@amParallelPosF		= 0.5
		amres_sgm@amOrthogonalPosF		= 0.5
		
		if(Vars_3d(0).eq."PRECT")then
			amres@amJust					= "TopRight"
			amres@amParallelPosF			= 0.5
			amres@amOrthogonalPosF			= -0.5
		end if
		
		annoid1                         = gsn_add_annotation(plot(1),lbid1,amres)
		annoid2                         = gsn_add_annotation(plot_diff(1),lbid2,amres)
		if(sgm_diff_plot)then
		annoid3                         = gsn_add_annotation(plot_diff_sgm(0),lbid3,amres_sgm)
		end if
		
		pres							= True
		pres@gsnPanelYWhiteSpacePercent = 5.0
		pres@gsnPanelXWhiteSpacePercent = 5.0
		pres@gsnPanelLeft				= 0.05
		pres@gsnPanelRight				= 0.95
		
		if(sgm_diff_plot)then
			pres_sgm						= pres
			pres_sgm@amJust					= "TopRight"
			
			tmp_plots						= (/plot_diff_sgm(0)/)
			
			if(dimsizes(tmp_plots).gt.1)then
			pres_sgm@gsnPanelFigureStrings	= (/"a)","b)"/)
			end if
			
			;gsn_panel(wks,(/plot_diff_sgm(0),plot_diff_sgm(1)/),(/1,2/),pres_sgm)
			gsn_panel(wks,(/plot_diff_sgm(0)/),(/1,1/),pres_sgm)
			exit
		end if
		
		;pres@gsnPanelRowSpec = True
		
		;gsn_panel(wks,plot_diff(1:),(/1,3/),pres)
		;gsn_panel(wks,plot(1:),(/1,3/),pres)
		
		if(dimsizes(Vars_3d).le.8)then
			if(mod(dimsizes(Vars_3d),2).ne.0)then
				gsn_panel(wks,plot_diff,(/dimsizes(Vars_3d)/2+1,2/),pres)
				gsn_panel(wks,plot,(/dimsizes(Vars_3d)/2+1,2/),pres)
			else
				gsn_panel(wks,plot_diff,(/dimsizes(Vars_3d)/2,2/),pres)
				gsn_panel(wks,plot,(/dimsizes(Vars_3d)/2,2/),pres)
			end if
		else
			if(mod(dimsizes(Vars_3d),3).ne.0)then
				gsn_panel(wks,plot_diff,(/dimsizes(Vars_3d)/3+1,3/),pres)
				gsn_panel(wks,plot,(/dimsizes(Vars_3d)/3+1,3/),pres)
			else
				gsn_panel(wks,plot_diff,(/dimsizes(Vars_3d)/3,3/),pres)
				gsn_panel(wks,plot,(/dimsizes(Vars_3d)/3,3/),pres)
			end if
		end if
		
		
		delete([/count_1, count_2/])
		exit
end if ; ALPIsub_bin_1d_plot

end if ; ALPI_revised


end

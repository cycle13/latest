
load "./functions_regrid.ncl"
load "./functions_getvar.ncl"
load "./functions_calculate.ncl"

load "./functions_matrix.ncl"

begin

;***************file 1********************

horizon_interp			= False ; used for different input data resolution
resolution				= "2deg"

ALPI_revised			= True ; overall control
ALPI_readout_data		= True
ALPI_mask_ocn			= True
ALPIsub_revised			= True ; from LTS to EIS

ALPIsub_mask_plot		= False
ALPIsub_bin_2d_plot		= False
ALPIsub_bin_1d_plot		= True 
sgm_diff_plot			= True	; plot sgm_diff for sgm_tota, sgm_shal, sgm_turb in one figure by overlay
sgm_shal_4terms			= False  ; term1, term2 and term3 for variance of shallow convection in one figure
								; sgm_shal_4terms should work with sgm_diff_plot = True.
if(sgm_shal_4terms)then
	sgm_diff_plot		= True
end if

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


;----- ALPI_revised
if(ALPI_revised)then
	pick_var_4d			= (/"CLOUD","CLDLIQ","RELHUM","sgm_tota","N1"/)
	pick_var_3d			= (/\
							"netCRF","SWCF","CLDLOW","TGCLDLWP","sgm_tota","N1",\ ; normal
							"sgm_tota","sgm_shal","sgm_turb",\	; sgm_diff_plot
							"sgm_shal","term1_shal","term2_shal","term3_shal"\ ; sgm_shal_4terms
							/)

	use_idx_4d		= new(dimsizes(pick_var_4d),integer)
	use_idx_3d		= new(dimsizes(pick_var_3d),integer)

	do ivar=0,dimsizes(pick_var_4d)-1
		use_idx_4d(ivar)	= ind(varsA.eq.pick_var_4d(ivar))
		print(varsA(ivar))
	end do ; do ivar=
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

	Vars_4d					= varsA(use_idx_4d)
	Vars_4d_levs			= varsA_levs(use_idx_4d,:)
	Vars_4d_dlevs_1			= varsA_dlevs_1(use_idx_4d,:)
	Vars_4d_dlevs_2			= varsA_dlevs_2(use_idx_4d,:)
	Vars_4d_dlevs_sp		= varsA_dlevs_sp(use_idx_4d,:)
	Vars_4d_units			= varsA_units(use_idx_4d)
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

if(ALPI_readout_data)then

; -------------------------read Vars_4D from monthly data set 

yearS				= 1980
yearE				= 1988
ntime				= (yearE-yearS+1)*12

vars_4d				= new((/dimsizes(Vars_4d),dimsizes(tag),dimsizes(case_nor),ntime,nlev,nlat,nlon/),float) ; (/time,lev,lat,lon/)
vars_4d!0			= "var"
vars_4d&var			= ispan(0,dimsizes(Vars_4d)-1,1)
vars_4d!1		 	= "tag"
vars_4d&tag			= fspan(0,dimsizes(tag)-1,dimsizes(tag))
vars_4d!2		 	= "case"
vars_4d&case		= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))

vars_3d				= new((/dimsizes(Vars_3d),dimsizes(tag),dimsizes(case_nor),ntime,nlat,nlon/),float) ; (/time,lat,lon/)
vars_3d!0			= "var"
vars_3d&var			= ispan(0,dimsizes(Vars_3d)-1,1)
vars_3d!1		 	= "tag"
vars_3d&tag			= fspan(0,dimsizes(tag)-1,dimsizes(tag))
vars_3d!2		 	= "case"
vars_3d&case		= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))


vars_mid_4d				= new((/dimsizes(Vars_mid_4d),dimsizes(tag),dimsizes(case_nor),ntime,nlev,nlat,nlon/),float) ; (/time,lev,lat,lon/)
vars_mid_4d!0			= "var"
vars_mid_4d&var			= ispan(0,dimsizes(Vars_mid_4d)-1,1)
vars_mid_4d!1		 	= "tag"
vars_mid_4d&tag			= fspan(0,dimsizes(tag)-1,dimsizes(tag))
vars_mid_4d!2		 	= "case"
vars_mid_4d&case		= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))

vars_mid_3d				= new((/dimsizes(Vars_mid_3d),dimsizes(tag),dimsizes(case_nor),ntime,nlat,nlon/),float) ; (/time,lat,lon/)
vars_mid_3d!0			= "var"
vars_mid_3d&var			= ispan(0,dimsizes(Vars_mid_3d)-1,1)
vars_mid_3d!1		 	= "tag"
vars_mid_3d&tag			= fspan(0,dimsizes(tag)-1,dimsizes(tag))
vars_mid_3d!2		 	= "case"
vars_mid_3d&case		= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))


do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1

fdir 				= "/home/lyl/WORK4/cesm1_2_1/archive/"+case(itag,icase)+"/atm/hist/"

; --- Get list of files from staryear to endyear
years                        = ispan(yearS,yearE,1)

; ------ get monthly data from specific start and end years --------- starts
do iyr=0,dimsizes(years)-1
yr_4d       				= sprinti("%04d",years(iyr))

if(iyr.eq.0)then
    tmpf         			= systemfunc("ls "+fdir+case(itag,icase)+".cam.h0.*"+yr_4d+"*.nc")
    tmpf1        			= tmpf
	delete(tmpf)
	if(iyr.eq.dimsizes(years)-1)then
		alls				= tmpf1
	end if
else
    tmpf         			= systemfunc("ls "+fdir+case(itag,icase)+".cam.h0.*"+yr_4d+"*.nc")
    tmpf2        			= array_append_record(tmpf1,tmpf,0)
	delete(tmpf)
    delete(tmpf1)
    tmpf1        			= tmpf2
    delete(tmpf2)
    if(iyr.eq.dimsizes(years)-1)then
        alls    			= tmpf1
		delete(tmpf1)
    end if
end if
end do ; do iyr=

all_files       			= alls
delete(alls)

print(all_files)
fall      					= addfiles (all_files, "r")   ; note the "s" of addfile
ListSetType (fall, "cat")             ; concatenate or "merge" (default)
printVarSummary(fall)

; get lev
lev							= fall[0]->lev
printVarSummary(lev)

tmp							= fall[0]->OCNFRAC(0,:,:)
if(horizon_interp)then
	OCNFRAC					= linint2_Wrap(tmp&lon,tmp&lat,tmp,True,klon,klat,0)
else
	OCNFRAC					= tmp
end if
delete(tmp)

; --- read vars_4d
do ivar=0,dimsizes(Vars_4d)-1
print(Vars_4d(ivar))
	if(isfilevar(fall[0],Vars_4d(ivar)))then
		tmp    						= fall[:]->$Vars_4d(ivar)$(:,:,:,:)                    ; note syntax [:]
		tmp@_FillValue  			= default_fillvalue("float")
		print("read "+Vars_4d(ivar)+" well;")
		if(horizon_interp)then
			vars_4d(ivar,itag,icase,:,:,:,:)		= linint2_Wrap(tmp&lon,tmp&lat,tmp,True,klon,klat,0)
			print("regrid "+Vars_4d(ivar)+" well;")
		else
			vars_4d(ivar,itag,icase,:,:,:,:)		= tmp
		end if
		if(Vars_4d(ivar).eq."sgm_tota" .or. Vars_4d(ivar).eq."sgm_shal" .or. Vars_4d(ivar).eq."sgm_turb")then
			vars_4d(ivar,itag,icase,:,:,:,:)	= abs(vars_4d(ivar,itag,icase,:,:,:,:))*1e6
		end if
;		if(Vars_4d(ivar).eq."CLOUD")then
;			vars_4d(ivar,itag,icase,:,:,:,:)	= vars_4d(ivar,itag,icase,:,:,:,:)*100.
;		end if
		if(Vars_4d(ivar).eq."CLDLIQ")then
			vars_4d(ivar,itag,icase,:,:,:,:)		= vars_4d(ivar,itag,icase,:,:,:,:)*1000.
		end if
		delete(tmp)
	else
		if(.not.ismissing(ind(vars_pdf.eq.Vars_4d(ivar))) .and. icase.eq.0)then ; make not existing vars in park-rh as missing value.
			print(Vars_4d(ivar))
			vars_4d(ivar,itag,icase,:,:,:,:)	= vars_4d@_FillValue
		end if
	end if
end do ; do ivar=

; --- read vars_3d
do ivar=0,dimsizes(Vars_3d)-1
	print(Vars_3d(ivar))
	if(isfilevar(fall[0],Vars_3d(ivar))) then
		tmp									= fall[:]->$Vars_3d(ivar)$
		tmp@_FillValue  					= default_fillvalue("float")
		if(dimsizes(dimsizes(tmp)).eq.4)then ; figure out variables having vertical level coordinate
			tmp1							= fall[:]->$Vars_3d(ivar)$(:,{spec_lev},:,:)
			print("read "+Vars_3d(ivar)+" well;")
			if(horizon_interp)then
				vars_3d(ivar,itag,icase,:,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
			else
				vars_3d(ivar,itag,icase,:,:,:)	= tmp1
			end if
			delete(tmp1)
		else
			tmp1						 	= fall[:]->$Vars_3d(ivar)$
			print("read "+Vars_3d(ivar)+" well;")
			if(horizon_interp)then
				vars_3d(ivar,itag,icase,:,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
			else
				vars_3d(ivar,itag,icase,:,:,:)	= tmp1
			end if
			delete(tmp1)
		end if
		delete(tmp)
		if(Vars_3d(ivar).eq."TGCLDLWP")then
			vars_3d(ivar,itag,icase,:,:,:)	= vars_3d(ivar,itag,icase,:,:,:)*1e3
		end if
		if(Vars_3d(ivar).eq."sgm_tota" .or. Vars_3d(ivar).eq."sgm_shal" .or. Vars_3d(ivar).eq."sgm_turb")then
			vars_3d(ivar,itag,icase,:,:,:)	= abs(vars_3d(ivar,itag,icase,:,:,:))*1e6
		end if
	else
		if(.not.ismissing(ind(vars_pdf.eq.Vars_3d(ivar))) .and. icase.eq.0)then ; make not existing vars in park-rh as missing value.
			print(Vars_3d(ivar))
			vars_3d(ivar,itag,icase,:,:,:)	= vars_4d@_FillValue
		else
			if(Vars_3d(ivar).eq."netCRF")then
				SWCF						= fall[:]->SWCF(:,:,:)
				LWCF						= fall[:]->LWCF(:,:,:)
				tmp1						= (/SWCF+LWCF/)
				copy_VarCoords(SWCF,tmp1)
				if(horizon_interp)then
					vars_3d(ivar,itag,icase,:,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
				else
					vars_3d(ivar,itag,icase,:,:,:)   = tmp1
				end if
				delete(tmp1)
				delete([/SWCF,LWCF/])
			end if
			if(Vars_3d(ivar).eq."term1_shal")then
				umf_shal_tmp			= fall[:]->umf_shal(:,{spec_lev},:,:)
				qwu_qw_tmp				= fall[:]->qwu_qw(:,{spec_lev},:,:)
				clddep2_tmp				= fall[:]->clddep2(:,:,:)
				wstarPBL_tmp			= fall[:]->wstarPBL(:,:,:)
				tmp1					= -0.5*umf_shal_tmp*qwu_qw_tmp*clddep2_tmp/where(wstarPBL_tmp.ne.0,wstarPBL_tmp,vars_4d@_FillValue)
				tmp1					= tmp1*1e6
				copy_VarCoords(umf_shal_tmp,tmp1)
				if(horizon_interp)then
					vars_3d(ivar,itag,icase,:,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
				else
					vars_3d(ivar,itag,icase,:,:,:)	= tmp1
				end if
				delete(tmp1)
				delete([/umf_shal_tmp,qwu_qw_tmp,clddep2_tmp,wstarPBL_tmp/])
			end if
			if(Vars_3d(ivar).eq."term2_shal")then
				umf_shal_tmp			= fall[:]->umf_shal(:,{spec_lev},:,:)
				thlu_thl_tmp			= fall[:]->thlu_thl(:,{spec_lev},:,:)
				clddep2_tmp				= fall[:]->clddep2(:,:,:)
				wstarPBL_tmp			= fall[:]->wstarPBL(:,:,:)
				tmp1					= -0.5*umf_shal_tmp*thlu_thl_tmp*clddep2_tmp/where(wstarPBL_tmp.ne.0,wstarPBL_tmp,vars_4d@_FillValue)
				tmp1					= tmp1*1e6
				copy_VarCoords(umf_shal_tmp,tmp1)
				if(horizon_interp)then
					vars_3d(ivar,itag,icase,:,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
				else
					vars_3d(ivar,itag,icase,:,:,:)	= tmp1
				end if 
				delete(tmp1)
				delete([/umf_shal_tmp,thlu_thl_tmp,clddep2_tmp,wstarPBL_tmp/])
			end if
			if(Vars_3d(ivar).eq."term3_shal")then
				umf_shal_tmp			= fall[:]->umf_shal(:,{spec_lev},:,:)
				qwu_thlu_tmp			= fall[:]->qwu_thlu(:,{spec_lev},:,:)
				clddep2_tmp				= fall[:]->clddep2(:,:,:)
				wstarPBL_tmp			= fall[:]->wstarPBL(:,:,:)
				tmp1					= -0.5*umf_shal_tmp*qwu_thlu_tmp*clddep2_tmp/where(wstarPBL_tmp.ne.0,wstarPBL_tmp,vars_4d@_FillValue)
				tmp1					= tmp1*1e6
				copy_VarCoords(umf_shal_tmp,tmp1)
				if(horizon_interp)then
					vars_3d(ivar,itag,icase,:,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
				else
					vars_3d(ivar,itag,icase,:,:,:)	= tmp1
				end if
				delete(tmp1)
				delete([/umf_shal_tmp,qwu_thlu_tmp,clddep2_tmp,wstarPBL_tmp/])
			end if
		end if
	end if
end do ; do ivar=

; --- read vars_mid_4d
do ivar=0,dimsizes(Vars_mid_4d)-1
print(Vars_mid_4d(ivar))
	if(isfilevar(fall[0],Vars_mid_4d(ivar)))then
		tmp    						= fall[:]->$Vars_mid_4d(ivar)$(:,:,:,:)                    ; note syntax [:]
		tmp@_FillValue  			= default_fillvalue("float")
		print("read "+Vars_mid_4d(ivar)+" well;")
		if(horizon_interp)then
			vars_mid_4d(ivar,itag,icase,:,:,:,:)		= linint2_Wrap(tmp&lon,tmp&lat,tmp,True,klon,klat,0)
			print("regrid "+Vars_mid_4d(ivar)+" well;")
		else
			vars_mid_4d(ivar,itag,icase,:,:,:,:)		= tmp
		end if
		delete(tmp)
	end if
end do ; do ivar=

; --- read vars_mid_3d
do ivar=0,dimsizes(Vars_mid_3d)-1
print(Vars_mid_3d(ivar))
	if(isfilevar(fall[0],Vars_mid_3d(ivar)))then
		tmp    						= fall[:]->$Vars_mid_3d(ivar)$(:,:,:)                    ; (time,lat,lon)
		tmp@_FillValue  			= default_fillvalue("float")
		print("read "+Vars_mid_3d(ivar)+" well;")
		if(horizon_interp)then
			vars_mid_3d(ivar,itag,icase,:,:,:)		= linint2_Wrap(tmp&lon,tmp&lat,tmp,True,klon,klat,0)
			print("regrid "+Vars_mid_3d(ivar)+" well;")
		else
			vars_mid_3d(ivar,itag,icase,:,:,:)		= tmp
		end if
		delete(tmp)
	else
		if(Vars_mid_3d(ivar).eq."PRECT")then
			precl_tmp				= fall[:]->PRECL(:,:,:)
			precc_tmp				= fall[:]->PRECC(:,:,:)
			tmp1					= precl_tmp+precc_tmp
			copy_VarCoords(precl_tmp,tmp1)
			if(horizon_interp)then
				vars_mid_3d(ivar,itag,icase,:,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
				print("regrid "+Vars_mid_3d(ivar)+" well;")
			else
				vars_mid_3d(ivar,itag,icase,:,:,:)   = tmp1
			end if
			delete(tmp1)
			vars_mid_3d(ivar,itag,icase,:,:,:)	= vars_mid_3d(ivar,itag,icase,:,:,:)*8.64e7
			delete(precl_tmp)
			delete(precc_tmp)
		end if
	end if
end do ; do ivar=


end do ; do icase=
end do ; do itag=

printVarSummary(vars_4d)
printMinMax(vars_4d,True)
printVarSummary(vars_3d)
printMinMax(vars_3d,True)
printVarSummary(vars_mid_4d)
printMinMax(vars_mid_4d,True)
printVarSummary(vars_mid_3d)
printMinMax(vars_mid_3d,True)

;=================================================
; ALPI IS DEFINE BY LTS AND PRECIPITATION.
;if(ALPI_bin_cal)then

;vars_mid_4d: [var | 4] x [tag | 2] x [case | 2] x [time | 12] x [lev | 30] x [lat | 96] x [lon | 144]

; ----------- 1.1 calculate LTS
LTS					= new((/dimsizes(tag),dimsizes(case_nor),ntime,nlat,nlon/),float) ; unit is degreeC
LTS!0		 		= "tag"
LTS&tag				= fspan(0,dimsizes(tag)-1,dimsizes(tag))
LTS!1		 		= "case"
LTS&case			= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))

EIS					= LTS

do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1
do itime=0,ntime-1
result							= cal_eis(vars_mid_4d(0,itag,icase,itime,:,:,:),vars_mid_4d(1,itag,icase,itime,:,:,:),lev)
LTS(itag,icase,itime,:,:)		= tofloat(result[0])
EIS(itag,icase,itime,:,:)		= tofloat(result[1])
end do ; do itime=
end do ; do icase=
end do ; do itag=

copy_VarCoords(vars_mid_4d(var|0,tag|:,case|:,time|:,lev|0,lat|:,lon|:),LTS)
copy_VarCoords(vars_mid_4d(var|0,tag|:,case|:,time|:,lev|0,lat|:,lon|:),EIS)

printVarSummary(LTS)
printMinMax(LTS,True)
printVarSummary(EIS)
printMinMax(EIS,True)

; ---------- 1.2 get precipitation
;--- separate PRECT from vars_mid_3d to avoid the plotting difficulty if I dont want to plot PRECT.
do ivar=0,dimsizes(Vars_mid_3d)-1
	if(Vars_mid_3d(ivar).eq."PRECT")then
		print(Vars_mid_3d(ivar))
		prect					= vars_mid_3d(var|ivar,tag|:,case|:,time|:,lat|:,lon|:)
		break
	end if
end do ; do ivar=
printVarSummary(prect)
printMinMax(prect,True)

; --- get OCNFRAC (time,lat,lon)
OCNFRAC_3d	= conform_dims(dimsizes(LTS(tag|0,case|0,time|:,lat|:,lon|:)),OCNFRAC,(/1,2/))
copy_VarCoords(LTS(tag|0,case|0,time|:,lat|:,lon|:),OCNFRAC_3d)
printVarSummary(OCNFRAC_3d)
printMinMax(OCNFRAC_3d,True)

if(ALPI_mask_ocn)then
	; addition: mask land region
	LTS_ocn			= LTS
	EIS_ocn			= EIS
	prect_ocn		= prect
	do itag=0,dimsizes(tag)-1
	do icase=0,dimsizes(case_nor)-1
		LTS_ocn(itag,icase,:,:,:)		= where(OCNFRAC_3d.gt.0.5, LTS(itag,icase,:,:,:),vars_4d@_FillValue)
		EIS_ocn(itag,icase,:,:,:)		= where(OCNFRAC_3d.gt.0.5, EIS(itag,icase,:,:,:),vars_4d@_FillValue)
		prect_ocn(itag,icase,:,:,:)		= where(OCNFRAC_3d.gt.0.5, prect(itag,icase,:,:,:),vars_4d@_FillValue)
	end do ; do icase=
	end do ; do itag=
	copy_VarCoords(LTS,LTS_ocn)
	copy_VarCoords(EIS,EIS_ocn)
	copy_VarCoords(prect,prect_ocn)
else
	LTS_ocn			= LTS
	EIS_ocn			= EIS 
	prect_ocn		= prect
end if ; ALPI_mask_ocn

LTS_tropics	 	= LTS_ocn(:,:,:,{latS_trop:latE_trop},:)
EIS_tropics	 	= EIS_ocn(:,:,:,{latS_trop:latE_trop},:)
prect_tropics	 = prect_ocn(:,:,:,{latS_trop:latE_trop},:)

printVarSummary(LTS_tropics)
printVarSummary(EIS_tropics)
printVarSummary(prect_tropics)

; define ALPI array
ALPI				= LTS_tropics

; --------- 1.3 get LTS_distance and Prec_distance
do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1
do itime=0,ntime-1
	;---1.3.1.1 find LTS_max and LTS_min
	LTS_max				= max(LTS_tropics(itag,icase,itime,:,:))
	LTS_min				= min(LTS_tropics(itag,icase,itime,:,:))
	;---1.3.1.2 find EIS_max and EIS_min
	EIS_max				= max(EIS_tropics(itag,icase,itime,:,:))
	EIS_min				= min(EIS_tropics(itag,icase,itime,:,:))
	;---1.3.2 find prec_max
	prect_max			= max(prect_tropics(itag,icase,itime,:,:))
	prect_min			= min(prect_tropics(itag,icase,itime,:,:))
	
	;---1.3.3 LTS_distance
	LTS_max_2D				= conform_dims(dimsizes(LTS_tropics(itag,icase,itime,:,:)),LTS_max,-1)
	copy_VarCoords(LTS_tropics(itag,icase,itime,:,:),LTS_max_2D)
	LTS_min_2D				= conform_dims(dimsizes(LTS_tropics(itag,icase,itime,:,:)),LTS_min,-1)
	copy_VarCoords(LTS_tropics(itag,icase,itime,:,:),LTS_min_2D)
	LTS_distance			= (LTS_max_2D-LTS_tropics(itag,icase,itime,:,:))/(LTS_max_2D-LTS_min_2D)+0.1
	printMinMax(LTS_distance,True)
	
	;---1.3.3.1 EIS_distance
	EIS_max_2D				= conform_dims(dimsizes(EIS_tropics(itag,icase,itime,:,:)),EIS_max,-1)
	copy_VarCoords(EIS_tropics(itag,icase,itime,:,:),EIS_max_2D)
	EIS_min_2D				= conform_dims(dimsizes(EIS_tropics(itag,icase,itime,:,:)),EIS_min,-1)
	copy_VarCoords(EIS_tropics(itag,icase,itime,:,:),EIS_min_2D)
	EIS_distance			= (EIS_max_2D-EIS_tropics(itag,icase,itime,:,:))/(EIS_max_2D-EIS_min_2D)+0.1
	printMinMax(EIS_distance,True)
	
	;---1.3.4 Prec_distance
	prect_max_2D				= conform_dims(dimsizes(prect_tropics(itag,icase,itime,:,:)),prect_max,-1)
	copy_VarCoords(prect_tropics(itag,icase,itime,:,:),prect_max_2D)
	prect_min_2D				= conform_dims(dimsizes(prect_tropics(itag,icase,itime,:,:)),prect_min,-1)
	copy_VarCoords(prect_tropics(itag,icase,itime,:,:),prect_min_2D)
	prect_distance				= (prect_max_2D-prect_tropics(itag,icase,itime,:,:))/(prect_max_2D-prect_min_2D)+0.1
	printMinMax(prect_distance,True)
	
	;------------ 1.4 get ALPI
	if(ALPIsub_revised)then
		ALPI(itag,icase,itime,:,:)		= atan2(prect_distance,EIS_distance)
	else
		ALPI(itag,icase,itime,:,:)		= atan2(prect_distance,LTS_distance)
	end if

end do ; do itime=
end do ; do icase=
end do ; do itag=

r2d							= 45.0/atan(1.0)     ; conversion factor (radians to degrees)
ALPI						= ALPI*r2d

printVarSummary(ALPI) ; (tag, case, lat, lon)
printMinMax(ALPI,True)

; =========== bin based on ALPI percentile
; ----- 1.1 sorting based on ALPI index: increasing

full_length					= dimsizes(ndtooned(ALPI(0,0,0,{latS_trop:latE_trop},:)))
nlength						= num(.not.ismissing(ndtooned(ALPI(0,0,0,{latS_trop:latE_trop},:))))
ALPI_bins					= new((/dimsizes(tag),dimsizes(case_nor),ntime,nlev,nbins/),float)
vars_4d_bin					= new((/dimsizes(Vars_4d),dimsizes(tag),dimsizes(case_nor),ntime,nlev,full_length,nbins/),float)
vars_4d_bin_mean			= new((/dimsizes(Vars_4d),dimsizes(tag),dimsizes(case_nor),ntime,nlev,nbins/),float)
vars_4d_bin_std				= new((/dimsizes(Vars_4d),dimsizes(tag),dimsizes(case_nor),ntime,nlev,nbins/),float)
vars_3d_bin					= new((/dimsizes(Vars_3d),dimsizes(tag),dimsizes(case_nor),ntime,full_length,nbins/),float)
vars_3d_bin_mean			= new((/dimsizes(Vars_3d),dimsizes(tag),dimsizes(case_nor),ntime,nbins/),float)
vars_3d_bin_std				= new((/dimsizes(Vars_3d),dimsizes(tag),dimsizes(case_nor),ntime,nbins/),float)


vars_4d_bin!0				= "var"
vars_4d_bin&var				= ispan(0,dimsizes(Vars_4d)-1,1)
vars_4d_bin!1				= "tag"
vars_4d_bin&tag				= ispan(0,dimsizes(tag)-1,1)
vars_4d_bin!2				= "case"
vars_4d_bin&case			= ispan(0,dimsizes(case_nor)-1,1)
vars_4d_bin!3				= "time"
vars_4d_bin&time			= vars_4d&time
vars_4d_bin!4				= "lev"
vars_4d_bin&lev				= vars_4d&lev
vars_4d_bin!6				= "bin"
vars_4d_bin&bin				= ispan(1,nbins,1)

copy_VarCoords(vars_4d_bin(:,:,:,:,:,0,:),vars_4d_bin_mean)
copy_VarCoords(vars_4d_bin(:,:,:,:,:,0,:),vars_4d_bin_std)

vars_3d_bin!0				= "var"
vars_3d_bin&var				= ispan(0,dimsizes(Vars_3d)-1,1)

copy_VarCoords(vars_4d_bin(0,:,:,:,0,:,:),vars_3d_bin(0,:,:,:,:,:))
copy_VarCoords(vars_3d_bin(:,:,:,:,0,:),vars_3d_bin_mean)
copy_VarCoords(vars_3d_bin(:,:,:,:,0,:),vars_3d_bin_std)


ALPI_mask					= new((/dimsizes(tag),dimsizes(case_nor),ntime,full_length/),float)

print(nlength)

do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1
do itime=0,ntime-1
	ALPI_1D						= ndtooned(ALPI(itag,icase,itime,{latS_trop:latE_trop},:))
	ALPI_1D_sorted				= ALPI_1D
	ip							= dim_pqsort(ALPI_1D_sorted,2) ; sorting and change the raw data
	printVarSummary(ALPI_1D_sorted)

	; ----for Vars_3d ----- 
	do ivar=0,dimsizes(Vars_3d)-1
		vars_3d_1D					= ndtooned(vars_3d(ivar,itag,icase,itime,{latS_trop:latE_trop},:))

		lat							= vars_3d&lat
		clat                		= tofloat(cos(0.01745329*lat))
		clat!0              		= "lat"
		clat&lat            		= vars_3d&lat
		clat_2D						= conform_dims(dimsizes(vars_3d(ivar,itag,icase,itime,:,:)),clat,0)
		copy_VarCoords(vars_3d(ivar,itag,icase,itime,:,:),clat_2D)
		clat_2D_1D					= ndtooned(clat_2D({latS_trop:latE_trop},:))

		do ibin=0,nbins-1
			print("itag="+itag+"	icase="+icase+"	itime="+itime+"	ivar="+ivar+"	ibin="+ibin)
			tmp1		= ALPI_1D_sorted(toint(nlength/nbins)*ibin)
			; avoid missing some points at the end
			; notion: pls dont write "tmp2 = ALPI_1D(dimsizes(ALPI_1D)-1)" when missing value exists.
			if(ibin.eq.(nbins-1))then
				tmp2		= ALPI_1D_sorted(nlength-1)
			else
				tmp2		= ALPI_1D_sorted(toint(nlength/nbins)*(ibin+1))
			end if
			vars_3d_bin(ivar,itag,icase,itime,:,ibin)	= where(ALPI_1D_sorted.gt.tmp1 .and. ALPI_1D_sorted.lt.tmp2, vars_3d_1D(ip), vars_3d@_FillValue)
			vars_3d_bin_mean(ivar,itag,icase,itime,ibin)		= dim_avg_wgt_n_Wrap(vars_3d_bin(ivar,itag,icase,itime,:,ibin),clat_2D_1D(ip),1,0) 	; the weighting effect is small. 2018.12.18
			vars_3d_bin_std(ivar,itag,icase,itime,ibin)		= dim_stddev(vars_3d_bin(ivar,itag,icase,itime,:,ibin))
	
	; ========= very brillient idea suggested by Jishi! in order to reverse to the Lat-Lon plot, we need to keep the ip index always.
			ALPI_mask(itag,icase,itime,ip(ind(ALPI_1D_sorted.gt.tmp1 .and. ALPI_1D_sorted.lt.tmp2)))	=  ALPI_mask_val(ibin) 
	; ========
		end do ; do ibin=
		delete([/tmp1,tmp2/])
	end do ; do ivar=
	
	; ----for Vars_4d ----- 
	do ivar=0,dimsizes(Vars_4d)-1
		do ilev=0,nlev-1
			vars_4d_1D					= ndtooned(vars_4d(ivar,itag,icase,itime,ilev,{latS_trop:latE_trop},:))
			printVarSummary(vars_4d_1D)
			lat							= vars_4d&lat
			clat                		= tofloat(cos(0.01745329*lat))
			clat!0              		= "lat"
			clat&lat            		= vars_4d&lat
			clat_2D						= conform_dims(dimsizes(vars_4d(ivar,itag,icase,itime,ilev,:,:)),clat,0)
			copy_VarCoords(vars_4d(ivar,itag,icase,itime,ilev,:,:),clat_2D)
			clat_2D_1D					= ndtooned(clat_2D({latS_trop:latE_trop},:))
			;	print(ALPI_1D_sorted+"		"+clat_2D_1D(ip))
	
			do ibin=0,nbins-1
				print("itag="+itag+"	icase="+icase+"	itime="+itime+"	ilev="+ilev+"	ivar="+ivar+"	ibin="+ibin)
				tmp1		= ALPI_1D_sorted(toint(nlength/nbins)*ibin)
				; avoid missing some points at the end
				; notion: pls dont write "tmp2 = ALPI_1D(dimsizes(ALPI_1D)-1)" when missing value exists.
				if(ibin.eq.(nbins-1))then 
					tmp2		= ALPI_1D_sorted(nlength-1)
				else
					tmp2		= ALPI_1D_sorted(toint(nlength/nbins)*(ibin+1))
				end if
				vars_4d_bin(ivar,itag,icase,itime,ilev,:,ibin)	= \
					where(ALPI_1D_sorted.gt.tmp1 .and. ALPI_1D_sorted.lt.tmp2, vars_4d_1D(ip), vars_4d@_FillValue)
				;print(num(.not.ismissing(vars_4d_bin(ivar,itag,icase,itime,ilev,:,ibin)))) ; 1366
	
				vars_4d_bin_mean(ivar,itag,icase,itime,ilev,ibin)	=\
					dim_avg_wgt_n_Wrap(vars_4d_bin(ivar,itag,icase,itime,ilev,:,ibin),clat_2D_1D(ip),1,0) 	; the weighting effect is small. 2018.12.18
				vars_4d_bin_std(ivar,itag,icase,itime,ilev,ibin)	=\
					dim_stddev(vars_4d_bin(ivar,itag,icase,itime,ilev,:,ibin))
	
				; ======== very brillient idea suggested by Jishi! in order to reverse to the Lat-Lon plot, we need to keep the ip index always.
				ALPI_mask(itag,icase,itime,ip(ind(ALPI_1D_sorted.gt.tmp1 .and. ALPI_1D_sorted.lt.tmp2)))	=  ALPI_mask_val(ibin) 
				; =======
			end do ; do ibin=
	
			;if(ilev.eq.21)then
			;	; -------- used to test the binned variables
			;	wks						= gsn_open_wks(plot_type,"./figure/ALPI_binned_variable_test")
			;	
			;	dims					= dimsizes(vars_4d_bin(ivar,itag,icase,itime,ilev,:,:))
			;	tmps					= transpose(vars_4d_bin(ivar,itag,icase,itime,ilev,:,:))
			;	
			;	res						= True
			;	res@vpWidthF			= 0.7
			;	res@vpHeightF			= 0.4
			;	res@xyLineColors		= ispan(2,250,250/dims(1))
			;	res@tiYAxisString		= Vars_4d(ivar)
			;	res@tiXAxisString		= "ndtooned Vars_4d [Lat-Lon]"
			;	
			;	plot					= gsn_csm_y(wks,tmps,res)
			;	delete([/plot,tmps,dims/])
			;	; ---------- 
			;end if
	
	delete([/tmp1,tmp2/])
	end do ; do ilev=
end do ; do ivar=

end do ; do itime=
end do ; do icase=
end do ; do itag=


printVarSummary(vars_4d_bin)
printMinMax(vars_4d_bin,True)
printVarSummary(vars_4d_bin_mean)
printMinMax(vars_4d_bin_mean,True)
printVarSummary(vars_3d_bin_mean)
printMinMax(vars_3d_bin_mean,True)

; ----------------- try to output the binned variables for convenience -----------------------------------

outdir						= "/home/share3/lyl/work3/qinyi/mid-data/cfmip/"
outfile						= "ALPI_binned_vars.nc"

Vars_out_4d					= Vars_4d
Vars_out_3d					= Vars_3d

system("/bin/rm -f "+outdir+outfile)
fout						= addfile(outdir+outfile,"c")

; ---- output these two variables to faciliate the ALPIsub_mask_plot
fout->ALPI					= ALPI
fout->ALPI_mask				= ALPI_mask

do itag=0,dimsizes(tag)-1
	do ivar=0,dimsizes(Vars_4d)-1
		Vars_out_4d(ivar)			= Vars_4d(ivar)+"_"+tag(itag)+"_4d_mean"
		print(Vars_out_4d(ivar))
		fout->$Vars_out_4d(ivar)$		= vars_4d_bin_mean(ivar,itag,:,:,:,:)

		Vars_out_4d(ivar)			= Vars_4d(ivar)+"_"+tag(itag)+"_4d_std"
		print(Vars_out_4d(ivar))
		fout->$Vars_out_4d(ivar)$		= vars_4d_bin_std(ivar,itag,:,:,:,:)
	end do ; do ivar=

	do ivar=0,dimsizes(Vars_3d)-1
		Vars_out_3d(ivar)			= Vars_3d(ivar)+"_"+tag(itag)+"_3d_mean"
		print(Vars_out_3d(ivar))
		fout->$Vars_out_3d(ivar)$		= vars_3d_bin_mean(ivar,itag,:,:,:)

		Vars_out_3d(ivar)			= Vars_3d(ivar)+"_"+tag(itag)+"_3d_std"
		print(Vars_out_3d(ivar))
		fout->$Vars_out_3d(ivar)$		= vars_3d_bin_std(ivar,itag,:,:,:)
	end do ; do ivar=
end do ; do itag=

exit
end if ; ALPI_readout_data

end if ; ALPI_revised


end

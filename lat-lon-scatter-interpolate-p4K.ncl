
load "./functions_regrid.ncl"
load "./functions_getvar.ncl"
load "./functions_calculate.ncl"

begin

;***************file 1********************

;fdir 				= "/home/lyl/WORK4/cesm1_2_1/archive/"
;fdir_append			= "/atm/hist/"

pattern_correlation		= False
regime_contribution 	= False
regime_mask_plot		= False ; this one depends on the regime_contribution
lat_lon_plot			= False
Sherwood_2014			= False

ALPI_bin_cal			= True;True 
ALPIsub_revised			= True;True ; from LTS to EIS
ALPIsub_mask_plot		= False;False
ALPIsub_bin_plot		= True;True 

scatter_p2p				= False ; True ; -- some arrays are too small. only for two variables now.
scatter_sub_ocn			= False ; True	; --- mask ocean or not
scatter_sub_one_var		= False ; False	; --- one variable in X and Y axis -- only two variables in use_idx
scatter_sub_sgm_climo	= False ; False ; --- climo: sgm related variables in Y and sgm in X axis
scatter_sub_sgm_diff	= True ; True ; --- p4K-cntl: sgm related variables in Y and sgm in X axis

sherwood_index		= (/"SS","DD","LTMI"/)

plot_type			= "pdf"

latS_trop			= -30.
latE_trop			= 30.
lonS_scatter		= 220.
lonE_scatter		= 270.

select_idx			= 0	; choose how many regimes would be plotted 
scale_factor		= 3.	; use to amplify diff to make them able to plot in non-diff plot
regime_mask_val		= (/0,1,2,3,4,5,6/)
LTS_criteria		= 17.

nbins				= 10
;ALPI_mask_val		= fspan(0,nbins-1,nbins)
ALPI_mask_val		= fspan(0,nbins-1,nbins)+0.5 ;-- the addition of 0.5 is to make sure that 
												 ;-- the value is in the center of cnLevels rather 
												 ;-- than the edge. After this change, RasterFill
												 ;-- is not needed and masking plot would be right as well.

fdir				= "/home/share3/lyl/work3/qinyi/mid-data/cfmip/"

;case_nor 				= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_MG15_amip-p4K"/)
;case_p4K		 		= (/"FAMIPC5_f09f09_mac2_amip_outSGM","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM"/)
;case_nor 				= (/"FAMIPC5_f09f09_MG15_amip","FAMIPC5_f09f09_mac2_amip_outSGM"/)
;case_p4K		 		= (/"FAMIPC5_f09f09_MG15_amip-p4K","FAMIPC5_f09f09_mac2_amip-p4K_debug_outSGM"/)

case_nor 				= (/"FAMIPC5_f19f19_MG15_amip","FAMIPC5_f19f19_mac2_amip_outSGM"/)
case_p4K		 		= (/"FAMIPC5_f19f19_MG15_amip-p4K","FAMIPC5_f19f19_mac2_amip-p4K_debug_outSGM"/)


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

sgm_shal_levs		= fspan(0.,0.84,15)
sgm_shal_dlevs_1	= fspan(-0.28,0.28,15)
sgm_shal_dlevs_2	= sgm_tota_dlevs_1
sgm_shal_dlevs_sp	= (/-1, -0.5, 0.5, 1/)

sgm_turb_levs		= fspan(0.,0.84,15)
sgm_turb_dlevs_1	= fspan(-0.28,0.28,15)
sgm_turb_dlevs_2	= sgm_tota_dlevs_1
sgm_turb_dlevs_sp	= (/-1, -0.5, 0.5, 1/)

N1_levs				= fspan(0,0.84,15)
N1_dlevs_1			= fspan(-0.28,0.28,15)
N1_dlevs_2			= N1_dlevs_1
N1_dlevs_sp			= (/-1, -0.5, 0., 1/)

netCRF_levs			= (/-80,-75,-70,-65,-60,-55,-50,-45,-40,-30,-20,-10,0,10,20/)
netCRF_dlevs_1		= (/-50,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,50/)/2.
netCRF_dlevs_2		= netCRF_dlevs_1
netCRF_dlevs_sp		= (/-20,-10,10,20/)

RELHUM_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)*100.
RELHUM_dlevs_1		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)*100.
RELHUM_dlevs_2		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)*100.
RELHUM_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)*100.

FREQSH_levs			= (/.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.8,.9/)
FREQSH_dlevs_1		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
FREQSH_dlevs_2		= (/-0.25,-0.2,-0.16,-0.12,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.12,0.16,0.20,0.25/)
FREQSH_dlevs_sp		= (/-0.1,-0.05,0.05,0.1/)

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


Vars_mid	= (/"T", "Q", "OMEGA", "RELHUM"/)

;				   0            1          2           3           4            5          6	    7      8       9		10	      11	  12 	 13        14         15           16        17         18      19			20		   21         22      23     24       25         26       27      28        29      30     31       32         33          34			35               36				    37				   38			   39			   40			    41					42	            43       44        45           46           47        48
varsA		= (/"PRECT",	"PRECC",	"PRECL",	"SWCF",		"LWCF",		"CLDLOW",	"CLDHGH", "LTS", "EIS", "ECTEI","TGCLDLWP","sgm_tota","N1", "netCRF","sgm_shal","sgm_turb", "RELHUM", "FREQSH", "CMFMC","umf_shal", "clddep2", "wstarPBL", "lengi", "shi", "aa_gp",	"bb_gp", "beta_gp", "dqwdz", "dthldz", "qw2", "thl2", "qwthl", "qwu_qw", "thlu_thl", "qwu_thlu", "qtu_shal", "vert_qwflx_shal", "vert_thlflx_shal","vert_qwthlflx_shal","tau_shal","vert_qwflx_turb","vert_thlflx_turb","vert_qwthlflx_turb","lengi2","dqwdz_2","dthldz_2","dqwdz_dthldz","aa2_gp", "bb2_gp"/)
varsA_levs	= (/PRECT_levs, PRECT_levs,     PRECT_levs,     SWCF_levs,      LWCF_levs, 	CLDLOW_levs,	CLDHGH_levs, LTS_levs, EIS_levs, ECTEI_levs,TGCLDLWP_levs,sgm_tota_levs, N1_levs, netCRF_levs,sgm_shal_levs, sgm_turb_levs, RELHUM_levs, FREQSH_levs, CMFMC_levs, umf_shal_levs, clddep2_levs, wstarPBL_levs, lengi_levs, shi_levs, aa_gp_levs, bb_gp_levs, beta_gp_levs, dqwdz_levs, dthldz_levs, qw2_levs, thl2_levs, qwthl_levs, qwu_qw_levs, thlu_thl_levs, qwu_thlu_levs, qtu_shal_levs,vert_qwflx_shal_levs, vert_thlflx_shal_levs, vert_qwthlflx_shal_levs, tau_shal_levs, vert_qwflx_turb_levs, vert_thlflx_turb_levs, vert_qwthlflx_turb_levs, lengi2_levs,dqwdz_2_levs, dthldz_2_levs, dqwdz_dthldz_levs, aa2_gp_levs, bb2_gp_levs/)
varsA_dlevs_1	= (/PRECT_dlevs_1,PRECT_dlevs_1,    PRECT_dlevs_1,    SWCF_dlevs_1,     LWCF_dlevs_1,	CLDLOW_dlevs_1,	CLDHGH_dlevs_1, LTS_dlevs_1, EIS_dlevs_1,	ECTEI_dlevs_1,TGCLDLWP_dlevs_1, sgm_tota_dlevs_1, N1_dlevs_1, netCRF_dlevs_1, sgm_shal_dlevs_1, sgm_turb_dlevs_1, RELHUM_dlevs_1, FREQSH_dlevs_1, CMFMC_dlevs_1, umf_shal_dlevs_1, clddep2_dlevs_1, wstarPBL_dlevs_1, lengi_dlevs_1, shi_dlevs_1, aa_gp_dlevs_1, bb_gp_dlevs_1, beta_gp_dlevs_1, dqwdz_dlevs_1, dthldz_dlevs_1, qw2_dlevs_1, thl2_dlevs_1, qwthl_dlevs_1, qwu_qw_dlevs_1, thlu_thl_dlevs_1, qwu_thlu_dlevs_1, qtu_shal_dlevs_1,vert_qwflx_shal_dlevs_1, vert_thlflx_shal_dlevs_1, vert_qwthlflx_shal_dlevs_1, tau_shal_dlevs_1, vert_qwflx_turb_dlevs_1, vert_thlflx_turb_dlevs_1, vert_qwthlflx_turb_dlevs_1, lengi2_dlevs_1, dqwdz_2_dlevs_1, dthldz_2_dlevs_1, dqwdz_dthldz_dlevs_1, aa2_gp_dlevs_1, bb2_gp_dlevs_1/)
varsA_dlevs_2	= (/PRECT_dlevs_2,PRECT_dlevs_2,    PRECT_dlevs_2,    SWCF_dlevs_2,     LWCF_dlevs_2,	CLDLOW_dlevs_2,	CLDHGH_dlevs_2, LTS_dlevs_2, EIS_dlevs_2,	ECTEI_dlevs_2,TGCLDLWP_dlevs_2, sgm_tota_dlevs_2, N1_dlevs_2, netCRF_dlevs_2, sgm_shal_dlevs_2, sgm_turb_dlevs_2, RELHUM_dlevs_2, FREQSH_dlevs_2, CMFMC_dlevs_2, umf_shal_dlevs_2, clddep2_dlevs_2, wstarPBL_dlevs_2, lengi_dlevs_2, shi_dlevs_2, aa_gp_dlevs_2, bb_gp_dlevs_2, beta_gp_dlevs_2, dqwdz_dlevs_2, dthldz_dlevs_2, qw2_dlevs_2, thl2_dlevs_2, qwthl_dlevs_2, qwu_qw_dlevs_2, thlu_thl_dlevs_2, qwu_thlu_dlevs_2, qtu_shal_dlevs_2,vert_qwflx_shal_dlevs_2, vert_thlflx_shal_dlevs_2, vert_qwthlflx_shal_dlevs_2, tau_shal_dlevs_2, vert_qwflx_turb_dlevs_2, vert_thlflx_turb_dlevs_2, vert_qwthlflx_turb_dlevs_2,lengi2_dlevs_2, dqwdz_2_dlevs_2, dthldz_2_dlevs_2, dqwdz_dthldz_dlevs_2, aa2_gp_dlevs_2, bb2_gp_dlevs_2/)
varsA_dlevs_sp	= (/PRECT_dlevs_sp,PRECT_dlevs_sp,    PRECT_dlevs_sp,    SWCF_dlevs_sp,     LWCF_dlevs_sp,	CLDLOW_dlevs_sp,	CLDHGH_dlevs_sp, LTS_dlevs_sp, EIS_dlevs_sp, ECTEI_dlevs_sp,TGCLDLWP_dlevs_sp, sgm_tota_dlevs_sp, N1_dlevs_sp, netCRF_dlevs_sp, sgm_shal_dlevs_sp, sgm_turb_dlevs_sp, RELHUM_dlevs_sp, FREQSH_dlevs_sp, CMFMC_dlevs_sp, umf_shal_dlevs_sp, clddep2_dlevs_sp, wstarPBL_dlevs_sp, lengi_dlevs_sp, shi_dlevs_sp, aa_gp_dlevs_sp, bb_gp_dlevs_sp, beta_gp_dlevs_sp, dqwdz_dlevs_sp, dthldz_dlevs_sp, qw2_dlevs_sp, thl2_dlevs_sp, qwthl_dlevs_sp, qwu_qw_dlevs_sp, thlu_thl_dlevs_sp, qwu_thlu_dlevs_sp, qtu_shal_dlevs_sp, vert_qwflx_shal_dlevs_sp, vert_thlflx_shal_dlevs_sp, vert_qwthlflx_shal_dlevs_sp, tau_shal_dlevs_sp, vert_qwflx_turb_dlevs_sp, vert_thlflx_turb_dlevs_sp, vert_qwthlflx_turb_dlevs_sp, lengi2_dlevs_sp, dqwdz_2_dlevs_sp, dthldz_2_dlevs_sp, dqwdz_dthldz_dlevs_sp, aa2_gp_dlevs_sp, bb2_gp_dlevs_sp/)

varsA_units	= (/"mm/day",   "mm/day",       "mm/day",       "W/m~S~2~N~",   "W/m~S~2~N~",	"fraction",	"fraction", "K",	"K",	"K", "g/m2", "-", "fraction", "W/m~S~2~N~", "-", "-", "%", "fraction","kg/m2/s", "kg/m2/s","m", "m/s", "m", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-","-","-","-","-","-","-","-","-","-","-","-","-","-"/)
plotsA_strings	= (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)","(m)","(n)","(o)","(p)","(q)","(r)","(s)","(t)","(u)","(v)","(w)","(x)","(y)","(z)"/)


; --- scatter
;use_idx			= (/14,17,19,20,21,24,25,27,28,32,33,34,35/)		; for scatter plot -- sgm_shal
;use_idx			= (/14,17,36,37,38,39/)		; for scatter plot -- sgm_shal
;use_idx			= (/15,22,23,24,25,27,28,29,30,31/)		; for scatter plot -- sgm_turb
;use_idx			= (/15,40,41,42,43,44,45,46/)

;use_idx				= (/3,5/)	; --- one variable in X and Y axis

; --- ALPI
;use_idx				= (/0,13,3,4,5,10,11,12/)	; for ALPI
;use_idx				= (/0,11,12,14,15/)	; for ALPI -- only for sgm_tota and N1, and sgm_turb, sgm_shal
;use_idx				= (/0,14,17,18,19,20,21,24,25,32,33,34,47,48,36,37,38,39/)	 	; sgm_shal related
use_idx				= (/0,15,22,23,24,25,27,28,47,48,40,41,42,43,44,45,46/)			; sgm_turb related 29,30,31 == 40,41,42


;use_idx			= (/3,4,13,5,6,10/)
;use_idx			= (/9,7,8/)
;use_idx			= (/13,3,5,10,11,12/)
;use_idx			= (/13,3,4,5,6,10/)
;use_idx				= (/0,5,10,16/)
;use_idx				= (/0,16,17,18/)
;use_idx				= (/0,16,17,18,19,20,21,22,23,24/)

; --- for variables only in Gauss-PDF
use_idx_pdf			= (/11,12,14,15,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48/)

vars				= varsA(use_idx)
vars_levs			= varsA_levs(use_idx,:)
vars_dlevs_1		= varsA_dlevs_1(use_idx,:)
vars_dlevs_2		= varsA_dlevs_2(use_idx,:)
vars_dlevs_sp		= varsA_dlevs_sp(use_idx,:)
vars_units			= varsA_units(use_idx)
plots_strings		= plotsA_strings

vars_pdf			= varsA(use_idx_pdf)


nlat				= 73		; number of output latitudes
nlon				= 144		; number of output longitudes
nlev				= 30

klon = ispan(0,nlon-1,1)*2.5	; to 2.5 x 2.5
klat = ispan(0,nlat-1,1)*2.5-90.

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
	tmp									= ftmp1->OCNFRAC(0,:,:)
	OCNFRAC								= linint2_Wrap(tmp&lon,tmp&lat,tmp,True,klon,klat,0)
	delete(tmp)

	do ivar=0,dimsizes(Vars_mid)-1
		tmp								= ftmp1->$Vars_mid(ivar)$(0,:,:,:)
		vars_mid(ivar,itag,icase,:,:,:)		= linint2_Wrap(tmp&lon,tmp&lat,tmp,True,klon,klat,0)
		delete(tmp)
	end do ; do ivar=

	do ivar=0,dimsizes(vars)-1
		if(isfilevar(ftmp1,vars(ivar))) then
			tmp									= ftmp1->$vars(ivar)$
			if(dimsizes(dimsizes(tmp)).eq.4)then ; figure out variables having vertical level coordinate
				tmp1							= ftmp1->$vars(ivar)$(0,{887},:,:)
				vars_all(ivar,itag,icase,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
				delete(tmp1)
			else
				tmp1						 	= ftmp1->$vars(ivar)$(0,:,:)
				vars_all(ivar,itag,icase,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
				delete(tmp1)
			end if

			if(vars(ivar).eq."TGCLDLWP")then
				vars_all(ivar,itag,icase,:,:)	= vars_all(ivar,itag,icase,:,:)*1e3
			end if
			if(vars(ivar).eq."sgm_tota" .or. vars(ivar).eq."sgm_shal" .or. vars(ivar).eq."sgm_turb")then
;				vars_all(ivar,itag,icase,:,:)	= sqrt(abs(vars_all(ivar,itag,icase,:,:)))*1e3
				vars_all(ivar,itag,icase,:,:)	= abs(vars_all(ivar,itag,icase,:,:))*1e6
			end if
			delete(tmp)
		else
			if(.not.ismissing(ind(vars_pdf.eq.vars(ivar))) .and. icase.eq.0)then
				print(vars(ivar))
				vars_all(ivar,itag,icase,:,:)	= vars_all@_FillValue
			else
				if(vars(ivar).eq."lengi2")then
					lengi_tmp			= ftmp1->lengi(0,{887},:,:)
					tmp1					= lengi_tmp^2
					copy_VarCoords(lengi_tmp,tmp1)
					vars_all(ivar,itag,icase,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
					delete(tmp1)
					delete(lengi_tmp)
				end if
				if(vars(ivar).eq."aa2_gp")then
					aa_gp_tmp			= ftmp1->aa_gp(0,{887},:,:)
					tmp1					= aa_gp_tmp^2
					copy_VarCoords(aa_gp_tmp,tmp1)
					vars_all(ivar,itag,icase,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
					delete(tmp1)
					delete(aa_gp_tmp)
				end if
				if(vars(ivar).eq."bb2_gp")then
					bb_gp_tmp			= ftmp1->bb_gp(0,{887},:,:)
					tmp1					= bb_gp_tmp^2
					copy_VarCoords(bb_gp_tmp,tmp1)
					vars_all(ivar,itag,icase,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
					delete(tmp1)
					delete(bb_gp_tmp)
				end if
				if(vars(ivar).eq."dqwdz_2")then
					dqwdz_tmp			= ftmp1->dqwdz(0,{887},:,:)
					tmp1					= dqwdz_tmp^2
					copy_VarCoords(dqwdz_tmp,tmp1)
					vars_all(ivar,itag,icase,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
					delete(tmp1)
					delete(dqwdz_tmp)
				end if
				if(vars(ivar).eq."dthldz_2")then
					dthldz_tmp			= ftmp1->dthldz(0,{887},:,:)
					tmp1					= dthldz_tmp^2
					copy_VarCoords(dthldz_tmp,tmp1)
					vars_all(ivar,itag,icase,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
					delete(tmp1)
					delete(dthldz_tmp)
				end if
				if(vars(ivar).eq."dqwdz_dthldz")then
					dthldz_tmp			= ftmp1->dthldz(0,{887},:,:)
					dqwdz_tmp			= ftmp1->dqwdz(0,{887},:,:)
					tmp1					= dthldz_tmp*dqwdz_tmp
					copy_VarCoords(dthldz_tmp,tmp1)
					vars_all(ivar,itag,icase,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
					delete(tmp1)
					delete([/dthldz_tmp,dqwdz_tmp/])
				end if

				if(vars(ivar).eq."vert_qwflx_turb")then
					tmp1					= ftmp1->qw2(0,{887},:,:)
					vars_all(ivar,itag,icase,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
					delete(tmp1)
				end if
				if(vars(ivar).eq."vert_thlflx_turb")then
					tmp1					= ftmp1->thl2(0,{887},:,:)
					vars_all(ivar,itag,icase,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
					delete(tmp1)
				end if
				if(vars(ivar).eq."vert_qwthlflx_turb")then
					tmp1					= ftmp1->qwthl(0,{887},:,:)
					tmp1					= tmp1*(-1.)
					vars_all(ivar,itag,icase,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
					delete(tmp1)
				end if
				if(vars(ivar).eq."vert_qwflx_shal")then
					umf_shal_tmp			= ftmp1->umf_shal(0,{887},:,:)
					qwu_qw_tmp				= ftmp1->qwu_qw(0,{887},:,:)
					tmp1					= -2.*umf_shal_tmp*qwu_qw_tmp
					copy_VarCoords(umf_shal_tmp,tmp1)
					vars_all(ivar,itag,icase,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
					delete(tmp1)
					delete(umf_shal_tmp)
					delete(qwu_qw_tmp)
				end if
				if(vars(ivar).eq."vert_thlflx_shal")then
					umf_shal_tmp			= ftmp1->umf_shal(0,{887},:,:)
					thlu_thl_tmp			= ftmp1->thlu_thl(0,{887},:,:)
					tmp1					= -2.*umf_shal_tmp*thlu_thl_tmp
					copy_VarCoords(umf_shal_tmp,tmp1)
					vars_all(ivar,itag,icase,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
					delete(tmp1)
					delete(umf_shal_tmp)
					delete(thlu_thl_tmp)
				end if
				if(vars(ivar).eq."vert_qwthlflx_shal")then
					umf_shal_tmp			= ftmp1->umf_shal(0,{887},:,:)
					qwu_thlu_tmp			= ftmp1->qwu_thlu(0,{887},:,:)
					tmp1					= umf_shal_tmp*qwu_thlu_tmp
					copy_VarCoords(umf_shal_tmp,tmp1)
					vars_all(ivar,itag,icase,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
					delete(tmp1)
					delete(umf_shal_tmp)
					delete(qwu_thlu_tmp)
				end if
				if(vars(ivar).eq."tau_shal")then
					clddep2_tmp			= ftmp1->clddep2(0,:,:)
					wstarPBL_tmp		= ftmp1->wstarPBL(0,:,:)
					tmp1					= clddep2_tmp/where(wstarPBL_tmp.ne.0,wstarPBL_tmp,vars_all@_FillValue)
					copy_VarCoords(clddep2_tmp,tmp1)
					vars_all(ivar,itag,icase,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
					delete(tmp1)
					delete([/clddep2_tmp,wstarPBL_tmp/])
				end if
			end if ; .not. ismissing

			if(vars(ivar).eq."PRECT")then
				precl_tmp				= ftmp1->PRECL(0,:,:)
				precc_tmp				= ftmp1->PRECC(0,:,:)
				tmp1					= precl_tmp+precc_tmp
				copy_VarCoords(precl_tmp,tmp1)
				vars_all(ivar,itag,icase,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
				delete(tmp1)
				vars_all(ivar,itag,icase,:,:)	= vars_all(ivar,itag,icase,:,:)*8.64e7
				delete(precl_tmp)
				delete(precc_tmp)
			end if
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
;			if(vars(ivar).eq."sgm_tota" .or. vars(ivar).eq."N1" .or. vars(ivar).eq."sgm_shal" .or. vars(ivar).eq."sgm_turb".or.\
;			   vars(ivar).eq."wstarPBL" .or. vars(ivar).eq."clddep2" .or. vars(ivar).eq."umf_shal" .or. vars(ivar).eq."lengi" .or. \
;			   vars(ivar).eq."shi" .or. vars(ivar).eq."aa_gp") then
;			if(.not.ismissing(ind(vars_pdf.eq.vars(ivar))))then
;				print(vars(ivar))
;				vars_all(ivar,itag,icase,:,:)	= vars_all@_FillValue
;			end if
			if(vars(ivar).eq."netCRF")then
				SWCF						= ftmp1->SWCF(0,:,:)
				LWCF						= ftmp1->LWCF(0,:,:)
				tmp1						= (/SWCF+LWCF/)
				copy_VarCoords(SWCF,tmp1)
				vars_all(ivar,itag,icase,:,:)	= linint2_Wrap(tmp1&lon,tmp1&lat,tmp1,True,klon,klat,0)
				delete(tmp1)
				delete([/SWCF,LWCF/])
			end if

		end if
		if(icase.eq.1)then
		print(vars(ivar))
		printMinMax(vars_all(ivar,itag,icase,:,:),True)
		end if
	end do ; do ivar=

end do ; do icase=
end do ; do itag=

; change dimension name
vars_all!3		= "lat"
vars_all!4		= "lon"

vars_mid!4		= "lat"
vars_mid!5		= "lon"

printVarSummary(vars_all)
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

delete([/res, res1, res2, tmp, tmp_xaxis, legends, lgres, amres, lbid1, annoid1/])
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
delete(tmp)

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


lon					= vars_all&lon
lat                 = vars_all&lat
nlat				= dimsizes(lat)
nlon				= dimsizes(lon)
clat                = tofloat(cos(0.01745329*lat))
clat!0              = "lat"
clat&lat            = vars_all&lat

; ################# scatter_p2p ########################
if(scatter_p2p)then

; ------ set the plotting variables

vars_all_trop_ocn			= vars_all
vars_all_tagdiff_trop_ocn	= vars_all_tagdiff
vars_all_casediff_trop_ocn	= vars_all_casediff
vars_all_diff2_trop_ocn		= vars_all_diff2

do ivar=0,dimsizes(vars)-1
do icase=0,dimsizes(case_nor)-1
do itag=0,dimsizes(tag)-1

omega500					= vars_mid(2,itag,icase,{500},:,:)
omega500					= omega500*864.
lat2d						= conform_dims(dimsizes(omega500),lat,0)
lon2d						= conform_dims(dimsizes(omega500),lon,1)

if(scatter_sub_ocn)then
; mask land region
	crit1					= (OCNFRAC.gt.0.8).and.(lat2d.gt.latS_trop.and.lat2d.lt.latE_trop).and.(omega500.gt.0)
	vars_all_trop_ocn(ivar,itag,icase,:,:)					= where(crit1, vars_all(ivar,itag,icase,:,:),vars_all@_FillValue)
	vars_all_tagdiff_trop_ocn(ivar,icase,:,:)				= where(crit1, vars_all_tagdiff(ivar,icase,:,:),vars_all@_FillValue)
	vars_all_casediff_trop_ocn(ivar,itag,:,:)				= where(crit1, vars_all_casediff(ivar,itag,:,:),vars_all@_FillValue)
	vars_all_diff2_trop_ocn(ivar,:,:)						= where(crit1, vars_all_diff2(ivar,:,:),vars_all@_FillValue)
else
	crit1					= (lat2d.gt.latS_trop.and.lat2d.lt.latE_trop).and.(omega500.gt.0)
	vars_all_trop_ocn(ivar,itag,icase,:,:)					= where(crit1,vars_all(ivar,itag,icase,:,:),vars_all@_FillValue)
	vars_all_tagdiff_trop_ocn(ivar,icase,:,:)				= where(crit1,vars_all_tagdiff(ivar,icase,:,:),vars_all@_FillValue)
	vars_all_casediff_trop_ocn(ivar,itag,:,:)				= where(crit1,vars_all_casediff(ivar,itag,:,:),vars_all@_FillValue)
	vars_all_diff2_trop_ocn(ivar,:,:)						= where(crit1,vars_all_diff2(ivar,:,:),vars_all@_FillValue)
end if

end do ; do itag=
end do ; do icase=
end do ; do ivar=

pltarry					 	= new((/dimsizes(vars),dimsizes(case_nor),2,nlat*nlon/),float)
pltx						= new((/dimsizes(vars),dimsizes(case_nor),nlat*nlon/),float)
rcarry						= new((/dimsizes(vars),2/),float)


do ivar=0,dimsizes(vars)-1
	if(vars(0).eq."sgm_shal" .or. vars(0).eq."sgm_turb")then
	start_idx=1
	else
	start_idx=0
	end if

do icase=start_idx,dimsizes(case_nor)-1	; only for Gauss-PDF

if(scatter_sub_one_var)then
; --- for cntl and p4K-cntl plot for the same variable
x							= ndtooned(vars_all_trop_ocn(ivar,0,icase,::3,::3))
y							= ndtooned(vars_all_casediff_trop_ocn(ivar,0,::3,::3))
;y							= ndtooned(vars_all_tagdiff_trop_ocn(ivar,icase,::3,::3))
end if

; --- for Gauss-PDF, different variables and try to capture their correlation
; --- x: sgm_shal or sgm_turb; y: any other related variables for sgm_shal and sgm_turb

if(scatter_sub_sgm_climo)then
; -- for climatology
x							= ndtooned(vars_all_trop_ocn(0,0,icase,::3,::3))
y							= ndtooned(vars_all_trop_ocn(ivar,0,icase,::3,::3))
end if

if(scatter_sub_sgm_diff)then
; -- for difference between p4K and cntl
x							= ndtooned(vars_all_tagdiff_trop_ocn(0,icase,::3,::3))
y							= ndtooned(vars_all_tagdiff_trop_ocn(ivar,icase,::3,::3))
end if

printVarSummary(x)
printVarSummary(y)

; -- make the match between x and y by eliminating any data in each one array
ytmp						= where(x.ne.vars_all@_FillValue,y,vars_all@_FillValue)
xtmp						= where(ytmp.ne.vars_all@_FillValue,x,vars_all@_FillValue)
printVarSummary(ytmp)
printVarSummary(xtmp)

; --- sorting is used to eliminate those missing values.
ii							= dim_pqsort_n(xtmp,1,0)
xx							= xtmp(ii)
yy							= ytmp(ii)
effnum						= (num(.not.ismissing(xx)))-1 ; -- number of no-missing values
print(effnum)

; linear regression
rc							= regline_stats(xx(:effnum),yy(:effnum))
printVarSummary(rc)

nx							= dimsizes(xx(:effnum))
print(nx)

pltarry(ivar,icase,0,:effnum)				= yy(:effnum)
pltarry(ivar,icase,1,:effnum)				= rc@Yest

pltx(ivar,icase,:effnum)				= xx(:effnum)
rcarry(ivar,icase)							= rc@r


end do ; do icase=
end do ; do ivar=


if(scatter_sub_sgm_climo)then
wks								= gsn_open_wks(plot_type,"./figure/scatter_sgm_climo_"+vars(0))
else
	if(scatter_sub_sgm_diff)then
		wks								= gsn_open_wks(plot_type,"./figure/scatter_sgm_diff_"+vars(0))
	else
		wks								= gsn_open_wks(plot_type,"./figure/scatter")
	end if
end if


plot							= new((/dimsizes(vars)/),graphic)
txid							= new((/dimsizes(vars)/),graphic)
annoid							= new((/dimsizes(vars)/),graphic)
if(start_idx.eq.0)then
colors							= (/"grey","lightskyblue"/)
else
colors							= (/"blue","lightskyblue"/)
end if


count	= 0

;do ivar=0,dimsizes(vars)-1
do ivar=start_idx,dimsizes(vars)-1

res								= True
res@gsnDraw						= False
res@gsnFrame					= False
res@xyMarkLineModes 			= "Markers"
res@xyMarkers					= 16
res@xyMarkerSizeF				= 0.01
res@tmLabelAutoStride			= True

res@xyMarkerColor				= colors(0)
res@xyLineColor					= colors(0)
res@xyLineThicknessF			= 3.0

res@gsnStringFontHeightF		= 0.03
if(scatter_sub_one_var)then
res@tiXAxisString				= "CNTL_"+vars(ivar)
res@tiYAxisString				= "P4K-CNTL_"+vars(ivar)
end if
if(scatter_sub_sgm_climo)then
res@tiXAxisString				= "climo_"+vars(0)
res@tiYAxisString				= "climo_"+vars(ivar)
end if
if(scatter_sub_sgm_diff)then
res@tiXAxisString				= "P4K-CNTL_"+vars(0)
res@tiYAxisString				= "P4K-CNTL_"+vars(ivar)
end if

res@gsnLeftString				= plots_strings(count)
res@gsnYRefLine					= 0
res@gsnXRefLine					= 0

res@tiYAxisFontHeightF			= 0.03
res@tiXAxisFontHeightF			= 0.03
res@tmXBLabelFontHeightF		= 0.03
res@tmYLLabelFontHeightF		= 0.03

res@gsnStringFontHeightF	= 0.03

res@tmYLFormat				= "0*+^e"
res@tmYLMinorOn				= False
res@tmXBMinorOn				= False
;res@tmXBPrecision			= 6


;plot(ivar)				= gsn_csm_xy(wks,pltx(ivar,0,:),pltarry(ivar,0,:,:),res)
plot(ivar)				= gsn_csm_xy(wks,pltx(ivar,start_idx,:),pltarry(ivar,start_idx,:,:),res)

count	= count + 1

if(start_idx.eq.0)then
do icase=1,dimsizes(case_nor)-1
	res@xyMarkerColor		= colors(icase)
	res@xyLineColor			= colors(icase)
	res@xyLineThicknessF	= 3.0

	ovplot					= gsn_csm_xy(wks,pltx(ivar,icase,:),pltarry(ivar,icase,:,:),res)
	overlay(plot(ivar),ovplot)
end do ; do icase=
end if

txres							= True
txres@txPerimOn					= True
txres@txFontHeightF				= 0.03
;txres@txBackgroundFillColor		= colors(ivar)

if(start_idx.eq.0)then
txid(ivar)				= gsn_create_text(wks,"r="+sprintf("%6.3f",rcarry(ivar,0))+"; "+sprintf("%6.3f",rcarry(ivar,1)),txres)
else
txid(ivar)				= gsn_create_text(wks,"r="+sprintf("%6.3f",rcarry(ivar,1)),txres)
end if

amres							= True
amres@amJust					= "TopLeft"
amres@amParallelPosF			= -0.5
amres@amOrthogonalPosF			= -0.5

annoid(ivar)				= gsn_add_annotation(plot(ivar),txid(ivar), amres)


end do ; do ivar=

resP								= True
resP@gsnPanelXWhiteSpacePercent		= 5.
resP@gsnPanelYWhiteSpacePercent		= 5.
resP@gsnPanelLeft					= 0.05
resP@gsnPanelRight					= 0.95

if(start_idx.eq.0)then
gsn_panel(wks,ndtooned(plot),(/dimsizes(ndtooned(plot))/2,2/),resP)
else
	if(mod(dimsizes(vars),4).ne.0)then
		gsn_panel(wks,ndtooned(plot),(/dimsizes(ndtooned(plot))/4+1,4/),resP)
	else
		gsn_panel(wks,ndtooned(plot),(/dimsizes(ndtooned(plot))/4,4/),resP)
	end if
end if

delete(count)
exit

end if ; scatter_p2p



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

res@lbLabelBarOn			= False
res@lbLabelPosition  		= "Center"           ; label position
res@lbLabelAlignment 		= "BoxCenters"       ; label orientation
;res@lbLabelStrings   		= ""+ res@cnLevels   ; trick 
res@lbLabelStrings   		= ""+ (/ 0,1,2,3,4,5,6/)

res@pmLabelBarHeightF		=  0.100
res@pmLabelBarWidthF 		=  0.80

res@mpCenterLonF			= 180

res@gsnStringFontHeightF	= 0.025
res@tmXBLabelFontHeightF  	= 0.018
res@tmXTLabelFontHeightF  	= 0.018
res@tmYLLabelFontHeightF  	= 0.018

res@tmXBMinorOn				= False
res@tmXTMinorOn				= False
res@tmYLMinorOn				= False
res@tmYRMinorOn				= False

;res@gsnLeftString			= "Selected regimes"

res1						= res
delete([/res1@cnLevels, res1@cnFillColors, res1@lbLabelStrings/])
;--------------------------- These settings are for regimes_mask_s -----------
; each regime is colored in single plot rather than all in one plot.
res1@cnLevels				= (/1/)
res1@cnFillColors			= (/"lightcyan2","white"/)
res1@lbLabelStrings   		= ""+ (/ 0,1/)

plot						= new((/dimsizes(tag),dimsizes(case_nor)/),graphic)
do ivar=0,0;dimsizes(vars)-1
count						= 0
do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1
print(count)
res@gsnLeftString			= plots_strings(count)+" Selected regimes"
res@gsnRightString = tag(itag)+"_"+case_name(icase)
plot(itag,icase)		= gsn_csm_contour_map_ce(wks,data2regimes_mask(ivar,itag,icase,:,:),res)

count						= count+1
;do ir=0,dimsizes(regimes_name)-1
;plot1		= gsn_csm_contour_map_ce(wks,data2regimes_mask_s(ir,ivar,itag,icase,:,:),res1)
;end do 

end do
end do
end do

pres								= True
pres@gsnPanelLabelBar				= True
pres@pmLabelBarWidthF				= 0.8
pres@lbLabelFontHeightF				= 0.015
pres@lbLabelPosition  				= "Center"           ; label position
pres@lbLabelAlignment 				= "BoxCenters"       ; label orientation
pres@lbLabelStrings   				= ""+ (/ 0,1,2,3,4,5,6/)
pres@gsnPanelYWhiteSpacePercent		= 5.0
pres@gsnPanelXWhiteSpacePercent		= 5.0


gsn_panel(wks,ndtooned(plot),(/2,2/),pres)
delete(count)

exit
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
sres@vpHeightF             = 0.4 ; -- changed from 0.5 to 0.4; Dec 21, 2018

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

sres@gsnStringFontHeightF	= 0.03
sres@gsnRightString        = ""
;sres@tiYAxisString         = "W/m2"

sres@tmYLMinorOn			= False
sres@tmYRMinorOn			= False

sres@gsnStringFontHeightF	= 0.03

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
count			= 0
do ivar=0,dimsizes(vars)-1
  if(.not.ismissing(ind(vars_pdf.eq.vars(ivar))))then ; for vars_pdf 
	sres@tiYAxisString		= vars(ivar)
	sres2@tiYAxisString		= vars(ivar)
	sres3@tiYAxisString		= vars(ivar)
  else
	sres@tiYAxisString		= vars(ivar)+"("+vars_units(ivar)+")"
	sres2@tiYAxisString		= vars(ivar)+"("+vars_units(ivar)+")"
	sres3@tiYAxisString		= vars(ivar)+"("+vars_units(ivar)+")"
	end if

	sres@trYMinF			= -1.*max((/max(regimes_tagdiff(select_idx:,ivar,:)),abs(min(regimes_tagdiff(select_idx:,ivar,:))),max(scale_factor*regimes_diff2(select_idx,ivar)),abs(min(scale_factor*regimes_diff2(select_idx,ivar)))/))*1.2
	sres@trYMaxF			= 1.*max((/max(regimes_tagdiff(select_idx:,ivar,:)),abs(min(regimes_tagdiff(select_idx:,ivar,:))),max(scale_factor*regimes_diff2(select_idx,ivar)),abs(min(scale_factor*regimes_diff2(select_idx,ivar)))/))*1.2
	sres@gsnLeftString		= plots_strings(count)+" un-weighted"

;	sres2@tiYAxisString		= vars(ivar)+"("+vars_units(ivar)+")"
	sres2@trYMinF			= -1.*max((/max(regimes_tagdiff_areawgt(select_idx:,ivar,:)),abs(min(regimes_tagdiff_areawgt(select_idx:,ivar,:))),max(scale_factor*regimes_diff2_areawgt(select_idx:,ivar)),abs(min(scale_factor*regimes_diff2_areawgt(select_idx:,ivar)))/))*1.2
	sres2@trYMaxF			= 1.*max((/max(regimes_tagdiff_areawgt(select_idx:,ivar,:)),abs(min(regimes_tagdiff_areawgt(select_idx:,ivar,:))),max(scale_factor*regimes_diff2_areawgt(select_idx:,ivar)),abs(min(scale_factor*regimes_diff2_areawgt(select_idx:,ivar)))/))*1.2

;	sres3@tiYAxisString		= vars(ivar)+"("+vars_units(ivar)+")"
	sres3@trYMinF			= -1.*max((/max(regimes_tagdiff_areawgt_1(select_idx:,ivar,:)),abs(min(regimes_tagdiff_areawgt_1(select_idx:,ivar,:))),max(scale_factor*regimes_diff2_areawgt_1(select_idx:,ivar)),abs(min(scale_factor*regimes_diff2_areawgt_1(select_idx:,ivar)))/))*1.2
	sres3@trYMaxF			= 1.*max((/max(regimes_tagdiff_areawgt_1(select_idx:,ivar,:)),abs(min(regimes_tagdiff_areawgt_1(select_idx:,ivar,:))),max(scale_factor*regimes_diff2_areawgt_1(select_idx:,ivar)),abs(min(scale_factor*regimes_diff2_areawgt_1(select_idx:,ivar)))/))*1.2

	sres2@gsnLeftString		= plots_strings(count)+" area-weighted"
	sres3@gsnLeftString		= plots_strings(count)+" area-weighted-but-same"

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

;	if(vars(ivar).ne."sgm_tota" .and. vars(ivar).ne."N1" .and. vars(ivar).ne."sgm_shal" .and. vars(ivar).and."sgm_turb".and.\
;	   vars(ivar).ne."wstarPBL" .and. vars(ivar).ne."clddep2" .and. vars(ivar).ne."umf_shal" .and. vars(ivar).ne."shi" .and.\
;	   vars(ivar).ne. "lengi" .and. vars(ivar).ne."aa_gp") then
  if(ismissing(ind(vars_pdf.eq.vars(ivar))))then
	lnres@gsLineColor       = colors(0)
	dum0(ivar,0)			= gsn_add_polyline(wks,plots(ivar,0),fspan(1-2*barwidth/2.,nregimes_used-2*barwidth/2,nregimes_used),regimes_tagdiff(select_idx:,ivar,0),lnres)
	dum1(ivar,0)			= gsn_add_polyline(wks,plots_wgt(ivar,0),fspan(1-2*barwidth/2.,nregimes_used-2*barwidth/2,nregimes_used),regimes_tagdiff_areawgt(select_idx:,ivar,0),lnres)
	dum2(ivar,0)			= gsn_add_polyline(wks,plots_wgt_1(ivar,0),fspan(1-2*barwidth/2.,nregimes_used-2*barwidth/2,nregimes_used),regimes_tagdiff_areawgt_1(select_idx:,ivar,0),lnres)

	lnres@gsLineColor		= colors(2)
	dum0(ivar,2) 			= gsn_add_polyline(wks,plots(ivar,2),fspan(1+2*barwidth/2,nregimes_used+2*barwidth/2,nregimes_used),regimes_diff2(select_idx:,ivar)*scale_factor,lnres)
	dum1(ivar,2) 			= gsn_add_polyline(wks,plots_wgt(ivar,2),fspan(1+2*barwidth/2,nregimes_used+2*barwidth/2,nregimes_used),regimes_diff2_areawgt(select_idx:,ivar)*scale_factor,lnres)
	dum2(ivar,2) 			= gsn_add_polyline(wks,plots_wgt_1(ivar,2),fspan(1+2*barwidth/2,nregimes_used+2*barwidth/2,nregimes_used),regimes_diff2_areawgt_1(select_idx:,ivar)*scale_factor,lnres)
	end if

count = count+1
; Overlay the last three plots on the first one.
       do j=1,dimsizes(case_nor)-1+1
         overlay(plots(ivar,0),plots(ivar,j))
		 overlay(plots_wgt(ivar,0),plots_wgt(ivar,j))
		 overlay(plots_wgt_1(ivar,0),plots_wgt_1(ivar,j))
       end do
end do

pres					= True
pres@gsnMaximize    	= True   ; Maximize in frame
;pres@gsnFrame       	= False  ; Don't advance frame.
;pres@gsnPanelBottom 	= 0.13   ; Leave room for labelbar

pres@gsnPanelLeft		= 0.05
pres@gsnPanelRight		= 0.95
pres@gsnPanelXWhiteSpacePercent		= 3.0
pres@gsnPanelYWhiteSpacePercent		= 3.0

; Create some labelbars for a legend
lbres                    = True          ; labelbar only resources
lbres@lbAutoManage		= True
lbres@vpWidthF           = 0.25          ; labelbar width
lbres@vpHeightF          = 0.10           ; labelbar height
lbres@lbBoxMajorExtentF  = 0.50          ; puts space between color boxes
lbres@lbMonoFillPattern  = True          ; Solid fill pattern
lbres@lbLabelFontHeightF = 0.10        ; font height. default is small
lbres@lbFillColors			= colors

lbres@lbLabelJust        = "TopLeft"  ; left justify labels
;lbres@lbPerimOn          = True

;if(vars(0).eq."sgm_tota")then
;lbres@vpHeightF				= 0.05
;lbres@lbFillColor			= colors(1)
;labels = (/"Gauss-PDF"/)
;else
labels = (/"Park-RH", "Gauss-PDF","PDF-RH"/)
;end if

lbid						= gsn_create_labelbar(wks,dimsizes(labels),labels,lbres)
lbid_wgt					= gsn_create_labelbar(wks,dimsizes(labels),labels,lbres)
lbid_wgt_1					= gsn_create_labelbar(wks,dimsizes(labels),labels,lbres)

amres						= True
amres@amJust				= "BottomLeft"
amres@amParallelPosF		= -0.5
amres@amOrthogonalPosF		= 0.5
annoid						= gsn_add_annotation(plots(0,0),lbid,amres)
annoid_wgt					= gsn_add_annotation(plots_wgt(0,0),lbid_wgt,amres)
annoid_wgt_1				= gsn_add_annotation(plots_wgt_1(0,0),lbid_wgt_1,amres)

; Panel the plots.
;First plot
gsn_panel(wks,plots(:,0),(/dimsizes(plots(:,0))/2,2/),pres)
;drawNDCGrid(wks)

; Second plot
gsn_panel(wks,plots_wgt(:,0),(/dimsizes(plots_wgt(:,0))/2,2/),pres)

; Third plot
gsn_panel(wks,plots_wgt_1(:,0),(/dimsizes(plots_wgt_1(:,0))/2,2/),pres)

delete(count)
exit

end if ; regime_contribution

;=================================================
; ALPI IS DEFINE BY LTS AND PRECIPITATION.
if(ALPI_bin_cal)then

; ----------- 1.1 calculate LTS
LTS					= new((/dimsizes(tag),dimsizes(case_nor),nlat,nlon/),float) ; unit is degreeC
LTS!0		 		= "tag"
LTS&tag				= fspan(0,dimsizes(tag)-1,dimsizes(tag))
LTS!1		 		= "case"
LTS&case			= fspan(0,dimsizes(case_nor)-1,dimsizes(case_nor))

EIS					= LTS

do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1
result					= cal_eis(vars_mid(0,itag,icase,:,:,:),vars_mid(1,itag,icase,:,:,:),lev)
LTS(itag,icase,:,:)		= tofloat(result[0])
EIS(itag,icase,:,:)		= tofloat(result[1])
end do ; do icase=
end do ; do itag=
copy_VarCoords(vars_all(0,:,:,:,:),LTS)
copy_VarCoords(vars_all(0,:,:,:,:),EIS)
printVarSummary(LTS)
printMinMax(LTS,True)
printVarSummary(EIS)
printMinMax(EIS,True)

; ---------- 1.2 get precipitation
do ivar=0,dimsizes(vars)-1
	if(vars(ivar).eq."PRECT")then
		print(vars(ivar))
		prect					= vars_all(ivar,:,:,:,:)
		break
	end if
end do ; do ivar=
printVarSummary(prect)
printMinMax(prect,True)

mask_ocn=True
if(mask_ocn)then
	; addition: mask land region
	LTS_ocn			= LTS
	EIS_ocn			= EIS
	prect_ocn		= prect
	do itag=0,dimsizes(tag)-1
	do icase=0,dimsizes(case_nor)-1
		LTS_ocn(itag,icase,:,:)		= where(OCNFRAC.gt.0.5, LTS(itag,icase,:,:),vars_all@_FillValue)
		EIS_ocn(itag,icase,:,:)		= where(OCNFRAC.gt.0.5, EIS(itag,icase,:,:),vars_all@_FillValue)
		prect_ocn(itag,icase,:,:)	= where(OCNFRAC.gt.0.5, prect(itag,icase,:,:),vars_all@_FillValue)
	end do ; do icase=
	end do ; do itag=
	copy_VarCoords(LTS,LTS_ocn)
	copy_VarCoords(EIS,EIS_ocn)
	copy_VarCoords(prect,prect_ocn)
else
	LTS_ocn			= LTS
	EIS_ocn			= EIS 
	prect_ocn		= prect
end if

LTS_tropics	 = LTS_ocn(:,:,{latS_trop:latE_trop},:)
EIS_tropics	 = EIS_ocn(:,:,{latS_trop:latE_trop},:)
prect_tropics	 = prect_ocn(:,:,{latS_trop:latE_trop},:)

printVarSummary(LTS_tropics)
printVarSummary(EIS_tropics)
printVarSummary(prect_tropics)

; define ALPI array
ALPI				= LTS_tropics

; --------- 1.3 get LTS_distance and Prec_distance
do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1
;---1.3.1.1 find LTS_max and LTS_min
LTS_max				= max(LTS_tropics(itag,icase,:,:))
LTS_min				= min(LTS_tropics(itag,icase,:,:))
;---1.3.1.2 find EIS_max and EIS_min
EIS_max				= max(EIS_tropics(itag,icase,:,:))
EIS_min				= min(EIS_tropics(itag,icase,:,:))
;---1.3.2 find prec_max
prect_max			= max(prect_tropics(itag,icase,:,:))
prect_min			= min(prect_tropics(itag,icase,:,:))

;---1.3.3 LTS_distance
LTS_max_2D				= conform_dims(dimsizes(LTS_tropics(itag,icase,:,:)),LTS_max,-1)
copy_VarCoords(LTS_tropics(itag,icase,:,:),LTS_max_2D)
LTS_min_2D				= conform_dims(dimsizes(LTS_tropics(itag,icase,:,:)),LTS_min,-1)
copy_VarCoords(LTS_tropics(itag,icase,:,:),LTS_min_2D)
LTS_distance			= (LTS_max_2D-LTS_tropics(itag,icase,:,:))/(LTS_max_2D-LTS_min_2D)+0.1
printMinMax(LTS_distance,True)

;---1.3.3.1 EIS_distance
EIS_max_2D				= conform_dims(dimsizes(EIS_tropics(itag,icase,:,:)),EIS_max,-1)
copy_VarCoords(EIS_tropics(itag,icase,:,:),EIS_max_2D)
EIS_min_2D				= conform_dims(dimsizes(EIS_tropics(itag,icase,:,:)),EIS_min,-1)
copy_VarCoords(EIS_tropics(itag,icase,:,:),EIS_min_2D)
EIS_distance			= (EIS_max_2D-EIS_tropics(itag,icase,:,:))/(EIS_max_2D-EIS_min_2D)+0.1
printMinMax(EIS_distance,True)

;---1.3.4 Prec_distance
prect_max_2D				= conform_dims(dimsizes(prect_tropics(itag,icase,:,:)),prect_max,-1)
copy_VarCoords(prect_tropics(itag,icase,:,:),prect_max_2D)
prect_min_2D				= conform_dims(dimsizes(prect_tropics(itag,icase,:,:)),prect_min,-1)
copy_VarCoords(prect_tropics(itag,icase,:,:),prect_min_2D)
prect_distance				= (prect_max_2D-prect_tropics(itag,icase,:,:))/(prect_max_2D-prect_min_2D)+0.1
printMinMax(prect_distance,True)

;------------ 1.4 get ALPI
if(ALPIsub_revised)then
ALPI(itag,icase,:,:)		= atan2(prect_distance,EIS_distance)
else
ALPI(itag,icase,:,:)		= atan2(prect_distance,LTS_distance)
end if

end do ; do icase=
end do ; do itag=

r2d							= 45.0/atan(1.0)     ; conversion factor (radians to degrees)
ALPI						= ALPI*r2d

printVarSummary(ALPI) ; (tag, case, lat, lon)
printMinMax(ALPI,True)

; =========== bin based on ALPI percentile
; ----- 1.1 sorting based on ALPI index: increasing

full_length					= dimsizes(ndtooned(ALPI(0,0,{latS_trop:latE_trop},:)))
nlength						= num(.not.ismissing(ndtooned(ALPI(0,0,{latS_trop:latE_trop},:))))
ALPI_bins					= new((/dimsizes(tag),dimsizes(case_nor),nbins/),float)
vars_all_bin				= new((/dimsizes(vars),dimsizes(tag),dimsizes(case_nor),full_length,nbins/),float)
vars_all_bin_mean			= new((/dimsizes(vars),dimsizes(tag),dimsizes(case_nor),nbins/),float)
vars_all_bin_std			= new((/dimsizes(vars),dimsizes(tag),dimsizes(case_nor),nbins/),float)

ALPI_mask					= new((/dimsizes(tag),dimsizes(case_nor),full_length/),float)

print(nlength)

do itag=0,dimsizes(tag)-1
do icase=0,dimsizes(case_nor)-1
	ALPI_1D						= ndtooned(ALPI(itag,icase,{latS_trop:latE_trop},:))
	ALPI_1D_sorted				= ALPI_1D
	ip							= dim_pqsort(ALPI_1D_sorted,2) ; sorting and change the raw data
	printVarSummary(ALPI_1D_sorted)

do ivar=0,dimsizes(vars)-1
	vars_all_1D					= ndtooned(vars_all(ivar,itag,icase,{latS_trop:latE_trop},:))
	printVarSummary(vars_all_1D)

	lat							= vars_all&lat
	clat                		= tofloat(cos(0.01745329*lat))
	clat!0              		= "lat"
	clat&lat            		= vars_all&lat
	clat_2D						= conform_dims(dimsizes(vars_all(ivar,itag,icase,:,:)),clat,0)
	copy_VarCoords(vars_all(ivar,itag,icase,:,:),clat_2D)
	clat_2D_1D					= ndtooned(clat_2D({latS_trop:latE_trop},:))

;	print(ALPI_1D_sorted+"		"+clat_2D_1D(ip))

	do ibin=0,nbins-1
		tmp1		= ALPI_1D_sorted(toint(nlength/nbins)*ibin)
		if(ibin.eq.(nbins-1))then ; avoid missing some points at the end --- notion: pls dont write "tmp2 = ALPI_1D(dimsizes(ALPI_1D)-1)" when missing value exists.
			tmp2		= ALPI_1D_sorted(nlength-1)
		else
			tmp2		= ALPI_1D_sorted(toint(nlength/nbins)*(ibin+1))
		end if
		vars_all_bin(ivar,itag,icase,:,ibin)		= where(ALPI_1D_sorted.gt.tmp1 .and. ALPI_1D_sorted.lt.tmp2, vars_all_1D(ip), vars_all@_FillValue)
;		print(num(.not.ismissing(vars_all_bin(ivar,itag,icase,:,ibin)))) ; 1366

;		vars_all_bin_mean(ivar,itag,icase,ibin)		= dim_avg(vars_all_bin(ivar,itag,icase,:,ibin))
		vars_all_bin_mean(ivar,itag,icase,ibin)		= dim_avg_wgt_n_Wrap(vars_all_bin(ivar,itag,icase,:,ibin),clat_2D_1D(ip),1,0) 	; the weighting effect is small. 2018.12.18
		vars_all_bin_std(ivar,itag,icase,ibin)		= dim_stddev(vars_all_bin(ivar,itag,icase,:,ibin))

; ========= very brillient idea suggested by Jishi! in order to reverse to the Lat-Lon plot, we need to keep the ip index always.
		ALPI_mask(itag,icase,ip(ind(ALPI_1D_sorted.gt.tmp1 .and. ALPI_1D_sorted.lt.tmp2)))	=  ALPI_mask_val(ibin) 
;		ALPI_mask(itag,icase,:)			= where(ALPI_1D.gt.tmp1 .and. ALPI_1D.lt.tmp2, ALPI_mask_val(ibin), ALPI_mask(itag,icase,:))  ; --- cannot work.
; ========
	end do ; do ibin=

;; -------- used to test the binned variables
;wks						= gsn_open_wks(plot_type,"./figure/ALPI_binned_variable_test")
;
;dims					= dimsizes(vars_all_bin(ivar,itag,icase,:,:))
;;tmps					= reshape(vars_all_bin(ivar,itag,icase,:,:),(/dims(1),dims(0)/)) ;--- reshape cannot work here, because it does not change the storage sequency of data actually.
;tmps					= transpose(vars_all_bin(ivar,itag,icase,:,:))
;
;res						= True
;res@vpWidthF			= 0.7
;res@vpHeightF			= 0.4
;res@xyLineColors		= ispan(2,250,250/dims(1))
;res@tiYAxisString		= vars(ivar)
;res@tiXAxisString		= "ndtooned vars [Lat-Lon]"
;
;plot					= gsn_csm_y(wks,tmps,res)
;exit
;delete([/plot,tmps,dims/])
;; ---------- 

delete([/tmp1,tmp2/])
end do ; do ivar=
end do ; do icase=
end do ; do itag=


printVarSummary(vars_all_bin)
printMinMax(vars_all_bin,True)
printVarSummary(vars_all_bin_mean)
printMinMax(ALPI_mask,True)

; get mask variable to verify the correctness of ALPI bin
dims			= dimsizes(ALPI(0,0,:,:))
; reshape into Lat-Lon dimension variable
ALPI_mask_latlon								= reshape(ALPI_mask,(/dimsizes(tag),dimsizes(case_nor),dims(0),dims(1)/))
copy_VarCoords(ALPI,ALPI_mask_latlon(:,:,:,:))
printVarSummary(ALPI_mask_latlon)
; get the difference to identify the distinction of Gauss-PDF and Park-RH
ALPI_mask_latlon_diff							= (/ALPI_mask_latlon(:,1,:,:)-ALPI_mask_latlon(:,0,:,:)/)
copy_VarCoords(ALPI_mask_latlon(:,0,:,:),ALPI_mask_latlon_diff)
ALPI_mask_latlon_diff2							= (/ALPI_mask_latlon(1,:,:,:)-ALPI_mask_latlon(0,:,:,:)/)
copy_VarCoords(ALPI_mask_latlon(0,:,:,:),ALPI_mask_latlon_diff2)

if(ALPIsub_mask_plot)then

; ======== the following part is used to verify whether the ALPI index sorting is right.
; ======== ---- yes, plots indicate that all regimes are right.

;wks			= gsn_open_wks(plot_type,"./figure/ALPI_mask_plot_nbins_"+nbins+"_Raster")
wks			= gsn_open_wks(plot_type,"./figure/ALPI_mask_plot_nbins_"+nbins)

gsn_define_colormap(wks,"NCV_manga")
;gsn_define_colormap(wks,"NCV_blue_red")

res							= True

res@gsnFrame				= False
res@gsnDraw					= False

res@cnFillOn				= True
res@cnLinesOn				= False
res@cnLineLabelsOn			= False
;res@cnFillMode       		= "RasterFill"
res@cnLevelSelectionMode	= "ExplicitLevels"
res@cnLevels         		= ispan(1,nbins-1,1);ALPI_mask_val(1:)
res@cnFillColors			= ispan(2,250,250/dimsizes(ALPI_mask_val)) ;(/"lightcyan2","lightgoldenrod2","lightyellow2","lightgrey","lightskyblue2","lightseagreen","lightpink2"/)

res@lbLabelBarOn			= False
res@lbLabelPosition  		= "Center"           ; label position
res@lbLabelAlignment 		= "BoxCenters"       ; label orientation
;res@lbLabelStrings   		= ""+ res@cnLevels   ; trick 
res@lbLabelStrings   		= ""+ (ALPI_mask_val(:)-0.5) ; --- slightly tricky here: to make the labeled value look clear.

fontheight					= 0.025
res@lbLabelFontHeightF		= fontheight

res@pmLabelBarHeightF		=  0.100
res@pmLabelBarWidthF 		=  0.80

res@mpCenterLonF			= 180
res@mpMinLatF				= -30.
res@mpMaxLatF				= 30.

res@tmXBLabelStride			= 2.
res@gsnStringFontHeightF	= fontheight
res@gsnLeftString			= "ALPI regimes"
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
plot_diff(itag)				= gsn_csm_contour_map_ce(wks, ALPI_mask_latlon_diff(itag,:,:),resd)
count_1 = count_1 + 1

do icase=0,dimsizes(case_nor)-1
res@gsnLeftString	= plots_strings(count_2)
res@gsnRightString = tag(itag)+"_"+case_name(icase)
plot(itag,icase)		= gsn_csm_contour_map_ce(wks,ALPI_mask_latlon(itag,icase,:,:),res)
count_2	= count_2 + 1
end do ; do itag=
end do ; do icase=

do icase=0,dimsizes(case_nor)-1
resd@gsnLeftString			= plots_strings(count_3)
resd@gsnRightString			= "p4K-cntl_"+case_name(icase)
plot_diff2(icase)			= gsn_csm_contour_map_ce(wks, ALPI_mask_latlon_diff2(icase,:,:),resd)
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

gsn_panel(wks,ndtooned(plot),(/2,2/),pres)
gsn_panel(wks,plot_diff,(/2,1/),presd)
gsn_panel(wks,plot_diff2,(/2,1/),presd)
delete([/count_1, count_2, count_3 /])
exit
end if ; ALPI_mask_plot

; ----------------------------------------------------
; -------- ALPI_binned_vars
; ----------------------------------------------------
if(ALPIsub_bin_plot)then
if(ALPIsub_revised)then
wks  						= gsn_open_wks (plot_type,"./figure/ALPI_binned_vars_nbins_"+nbins+"_revised_ALPI"+"_"+vars(dimsizes(vars)-1))
else
wks  						= gsn_open_wks (plot_type,"./figure/ALPI_binned_vars_nbins_"+nbins)
end if

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
;res@tmXBLabels				= fspan(0,nbins,nbins+1)*100/nbins
res@tmXBLabels				= fspan(0,nbins,nbins+1)*100/nbins
res@trXMaxF					= (nbins)

res@tiYAxisFontHeightF			= 0.03
res@tiXAxisFontHeightF			= 0.03
res@tmXBLabelFontHeightF		= 0.03
res@tmYLLabelFontHeightF		= 0.03

res@gsnStringFontHeightF	= 0.03

res@tmYLFormat				= "0*+^e"
res@tmYLMinorOn				= False
res@tmXBPrecision			= 6

res_diff					= res

plot						= new(dimsizes(vars),graphic)
error_bar					= new((/dimsizes(vars),dimsizes(tag)*dimsizes(case_nor),nbins/),graphic)
upp_bound					= new((/dimsizes(vars),dimsizes(tag)*dimsizes(case_nor),nbins/),graphic)
dwn_bound					= new((/dimsizes(vars),dimsizes(tag)*dimsizes(case_nor),nbins/),graphic)
centers						= new((/dimsizes(vars),dimsizes(tag)*dimsizes(case_nor),nbins/),graphic)
lines						= new((/dimsizes(vars),dimsizes(tag)*dimsizes(case_nor)/),graphic)
polygons					= new((/dimsizes(vars),dimsizes(tag)*dimsizes(case_nor)/),graphic)

plot_diff					= new(dimsizes(vars),graphic)
lines_diff					= new((/dimsizes(vars),dimsizes(case_nor)+1/),graphic)
lines_ref0					= new((/dimsizes(vars),dimsizes(case_nor)+1/),graphic)
centers_diff				= new((/dimsizes(vars),dimsizes(case_nor)+1,nbins/),graphic)

plot_diff_sgm					= new(dimsizes(vars),graphic)
lines_diff_sgm					= new((/dimsizes(vars),dimsizes(case_nor)+1/),graphic)
lines_ref0_sgm					= new((/dimsizes(vars),dimsizes(case_nor)+1/),graphic)
centers_diff_sgm				= new((/dimsizes(vars),dimsizes(case_nor)+1,nbins/),graphic)



count_1		= 0
count_2		= 0

;do ivar=0,dimsizes(vars)-1
if(vars(1).eq."sgm_tota")then
start_idx	= 1
else
start_idx	= 0
end if

do ivar=start_idx,dimsizes(vars)-1 ; don't draw PRECT

res@gsnLeftString			= plots_strings(count_1)
res@tiYAxisString			= vars(ivar)+" ("+vars_units(ivar)+")"
tmp							= reshape(vars_all_bin_mean(ivar,:,:,:),(/dimsizes(tag)*dimsizes(case_nor),nbins/))
tmp2						= reshape(vars_all_bin_std(ivar,:,:,:),(/dimsizes(tag)*dimsizes(case_nor),nbins/))

tmp_diff1					= (/vars_all_bin_mean(ivar,1,:,:)-vars_all_bin_mean(ivar,0,:,:)/) ; difference between p4K and cntl
tmp_diff2					= (/tmp_diff1(1,:)-tmp_diff1(0,:)/) ; difference between PDF_p4K-cntl and RH_p4K-cntl
tmp_diff_final				= new((/dimsizes(case_nor)+1,nbins/),float)
tmp_diff_final(0,:)			= tmp_diff1(0,:)
tmp_diff_final(1,:)			= tmp_diff1(1,:)
tmp_diff_final(2,:)			= tmp_diff2(:)

; --- make the sgm related variables share the same Y axis range
if(start_idx.eq.1 .and. (vars(ivar).eq."sgm_tota" .or. vars(ivar).eq."sgm_shal" .or. vars(ivar).eq."sgm_turb"))then
sgm_all_max					= max(vars_all_bin_mean((/1,3,4/),:,:,:)+vars_all_bin_std((/1,3,4/),:,:,:))
sgm_all_min					= min(vars_all_bin_mean((/1,3,4/),:,:,:)-vars_all_bin_std((/1,3,4/),:,:,:))
res@trYMaxF					= sgm_all_max
res@trYMinF					= sgm_all_min
else
res@trYMaxF					= max(tmp+tmp2)
res@trYMinF					= min(tmp-tmp2)
end if

;plot(ivar) 					= gsn_csm_y(wks, tmp, res)
if(start_idx.eq.1)then
plot(ivar-1) 					= gsn_csm_y(wks, tmp, res)
else
plot(ivar) 					= gsn_csm_y(wks, tmp, res)
end if
count_1 = count_1 + 1

;;;;;;;;;;;;;;;;;
res_diff@gsnLeftString			= plots_strings(count_2)
res_diff@tiYAxisString			= vars(ivar)+" ("+vars_units(ivar)+")"

; --- make the sgm related variables share the same Y axis range
if(start_idx.eq.1)then
	if(vars(ivar).eq."sgm_tota" .or. vars(ivar).eq."sgm_shal" .or. vars(ivar).eq."sgm_turb")then
	sgm_all_max							= max(vars_all_bin_mean((/1,3,4/),1,:,:)-vars_all_bin_mean((/1,3,4/),0,:,:))
	sgm_all_min							= min(vars_all_bin_mean((/1,3,4/),1,:,:)-vars_all_bin_mean((/1,3,4/),0,:,:))
	res_diff@trYMaxF					= sgm_all_max
	res_diff@trYMinF					= sgm_all_min
	else
	res_diff@trYMaxF					= max(vars_all_bin_mean(ivar,1,:,:)-vars_all_bin_mean(ivar,0,:,:))
	res_diff@trYMinF					= min(vars_all_bin_mean(ivar,1,:,:)-vars_all_bin_mean(ivar,0,:,:))
	end if
end if

; resource for sgm variables, because I will overlay them on one plot
res_diff_sgm					= res_diff
if(start_idx.eq.1 .and. (ivar-1).eq.0)then
res_diff_sgm@gsnLeftString		= ""
res_diff_sgm@tiYAxisString		= "SGM"
else
res_diff_sgm@gsnLeftString		= ""
end if

;plot_diff(ivar)				= gsn_csm_y(wks, tmp_diff_final, res_diff)
if(start_idx.eq.1)then
plot_diff(ivar-1)				= gsn_csm_y(wks, tmp_diff_final, res_diff)
plot_diff_sgm(ivar-1)			= gsn_csm_y(wks, tmp_diff_final, res_diff_sgm)
else
plot_diff(ivar)					= gsn_csm_y(wks, tmp_diff_final, res_diff)
;plot_diff_sgm(ivar)				= gsn_csm_y(wks, tmp_diff_final, res_diff_sgm)
end if
count_2	= count_2 + 1


;-------- add error bar ---------
polyres						= True
polyres@gsMarkerIndex		= 1
polyres@gsMarkerSizeF		= 0.02

polyres_sgm					= polyres

colors						= (/"red","red","blue","blue","black","black"/)
dashes						= (/0,12,0,12,0,12/)
thickes						= (/2,2,2,2,4,4/)
									;PRECT sgm_tota, N1, sgm_shal, sgm_turb
colors_sgm						= (/"red","red","red","blue","cyan","black","black"/)
dashes_sgm						= (/0,12,0,12,0,12/)
thickes_sgm						= (/2,2,2,2,4,4/)

;polyres@gsLineColor			= (/"red","blue","lightpink","cyan"/)

; loop through the points
do idx1=0,dimsizes(tag)*dimsizes(case_nor)-1
;if((vars(ivar).eq."sgm_tota" .or. vars(ivar).eq."N1" .or. \
;    vars(ivar).eq."sgm_shal" .or. vars(ivar).eq."sgm_turb" .or. \
;	vars(ivar).eq."wstarPBL" .or. vars(ivar).eq."umf_shal" .or. \
;	vars(ivar).eq."clddep2" .or. vars(ivar).eq."lengi" .or. \
;	vars(ivar).eq."shi" .or. vars(ivar).eq."aa_gp") .and. \
;	(idx1.eq.0 .or. idx1.eq.2))then ; -- avoid missing value for sgm_tota and N1 to be plotted
	  if((.not.ismissing(ind(vars_pdf.eq.vars(ivar)))) .and. (idx1.eq.0 .or. idx1.eq.2))then
	continue
else
polyres@gsLineColor			= colors(idx1)
polyres@gsLineDashPattern	= dashes(idx1)
polyres@gsLineThicknessF	= thickes(idx1)

;lines(ivar,idx1)				= gsn_add_polyline(wks,plot(ivar),ispan(0,nbins-1,1)+0.5+idx1*0.1,tmp(idx1,:),polyres)
;lines(ivar,idx1)				= gsn_add_polyline(wks,plot(ivar),ispan(0,nbins-1,1)+0.5,tmp(idx1,:),polyres)
if(start_idx.eq.1)then
lines(ivar-1,idx1)				= gsn_add_polyline(wks,plot(ivar-1),ispan(0,nbins-1,1)+0.5,tmp(idx1,:),polyres)
else
lines(ivar,idx1)				= gsn_add_polyline(wks,plot(ivar),ispan(0,nbins-1,1)+0.5,tmp(idx1,:),polyres)
end if

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

do idx2=0,nbins-1
;	centers(ivar,idx1,idx2)		= gsn_add_polymarker(wks,plot(ivar),idx2+0.5+idx1*0.1,tmp(idx1,idx2),polyres)
;	error_bar(ivar,idx1,idx2)	= gsn_add_polyline(wks,plot(ivar),(/idx2,idx2/)+0.5+idx1*0.1,(/tmp(idx1,idx2)+tmp2(idx1,idx2),tmp(idx1,idx2)-tmp2(idx1,idx2)/),polyres)
;	upp_bound(ivar,idx1,idx2)	= gsn_add_polyline(wks,plot(ivar),(/idx2-0.1,idx2+0.1/)+0.5+idx1*0.1,(/tmp(idx1,idx2)+tmp2(idx1,idx2),tmp(idx1,idx2)+tmp2(idx1,idx2)/),polyres)
;	dwn_bound(ivar,idx1,idx2)	= gsn_add_polyline(wks,plot(ivar),(/idx2-0.1,idx2+0.1/)+0.5+idx1*0.1,(/tmp(idx1,idx2)-tmp2(idx1,idx2),tmp(idx1,idx2)-tmp2(idx1,idx2)/),polyres)

if(start_idx.eq.1)then
	error_bar(ivar-1,idx1,idx2)	= gsn_add_polyline(wks,plot(ivar-1),(/idx2,idx2/)+0.5,(/tmp(idx1,idx2)+tmp2(idx1,idx2),tmp(idx1,idx2)-tmp2(idx1,idx2)/),polyres)
	upp_bound(ivar-1,idx1,idx2)	= gsn_add_polyline(wks,plot(ivar-1),(/idx2-0.1,idx2+0.1/)+0.5,(/tmp(idx1,idx2)+tmp2(idx1,idx2),tmp(idx1,idx2)+tmp2(idx1,idx2)/),polyres)
	dwn_bound(ivar-1,idx1,idx2)	= gsn_add_polyline(wks,plot(ivar-1),(/idx2-0.1,idx2+0.1/)+0.5,(/tmp(idx1,idx2)-tmp2(idx1,idx2),tmp(idx1,idx2)-tmp2(idx1,idx2)/),polyres)
else
	error_bar(ivar,idx1,idx2)	= gsn_add_polyline(wks,plot(ivar),(/idx2,idx2/)+0.5,(/tmp(idx1,idx2)+tmp2(idx1,idx2),tmp(idx1,idx2)-tmp2(idx1,idx2)/),polyres)
	upp_bound(ivar,idx1,idx2)	= gsn_add_polyline(wks,plot(ivar),(/idx2-0.1,idx2+0.1/)+0.5,(/tmp(idx1,idx2)+tmp2(idx1,idx2),tmp(idx1,idx2)+tmp2(idx1,idx2)/),polyres)
	dwn_bound(ivar,idx1,idx2)	= gsn_add_polyline(wks,plot(ivar),(/idx2-0.1,idx2+0.1/)+0.5,(/tmp(idx1,idx2)-tmp2(idx1,idx2),tmp(idx1,idx2)-tmp2(idx1,idx2)/),polyres)
end if


end do ; do idx2=

end if ; vars(ivar).eq."sgm_tota"...
end do ; do idx1=

do idx1=0,dimsizes(case_nor)
;if((vars(ivar).eq."sgm_tota" .or. vars(ivar).eq."N1" .or. \
;    vars(ivar).eq."sgm_shal" .or. vars(ivar).eq."sgm_turb" .or. \
;    vars(ivar).eq."umf_shal" .or. vars(ivar).eq."wstarPBL" .or. vars(ivar).eq."clddep2".or.\
;	vars(ivar).eq."lengi" .or. vars(ivar).eq."shi" .or. vars(ivar).eq."aa_gp") .and. (idx1.eq.0 .or. idx1.eq.2))then
 if((.not.ismissing(ind(vars_pdf.eq.vars(ivar)))).and.(idx1.eq.0 .or. idx1.eq.2))then
	continue
else
	polyres@gsLineColor					= colors(idx1*2)
	polyres@gsLineDashPattern 			= dashes(idx1*2)
	polyres@gsLineThicknessF			= thickes(idx1*2)

	if(start_idx.eq.1)then
	polyres_sgm@gsLineColor					= colors_sgm(ivar)
	polyres_sgm@gsLineDashPattern 			= dashes_sgm(idx1*2)
	polyres_sgm@gsLineThicknessF			= thickes_sgm(idx1*2)
	end if

;	lines_diff(ivar,idx1)				= gsn_add_polyline(wks,plot_diff(ivar),ispan(0,nbins-1,1)+0.5+idx1*0.1,tmp_diff_final(idx1,:),polyres)
if(start_idx.eq.1)then
	lines_diff(ivar-1,idx1)				= gsn_add_polyline(wks,plot_diff(ivar-1),ispan(0,nbins-1,1)+0.5,tmp_diff_final(idx1,:),polyres)
	lines_diff_sgm(ivar-1,idx1)			= gsn_add_polyline(wks,plot_diff_sgm(ivar-1),ispan(0,nbins-1,1)+0.5,tmp_diff_final(idx1,:),polyres_sgm)
else
	lines_diff(ivar,idx1)				= gsn_add_polyline(wks,plot_diff(ivar),ispan(0,nbins-1,1)+0.5,tmp_diff_final(idx1,:),polyres)
;	lines_diff_sgm(ivar,idx1)				= gsn_add_polyline(wks,plot_diff_sgm(ivar),ispan(0,nbins-1,1)+0.5,tmp_diff_final(idx1,:),polyres_sgm)
end if

	do idx2=0,nbins-1
		if(start_idx.eq.1)then
				centers_diff(ivar-1,idx1,idx2)	= gsn_add_polymarker(wks,plot_diff(ivar-1),idx2+0.5,tmp_diff_final(idx1,idx2),polyres)
				centers_diff_sgm(ivar-1,idx1,idx2)	= gsn_add_polymarker(wks,plot_diff_sgm(ivar-1),idx2+0.5,tmp_diff_final(idx1,idx2),polyres_sgm)
		else
				centers_diff(ivar,idx1,idx2)	= gsn_add_polymarker(wks,plot_diff(ivar),idx2+0.5,tmp_diff_final(idx1,idx2),polyres)
;				centers_diff_sgm(ivar,idx1,idx2)	= gsn_add_polymarker(wks,plot_diff_sgm(ivar),idx2+0.5,tmp_diff_final(idx1,idx2),polyres_sgm)
		end if
	end do ; do idx2=

	polyres@gsLineColor					= "grey"
	polyres@gsLineDashPattern			= 14
	
	polyres_sgm@gsLineColor					= "grey"
	polyres_sgm@gsLineDashPattern			= 14

	if(start_idx.eq.1)then
	lines_ref0(ivar-1,idx1)				= gsn_add_polyline(wks,plot_diff(ivar-1),ispan(0,nbins,1),fspan(0,0,nbins+1),polyres)
	lines_ref0_sgm(ivar-1,idx1)				= gsn_add_polyline(wks,plot_diff_sgm(ivar-1),ispan(0,nbins,1),fspan(0,0,nbins+1),polyres_sgm)
	else
	lines_ref0(ivar,idx1)				= gsn_add_polyline(wks,plot_diff(ivar),ispan(0,nbins,1),fspan(0,0,nbins+1),polyres)
;	lines_ref0_sgm(ivar,idx1)				= gsn_add_polyline(wks,plot_diff_sgm(ivar),ispan(0,nbins,1),fspan(0,0,nbins+1),polyres_sgm)
	end if
end if ; vars(ivar).eq."sgm_tota"...
end do ; do idx1=

end do ; do ivar=

if(start_idx.eq.1)then
overlay(plot_diff_sgm(0),plot_diff_sgm(2))
overlay(plot_diff_sgm(0),plot_diff_sgm(3))
end if

; ------------------------Add legends
;legends                         = vars
legends							= (/"RH_cntl","PDF_cntl","RH_p4K","PDF_p4K"/)
legends_diff1					= (/"RH_p4K-cntl", "PDF_p4K-cntl","PDF-RH_p4K-cntl"/)
legends_sgm						= (/"SGM_tota", "SGM_shal", "SGM_turb"/)

lgres							= True

lgres@lgLabelFontHeightF        = 0.10
lgres@vpWidthF                  = 0.25
lgres@vpHeightF                 = 0.15
lgres@lgPerimOn                 = False

lgres1							= lgres
lgres2							= lgres
lgres_sgm						= lgres

lgres1@lgLineColors				= colors
lgres1@lgLineThicknessF			= res@xyLineThicknessF
lgres1@lgDashIndexes			= dashes

lgres2@lgLineColors				= colors(::2)
lgres2@lgLineThicknessF			= res@xyLineThicknessF
lgres2@lgDashIndexes			= dashes(::2)

lgres_sgm@lgLineColors				= colors_sgm((/1,3,4/))
lgres_sgm@lgLineThicknessF			= res@xyLineThicknessF
lgres_sgm@lgDashIndexes				= dashes_sgm(::2)

lbid1                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres1)
lbid2                           = gsn_create_legend(wks,dimsizes(legends_diff1),legends_diff1,lgres2)
if(start_idx.eq.1)then
lbid3                           = gsn_create_legend(wks,dimsizes(legends_sgm),legends_sgm,lgres_sgm)
end if

amres                           = True
amres@amJust					= "BottomLeft"
amres@amParallelPosF            = -0.5
amres@amOrthogonalPosF          = 0.5

amres_sgm						= amres
amres_sgm@amJust				= "BottomRight"
amres_sgm@amParallelPosF		= 0.5
amres_sgm@amOrthogonalPosF		= 0.5

if(vars(0).eq."PRECT")then
	amres@amJust					= "TopRight"
	amres@amParallelPosF			= 0.5
	amres@amOrthogonalPosF			= -0.5
end if

annoid1                         = gsn_add_annotation(plot(0),lbid1,amres)
annoid2                         = gsn_add_annotation(plot_diff(0),lbid2,amres)
if(start_idx.eq.1)then
annoid3                         = gsn_add_annotation(plot_diff_sgm(0),lbid3,amres_sgm)
end if

pres							= True
pres@gsnPanelYWhiteSpacePercent = 5.0
pres@gsnPanelXWhiteSpacePercent = 5.0
pres@gsnPanelLeft				= 0.05
pres@gsnPanelRight				= 0.95

if(start_idx.eq.1)then
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

;gsn_panel(wks,plot_diff,(/dimsizes(vars)/2,2/),pres)
;gsn_panel(wks,plot,(/dimsizes(vars)/2,2/),pres)
if(dimsizes(vars).le.6)then
if(mod(dimsizes(vars),2).ne.0)then
	gsn_panel(wks,plot_diff,(/dimsizes(vars)/2+1,2/),pres)
	gsn_panel(wks,plot,(/dimsizes(vars)/2+1,2/),pres)
else
	gsn_panel(wks,plot_diff,(/dimsizes(vars)/2,2/),pres)
	gsn_panel(wks,plot,(/dimsizes(vars)/2,2/),pres)
end if
else
if(mod(dimsizes(vars),4).ne.0)then
	gsn_panel(wks,plot_diff,(/dimsizes(vars)/4+1,4/),pres)
	gsn_panel(wks,plot,(/dimsizes(vars)/4+1,4/),pres)
else
	gsn_panel(wks,plot_diff,(/dimsizes(vars)/4,4/),pres)
	gsn_panel(wks,plot,(/dimsizes(vars)/4,4/),pres)
end if

end if


delete([/count_1, count_2/])
exit
end if ; ALPI_bin_plot

end if ; ALPI_bin_cal

;************************create plot***************
if(lat_lon_plot)then
do ivar=0,dimsizes(vars)-1

;wks 					= gsn_open_wks(plot_type,"figure/lat-lon-"+vars(ivar)+"-HammerProj")
wks 					= gsn_open_wks(plot_type,"figure/lat-lon-"+vars(ivar)+"-p4K-"+vars(dimsizes(vars)-1))

plot 					= new(9,graphic)
plot1 					= new(9,graphic)

;gsn_define_colormap(wks,"BlueDarkRed18")

cmap1 					= read_colormap_file("WhiteBlueGreenYellowRed")
cmap2 					= read_colormap_file("BlueRed")
;cmap2					= read_colormap_file("BlueWhiteOrangeRed")
;cmap2				 	= read_colormap_file("NCV_roullet")
;cmap2				 	= read_colormap_file("NCV_blue_red")

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
;resl@gsnLeftString				= ""
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


count	= 0

do it=0,8
if(it.eq.0.or.it.eq.1.or.it.eq.3.or.it.eq.4)then
	res@gsnLeftString			= plots_strings(count)+" "+strings4K(it)
;	plot(it)				= gsn_csm_contour_map(wks,var4K(it,:,:),res)
	plot(it)				= gsn_csm_contour_map(wks,smth9_Wrap(var4K(it,:,:),p,q,False),res)

;	plot_ov					= gsn_csm_contour(wks,var4K(it,:,:),resl)
;	overlay(plot(it),plot_ov)
else
	if(it.eq.6 .or. it.eq.7)then
		res3@gsnLeftString			= plots_strings(count)+" "+strings4K(it)
		plot(it)					= gsn_csm_contour_map(wks,smth9_Wrap(var4K(it,:,:),p,q,False),res3)
	else
		res2@gsnLeftString			= plots_strings(count)+" "+strings4K(it)
;		if(it.eq.8)then
;			res2@gsnCenterString	= sprintf("%5.2f",ra(ivar))
;		else
;			res2@gsnCenterString	= ""
;		end if
		plot(it)				= gsn_csm_contour_map(wks,smth9_Wrap(var4K(it,:,:),p,q,False),res2)
	end if
end if
count = count + 1
end do ; do it=


resP = True
resP@gsnPanelMainString = vars(ivar)

resP@gsnPanelXWhiteSpacePercent		= 2.5
resP@gsnPanelYWhiteSpacePercent		= 5.

;if(vars(ivar).eq."SWCF".or.vars(ivar).eq."CLDLOW")then
gsn_panel(wks,plot,(/3,3/),resP)
delete(count)
;gsn_panel(wks,plot1,(/3,3/),resP)
;end if

end do ; do ivar=

end if ; lat_lon_plot
end

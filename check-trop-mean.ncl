
casename=(/\
		"FAMIPC5_f19f19_MG15_amip", "FAMIPC5_f19f19_MG15_amip-p4K"\
		,"FAMIPC5_f19f19_mac2_amip_outSGM", "FAMIPC5_f19f19_mac2_amip-p4K_debug_outSGM"\
		/)

datadir="/home/share3/lyl/work3/qinyi/mid-data/cfmip/"

do icase=0,dimsizes(casename)-1
	f1		= addfile(datadir+casename(icase)+".ANN.climo.tropmean.nc","r")
	tmp		= f1->CLDTOT_ISCCP
	print("casename="+casename(icase)+"		CLDTOT_ISCCP="+tmp)

end do ; do icase=

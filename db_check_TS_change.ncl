
; 2018-12-10
; enlighted by Webb et al. (2015), in which indicated the temperature change is the near surface temperature change \delta_TS
; so I double check whether the TS indeed changes much between two simulations.
; shell script was run to get the global mean value of TS and their difference between p4K and cntl experiments.
; here, just simply print these values on the screen.

fdir	= "/home/lyl/WORK3/qinyi/mid-data/cfmip/"
casename=(/"MG10_diff_TS_globmean.nc","MG15_diff_TS_globmean.nc","Mac2_MG15_diff_TS_globmean.nc","Mac2_MG10_debugoff_diff_TS_globmean.nc","Mac2_MG10_debugon_diff_TS_globmean.nc"/)

;vars	= new((/dimsizes(casename)/),float)
;do icase=0,dimsizes(casename)-1

f1			= addfiles(fdir+casename,"r")

vars		= f1[:]->TS

;end do ; do icase=

printVarSummary(vars)
print(casename+"		"+vars)



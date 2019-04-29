function read_data_monthly(fdir:string,fappend:string,case:string,var:string,landmask:logical,ocnmask:logical \
						   , yearS:numeric,yearE:numeric \
						   )

begin 
years		= ispan(yearS,yearE,1)

; ------ get monthly data from specific start and end years --------- starts
do iyr=0,dimsizes(years)-1
yr_4d		= sprinti("%04d",years(iyr))

if(iyr.eq.0)then
	tmp			= systemfunc("ls "+fdir+case+fappend+"*"+yr_4d+"*.nc")
	tmp1		= tmp
else
	tmp			= systemfunc("ls "+fdir+case+fappend+"*"+yr_4d+"*.nc")
	tmp2		= array_append_record(tmp1,tmp,0)
	delete(tmp1)
	tmp1		= tmp2
	delete(tmp2)
	if(iyr.eq.dimsizes(years)-1)then
		alls	= tmp1
	end if
end if
end do ; do iyr=
;print(alls)
; ------ get monthly data from specific start and end years -------- ends

;print(fdir+case+fappend)
;all_files = systemfunc ("ls "+fdir+case+fappend+"*.nc | grep '198[01234]-..\.nc'")
;all_files = systemfunc ("ls "+fdir+case+fappend+"* | grep '000[234]-..\.nc'")
;print(all_files)
;exit
all_files		= alls
delete(alls)

     fall      = addfiles (all_files, "r")   ; note the "s" of addfile
  ;========================
  ; choose how files are combined and read in variable across files
  ;========================  
     ListSetType (fall, "cat")             ; concatenate or "merge" (default)

     xx_tmp    = fall[:]->$var$                    ; note syntax [:]

	 xx_tmp@_FillValue 	= default_fillvalue("float")
	 print(xx_tmp@_FillValue)

if(landmask)then
	  ; 2018.09.07
	  ;=======================
	  ; mask land grid
	  ;=======================
	  ocnfrac_tmp	= fall[:]->OCNFRAC
	  if(dimsizes(dimsizes(xx_tmp)).eq.4)then
	  ocnfrac	= conform(xx_tmp,ocnfrac_tmp,(/0,2,3/))
	  copy_VarCoords(xx_tmp,ocnfrac)
	  else
	  ocnfrac	= ocnfrac_tmp
	  end if
	;  printVarSummary(ocnfrac)
	
	  xx		= where(ocnfrac.gt.0.5,xx_tmp,xx_tmp@_FillValue)
	  copy_VarCoords(xx_tmp,xx)
	  printVarSummary(xx)
	
	  delete([/ocnfrac_tmp,ocnfrac,xx_tmp/])
else
	if(ocnmask)then
	  ; 2018.09.09
	  ;=======================
	  ; mask ocn grid
	  ;=======================
	  landfrac_tmp	= fall[:]->LANDFRAC
	  if(dimsizes(dimsizes(xx_tmp)).eq.4)then
	  landfrac	= conform(xx_tmp,landfrac_tmp,(/0,2,3/))
	  copy_VarCoords(xx_tmp,landfrac)
	  else
	  landfrac	= landfrac_tmp
	  end if
	
	  xx		= where(landfrac.gt.0.5,xx_tmp,xx_tmp@_FillValue)
	  copy_VarCoords(xx_tmp,xx)
	  printVarSummary(xx)
	
	  delete([/landfrac_tmp,landfrac,xx_tmp/])
	else
	  xx		= xx_tmp
	  delete(xx_tmp)
	end if
end if

  return(xx)
  delete(all_files)

end



function read_omega(fdir:string,case:string)

begin 
  ;========================
  ; get list of all files and open as "one big file"
  ;========================                             
;     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23]-..\.nc'")
;     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '....-..\.nc'")
     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23456789]-..\.nc'")

     fall      = addfiles (all_files, "r")   ; note the "s" of addfile
  ;========================
  ; choose how files are combined and read in variable across files
  ;========================  
     ListSetType (fall, "cat")             ; concatenate or "merge" (default)

     xx    = fall[:]->OMEGA                     ; note syntax [:]
	 
	 return(xx)
end


function read_SWCF(fdir:string,case:string)

begin 
  ;========================
  ; get list of all files and open as "one big file"
  ;========================                             
;     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23]-..\.nc'")
;     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '....-..\.nc'")
     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23456789]-..\.nc'")


     fall      = addfiles (all_files, "r")   ; note the "s" of addfile
  ;========================
  ; choose how files are combined and read in variable across files
  ;========================  
     ListSetType (fall, "cat")             ; concatenate or "merge" (default)

     xx    = fall[:]->SWCF                     ; note syntax [:]
	 
	 return(xx)

end

function read_LWCF(fdir:string,case:string)

begin 
  ;========================
  ; get list of all files and open as "one big file"
  ;========================                             
;     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23]-..\.nc'")
;all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '....-..\.nc'")
     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23456789]-..\.nc'")

     fall      = addfiles (all_files, "r")   ; note the "s" of addfile
  ;========================
  ; choose how files are combined and read in variable across files
  ;========================  
     ListSetType (fall, "cat")             ; concatenate or "merge" (default)

     xx    = fall[:]->LWCF                     ; note syntax [:]
	 
	 return(xx)

end

function read_CLDLOW(fdir:string,case:string)

begin 
  ;========================
  ; get list of all files and open as "one big file"
  ;========================                             
;     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23]-..\.nc'")
;all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '....-..\.nc'")
     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23456789]-..\.nc'")

     fall      = addfiles (all_files, "r")   ; note the "s" of addfile
  ;========================
  ; choose how files are combined and read in variable across files
  ;========================  
     ListSetType (fall, "cat")             ; concatenate or "merge" (default)

     xx    = fall[:]->CLDLOW                     ; note syntax [:]
	 
	 return(xx)

end

function read_CLDMED(fdir:string,case:string)

begin 
  ;========================
  ; get list of all files and open as "one big file"
  ;========================                             
;     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23]-..\.nc'")
;all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '....-..\.nc'")
     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23456789]-..\.nc'")

     fall      = addfiles (all_files, "r")   ; note the "s" of addfile
  ;========================
  ; choose how files are combined and read in variable across files
  ;========================  
     ListSetType (fall, "cat")             ; concatenate or "merge" (default)

     xx    = fall[:]->CLDMED                     ; note syntax [:]
	 
	 return(xx)

end


function read_CLDHGH(fdir:string,case:string)

begin 
  ;========================
  ; get list of all files and open as "one big file"
  ;========================                             
;     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23]-..\.nc'")
;all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '....-..\.nc'")
     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23456789]-..\.nc'")

     fall      = addfiles (all_files, "r")   ; note the "s" of addfile
  ;========================
  ; choose how files are combined and read in variable across files
  ;========================  
     ListSetType (fall, "cat")             ; concatenate or "merge" (default)

     xx    = fall[:]->CLDHGH                     ; note syntax [:]
	 
	 return(xx)

end

function read_RELHUM(fdir:string,case:string)

begin 
  ;========================
  ; get list of all files and open as "one big file"
  ;========================                             
;     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23]-..\.nc'")
;all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '....-..\.nc'")
     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23456789]-..\.nc'")

     fall      = addfiles (all_files, "r")   ; note the "s" of addfile
  ;========================
  ; choose how files are combined and read in variable across files
  ;========================  
     ListSetType (fall, "cat")             ; concatenate or "merge" (default)

     xx    = fall[:]->RELHUM                     ; note syntax [:]
	 
	 return(xx)

end

function read_T(fdir:string,case:string)

begin 
  ;========================
  ; get list of all files and open as "one big file"
  ;========================                             
;     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23]-..\.nc'")
;all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '....-..\.nc'")
     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23456789]-..\.nc'")

     fall      = addfiles (all_files, "r")   ; note the "s" of addfile
  ;========================
  ; choose how files are combined and read in variable across files
  ;========================  
     ListSetType (fall, "cat")             ; concatenate or "merge" (default)

     xx    = fall[:]->T                     ; note syntax [:]
	 
	 return(xx)

end

function read_PBLH(fdir:string,case:string)

begin 
  ;========================
  ; get list of all files and open as "one big file"
  ;========================                             
;     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23]-..\.nc'")
;all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '....-..\.nc'")
     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23456789]-..\.nc'")

     fall      = addfiles (all_files, "r")   ; note the "s" of addfile
  ;========================
  ; choose how files are combined and read in variable across files
  ;========================  
     ListSetType (fall, "cat")             ; concatenate or "merge" (default)

     xx    = fall[:]->PBLH                     ; note syntax [:]
	 
	 return(xx)

end



function read_level(fdir:string,case:string)

begin 
  ;========================
  ; get list of all files and open as "one big file"
  ;========================                             
     all_files = systemfunc ("ls "+fdir+case+"/atm/hist/* | grep '000[23]-..\.nc'")
     fall      = addfiles (all_files, "r")   ; note the "s" of addfile
  ;========================
  ; choose how files are combined and read in variable across files
  ;========================  
     ListSetType (fall, "cat")             ; concatenate or "merge" (default)

     xx    = fall[:]->lev                     ; note syntax [:]
	 
	 return(xx)

end



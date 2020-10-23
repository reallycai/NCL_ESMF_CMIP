; This file still has to be loaded manually
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"


;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
;NCL batch interpolation CMIP
;@author:reallycai
;@data:2020/10/23
;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
begin
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
;       Define parameters,You must set this part to run the program
;       including interpolation year, range, method, accuracy
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
;----------------------------------------------------------------------
;       set interpolation year         
;----------------------------------------------------------------------	
	yrStrt     = 1950
	yrLast     = 2014  
	month = (yrLast -yrStrt+1)*12      
;----------------------------------------------------------------------
;       set interpolation scope
;----------------------------------------------------------------------	 
	latS  = 5    
	latN  = 42    
	lonW  = 102  
	lonE  = 130
;----------------------------------------------------------------------
;       set interpolation method 
;       and accuracy
;----------------------------------------------------------------------
	InterpMethods = (/"bilinear","patch","conserve","neareststod"/)
	Method = InterpMethods(0) 
	GridType = 0.5
;----------------------------------------------------------------------
;       CMIP name of variable
;	  including lat,lon,varname
;----------------------------------------------------------------------
	coordlat=(/"latitude","lat","nav_lat"/)        
	coordlon=(/"longitude","lon","nav_lon"/)                
	var_name=(/"sst","tos"/)
;----------------------------------------------------------------------
;       Sets the file read/write location and name
;----------------------------------------------------------------------
	path = "/home/jack/Downloads/cmipdata/sstdata/"	
	pathout = "/home/jack/Downloads/cmipdata/regrid/"             
	outputname = "all.nc"
	
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
;           pretreatment for batch
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	print(" -----------------------------------------------------------------")
	print("|	This program completes CMIP spatiotemporal extraction 	 |")
	print("|	          and grid interpolation in batches        	 |")
	print("| 	        Interpolation years:"+yrStrt+"-"+yrLast+"        	         |")	
	print("|               Interpolation scope:"+latS+"-"+latN+"N,"+lonW+"-"+lonE+"E               |")
	print("| 	        Interpolation method:"+Method+"        	         |")	
	print(" -----------------------------------------------------------------")

	starttime = toint(systemfunc("date +%s"))	

	setfileoption("nc","Format","NetCDF4")
	ncdf = addfile(pathout+outputname ,"c")
	p = systemfunc("ls "+path)
	pn = dimsizes(p)

	modname = systemfunc("ls "+ path +"|awk -F_ '{print $3}'| uniq ")
	modnum = dimsizes(modname)
	

	modcount = new(modnum+1,integer)
	modcount(0) = 0
	modcount(modnum) = 0
	
	do j= 1,modnum-1
		modcount(j) = toint(systemfunc("ls "+ path +"|grep _"+modname(j-1)+"_ |wc -l"))
	end do
	modcount = cumsum(modcount,0)
	j = 0
	s = 0
	
	srcGridName = "src_"+modname+"_SCRIP.nc"
	dstGridName = "dst_"+modname+"_SCRIP.nc"
	wgtFileName = "wgt_"+modname+".nc"

	;print(p)
	;print(modname)
	
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
;	Create a matrix to hold the interpolation results
;     including (time,models,lat,lon)
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	nlat = toint((latN-latS)/GridType+1)
	nlon = toint((lonE-lonW)/GridType+1)	
	newp = new((/month,modnum,nlat,nlon/),float,1e20)
	newp!0="time"
	newp!1="model"	
	newp!2="lat"
	newp!3="lon"
	
	newp&model  = ispan(1,modnum,1)
	newp&time = ispan(1,month,1)
	newp&lat  = fspan(latS,latN,nlat)
	newp&lon = fspan(lonW,lonE,nlon)
	newp&time@units = "months since 1950-1"
	newp&lat@units = "degree_nort"
	newp&lon@units = "degree_east"

	newp@remap= "ESMF"
	newp@long_name = "Sea Surface Temperature" 
	newp@units = "degC"
	newp@coordinates = "latitude longitude"

;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
;	main programm
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	do i=0,pn-1
		system("echo 3 > /proc/sys/vm/drop_caches")		
		;print("*************"+"No."+(i+1)+" looping"+"*************")
;----------------------------------------------------------------------
;       set interpolation scope
;----------------------------------------------------------------------
		if(i.eq.modcount(j))
			print("************No."+(j+1)+" model:"+modname(j)+" starts Interpolation************")
			dims = 1 
			f = addfile(path+p(i),"r")
;----------------------------------------------------------------------
;       set interpolation scope
;----------------------------------------------------------------------	
			do k=0,dimsizes(coordlat)-1
				if (isfilevar(f, coordlat(k))) then
					lat2d = f->$coordlat(k)$	
					lon2d = f->$coordlon(k)$			
					if(dimsizes(getvardims(lat2d)).eq.2)then			
						dims = 2
						break
					end if
				end if
  			end do
			if  (isfilevar(f, var_name(0)))then
				valu = var_name(0)
			else
				valu = var_name(1)
			end if
;----------------------------------------------------------------------
;   Write WRF grid description to SCRIP file
;----------------------------------------------------------------------
			temp = f->$valu$(0,:,:)

			Opt                = True
			Opt@ForceOverwrite = True
			Opt@Title          = "CMIP6 grid"
			Opt@Mask2D  = where(.not.ismissing(temp),1,0)
			;Opt@PrintTimings   = True

			if (dims.eq.2)
				curvilinear_to_SCRIP(srcGridName(j),lat2d,lon2d,Opt)
			else 
				rectilinear_to_SCRIP(srcGridName(j),lat2d,lon2d,Opt)
			end if
			delete(Opt)
;----------------------------------------------------------------------
; Convert  world grid to SCRIP file
;----------------------------------------------------------------------
			Opt                = True
			Opt@ForceOverwrite = True
			Opt@LLCorner       = (/latS,lonW /)
			Opt@URCorner       = (/ latN, lonE /) 
			Opt@Title          = "World Grid 1-degree Resolution"
			;Opt@PrintTimings   = True
			;Opt@Debug          = True
	
			latlon_to_SCRIP(dstGridName(j),"0.5deg",Opt)
			delete(Opt)
;----------------------------------------------------------------------
; generate a weight file for echo model
;----------------------------------------------------------------------	
			Opt                      = True
			Opt@ForceOverwrite       = True
			Opt@SrcESMF              = False
			Opt@DstESMF              = False
			Opt@InterpMethod         = Method(0)
			;Opt@PrintTimings         = True

			ESMF_regrid_gen_weights(srcGridName(j), dstGridName(j), wgtFileName(j), Opt)
			print(modname(j)+"'s interpolation weight file has been generated")
			delete(Opt)
			delete(lat2d)
			delete(lon2d)
			delete(temp)
			delete(f)

			Opt                = True
			;Opt@PrintTimings   = True
			;Opt@Debug          = True
			Opt@DstGridType = "rectilinear"

			j = j+1
		end if
;----------------------------------------------------------------------
; select the interpolation years from model
;----------------------------------------------------------------------
		f = addfile(path+p(i),"r")
		TIME       =  f->time
		mon       =  cd_calendar(TIME,-1)/100
		idmon     =  ind(mon.ge.yrStrt .and. mon.le.yrLast)
		moncount	    =   dimsizes(idmon)
		delete([/TIME,mon/])

		print("Load: "+p(i))
		print("****************this doc contains "+ moncount +" months****************")
		;print(mon)
		;print(idmon)
		;print(TIME)
;----------------------------------------------------------------------
; begin to interpolation 
;----------------------------------------------------------------------		
		if(moncount.lt.month)then
			if(dims.eq.2)then	
				print("****************(monthparting)Curvilinear: "+modname(j-1)+"****************")				
				do t=0,moncount-1
					sst = f->$valu$(idmon(t),:,:)
					tt =  ESMF_regrid_with_weights(sst,wgtFileName(j-1),Opt) 
					newp(t+s,j-1,:,:) =  conform_dims((/1,1,75,57/),tt,(/2,3/))	
					
				end do
				delete([/idmon,f,t/])	
				;print("！！！！第 "+(j+1)+" 个模式"+(t+s+1)+"个月处理完成！！！！")
			else		
				print("****************(monthparting)Rectilinear: "+modname(j-1)+"****************")				
				do t=0,moncount-1
					sst = f->$valu$(idmon(t),:,:)
					tt =  ESMF_regrid_with_weights(sst,wgtFileName(j-1),Opt) 
					newp(t+s,j-1,:,:) =  conform_dims((/1,1,75,57/),tt,(/2,3/))
				end do
				delete([/idmon,f,t/])
			end if	
			s = s+moncount
		
			if (s.eq.month)then
				print("！！！！！！！！NO."+j+" model has been finished！！！！！！！！")
				print("！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！")
				print("")	
				delete(sst)	
				delete(tt)
				s = 0
			end if
	
		else
			if(dims.eq.2)then
				print("****************(allmonth)Curvilinear: "+modname(j-1)+"****************")			
				do t=0,month-1
					sst = f->$valu$(idmon(t),:,:)
					tt =  ESMF_regrid_with_weights(sst,wgtFileName(j-1),Opt) 
					newp(t,j-1,:,:) =  conform_dims((/1,1,75,57/),tt,(/2,3/))	      
				end do
				print("！！！！！！！！NO."+j+" model has been finished！！！！！！！！")
				print("！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！")
				print("")
				delete([/idmon,t,f/])
				delete(sst)	
				delete(tt)
			else
				print("****************(allmonth)Rectilinear: "+modname(j-1)+"****************")				
				do t=0,month-1
					sst = f->$valu$(idmon(t),:,:)
					tt =  ESMF_regrid_with_weights(sst,wgtFileName(j-1),Opt) 
					newp(t,j-1,:,:) =  conform_dims((/1,1,75,57/),tt,(/2,3/))	      
				end do
				print("！！！！！！！！NO."+j+" model has been finished！！！！！！！！")
				print("！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！")
				print("")
				delete([/idmon,t,f/])
				delete(sst)
				delete(tt)	
			end if 
		end if
	end do 

;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
;	Output the results to NETCDF
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	fileAtt = True
	fileAtt@title = "CMIP interpolation file"
	fileAtt@usemothd = "Earth System Modeling Framework (ESMF) software"
	fileAtt@mydarling = "李淼"
	fileattdef(ncdf,fileAtt)

	dimNames = (/"time","model","lat","lon"/)
	dimSizes = (/-1,modnum,nlat,nlon/)	
	dimUnlim = (/True,False,False,False/)
	filedimdef(ncdf,dimNames,dimSizes,dimUnlim)

	filevardef(ncdf,"time",typeof(newp&time),getvardims(newp&time))
	filevardef(ncdf,"model",typeof(newp&model),getvardims(newp&model))
	filevardef(ncdf,"lat",typeof(newp&lat),getvardims(newp&lat))
	filevardef(ncdf,"lon",typeof(newp&lon),getvardims(newp&lon))
	filevardef(ncdf,"sst",typeof(newp),getvardims(newp))

	ncdf->sst = (/newp/)
	ncdf->model = (/newp&model/)
	ncdf->time = (/newp&time/)
	ncdf->lat = (/newp&lat/)
	ncdf->lon = (/newp&lon/)

	endtime = toint(systemfunc("date +%s"))
	print("You Have Finished Interpolation for All Models ")
	print("This programm ran for "+(endtime-starttime)+"s")
end 

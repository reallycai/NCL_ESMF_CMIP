; This file still has to be loaded manually
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"


;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
;                                  NCL batch interpolation CMIP
;@author:reallycai
;@data:2020/10/23
;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
undef ("getFre")
function getFre(path)
begin
	p = systemfunc("ls "+path)
	f = addfile(path+p(1),"r")
	fAtts = getfileatts(f) 
	do na=0,dimsizes(fAtts)-1
    		;print ("f@" + fAtts(na) + " = " + f@$fAtts(na)$)
		if (fAtts(na).eq."frequency")then
			Fre = f@$fAtts(na)$
			if (Fre.eq."mon")then
				summ = 12
			else if (Fre.eq."day")then
				summ = 365
			else if (Fre.eq."6hr")then
				summ = 1460
			else 
				summ = 2920
			end if
			end if
			end if
		end if
  	end do
	return summ
end

undef ("start")
procedure start(j,modVar,modFre,modname,modExp,modVL,dims)
begin
	print("*********************************************************")
	print("       No."+(j+1)+" model:"+modname(j)+" starts Interpolation      ")
	print("*******Variable: "+modVar(j)+"                         	  *******")
	print("*******Frequency: "+modFre(j)+"                         	  *******")
	print("*******Experiment_ID: "+modExp(j)+"                  *******")
	print("*******Variant_Label: "+modVL(j)+"                    *******")
	if(dims.eq.1)then
		print("*******Coordinate : Rectilinear                   *******")
	else
		print("*******Coordinate : Curvilinear                   *******")
	end if
	print("*********************************************************")
end

undef ("finish")
procedure finish(j)
begin
	print("！！！！！！！！NO."+j+" model has been finished！！！！！！！！")
	print("！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！")
	print("")
end



undef ("nclcmip")
procedure nclcmip(path,pathout,outputname,yrStrt,yrLast,latS,latN,lonW,lonE,Method,GridType)
begin
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
;       Some parameters that need to be defined in advance
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 	Fre = getFre(path)
	month = (yrLast -yrStrt+1)*Fre     
;----------------------------------------------------------------------
;       CMIP name of variable
;	  including lat,lon,varname
;----------------------------------------------------------------------
	coordlat=(/"latitude","lat","nav_lat"/)        
	coordlon=(/"longitude","lon","nav_lon"/)                

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
	
	cmd1 = "'ls' "+ path +" |awk -F_ '{print $1"
	cmd2 = integertochar(34)+"_"+integertochar(34)
	cmd3 = "$5}'| uniq"
	cmd = cmd1+cmd2+"$2"+cmd2+"$3"+cmd2+"$4"+cmd2+cmd3
	modVar =  systemfunc(cmd+"|awk -F_ '{print $1}'")
	modFre = systemfunc(cmd+"|awk -F_ '{print $2}'")
	modname = systemfunc(cmd+"|awk -F_ '{print $3}'")
	modExp = systemfunc(cmd+"|awk -F_ '{print $4}'")
	modVL = systemfunc(cmd+"|awk -F_ '{print $5}'")	
	modfirstlc = modVar+"_"+modFre+"_"+modname+"_"+modExp+"_"+modVL
	modnum = dimsizes(modname)
	
	print("model: "+modfirstlc)

	modcount = new(modnum+1,integer)
	modcount(0) = 0
	modcount(modnum) = 0
	

	do j= 1,modnum-1
		modcount(j) = toint(systemfunc("ls "+ path +"|grep "+modfirstlc(j-1)+" |wc -l"))
	end do

	modcount = cumsum(modcount,0)
	;print(modcount)
	j = 0
	s = 0
	
	srcGridName = "src_"+modfirstlc+"_SCRIP.nc"
	dstGridName = "dst_"+modfirstlc+"_SCRIP.nc"
	wgtFileName = "wgt_"+modfirstlc+".nc"

	;print(p)
	;print(modname)
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
;	Create a matrix to hold the interpolation results
;     including (time,models,lat,lon)
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	nlat = toint((latN-latS)/GridType+1)
	nlon = toint((lonE-lonW)/GridType+1)	
	newp = new((/month,modnum,nlat,nlon/),double,1e20)
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

			valu = modVar(j)
			start(j,modVar,modFre,modname,modExp,modVL,dims)
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
		if(moncount.eq.1)then
			print("****************this doc contains "+ (moncount-1) +" months****************")
			print("！！！！！！！No Interpolation Period In the File！！！！！！！")	
			print("")		
			delete([/idmon,f/])
			continue
		end if
		print("****************this doc has "+ moncount +" needed months****************")
		print("")
		;print(mon)
		;print(idmon)
		;print(TIME)
;----------------------------------------------------------------------
; begin to interpolation 
;----------------------------------------------------------------------		
		if(moncount.lt.month)then
			if(dims.eq.2)then	
				do t=0,moncount-1
					value = f->$valu$(idmon(t),:,:)
					tt =  ESMF_regrid_with_weights(value,wgtFileName(j-1),Opt) 
					newp(t+s,j-1,:,:) =  conform_dims((/1,1,75,57/),tt,(/2,3/))	
				end do
				delete([/idmon,f,t/])	
			else		
				do t=0,moncount-1
					value = f->$valu$(idmon(t),:,:)
					tt =  ESMF_regrid_with_weights(value,wgtFileName(j-1),Opt) 
					newp(t+s,j-1,:,:) =  conform_dims((/1,1,75,57/),tt,(/2,3/))
				end do
				delete([/idmon,f,t/])
			end if	
			
			s = s+moncount
		
			if (s.eq.month)then
				finish(j)
				delete(value)	
				delete(tt)
				s = 0
			end if
	
		
		else
			if(dims.eq.2)then
				do t=0,month-1
					value = f->$valu$(idmon(t),:,:)
					tt =  ESMF_regrid_with_weights(value,wgtFileName(j-1),Opt) 
					newp(t,j-1,:,:) =  conform_dims((/1,1,75,57/),tt,(/2,3/))	      
				end do
				finish(j)
				delete(value)	
				delete(tt)
				delete([/idmon,t,f/])
			else
				do t=0,month-1
					value = f->$valu$(idmon(t),:,:)
					tt =  ESMF_regrid_with_weights(value,wgtFileName(j-1),Opt) 
					newp(t,j-1,:,:) =  conform_dims((/1,1,75,57/),tt,(/2,3/))	      
				end do
				finish(j)
				delete(value)	
				delete(tt)
				delete([/idmon,t,f/])
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
	filevardef(ncdf,"value",typeof(newp),getvardims(newp))

	ncdf->value = (/newp/)
	ncdf->model = (/newp&model/)
	ncdf->time = (/newp&time/)
	ncdf->lat = (/newp&lat/)
	ncdf->lon = (/newp&lon/)

	endtime = toint(systemfunc("date +%s"))
	print("You Have Finished Interpolation for All Models ")
	print("This programm ran for "+(endtime-starttime)+"s")
end 

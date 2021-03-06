; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; This file still has to be loaded manually
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
;                                                                   NCL batch interpolation CMIP
;@author:reallycai
;@data:2020/11/04
;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
undef ("getFre")
function getFre(f,modVar)
local tempvar,f,fAtts,na,Fre,summ
begin
	fAtts = getfileatts(f) 
	tempvar = new(2,integer)	
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
	tempvar(0) = summ
	tempvar(1) = dimsizes(getfilevardims(f,modVar(0)))
	return tempvar	
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

undef ("costtime")
procedure costtime(starttime)
local var
begin
	endtime = toint(systemfunc("date +%s"))
	print("You Have Finished Interpolation for All Models ")
	print("This programm ran for "+(endtime-starttime)+"s")
end

undef ("nclcmip")
procedure nclcmip(path,pathout,outputname,yrStrt,yrLast,latS,latN,lonW,lonE,Method,GridType)
begin
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
;       Some parameters that need to be defined in advance
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
;----------------------------------------------------------------------
;       CMIP name of variable
;	  including lat,lon,varname
;----------------------------------------------------------------------
	coordlat=(/"latitude","lat","nav_lat"/)        
	coordlon=(/"longitude","lon","nav_lon"/)                
	nlat = toint((latN-latS)/GridType+1)
	nlon = toint((lonE-lonW)/GridType+1)
	setfileoption("nc","Format","NetCDF4")
	ncdf = addfile(pathout+outputname ,"c")
	starttime = toint(systemfunc("date +%s"))	

	p = systemfunc("ls "+path)
	pn = dimsizes(p)
	f = addfile(path+p(0),"r")
	
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
	;print("model: "+modfirstlc)
	print("****************There are "+modnum+" models in the path****************")
	modcount = new(modnum+1,integer)
	modcount(0) = 0
	modcount(modnum) = 0
	do j= 1,modnum-1
		modcount(j) = toint(systemfunc("ls "+ path +"|grep "+modfirstlc(j-1)+" |wc -l"))
	end do
	modcount = cumsum(modcount,0)
	;print(modcount)

 	tempvar = getFre(f,modVar)
	;print(tempvar)
	Duration = (yrLast -yrStrt+1)*tempvar(0)
	Dimension = tempvar(1)   
	
	if (Dimension.eq.3)
		print("****************This is a 3D interpolation task****************")			
		srcGridName = "src_"+modfirstlc+"_SCRIP.nc"
		dstGridName = "dst_"+modfirstlc+"_SCRIP.nc"
		wgtFileName = "wgt_"+modfirstlc+".nc"

		;print(p)
		;print(modname)
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
;	Create a matrix to hold the interpolation results
;     including (time,models,lat,lon)
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&	
		newp = new((/Duration,modnum,nlat,nlon/),float,1e20)
		var = f->$modVar(0)$
		newp!0="time"
		newp!1="model"	
		newp!2="lat"
		newp!3="lon"
	
		newp&model  = ispan(1,modnum,1)	
		newp&time = ispan(1,Duration,1)
		newp&lat  = fspan(latS,latN,nlat)
		newp&lon = fspan(lonW,lonE,nlon)
		newp&time@units = "months since 1950-1"
		newp&lat@units = "degree_nort"
		newp&lon@units = "degree_east"

		newp@remap= "ESMF"
		newp@long_name = var@long_name 
		newp@units = var@units
		newp@coordinates = "latitude longitude"
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
;	main programm
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		j = 0
		s = 0
		do i=0,pn-1
			system("echo 3 > /proc/sys/vm/drop_caches")		
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
				Opt@URCorner       = (/latN, lonE/) 
				Opt@Title          = "World Grid 1-degree Resolution"
				;Opt@PrintTimings   = True
				;Opt@Debug          = True
		
				latlon_to_SCRIP(dstGridName(j),tostring(GridType)+"deg",Opt)
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
			data       =  cd_calendar(TIME,-1)/100
			iddata      =  ind(data .ge.yrStrt .and. data .le.yrLast)
			datacount	    =   dimsizes(iddata )
			delete([/TIME,data/])
		
			print("Load: "+p(i))
			if(datacount.eq.1)then
				print("****************this doc contains "+ (datacount-1) +" periods****************")
				print("！！！！！！！No Interpolation Period In the File！！！！！！！")	
				print("")		
				delete([/iddata,f/])
				continue
			end if
			print("****************this doc has "+ datacount +" needed periods****************")
			print("")
			;print(data)
			;print(iddata)
			;print(TIME)
;----------------------------------------------------------------------
; begin to interpolation 
;----------------------------------------------------------------------		
			if(datacount.lt.Duration)then
				if(dims.eq.2)then	
					do t=0,datacount-1
						value = f->$valu$(iddata(t),:,:)
						tt =  ESMF_regrid_with_weights(value,wgtFileName(j-1),Opt) 
						newp(t+s,j-1,:,:) =  conform_dims((/1,1,nlat,nlon/),tt,(/2,3/))	
					end do
					delete([/iddata,f/])	
				else		
					do t=0,datacount-1
						value = f->$valu$(iddata(t),:,:)
						tt =  ESMF_regrid_with_weights(value,wgtFileName(j-1),Opt) 
						newp(t+s,j-1,:,:) =  conform_dims((/1,1,nlat,nlon/),tt,(/2,3/))
					end do
					delete([/iddata,f/])
				end if	
				s = s+datacount	
				if (s.eq.Duration)then
					finish(j)
					delete(value)	
					delete(tt)
					s = 0
				end if		
			else
				if(dims.eq.2)then
					do t=0,Duration-1
						value = f->$valu$(iddata(t),:,:)
						tt =  ESMF_regrid_with_weights(value,wgtFileName(j-1),Opt) 
						newp(t,j-1,:,:) =  conform_dims((/1,1,nlat,nlon/),tt,(/2,3/))	      
					end do
					finish(j)
					delete(value)	
					delete(tt)
					delete([/iddata,f/])
				else
					do t=0,Duration-1
						value = f->$valu$(iddata(t),:,:)
						tt =  ESMF_regrid_with_weights(value,wgtFileName(j-1),Opt) 
						newp(t,j-1,:,:) =  conform_dims((/1,1,nlat,nlon/),tt,(/2,3/))	      
					end do
					finish(j)
					delete(value)	
					delete(tt)
					delete([/iddata,f/])
				end if 
			end if		
 		end do 
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
;	Output the results to NETCDF
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		fileAtt = True
		fileAtt@title = "CMIP interpolation file"
		fileAtt@usemothd = "Earth System Modeling Framework (ESMF) software"
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
		costtime(starttime)
	else
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
;	4D
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&		
		print("****************This is a 4D interpolation task****************")		
		pleve = f->plev
		plevenum = dimsizes(pleve)
		print("****************The following layers to choose****************")
		print((/pleve/))	
		lev = new(plevenum+1,string)
		system("echo please choose a range of plev,enter go to begin")
		do pnum = 0,plevenum
			lev(pnum) = systemfunc("read lev; echo $lev")
			if (lev(pnum ).eq."go")
				break
			end if
		end do 
		
		laycount = new(pnum ,integer)
		do j=0,pnum -1
			laycount(j) = toint(lev(j))
		end do	
		
		srcGridName = new((/pnum ,modnum/),string)
		dstGridName = new((/pnum ,modnum/),string)
		wgtFileName  = new((/pnum ,modnum/),string)
		do j=0,pnum -1
			srcGridName(j,:) = "src_"+modfirstlc+"_"+laycount(j)+"_SCRIP.nc"
			dstGridName(j,:) = "dst_"+modfirstlc+"_"+laycount(j)+"_SCRIP.nc"
			wgtFileName(j,:) = "wgt_"+modfirstlc+"_"+laycount(j)+".nc"
		end do

		newp = new((/modnum,Duration,pnum ,nlat,nlon/),float,1e20)
		var = f->$modVar(0)$		
		newp!0="model"
		newp!1="time"
		newp!2="plev"
		newp!3="lat"
		newp!4="lon"
	
		newp&model  = ispan(1,modnum,1)
		newp&time = ispan(1,Duration,1)
		newp&plev = (/laycount/)
		newp&lat  = fspan(latS,latN,nlat)
		newp&lon = fspan(lonW,lonE,nlon)
		newp&time@units = "months since 1950-1"
		newp&plev@units = "Pa"
		newp&lat@units = "degree_nort"
		newp&lon@units = "degree_east"

		newp@remap= "ESMF"
		newp@long_name = var@long_name 
		newp@units = var@units
		newp@coordinates = "latitude longitude"
		
		j = 0
		s = 0
		do i=0,pn-1
			system("echo 3 > /proc/sys/vm/drop_caches")		
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
				do k = 0,pnum-1
					temp = f->$valu$(0,ind(laycount(k).eq.pleve),:,:)
					Opt                = True
					Opt@ForceOverwrite = True
					Opt@Title          = "CMIP6 grid"
					Opt@Mask2D  = where(.not.ismissing(temp),1,0)
					;Opt@PrintTimings   = True
		
					if (dims.eq.2)
						curvilinear_to_SCRIP(srcGridName(k,j),lat2d,lon2d,Opt)
					else 
						rectilinear_to_SCRIP(srcGridName(k,j),lat2d,lon2d,Opt)
					end if
					delete(Opt)
;----------------------------------------------------------------------
; Convert  world grid to SCRIP file
;----------------------------------------------------------------------
					Opt                = True
					Opt@ForceOverwrite = True
					Opt@LLCorner       = (/latS,lonW /)
					Opt@URCorner       = (/latN, lonE/) 
					Opt@Title          = "World Grid 1-degree Resolution"
					;Opt@PrintTimings   = True
					;Opt@Debug          = True
		
					latlon_to_SCRIP(dstGridName(k,j),tostring(GridType)+"deg",Opt)
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
		
					ESMF_regrid_gen_weights(srcGridName(k,j), dstGridName(k,j), wgtFileName(k,j), Opt)
					print(modname(j)+"_"+laycount(k)+"'s interpolation weight file has been generated")
					delete(Opt)
					delete(temp)
				end do
				j = j+1
				Opt                = True
				;Opt@PrintTimings   = True
				;Opt@Debug          = True
				Opt@DstGridType = "rectilinear"
			
				delete(lat2d)
				delete(lon2d)
				delete(f)
			end if

			f = addfile(path+p(i),"r")
			TIME       =  f->time
			data       =  cd_calendar(TIME,-1)/100
			;print(data)
			iddata     =  ind(data.ge.yrStrt .and.data.le.yrLast)
			;print(iddata)
			datacount	    =   dimsizes(iddata)
			delete([/TIME,data/])
			
			print("Load: "+p(i))
			if(datacount.eq.1)then
				print("****************this doc contains "+ (datacount-1) +" peiods****************")
				print("！！！！！！！No Interpolation Period In the File！！！！！！！")	
				print("")		
				delete([/iddata,f/])
				continue
			end if
			print("****************this doc has "+ datacount +" needed periods****************")
			print("")
	
			if(datacount.lt.Duration)then
				if(dims.eq.2)then	
					do k=0,pnum-1				
						do t=0,datacount-1
							value = f->$valu$(iddata(t),ind(laycount(k).eq.pleve),:,:)
							tt =  ESMF_regrid_with_weights(value,wgtFileName(k,j-1),Opt) 
							newp(j-1,t+s,k,:,:) =  conform_dims((/1,1,1,nlat,nlon/),tt,(/3,4/))	
						end do
						print("！！！！！！！！"+laycount(k)+"th layer has finished！！！！！！！！")
						delete(tt)
						delete(value)
					end do
					print("")
					delete([/iddata,f/])	
				else		
					do k=0,pnum-1					
						do t=0,datacount-1
							value = f->$valu$(iddata(t),ind(laycount(k).eq.pleve),:,:)
							tt =  ESMF_regrid_with_weights(value,wgtFileName(k,j-1),Opt) 
							newp(j-1,t+s,k,:,:)=  conform_dims((/1,1,1,nlat,nlon/),tt,(/3,4/))
						end do
						print("！！！！！！！！"+laycount(k)+"th layer has finished！！！！！！！！")
						delete(tt)	
						delete(value)					
					end do
					print("")
					delete([/iddata,f/])
				end if	
				s = s+datacount
				if (s.eq.Duration)then
					finish(j)
					s = 0
				end if
			else
				if(dims.eq.2)then
					do k=0,pnum-1
						do t=0,Duration-1
							value = f->$valu$(iddata(t),ind(laycount(k).eq.pleve),:,:)
							tt =  ESMF_regrid_with_weights(value,wgtFileName(k,j-1),Opt) 
							newp(j-1,t,k,:,:) =  conform_dims((/1,1,1,nlat,nlon/),tt,(/3,4/))      
						end do
						print("！！！！！！！！"+laycount(k)+"th layer has finished！！！！！！！！")
						delete(value)	
						delete(tt)
					end do
					finish(j)
					delete([/iddata,f/])
				else
					do k=0,pnum-1				
						do t=0,Duration-1
							value = f->$valu$(iddata(t),ind(laycount(k).eq.pleve),:,:)
							tt =  ESMF_regrid_with_weights(value,wgtFileName(k,j-1),Opt) 
							newp(j-1,t,k,:,:) = conform_dims((/1,1,1,nlat,nlon/),tt,(/3,4/))	      
						end do
						print("！！！！！！！！"+laycount(k)+"th layer has finished！！！！！！！！")
						delete(value)	
						delete(tt)
					end do
					finish(j)
					delete([/iddata,f/])
				end if 
			end if
		end do
			
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
;	Output the results to NETCDF
;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		fileAtt = True
		fileAtt@title = "CMIP interpolation file"
		fileAtt@usemothd = "Earth System Modeling Framework (ESMF) software"
		fileattdef(ncdf,fileAtt)
	
		dimNames = (/"model","time","plev","lat","lon"/)
		dimSizes = (/modnum,-1,pnum,nlat,nlon/)	
		dimUnlim = (/False,True,False,False,False/)
		filedimdef(ncdf,dimNames,dimSizes,dimUnlim)
	
		filevardef(ncdf,"time",typeof(newp&time),getvardims(newp&time))
		filevardef(ncdf,"model",typeof(newp&model),getvardims(newp&model))
		filevardef(ncdf,"plev",typeof(newp&plev),getvardims(newp&plev))
		filevardef(ncdf,"lat",typeof(newp&lat),getvardims(newp&lat))
		filevardef(ncdf,"lon",typeof(newp&lon),getvardims(newp&lon))
		filevardef(ncdf,"value",typeof(newp),getvardims(newp))
	
		ncdf->value = (/newp/)
		ncdf->plev = (/newp&plev/)
		ncdf->model = (/newp&model/)
		ncdf->time = (/newp&time/)
		ncdf->lat = (/newp&lat/)
		ncdf->lon = (/newp&lon/)
		costtime(starttime)
	end if
end 

; to save typing this little program uses often used settings
; warning: it does NOW know something about flag/bad data etc.
;(SH Feb 11 2000) Added own color routine and now has:
; aband = [1,2,3,4] plot bands 1-4
; adet  = [1,4,8] plot dets 1,4,8
; ascan = -1 only plot scan -1
; all=1 plot also the points flagged as no data otherwise do not plot
; these.
; Added option for frequency xaxis
;(SH Feb 21 2000)
; Added option for autoscaling
;(SH Feb 28 2000)
;Changed pl_select, to always return some aar. This prevents error
;when one does for example pl,aar,det=[1213]. Also a bug with a single
;second array(B) in pl_comparray was fixed. It would return 0 unless the
;first element was (coincidently) equal to the first element in the
;first array (A).
;(SH Mar  2 2000) Fixed mark with ylog plot.
;(SH Mar  3 2000) Start of error and plambda added, make bands, scans
;and dets outside of recur loop as to make pl,a,/s,/d make sense
;(SH Mar  6 2000) Fixed recursive call trashing axes titles.
;Added mstyle and mrange keywords.
;Improved the yaxis labels now use quantity and unit
;(SH Mar  7 2000) xontop now working
;(SH Mar  9 2000) added string values accepted to aband and ascan
; xaxis and yaxis keywords for plotting different things with still
; holding knowledges of band dets etc.
; added ramps keyword.
; updated colors from ramps,dets,bands,scans loop to still increase
; when oveplotting
;(SH Mar 13 2000) added ndet,nband,nscan keywords for NOT plotting
;specfic lines.
;(SH Mar 14 2000) Added +name option in title to be substituted by the
;sources name
; Fixed special truetype fonts being scrambled.
; (SH Apr  4 2000) added return option to simply get plotted data back
; as an aar
; (SH Dec 11 2000)
; Added options to print user,observer and proposal in title.
; Changed +name to +object for consistency
; (SH Dec 19 2000)
; Added pl_headkey to use instead of read_fits_key
;(SH May 30 2001)
; Added pl_define_aar to make it less dependent on IA and more quick

; compares to arrays (A,B) and gives back whether the elements
; in A matches any of those in B

FUNCTION pl_comparray,A,B
;Check for simple input
  IF n_elements(a) EQ 1 THEN return,total(a[0] EQ b)
  IF n_elements(b) EQ 1 THEN return,a EQ b[0]
; stretch the values in both directions
  AA = A#(b EQ b)
  BB = (a EQ a)#B
;compare the two matrices we just created and collapse the result
  return,total((AA EQ BB),2)
END

;Give back unique values
FUNCTION pl_uniq,a
  return,a[uniq(a,sort(a))]
END

;stupid thing for xlog plots
FUNCTION pl_xrange
  return,(!x.type EQ 0)*!x.crange+(!x.type EQ 1)*10d0^!x.crange
END

;stupid thing for ylog plots
FUNCTION pl_yrange
  return,(!y.type EQ 0)*!y.crange+(!y.type EQ 1)*10d0^!y.crange
END

;; Simple function to get the value of a fitsheader in a string
FUNCTION pl_headkey, h, keyword

  k = strupcase(keyword)
;; pad the keyword with some spaces
  k = strmid( k + '        ', 0, 8 )

; find the position of the keyword
  loc = 0
  loc = strpos(h,k,loc)
  IF loc EQ -1 THEN BEGIN
      return,''
  ENDIF

; Get the important part of this line (from the = to the /)
  result = strmid(h,loc+9,22)

  loc1 = strpos(result,"'")
  loc_temp = loc1
;  find the last '
  REPEAT BEGIN
      loc2 = loc_temp
      loc_temp = strpos(result,"'",loc2+1)
  ENDREP UNTIL loc_temp EQ -1

  IF (loc1 NE -1) AND (loc2 NE -1) THEN BEGIN
      result = strtrim(strmid(result,loc1+1,loc2-loc1-1),2)
  ENDIF

  return,result
END 

; Extract itk0 from header or else from first itk point
FUNCTION pl_itk0,a
  utcs = ''
;;  sts   = READ_FITS_KEY( a.header, 'EOHAUTCS', utcs, comment )
  utcs   = pl_headkey( a.header, 'EOHAUTCS')
  utcs = strcompress(utcs,/remove_all)
  IF strlen(utcs) EQ 11 THEN $
    itks = ut_to_itk(a.header,utcs) $
  ELSE itks=a.data[0].itk
  return,itks
END

;; Number as string without leading or trailing spaces
FUNCTION pl_n2s,in,format=format
  IF n_elements(format) EQ 0 THEN format=''
  return,strtrim(string(format=format,in),2)
END

; Function to take a float and return a string with that value and no
; leading or trailing spaces valus above 9999.99 or so small that they
; would be 0.00 are output in exp function. Integers are return as integers
FUNCTION pl_f2s,f,p2,deci=deci
;;Check params. with special attention for possible integer values
  CASE n_params() OF 
      0: return,''
      1: BEGIN
          IF n_elements(deci) EQ 0 THEN BEGIN
              ;; Is it integer?
              s = size(f)
              IF (total(s[s[0]+1] EQ [1,2,3,12,13,14,15]) NE 0) AND (abs(f) LT 1d5) THEN $
                return,pl_n2s(f)
              deci = 2
          ENDIF
      END
      2: deci=p2
  ENDCASE
  dec = floor(deci)
  IF dec EQ 0 THEN return,pl_n2s(round(f))
  f = double(f)
;; Now we get the magnitude of f, after rounding
;; therefore we add 0.5 after the last decimal place
  mag_f = floor(alog10(abs(f)+0.5d0*(1d-1)^double(dec)))

;; Now we have 3 possibilities:
;; 1) the number is too large -> make an exponent
;; 2) the number would round to 0 -> make an exponent
;; 3) make a float number
;; 1) and 2 have the same format so catch them in the same case
;; statement:
  IF (mag_f GE 5) OR (mag_f LT -dec) THEN BEGIN
      ;; How many characters needed for exponent (minimum is 2)?
      width_exp = (floor(alog10(abs(mag_f)))+1)>2
      ;; will look like (-)a.bbE-cc
      ;;            sign    a . bb  E - cc
      frm='(E'+pl_n2s((f LT 0)+1+1+dec+1+1+width_exp)+'.'+pl_n2s(dec)+'E'+pl_n2s(width_exp)+')' 
  ENDIF ELSE BEGIN
      ;;            sign     before   . after 
      frm='(F'+pl_n2s((f LT 0)+((mag_f>0)+1)+1+dec)+'.'+pl_n2s(dec)+')' 
  ENDELSE
  return,string(format=frm,f)
END

;Convert string names to sdir numbers allowed are:
;up/down/first/second/-1/+1/1/0
FUNCTION pl_scanname2sdir,names
  scans=['u','d','f','s','-','1','+','0']
  sdirs=[ -1,  1,  1, -1,  -1,  1, 1,  0]
;;Only look at the first letter u=up,d=down,f=first,s=second  
  names = strlowcase(strmid(names,0,1))
; build two matrices to compare
  Nnames  = n_elements(names)
  Nscans  = n_elements(scans)
  l = lindgen(Nnames,Nscans)
  AA = names[l MOD Nnames]
  BB = scans[l / Nnames]
;compare the two matrices we just created
  I = where(AA EQ BB)
  IF i[0] NE -1 THEN return,pl_uniq(sdirs[I / Nnames]) ELSE return,0
END
  
;Convert string names to line numbers allowed are:
;'1a','1b' etc and '1','2','3','4','4o' for the offband data
FUNCTION pl_bandname2line,names
  bands=['1a','1b','1d','1e','2a','2b','2c','3a','3c','3d','3e', '4','4a','4c','4d']
  lines=[   1,   2,   3,   4,   5,   6,   7,   9,  10,  11,  12,  13,  20,  21,  22]
  bands2= [     '1',      '2',         '3',          '4',        '4o']
  lines2=[[1,2,3,4],[5,6,7,7],[9,10,11,12],[13,13,13,13],[20,21,22,22]]

  names = strlowcase(names)
; build two matrices to compare
  Nnames  = n_elements(names)
  Nbands  = n_elements(bands)
  l = lindgen(Nnames,Nbands)
  AA = names[l MOD Nnames]
  BB = bands[l / Nnames]
;compare the two matrices we just created
  I = where(AA EQ BB)
  IF i[0] NE -1 THEN out = lines[I / Nnames] ELSE out = 0

;; This is to allow for names like aband='2'  
  Nbands = n_elements(bands2)
  l = lindgen(Nnames,Nbands)
  AA = names[l MOD Nnames]
  BB = bands2[l / Nnames]
;;compare the two matrices we just created
  I = where(AA EQ BB)
  IF i[0] NE -1 THEN BEGIN
      lines2 = lines2[*,I / Nnames]
      out = [out,reform(lines2,n_elements(lines2))]
  ENDIF
  return,pl_uniq(out)
END
  
FUNCTION pl_select,in,cond
  idx = where(cond)
  IF idx[0] NE -1 THEN BEGIN
      return,{type   : 'SAAR'        ,$ ; structure type
              header : in.header        ,$ ; the header string
              history: in.history       ,$ ; the header string
              data   : in.data[idx] } ; the given data  aar
  ENDIF ELSE BEGIN
      print,'PL_SELECT: Warning select-condition throws out ALL data. NO selection!'
      return,in
  ENDELSE
END

PRO pl_key,arr,n,colors
IF !p.charsize EQ 0d0 THEN $
  pheight=(convert_coord([!d.y_ch_size,1],/dev,/to_n))[0]*1.1 $
ELSE $
  pheight=(convert_coord([!d.y_ch_size,1],/dev,/to_n))[0]*!p.charsize*1.1

;;retrieve the normalized coords of the upperleft corner
xrange = pl_xrange()
yrange = pl_yrange()
x0 = (convert_coord([xrange[0],1],/to_n))[0]
y0 = (convert_coord([1,yrange[1]],/to_n))[1]

FOR i = 0,n_elements(arr)-1 DO BEGIN
xyouts,x0+0.05,y0-0.05-pheight*i, $
       string(format='(A,I5)',n,arr[i]),color=colors[i],/normal 
ENDFOR
END

FUNCTION pl_xtickvalues,range
  cmrange = 1d4/range
  rmin = min(cmrange,max=rmax)

  minstep = 10d0^floor(alog10(abs(rmin)))
  first = minstep*(ceil(rmin/minstep))
  nstep = floor((rmax-first)/minstep)
;; If we do not have a magnitude range than  
  IF nstep LT 4 THEN BEGIN
      step = 10d0^floor(alog10(rmax-rmin))
      first = step*(ceil(rmin/step))
      nstep = floor((rmax-first)/step)
; now make sure we have between 4-7 steps
      CASE 1 OF
          nstep LT 4: minstep = minstep/2d0
          nstep GT 7: minstep = minstep*2
          ELSE:
      ENDCASE
      first = minstep*(ceil(rmin/minstep))
      nstep = floor((rmax-first)/minstep)
      steps = indgen(nstep+1)*minstep+first
      return,1d4/steps
  ENDIF
  steps = indgen(nstep+1)*minstep+first
; Now we want to throw out things like 2100 2200 etc but not 2000 and
; 3000
  mag = 10d0^floor(alog10(steps))
  idx= where((steps MOD mag) EQ 0)
  return,1d4/steps[idx]
END

FUNCTION pl_xtickformat,axis,index,value,level
  COMMON COMMON_PLS, sh_plotcount,sh_plotflambda,sh_plotiras,sh_plotcm, $
    sh_plottime,sh_plotnu,sh_plotplambda

  CASE 1 OF 
      sh_plotnu: fct=3d2
      sh_plotcm: fct=1d4
      ELSE: fct=1d4
  ENDCASE

  return,pl_f2s(fct/value,deci=2*((sh_plotcm) OR (sh_plotnu)))
END

PRO pl_xontop

COMMON COMMON_PLS, sh_plotcount,sh_plotflambda,sh_plotiras,sh_plotcm, $
  sh_plottime,sh_plotnu,sh_plotplambda

IF !p.font EQ -1 THEN mu      = '!7l!X' ELSE mu   = '!Mm!X' 

CASE 1 OF 
    sh_plotcm: BEGIN
        xtit = 'Wavelength ['+mu+'m]'
        divisor = 1d4
    END
    sh_plotnu: BEGIN
        xtit = 'Wavelength ['+mu+'m]'
        divisor = 3d2
    END
    sh_plottime GT 0L: BEGIN
;;This makes no sense just exit
        return
    END
    ELSE: BEGIN
        xtit='Wavenumber [cm!U-1!N]'
        divisor = 1d4
    END
ENDCASE
values = pl_xtickvalues(pl_xrange())
;; The standard title gets too close to the labels so put it in
;; myself
axis,xaxis=1,xstyle=1,xtickformat='pl_xtickformat', $
     xticks=n_elements(values)-1, $
     xtickv=values,ymargin=0
label_pos = convert_coord((!x.crange[1]+!x.crange[0])/2d0,!y.crange[1], $
                          /data,/to_device)
char_height = !d.y_ch_size*((!p.charsize EQ 0)+!p.charsize*(!p.charsize NE 0))
xyouts,label_pos[0],label_pos[1]+char_height*1.75,xtit,align=0.5,/device
END
  
PRO pl_mark,mark,mrange=mrange,_extra=_extra
  IF n_elements(mrange) LT 2 THEN mrange = pl_yrange()
  FOR i = 0, n_elements(mark)-1 DO BEGIN
      oplot,[mark[i],mark[i]],mrange,_extra=_extra
  ENDFOR
END

;; To define a basic aar_structure
FUNCTION pl_define_aar,length=length

  return,{type   : 'SAAR'  ,$   ; structure type
          header : ''      ,$   ; the header string
          history: ''      ,$   ; history string
          data   : replicate({wave:0.0,flux:0.0,stdev:0.0,tint:0L,det:0L,itk:0L, $
                              line:0L,sdir:0L,flag:0L},length) } 
END

;; small routine to check for aar input
FUNCTION pl_is_aar,var
;; fail by default
  result = 0

  IF n_tags(var) GE 3 THEN BEGIN
      names = tag_names(var)
      IF (WHERE(names EQ 'TYPE'))[0] NE -1 THEN BEGIN
          IF strpos(var.type,'AAR') GE 0 THEN BEGIN
              result = 1
          ENDIF
      ENDIF
  ENDIF
  return, result
END


;list of options for pl
;oplot   = 1          Plot ontop of previous
;ll      = 1          Make axis xlog and ylog
;flambda = 1          Convert fluxes to f_lambda
;iras    = 1          Convert fluxes to lambda*f_lambda
;white   = 1          Make plot in default color (use with oplot)
;cm      = 1          X-axis in wavenumbers
;time    = 1          Plot in itk
;cfact   = f          Multiply all fluxes by f
;coffset = o          Add o to all fluxes
;dets    = 1          Plot detectors separately
;scans   = 1          Plot up/down separately
;bands   = 1          Plot aotbands separately
;nooff   = 1          Do no plot offband data  
;key     = 1          Make an legend to the detectors/scan/bands plotted
;xontop  = 1          Make double xaxis for wavenumber and micron
;ymin    = y          Make lower value of the y-axis fixed: y 
;mark    = [x1,x2,..] Put vertical line(s) in plot at value(s) of x1,x2,..
;mstyle  = i          Line style for markline
;mrange  = [y1,y2]    range to draw line in the marking
;xrange  = [x1,x2]    Make plot from x1 to x2
;prange  = [x1,x2]    Plot data in range from x1 to x2
;all     = 1          Show nodata too
;cutlow  = y          Cut out values below y (in Jy)
;adet    = [d1,d2,..] Show only dets d1,d2
;aband   = [b1,b2,..] Show only bands b1,b2
;ascan   = [s1]       Show only scan s1
;autoy   = 1/w        Automatically scale for a 'nice' plot, optional w
;                     for width of median to determine the peak
;error   = 1          Over plot errors
;plambda = p          Plot data with fluxe*wave*p
;xaxis   = 's'        make xaxis = aar.data.s
;yaxis   = 's'        make yaxis = aar.data.s
;ramps   = 1          plot data with different tint seperately 
;cbands  = 1          cut the extra data out
;jump    = j          dont connect parts with separation in x>j
;message = m          write m to the screen
;normal  = [x,y]|y    normalize to got through point x,y or the peak to y
;nrange  = [x1,x2]    do not plot data for range x1..x2
PRO pl,p1,p2,oplot=oplot,ll=ll,flambda=flambda,iras=iras,white=white,$
       cm=cm,time=time,cfact=cfact,coffset=coffset,dets=dets,scans=scans, $
       nu=nu,bands=bands,no_off=nooff,key=key,xontop=xontop,ymin=ymin, $
       mark=mark,xrange=xrange,prange=prange,recur=recur,all=all, $
       cutlow=cutlow,adet=adet,aband=aband,ascan=ascan,autoy=autoy, $
       ylog=ylog,error=error,plambda=plambda,mstyle=mstyle,mrange=mrange, $
       mthick=mthick,mcolor=mcolor,ramps=ramps,yaxis=yaxis,xaxis=xaxis, $
       ndet=ndet,nband=nband,nscan=nscan,noclip=noclip, $
       title=ttl,cshift=cshift,cbands=cbands,jump=jump,message=message, $
       normal=normal,return=return,nrange=nrange,_extra=_extra

  COMMON COMMON_PLS, sh_plotcount,sh_plotflambda,sh_plotiras,sh_plotcm, $
    sh_plottime,sh_plotnu,sh_plotplambda

; check for valid input which can be two aar an array of data or two
; simple arrays
  usage = 'Usage: pl,aar, pl,arr1 or pl,arr1,arr2'
  CASE n_params() OF 
      0: BEGIN
          print,usage
          return
      END
      1: BEGIN
          IF NOT pl_is_aar(p1) THEN BEGIN
              sp1 = size(p1)
              CASE sp1[0] OF
                  0: BEGIN
                      print,usage
                      return
                  END
                  1: BEGIN
                      aar = pl_define_aar(length=sp1[1])
                      aar.data.wave = indgen(sp1[1])
                      aar.data.flux = p1
                  END
                  ELSE: BEGIN
                      CASE sp1[1] OF
                          2: BEGIN
                              aar = pl_define_aar(length=sp1[2])
                              aar.data.wave = reform(p1[0,*],sp1[2])
                              aar.data.flux = reform(p1[1,*],sp1[2])
                          END
                          ELSE: BEGIN
                              aar = pl_define_aar(length=sp1[2])
                              aar.data.wave = reform(p1[0,*],sp1[2])
                              aar.data.flux = reform(p1[1,*],sp1[2])
                              aar.data.stdev= reform(p1[2,*],sp1[2])
                          END
                      ENDCASE
                  END
              ENDCASE
          ENDIF ELSE BEGIN
              aar = p1
          ENDELSE
      END
      2: BEGIN
          np1 = n_elements(p1)
          IF np1 EQ n_elements(p2) THEN BEGIN
              aar = pl_define_aar(length=np1)
              aar.data.wave = reform(p1,np1)
              aar.data.flux = reform(p2,np1)
          ENDIF ELSE BEGIN
              print,usage
              return
          ENDELSE
      END
  END

; Init  
; Make sure all common variables are initialised
  IF n_elements(sh_plotcount) EQ 0 THEN BEGIN
      sh_plotcount=0
      sh_plotflambda=0
      sh_plotiras=0
      sh_plotcm=0
      sh_plotnu=0
      sh_plottime=0L
      sh_plotplambda=0d0
  ENDIF

  oplot  = keyword_set(oplot  )
  time   = keyword_set(time   )
  cm     = keyword_set(cm     )
  nu     = keyword_set(nu     )
  iras   = keyword_set(iras   )
  flambda= keyword_set(flambda)
  clip   = NOT keyword_set(noclip)

;; Stupid things needed to plot ps fonts too  
  IF !p.font EQ -1 THEN BEGIN
      lambda  = '!7k!X'
      mu      = '!7l!X'
      times   = '!MX!X'
  ENDIF ELSE BEGIN 
      lambda= '!Ml!X'
      mu    = '!Mm!X'
      times = '!M´!X'
  ENDELSE

  cls = !p.color
  xstl = 1
;; standard symbol for plotting 10(histogram) or 1 for time plot  
  psm = 10-9*time
; axes labels defaults:
  xtxt = 'Wavelength ['+mu+'m]'
  yunit = '[Jy]'
  yquant= 'Flux density'
;; Additions to the command to issue
  suffix = ',yrange=[minf,maxf],xs=xstl,/ys,ps=psm,col=cls,xtitl=xtxt,'+ $
    'ytitl=yquant+" "+yunit,title=ttl,_extra=_extra'
  prefix = ''
;; largest float is smaller then infinity
  Inf = !values.f_infinity

; get plotcolors to make overplots in different colors
; Use the dynamic range of the device !p.color*(.1,1.0)  
; Might better use !d.table_size than !d.n_colors
  IF keyword_set(white) THEN BEGIN
      colors = make_array(450,value=!p.color)
  ENDIF ELSE BEGIN
      cols=[!p.color,!d.n_colors*(.1+0.06*[8,13,5,1,11,6,2,10,3,14,9,4,12,7])]
      colors = [cols,cols,cols,cols,cols,cols,cols,cols,cols,cols, $
                cols,cols,cols,cols,cols,cols,cols,cols,cols,cols, $
                cols,cols,cols,cols,cols,cols,cols,cols,cols,cols]
  ENDELSE 

;;Simple arrays make things quicker
  aardataline = aar.data.line
  aardatasdir = aar.data.sdir
  aardatadet  = aar.data.det
  aardatawave = aar.data.wave
  aardataflux = aar.data.flux
  aardatastdev= aar.data.stdev

; Only use the data needed and do the calculations but only do this
; when called by the user:
  IF NOT keyword_set(recur) THEN BEGIN
;; Is this oplot -> get the previous settings. 
      IF oplot THEN BEGIN
; Other units
; check if previous plot was f_lambda
          IF sh_plotflambda THEN flambda = 1
; check if previous plot was lambda*f_lambda
          IF sh_plotiras THEN iras = 1
; check if previous plot was cm-1
          IF sh_plotcm THEN cm = 1
; check if previous plot was nu
          IF sh_plotnu THEN nu = 1
; check if previous plot was time
          IF sh_plottime GT 0L THEN time = 1
; check if previous plot was plambda
          IF sh_plotplambda NE 0d0 THEN plambda = sh_plotplambda
      ENDIF ELSE BEGIN
          sh_plotcount = 0
      ENDELSE
      
      IF keyword_set(xaxis) THEN BEGIN
          foo  = execute('aar.data.wave = aar.data.'+xaxis)
          aardatawave = aar.data.wave
      ENDIF
      
      IF keyword_set(yaxis) THEN BEGIN
          foo  = execute('aar.data.flux = aar.data.'+yaxis)
          aardataflux = aar.data.flux
      ENDIF
      
; Determine the 0-point for the x-axis
      IF time THEN BEGIN
; We use the sh_plottime as the offset
          IF sh_plottime EQ 0L THEN BEGIN
              sh_plottime = pl_itk0(aar)
          ENDIF
          xtxt = "Time [s/24] (offset = "+pl_n2s(sh_plottime)+")"
      ENDIF ELSE BEGIN
          sh_plottime = 0L
      ENDELSE
      
;Index to keep good data
      igood = make_array(n_elements(aardatawave),value=1)
;; Throw out nodata
      igood = igood AND (keyword_set(all) OR ((aar.data.flag AND 16L) NE 16L))
;; Throw out masked data
      igood = igood AND ((aar.data.flag AND 2L^30) NE 2L^30)
      
; Do not use the offband data    
      igood = igood AND ((NOT keyword_set(nooff)) OR (aardataline LT 20))
      
; only want specific line(s)     
      IF n_elements(aband) NE 0 THEN BEGIN
;;Check if this a string(array)
          s = size(aband)
          IF s[s[0]+1] EQ 7 THEN aband = pl_bandname2line(aband)
          igood = igood AND pl_comparray(aardataline,aband)
      ENDIF 
      
; Exclude specific line(s)     
      IF n_elements(nband) NE 0 THEN BEGIN
;;Check if this a string(array)
          s = size(nband)
          IF s[s[0]+1] EQ 7 THEN nband = pl_bandname2line(nband)
          igood = igood AND (pl_comparray(aardataline,nband) EQ 0)
      ENDIF 
      
;Only plot specific detector(s)    
      IF n_elements(adet) NE 0 THEN BEGIN
          igood = igood AND pl_comparray(aardatadet,adet)
      ENDIF 
      
;Exclude specific detector(s)    
      IF n_elements(ndet) NE 0 THEN BEGIN
          igood = igood AND (pl_comparray(aardatadet,ndet) EQ 0)
      ENDIF 
      
;Only want specific scandir(s)    
      IF n_elements(ascan) NE 0 THEN BEGIN
          s = size(ascan)
          IF s[s[0]+1] EQ 7 THEN ascan = pl_scanname2sdir(ascan)
          igood = igood AND pl_comparray(aardatasdir,ascan)
      ENDIF 
      
;Exclude specific scandir(s)    
      IF n_elements(nscan) NE 0 THEN BEGIN
          s = size(nscan)
          IF s[s[0]+1] EQ 7 THEN nscan = pl_scanname2sdir(nscan)
          igood = igood AND (pl_comparray(aardatasdir,ascan) EQ 0)
      ENDIF 
      
;; We only take the data that is actually going to be plotted. There
;; are 4 ways by which data might be outside the plotting window (in
;; xrange).
;; 1) The xrange is specified.
;; 2) The prange is specified.
;; 3) This is oplot and data falls outside the !x.crange
;; 4) The values of !x.range are set
;; We take the narrowest range from this
      xmin = [-Inf]
      xmax = [ Inf]
      IF n_elements(xrange) GE 2 THEN BEGIN
          xmin = [xmin,xrange[0]<xrange[1]]
          xmax = [xmax,xrange[0]>xrange[1]]
      ENDIF
      IF n_elements(prange) GE 2 THEN BEGIN
          xmin = [xmin,prange[0]<prange[1]]
          xmax = [xmax,prange[0]>prange[1]]
      ENDIF
      IF oplot THEN BEGIN
          wrange = pl_xrange()
          xmin = [xmin,wrange[0]<wrange[1]]
          xmax = [xmax,wrange[0]>wrange[1]]
      ENDIF
      IF !x.range[0] NE !x.range[1] THEN BEGIN
          xmin = [xmin,!x.range[0]<!x.range[1]]
          xmax = [xmax,!x.range[0]>!x.range[1]]
      ENDIF
      
;; Take the tightest range if one of the 3 was set
      IF n_elements(xmin) GT 1 THEN BEGIN
          xmin = max(xmin)
          xmax = min(xmax)
; only take the needed data to speed up things
          IF xmin NE  xmax THEN BEGIN
              CASE 1 OF 
                  ;; Select all in wavelength
                  time: xlimits = [min(aardatawave),max(aardatawave)]
                  cm: xlimits = [1d4/xmax,1d4/xmin]
                  nu: xlimits = [3d2/xmax,3d2/xmin]
                  ELSE: xlimits = [xmin,xmax]
              END
              IF clip THEN begin
                  igood = igood AND (aardatawave GE xlimits[0]) AND $
                    (aardatawave LE xlimits[1])
                  IF time THEN BEGIN
                      tlimits = [xmin,xmax]
                      igood = igood AND $
                        ((aar.data.itk-sh_plottime) GE tlimits[0]) AND $
                        ((aar.data.itk-sh_plottime) LE tlimits[1])
                  ENDIF
              ENDIF
          ENDIF
      ENDIF
      
; Now throw out data below certain limit.     
      IF n_elements(cutlow) EQ 1 THEN BEGIN
          igood = igood AND (aardataflux GT cutlow)
      ENDIF
      
; cut out of band data
      IF n_elements(cbands) EQ 1 THEN BEGIN
          edges = [ [0,2000], $ ;unknown
                    [2.38,2.62],[2.6,3.02],[3.0,3.52],[3.5,4.1], $ ;1a,1b,1d,1e
                    [4.08,5.32],[5.3,7.05],[7.0,12.5], $ ;2a,2b,2c
                    [0,0], $    ;non existent
                    [12.45,16.55],[16.35,19.55],[19.5,27.55],[27.5,29.0], $ ;3a,3c,3d,3e
                    [28.9,45.2], $ ;4
                    [0,0],[0,0],[0,0],[0,0],[0,0],[0,0], $ ;nonexist
                    [12.45,16.55],[16.35,19.55],[19.5,27.55], $ ;4a,4c,3d
                    [0,0],[0,0], $
                    [45.0,200] ] ; Maybe LWS
          
          igood = igood AND ((aardatawave GT edges[0,aardataline]) AND $
                             (aardatawave LT edges[1,aardataline]))
      END
      
;remove the nrange data
      IF n_elements(nrange) EQ 2 THEN BEGIN
          igood = igood AND ((aardatawave LT nrange[0]) OR $
                             (aardatawave GT nrange[1]))
      ENDIF

      IF total(igood) EQ 0 THEN BEGIN
          print,'PL: Warning: no data in this wavelength range'
          return
      ENDIF ELSE BEGIN 
          aar = pl_select(aar,igood)
      END
      
;; Store this in the return value
      return = aar
      
      aardataline = aar.data.line
      aardatasdir = aar.data.sdir
      aardatawave = aar.data.wave
      aardataflux = aar.data.flux
      aardatadet  = aar.data.det
      aardatastdev= aar.data.stdev
      
; Now do the calculations requested
;; Normalize the spectra first if wanted
      IF keyword_set(normal) THEN BEGIN
          IF n_elements(normal) EQ 1 THEN BEGIN
              aardataflux = aardataflux/max(aardataflux)*normal
          ENDIF ELSE BEGIN
              idxwave = sort(aardatawave)
              aardataflux = aardataflux * normal[1] / $
                interpol(aardataflux[idxwave],aardatawave[idxwave],normal[0]) 
          ENDELSE
      ENDIF

;;Change the flux scales    
      IF flambda THEN BEGIN
          aardataflux = 3e-12*aardataflux/aardatawave^2d0*1d13
          aardatastdev= 3e-12*aardatastdev/aardatawave^2d0*1d13
          sh_plotflambda = 1
          yquant = 'Flux density'
          yunit = '[10!U-13!N W/m!U2!N/'+mu+'m]'
      ENDIF ELSE BEGIN
          sh_plotflambda = 0
      ENDELSE
      
      IF iras THEN BEGIN
          aardataflux = 3e-12*aardataflux/aardatawave*1d12
          aardatastdev= 3e-12*aardatastdev/aardatawave*1d12
          sh_plotiras = 1
          yquant = 'Flux'
          yunit = '[10!U-12!N W/m!U2!N]'
      ENDIF ELSE BEGIN
          sh_plotiras = 0
      ENDELSE
      
; check if factor is required
      IF keyword_set(cfact) THEN BEGIN
          aardataflux  = aardataflux*cfact 
          aardatastdev = aardatastdev*cfact 
          yquant = yquant+times+pl_f2s(cfact,deci=2) 
      ENDIF
      
; check if offset is required
      IF keyword_set(coffset) THEN BEGIN
          aardataflux  = aardataflux+coffset
          yquant=yquant+'+'+pl_f2s(coffset,deci=2)
      ENDIF
      
;;Do we want a p_lambda plot?    
      IF keyword_set(plambda) THEN BEGIN
          aardataflux = aardataflux*aardatawave^plambda
          aardatastdev= aardatastdev*aardatawave^plambda
          sh_plotplambda = plambda
          yunit = yunit+times+lambda+'!U'+pl_f2s(sh_plotplambda)+'!N'
      ENDIF ELSE BEGIN
          sh_plotplambda = 0d0
      ENDELSE
      
; Do requested modifications to the x-axis    
      IF cm THEN BEGIN
          aardatawave = 1d4/aardatawave
          sh_plotcm = 1
          xtxt = 'Wavenumber [cm!U-1!N]'
      ENDIF ELSE BEGIN
          sh_plotcm = 0
      ENDELSE
      
      IF nu THEN BEGIN
          aardatawave = 3d2/aardatawave
          sh_plotnu = 1
          xtxt = 'Frequency[THz]'
      ENDIF ELSE BEGIN
          sh_plotnu = 0
      ENDELSE
      
; check if shift is required
      IF keyword_set(cshift) THEN BEGIN
          aardatawave  = aardatawave+cshift
      ENDIF
      
      IF time THEN BEGIN
;now assign the time to the wave and sort everthing
          tt = aar.data.itk - sh_plottime
          idx = sort(tt)
          aardatawave  = tt[idx]
          aardataflux  = aardataflux[idx]
          aardatadet   = aardatadet[idx]
          aardatastdev = aardatastdev[idx]
          aardataline  = aardataline[idx] 
          aardatasdir  = aardatasdir[idx]
      ENDIF
      
      aar.data.wave  = aardatawave 
      aar.data.flux  = aardataflux     
      aar.data.det   = aardatadet      
      aar.data.stdev = aardatastdev    
      aar.data.line  = aardataline     
      aar.data.sdir  = aardatasdir 

      IF NOT oplot THEN BEGIN
          
          IF n_elements(ttl) NE 0 THEN BEGIN
;; this allows for calls like pl,aar,title='continuum subtracted
;; spectrum: +object +observer +proposal +user'
              in_keys =  ['object','observer','proposal','user','filename','ra','dec']
              head_keys= ['OBJECT','OBSERVER','EOHAPLID','USERNAME','FILENAME','ATTRA','ATTDEC']
              FOR i=0,n_elements(in_keys)-1 DO BEGIN
                  key_pos = strpos(strlowcase(ttl),'+'+in_keys[i])
                  IF key_pos NE -1 THEN BEGIN
                      aa = pl_headkey(aar.header,head_keys[i])
                      ttl = strmid(ttl,0,key_pos)+aa+ $
                        strmid(ttl,key_pos+strlen(in_keys[i])+1,strlen(ttl))
                  ENDIF
              ENDFOR
          ENDIF ELSE BEGIN
;;  Get title from object
              ttl = pl_headkey(aar.header,'OBJECT')
          ENDELSE
          
          IF keyword_set(xaxis) THEN BEGIN
              xtxt = xaxis
          ENDIF
          
          IF keyword_set(yaxis) THEN BEGIN
              yquant = yaxis
              yunit=''
          ENDIF
          
          IF n_elements(xrange) GE 2 THEN BEGIN
              suffix = suffix+',xrange=xrange[0:1]'
          ENDIF 
;;      ELSE BEGIN
;;        IF n_elements(prange) GE 2 THEN suffix = suffix+',xrange=prange[0:1]'
;;      ENDELSE
          
          IF keyword_set(ll) THEN BEGIN
              suffix = suffix+',/xl,/yl'
          ENDIF ELSE BEGIN
              IF keyword_set(ylog) THEN BEGIN
                  suffix = suffix+',/yl'
              ENDIF
          ENDELSE
          
;; A double axis plot
          IF keyword_set(xontop) THEN BEGIN
              xstl = 9
              suffix = suffix+',ymargin=[4,6]'
              ttl=ttl+'!C!C'
          ENDIF
          
; This is for the yaxis scaling first set default values for the min
; and max flux axis
          IF !y.range[0] EQ !y.range[1] THEN BEGIN
              minf = min(aardataflux,max=maxf)
;        plot,aardataflux,/ys
;        stop
          ENDIF ELSE BEGIN
              minf = !y.range[0]
              maxf = !y.range[1]
          ENDELSE
          
; Do we want autoscaling for the yaxis ? then override the standard
; values and use the yrange settings
          IF keyword_set(autoy) THEN BEGIN
              IF keyword_set(ylog) OR keyword_set(ll) THEN BEGIN
                  minf=.1>minf
              ENDIF ELSE BEGIN
                  minf=-5>minf
              ENDELSE
              width = n_elements(aardatawave)<(31*(autoy EQ 1) + autoy*(autoy NE 1))
              maxf=max(median(aardataflux,width))*1.10
          ENDIF
          
; Is there a minimum value for the yrange (used for negative values
; mostly) overrides the default and autovalues
          IF n_elements(ymin) EQ 1 THEN BEGIN
              minf = ymin
          ENDIF
      ENDIF    ;; Not oplot
  ENDIF ;;not recursive

;; Command to issue to make box
  plbox_command = 'plot,[aardatawave],[aardataflux]'+suffix+ $
    ',/nodata'

;; Do we want the detectors to be printed separately ?
  IF keyword_set(ramps) THEN BEGIN
      aardatatint = aar.data.tint
      psm=time                  ; 0 or 1 if this is time plot
      plotcount = sh_plotcount
      tint = pl_uniq(aardatatint)
      IF NOT oplot THEN BEGIN
          f = execute(plbox_command)
          IF keyword_set(xontop) THEN pl_xontop
          IF keyword_set(error) THEN oploterr,[aardatawave],[aardataflux],[aardatastdev],3
          IF n_elements(mark) GT 0 THEN pl_mark,mark,linestyle=mstyle,mrange=mrange,thick=mthick,color=mcolor
          IF keyword_set(key) THEN pl_key,tint,"Tint: ",colors[plotcount:*]
      ENDIF
      FOR i = n_elements(tint)-1,0,-1 DO BEGIN
          pl,pl_select(aar,aardatatint EQ tint[i]),ll=ll, $
             col=colors[plotcount+i],dets=dets,scans=scans, $
             bands=bands,ps=psm,jump=jump,message=message, $
             _extra=_extra,/oplot,/recur
      ENDFOR
      return
  ENDIF

;; Do we want the detectors to be printed separately ?
  IF keyword_set(dets) THEN BEGIN
      psm=time                  ; 0 or 1 if this is time plot
      plotcount = sh_plotcount
      det = pl_uniq(aardatadet)
      IF NOT oplot THEN BEGIN
          f = execute(plbox_command)
          IF keyword_set(xontop) THEN pl_xontop
          IF keyword_set(error) THEN oploterr,[aardatawave],[aardataflux],[aardatastdev],3
          IF n_elements(mark) GT 0 THEN pl_mark,mark,linestyle=mstyle,mrange=mrange,thick=mthick,color=mcolor
          IF keyword_set(key) THEN pl_key,det,"Det: ",colors[plotcount:*]
      ENDIF
      FOR i = 0,n_elements(det)-1 DO BEGIN
          pl,pl_select(aar,aardatadet EQ det[i]),ll=ll, $
             col=colors[plotcount+i],scans=scans,bands=bands, $
             ps=psm,jump=jump,message=message,_extra=_extra,/oplot,/recur
      ENDFOR
      return
  ENDIF

;; Do we want the scans to be printed separately ?
  IF keyword_set(scans) THEN BEGIN
      plotcount = sh_plotcount
      sdir = pl_uniq(aardatasdir)
      IF NOT oplot THEN BEGIN
          f = execute(plbox_command)
          IF keyword_set(xontop) THEN pl_xontop
          IF keyword_set(error) THEN oploterr,[aardatawave],[aardataflux],[aardatastdev],3
          IF n_elements(mark) GT 0 THEN pl_mark,mark,linestyle=mstyle,mrange=mrange,thick=mthick,color=mcolor
          IF keyword_set(key) THEN pl_key,sdir,"Scan: ",colors[plotcount:*]
      ENDIF
      FOR i = 0,n_elements(sdir)-1 DO BEGIN
          pl,pl_select(aar,aardatasdir EQ sdir[i]),ll=ll, $
             col=colors[plotcount+i],bands=bands,ps=psm,jump=jump, $
             message=message,_extra=_extra,/oplot,/recur
      ENDFOR
      return
  ENDIF

;; Do we want the lines to be printed separately ?
  IF keyword_set(bands) THEN BEGIN
      plotcount = sh_plotcount
      line = pl_uniq(aardataline)
      IF NOT oplot THEN BEGIN
          f = execute(plbox_command)
          IF keyword_set(xontop) THEN pl_xontop
          IF keyword_set(error) THEN oploterr,[aardatawave],[aardataflux],[aardatastdev],3
          IF n_elements(mark) GT 0 THEN pl_mark,mark,linestyle=mstyle,mrange=mrange,thick=mthick,color=mcolor
          IF keyword_set(key) THEN pl_key,line,"Line: ",colors[plotcount:*]
      ENDIF
      FOR i = 0,n_elements(line)-1 DO BEGIN
          pl,pl_select(aar,aardataline EQ line[i]),ll=ll, $
             col=colors[plotcount+i],ps=psm,jump=jump, $
             message=message,_extra=_extra,/oplot,/recur
      ENDFOR
      return
  ENDIF

;; Do we want the jump to be disconnected?
  IF keyword_set(jump) THEN BEGIN
      ;; split accoring to holes
      x_difference = abs(aardatawave - shift(aardatawave,1))
      where_jump = where(x_difference GT jump,n_jumps)
      ;; Only if we find two points which are marked as jump do we need
      ;; to plot the data separately
      CASE 1 OF 
          n_jumps GE 256: BEGIN
              print,'PL: Too many jumps detected,connecting points'
          END
          n_jumps GE 2: BEGIN
              plotcount = sh_plotcount
              wave_jump=aardatawave[[where_jump,n_elements(aardatawave)-1]]
              IF NOT oplot THEN BEGIN
                  f = execute(plbox_command)
                  IF keyword_set(xontop) THEN pl_xontop
                  IF keyword_set(error) THEN oploterr,[aardatawave],[aardataflux], $
                    [aardatastdev],3
                  IF n_elements(mark) GT 0 THEN pl_mark,mark,linestyle=mstyle,mrange=mrange,thick=mthick,color=mcolor
              ENDIF
              FOR i = 0,n_jumps-1 DO BEGIN
                  pl,pl_select(aar,(aardatawave GE wave_jump[i]) AND  $
                               (aardatawave LT wave_jump[i+1])), $
                     ll=ll,col=colors[plotcount],ps=psm, $
                     message=message,_extra=_extra,/oplot,/recur
              ENDFOR
              return
          END
          ELSE: BEGIN
              ;; No jumps detected
          END
      ENDCASE
  ENDIF

;; We want errors so first plot the frame then the errors and next do
;; an oplot of the data    
  IF keyword_set(error) THEN BEGIN
      IF NOT oplot THEN foo = execute('plot,[aardatawave],[aardataflux],/nodata'+suffix)
      oploterr,[aardatawave],[aardataflux],[aardatastdev],3
      oplot = 1
;    sh_plotcount = sh_plotcount+1
  ENDIF

  IF oplot THEN BEGIN
; Increase the plotcount
      cls = colors[sh_plotcount]
      prefix = 'o'
      suffix = ',ps=psm,col=cls,_extra=_extra'
  ENDIF

  foo = execute(prefix+'plot,[aardatawave],[aardataflux]'+suffix)
  IF keyword_set(message) THEN print,message
  sh_plotcount = sh_plotcount + 1

  IF keyword_set(xontop) THEN pl_xontop
  IF n_elements(mark) GT 0 THEN pl_mark,mark,linestyle=mstyle,mrange=mrange,thick=mthick,color=mcolor

END

; tin and tout should be specified as reals and not as integers
function ck_modbb,q,tin=tin,tout=tout,index=index,$
                  n0=n0,r0=r0,distance=distance,grainsize=grainsize,$
                  steps=steps

  if ((n_params() eq 0 ) or keyword_set(help)) then begin
  endif

; Check for valid input
  if not is_aar(q) then error,'F','No valid AAR structure specified!'
  if (n_elements(tin) eq 0) then error,'F','No T_in specified!'
  if (n_elements(tout) eq 0) then error,'F','No T_out specified!'
  if (n_elements(index) eq 0) then begin
    error,'W','No density index specified, value set to 2'
    index = 2
  endif
  if (n_elements(n0) eq 0) then error,'F','No zero density specified!'
  if (n_elements(r0) eq 0) then error,'F','No inner radius specified!'
  if (n_elements(distance) eq 0) then error,'F','No distance specified!'
  if (n_elements(grainsize) eq 0) then begin
    error,'W','No grain size specified, value set to 0.1 micron'
    grainsize = 0.1
  endif
  if (n_elements(steps) eq 0) then begin 
    error,'W','Number of steps not defined, value set to 10 steps'
    steps = 10
  endif
  
  d = distance  * 3.0857e18 ; correction from pc to cm
  a = grainsize  * 1e-4 ; correction from micron to cm 
  
  ret = q
  fnu = q
  ret.data.flux = 0
;  pl,ret,/xl,yra=[0,5e13]
 
  for i = 0,steps-1 do begin
    t = tin - i * (tin-tout)/steps
    power = (t/tin)^(2*index - 6)
    print,'t = ',t
    fnu = multiply(q,sh_bb(q,t,0))
 ;   pl,fnu,/opl
    fnu.data.flux = fnu.data.flux * power * ((tin-tout)/steps)
 ;   pl,fnu,/opl
    ret = sh_calcaar(ret,add=fnu)
 ;   pl,ret,/opl
  endfor
  
  extra = r0/D
  factor = 4 * !dpi * a * a * r0 * n0 * extra * extra / (3 - index)
  print,'factor = ',factor
  ret.data.flux = ret.data.flux * factor
  
  return,ret
end

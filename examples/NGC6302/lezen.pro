function lezen,tabelnaam,header=header

openr,lun,tabelnaam,/get_lun

dummy = string(80)

for j = 1,header do begin

readf,lun,dummy
endfor

tjop = fltarr(30000)
tjap = fltarr(30000)
i = 0
x = 0.0
y = 0.0
while (not eof(lun)) do begin
   readf,lun,x,y
   tjop(i) = x
   tjap(i) = y
   i = i + 1
endwhile
tabel = fltarr(2,i)
for k=0,i-1 do begin
   tabel(0,k) = tjop(k)
   tabel(1,k) = tjap(k)
endfor
close,lun
free_lun, lun
print,'Aantal elementen =',i
return,tabel
end




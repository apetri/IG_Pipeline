
fname = "ICs/inputspec_ics.txt"

spawn,"wc "+fname,result
result=strtrim(result, 1)
lines=long(result)
lines=lines(0)



en=fltarr(4, lines-1)
openr,1,fname
zstart=0.0
growthfactor=0.0
readf,1,zstart, growthfactor
readf,1,en
close,1

k_linear = en(0,*) 
d2_linear = en(1,*)

k_nonlinear = en(2,*)
d2_nonlinear = en(3,*)

plot, k_linear, d2_linear, /xlog, /ylog, yrange=[0.0001,1000]

oplot, k_nonlinear, d2_nonlinear, color=255



end



	function ran2(idum)
c	implicit real*8 (a-h,o-z)
	parameter(m=714025,ia=1366,ic=150889,rm=1./m)
	dimension ir(97)
	data iff /0/
	if(idum.lt.0.or.iff.eq.0) then
	iff=1
	idum=mod(ic-idum,m)
	do 11 j=1,97
	 idum=mod(ia*idum+ic,m)
	 ir(j)=idum
11	continue
	idum=mod(ia*idum+ic,m)
	iy=idum
	endif
	j=1+(97*iy)/m
	if(j.gt.97 .or. j.lt.1) pause
	iy=ir(j)
	ran2=iy*rm
	idum=mod(ia*idum+ic,m)
	ir(j)=idum
	return
	end


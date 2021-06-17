	double precision function runif(seed)
	implicit real*8 (a-h,o-z)

	runif=4.d0*seed*(1.d0-seed)
	seed=runif
	return
	end

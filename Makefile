FC              = gfortran
FCFLAGS         = -g -c -ffixed-line-length-132 -w -fPIC
FCOPT           = -O3
PROGRAMS        = Bhat
A_TARGET        = libBhat_f95.a
SO_TARGET       = libBhat_f95.so

.PHONY:  all a_lib so_lib

OBJS =  \
	hessian.o\
	newton.o\
	simplex.o\
	migrad.o\
	dqstep.o\
	mcmc.het-gam2.o\
	mcmc.het-lnorm2.o\
	berkson.o\
	global.o\
	plkhci.o\
	ftrafo.o\
	btrafo.o\
	ran2.o\
	randgs.o

OBJS_NOPT1 = axeb.o rs.o
OBJS_NOPT2 = alngam.o

all:	a_lib so_lib

a_lib:	$(A_TARGET)

so_lib: $(SO_TARGET)

%.o:	%.f
	$(FC) $(FCFLAGS) $(FCOPT) -c $<

$(OBJS_NOPT1): %.o: %.f
	$(FC) -c -O0 $(FCFLAGS) $<

$(OBJS_NOPT2): %.o: %.f90
	$(FC) -c $(FCFLAGS) $<


ran_lib.stamp:
	@echo building ranlib
	@cd ./ranlib && $(MAKE) && touch $@


$(A_TARGET): ran_lib.stamp $(OBJS) $(OBJS_NOPT1) $(OBJS_NOPT2) main.o
	#gfortran -g -c -w alngam.f90
	#gfortran -g -c -w -O0 axeb.f
	#gfortran -g -c -w -O0 rs.f
	@echo building archive
	$(AR) crs $@ $(filter-out ran_lib.stamp,$^) ./ranlib/*.o

$(SO_TARGET): ran_lib.stamp $(OBJS) $(OBJS_NOPT1) $(OBJS_NOPT2)
	@echo building shared object
	$(FC) -shared -o $@ $(filter-out ran_lib.stamp,$^) ./ranlib/*.o

clean:
	cd ./ranlib && $(MAKE) clean
	rm -f ran_lib.stamp
	rm -f $(OBJS) $(OBJS_NOPT1) $(OBJS_NOPT2) main.o
	rm -f $(A_TARGET) $(SO_TARGET)

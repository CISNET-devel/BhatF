#Makefile for building ranlib

FC = gfortran
FFLAGS = -g -c -O0 -w -fPIC

.PHONY: all clean

RANLIBOBJECTS = advnst.o  genmn.o   genunf.o  ignpoi.o ranf.o sexpo.o \
	genbet.o  genmul.o  getcgn.o  ignuin.o  sgamma.o \
	genchi.o  gennch.o  getsd.o   initgn.o  setall.o  snorm.o \
	genexp.o  gennf.o   ignbin.o  inrgcm.o  mltmod.o  setant.o  spofa.o \
	genf.o    gennor.o  ignlgi.o  lennob.o  phrtsd.o  setgmn.o \
	gengam.o  genprm.o  ignnbn.o  qrgnin.o  setsd.o  sdot.o

all:	$(RANLIBOBJECTS)

.o:	%.f
	$(FC) $(FFLAGS) $<

clean:
	rm -f $(RANLIBOBJECTS)

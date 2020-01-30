SHELL = /bin/sh
FC = pgf90
FFLAGS= -O -Mvect=sse -Mextend -Msecond_underscore -Bstatic_pgi -Mlfs
PROGRAM = txt2ioapi
LIBS =   -L/models/CMAQ/lib64/ioapi-3.1/Linux2_x86_64pg -lioapi \
         -L/models/CMAQ/lib/netcdf-3.6.0-p1/lib -lnetcdf
INCLUDE = -I/models/CMAQ/lib64/ioapi-3.1/ioapi/fixed_src \
          -I/models/CMAQ/lib/netcdf-3.6.0-p1/include \
          -I$(CURDIR)/INCLDIR

RAW = parsline.o \
      checkmem.o \
      cnvt_efs.o \
      cnvt_cantyp.o \
      cnvt_lai.o \
      cnvt_w126.o \
      cnvt_LDF.o \
      cnvt_arid.o \
      cnvt_nonarid.o \
      cnvt_fert.o \
      cnvt_nitrogen.o \
      cnvt_landtype.o \
      $(PROGRAM).o

.f.o:
	$(FC) $(FFLAGS) $(INCLUDE) -c -o $@ $<

.F.o:
	$(FC) $(FFLAGS) $(INCLUDE) -c -o $@ $<

$(PROGRAM):	$(RAW)
	$(FC) $(FFLAGS) $(INCLUDE) -o $(@) $(RAW) $(LIBS)

clean:
	rm -f $(PROGRAM)  *.o *.core

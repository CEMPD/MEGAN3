SHELL=/bin/sh
FC= pgf90
FFLAGS= -O -Mvect=sse -Mextend -Msecond_underscore -Mlfs -Bstatic_pgi
PROGRAM= megvea
LIBS =   -L/models/CMAQ/lib64/ioapi-3.1/Linux2_x86_64pg -lioapi \
         -L/models/CMAQ/lib/netcdf-3.6.0-p1/lib -lnetcdf
INCLUDE = -I/models/CMAQ/lib64/ioapi-3.1/ioapi/fixed_src \
          -I/models/CMAQ/lib/netcdf-3.6.0-p1/include \
          -I$(CURDIR)


OBJS = \
     checkmem.o \
     findlai.o  \
     $(PROGRAM).o

MODULES= *mod.f

#-----------------------------------------------------------------------
# line to allow file.mod files to be recognised by GNU make
%.o : %.mod
#-----------------------------------------------------------------------

.SUFFIXES : .f .f .o 

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE) $<
.F.o:
	$(FC) -c $(FFLAGS) $(INCLUDE) $<

#-----------------------------------------------------------------------

$(PROGRAM):	$(OBJS)
	$(FC) $(FFLAGS) $(INCLUDE) -o $(@) $(OBJS) $(LIBS)

#-----------------------------------------------------------------------

clean:
	rm -f $(PROGRAM) *.o *.mod *.core
#-----------------------------------------------------------------------



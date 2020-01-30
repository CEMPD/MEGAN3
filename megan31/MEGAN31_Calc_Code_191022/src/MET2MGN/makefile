PGF90 = /usr/local/pgi/linux86-64/bin/pgf90
FCFLGS = -O -Mextend -Mlfs -Msecond_underscore
FCFLGS2 = -O -Mbyteswapio -Mextend -Mlfs -Msecond_underscore
LNKFLGS = -Bstatic_pgi

OBJCTS = \
met2mgn.o \
readpar.o \
readmcip.o \
micromet.o \
caldate.o \
lcpgeo.o \
interp_lcp.o \
procmm5.o \
readmm5.o \
xtod.o \
vertmap.o \
juldate.o \
utmgeo.o


LIBS = -L/models/CMAQ/lib64/ioapi_3/Linux2_x86_64pg -lioapi \
       -L/models/CMAQ/lib64/netcdf-3.6.0-p1/lib -lnetcdf
INC = -I/models/CMAQ/lib64/ioapi_3/ioapi/fixed_src

met2mgn: $(OBJCTS)
	$(PGF90) $(LNKFLAGS) $(OBJCTS) $(LIBS) -o met2mgn
met2mgn.o:
	$(PGF90) -c $(FCFLGS2) $(INC) met2mgn.f
readpar.o:
	$(PGF90) -c $(FCFLGS) $(INC) readpar.f
readmcip.o:
	$(PGF90) -c $(FCFLGS2) $(INC) readmcip.f
micromet.o:
	$(PGF90) -c $(FCFLGS2) $(INC) micromet.f
caldate.o:
	$(PGF90) -c $(FCFLGS2) $(INC) caldate.f
lcpgeo.o:
	$(PGF90) -c $(FCFLGS2) ${INC} lcpgeo.f
interp_lcp.o:
	$(PGF90) -c $(FCFLGS2) ${INC} interp_lcp.f
procmm5.o:
	$(PGF90) -c $(FCFLGS2) ${INC} procmm5.f
readmm5.o:
	$(PGF90) -c $(FCFLGS2) ${INC} readmm5.f
xtod.o:
	$(PGF90) -c $(FCFLGS2) ${INC} xtod.f
vertmap.o:
	$(PGF90) -c $(FCFLGS2) ${INC} vertmap.f
juldate.o:
	$(PGF90) -c $(FCFLGS2) ${INC} juldate.f
utmgeo.o:
	$(PGF90) -c $(FCFLGS2) ${INC} utmgeo.f
clean:
	rm -rf *.o met2mgn


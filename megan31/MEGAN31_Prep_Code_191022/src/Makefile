.SUFFIXES: .o .f90

F90    = $(FC) 
#LIBS   = -L$(NETCDF_DIR)/lib -lnetcdf 
LIBS   = -L$(NETCDF_DIR)/lib -lnetcdf
INCLUDE_MODULES = -I$(NETCDF_DIR)/include
FFLAGS = -g -Bstatic_pgi

OBJS_PREPMEGAN4CMAQ_LAI = misc_definitions_module.o\
       constants_module.o\
       bio_types.o\
       area_mapper.o\
       prepmegan4cmaq_lai.o

OBJS_PREPMEGAN4CMAQ_CANTYPE = misc_definitions_module.o\
       constants_module.o\
       bio_types.o\
       area_mapper.o\
       prepmegan4cmaq_cantype.o

OBJS_PREPMEGAN4CMAQ_EF = misc_definitions_module.o\
       constants_module.o\
       bio_types.o\
       area_mapper.o\
       prepmegan4cmaq_ef.o

OBJS_PREPMEGAN4CMAQ_GRWFORM = misc_definitions_module.o\
       constants_module.o\
       bio_types.o\
       area_mapper.o\
       prepmegan4cmaq_grwform.o

OBJS_PREPMEGAN4CMAQ_ECOTYPE = misc_definitions_module.o\
       constants_module.o\
       bio_types.o\
       area_mapper.o\
       prepmegan4cmaq_ecotype.o

OBJS_PREPMEGAN4CMAQ_W126 = misc_definitions_module.o\
       constants_module.o\
       bio_types.o\
       area_mapper.o\
       prepmegan4cmaq_w126.o

OBJS_PREPMEGAN4CMAQ_FERT = misc_definitions_module.o\
       constants_module.o\
       bio_types.o\
       area_mapper.o\
       prepmegan4cmaq_fert.o

OBJS_PREPMEGAN4CMAQ_NON_ARID = misc_definitions_module.o\
       constants_module.o\
       bio_types.o\
       area_mapper.o\
       prepmegan4cmaq_non_arid.o

OBJS_PREPMEGAN4CMAQ_ARID = misc_definitions_module.o\
       constants_module.o\
       bio_types.o\
       area_mapper.o\
       prepmegan4cmaq_arid.o

OBJS_PREPMEGAN4CMAQ_LANDTYPE = misc_definitions_module.o\
       constants_module.o\
       bio_types.o\
       area_mapper.o\
       prepmegan4cmaq_landtype.o

OBJS_PREPMEGAN4CMAQ_NITROGEN = misc_definitions_module.o\
       constants_module.o\
       bio_types.o\
       area_mapper.o\
       prepmegan4cmaq_nitrogen.o

.f90.o:
		${F90} ${FFLAGS} -c ${INCLUDE_MODULES} $<

##dependencies

prepmegan4cmaq_lai.x: ${OBJS_PREPMEGAN4CMAQ_LAI}
		 ${F90} -o $@ ${FFLAGS} ${OBJS_PREPMEGAN4CMAQ_LAI} ${LIBS} 

prepmegan4cmaq_cantype.x: ${OBJS_PREPMEGAN4CMAQ_CANTYPE}
		 ${F90} -o $@ ${FFLAGS} ${OBJS_PREPMEGAN4CMAQ_CANTYPE} ${LIBS} 

prepmegan4cmaq_ef.x: ${OBJS_PREPMEGAN4CMAQ_EF}
		 ${F90} -o $@ ${FFLAGS} ${OBJS_PREPMEGAN4CMAQ_EF} ${LIBS} 

prepmegan4cmaq_grwform.x: ${OBJS_PREPMEGAN4CMAQ_GRWFORM}
		 ${F90} -o $@ ${FFLAGS} ${OBJS_PREPMEGAN4CMAQ_GRWFORM} ${LIBS} 

prepmegan4cmaq_ecotype.x: ${OBJS_PREPMEGAN4CMAQ_ECOTYPE}
		 ${F90} -o $@ ${FFLAGS} ${OBJS_PREPMEGAN4CMAQ_ECOTYPE} ${LIBS} 

prepmegan4cmaq_w126.x: ${OBJS_PREPMEGAN4CMAQ_W126}
		 ${F90} -o $@ ${FFLAGS} ${OBJS_PREPMEGAN4CMAQ_W126} ${LIBS} 

prepmegan4cmaq_fert.x: ${OBJS_PREPMEGAN4CMAQ_FERT}
		 ${F90} -o $@ ${FFLAGS} ${OBJS_PREPMEGAN4CMAQ_FERT} ${LIBS}

prepmegan4cmaq_non_arid.x: ${OBJS_PREPMEGAN4CMAQ_NON_ARID}
		 ${F90} -o $@ ${FFLAGS} ${OBJS_PREPMEGAN4CMAQ_NON_ARID} ${LIBS}

prepmegan4cmaq_arid.x: ${OBJS_PREPMEGAN4CMAQ_ARID}
		 ${F90} -o $@ ${FFLAGS} ${OBJS_PREPMEGAN4CMAQ_ARID} ${LIBS}

prepmegan4cmaq_landtype.x: ${OBJS_PREPMEGAN4CMAQ_LANDTYPE}
		 ${F90} -o $@ ${FFLAGS} ${OBJS_PREPMEGAN4CMAQ_LANDTYPE} ${LIBS}

prepmegan4cmaq_nitrogen.x: ${OBJS_PREPMEGAN4CMAQ_NITROGEN}
		 ${F90} -o $@ ${FFLAGS} ${OBJS_PREPMEGAN4CMAQ_NITROGEN} ${LIBS}

cleanup:
		rm -f *.o *.mod

clean:
		rm -f core *.o *.mod *~
		rm -f prepmegan4cmaq_lai.x prepmegan4cmaq_cantype.x prepmegan4cmaq_ef.x prepmegan4cmaq_grwform.x prepmegan4cmaq_ecotype.x prepmegan4cmaq_w126.x prepmegan4cmaq_nitrogen 

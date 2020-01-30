#! /bin/csh -f
########################################################################
source setcase.csh
## Directory setups
setenv PRJ USA 
setenv PROMPTFLAG N

# Program directory
setenv PROG   megsea
setenv EXEDIR /disk8/MEGAN3_1/source_code/MEGAN3_1/src/MEGSEA
setenv EXE    $EXEDIR/$PROG

# Input map data directory
setenv INPDIR $MGNINP

# Met data directory
setenv METDIR $MGNINP/MGNMET

# Output directory
setenv OUTDIR $MGNINT

# Log directory
setenv LOGDIR $MGNLOG/megsea
mkdir -p $LOGDIR

########################################################################

setenv EPISDATE 2013145 # Episode start date

foreach dom ( 36 12 )
set JD = $EPISDATE
while ($JD <= 2013198)
set GD = `yj2ymd $JD`
set m = `echo $GD | cut -c 5-6 | sed 's/^0//g'`
########################################################################
# Set up time and date to process
setenv SDATE $JD        #start date
setenv STIME 0
setenv RLENG 250000
setenv MONTH $m

@ JDnext = $JD + 1

########################################################################

########################################################################
# Set up for MEGAN
setenv RUN_MEGAN   Y       # Run megan?

# Grid definition
setenv GRIDDESC $cwd/GRIDDESC
setenv GDNAM3D tceq_${dom}km 

########### START OF SETTING BDSNP SOIL NO ###############
# New in MEGAN3.1
# YL (BDSNP=N) or Berekely Dalhousie (BDSNP=Y) for soil NO estimation?
setenv BDSNP_YN        Y      # [Y/N]: Y uses Berekely Dalhousie for soil NO calculation;
                              # N uses YL95 (same as older MEGAN versions);
                              # if Y is chosen, user needs to provide additional input data
setenv EPIC            N      # [Y/N]: Y uses EPIC model output for fertilizer;
                              # N uses default MEGAN provided fertilizer input
if ( $BDSNP_YN == Y ) then
setenv PX_VERSION      Y      # MCIP must be PX version when using BDSNP option

if ( $JD == $EPISDATE ) then
setenv INITIAL_RUN     Y      # Initial run if the model hasn't been run before;
                              # otherwise set to false and a restart file is needed
else
setenv INITIAL_RUN     N      # Initial run if the model hasn't been run before;
                              # otherwise set to false and a restart file is needed
# Restart file from running previous day
setenv SOILINSTATE $MGNINT/SOILINSTATE.$GDNAM3D.$JD.ncf
endif # INITIAL_RUN
# Additional input files for BDSNP algorithm
# climate file, nonarid
setenv CLIMAFILE $INPDIR/MAP/ARID_${GDNAM3D}.ncf
# climate file, arid
setenv CLIMNAFILE $INPDIR/MAP/NONARID_${GDNAM3D}.ncf
# biome type file
setenv LANDTYPEFILE $INPDIR/MAP/LANDTYPE_${GDNAM3D}.ncf
# Nitrogen deposition file
setenv NDEPFILE $INPDIR/MAP/NDEP_${GDNAM3D}.ncf
# fertilizer reservoir file
if ( $EPIC == N ) then
setenv FERTRESFILE $INPDIR/MAP/FERT_${GDNAM3D}.ncf
else
# if using EPIC model for fertilizer input
setenv EPICRESFILE $INPDIR/MAP/EPIC_$GDNAM3D.ncf
endif

# Restart file for next day
setenv SOILOUT $MGNINT/SOILINSTATE.$GDNAM3D.$JDnext.ncf
rm -f $SOILOUT

endif
######## END OF SETTINGS FOR BDSNP NO ##############

# CANTYP
setenv CANTYP $INPDIR/MAP/CT3_$GDNAM3D.ncf

# LAIS46
setenv LAIS46 $INPDIR/MAP/LAI3_$GDNAM3D.ncf

# MGNMET
setenv MGNMET $METDIR/MET.MEGAN.$GDNAM3D.rad45.${JD}.ncf

# Output
setenv MGNSEA $MGNINT/MGNSEA.$GDNAM3D.${SDATE}.ncf

########################################################################
## Run MEGAN
if ( $RUN_MEGAN == 'Y' ) then
   rm -f $MGNSEA
   $EXE | tee $LOGDIR/log.run.$PROG.$GDNAM3D.$SDATE.txt
endif

@ JD++
end  # End while JD

end # dom

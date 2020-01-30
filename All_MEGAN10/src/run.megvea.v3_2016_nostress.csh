#! /bin/csh -f
########################################################################
#source setcase.csh
## Directory setups
setenv PRJ USA 
setenv PROMPTFLAG N

# Program directory
setenv PROG   megvea
#setenv EXEDIR $MGNHOME/bin
setenv EXE     ./$PROG

# Input map data directory
setenv INPDIR ../inputs

# Output directory
setenv OUTDIR ../outputs

# Log directory
setenv LOGDIR ./LOGS
mkdir -p $LOGDIR

########################################################################

setenv Layers 5	         # canopy vertical layers, default is 5

########################################################################

foreach dom ( 12US1 )
foreach scen ( J3 ) #J4
set JD = 2016002
while ($JD <= 2016002)
########################################################################
# Set up time and date to process
setenv SDATE $JD        #start date
setenv STIME 0
setenv RLENG 250000

########################################################################

########################################################################
# Set up for MEGAN
setenv RUN_MEGAN   Y       # Run megan?

# Grid definition
setenv GRIDDESC $cwd/GRIDDESC
setenv GDNAM3D    ${dom} 

# LAIS46
setenv LAIS46 $INPDIR/LAI3_$GDNAM3D.ncf

# CANMET
setenv CANMET $INPDIR/CANMET.$GDNAM3D.${SDATE}.ncf

# DailyMET
setenv DailyMET $INPDIR/MET/DAYMET.$GDNAM3D.${SDATE}.ncf

# LDFILE
setenv LDFILE $INPDIR/LDF_$GDNAM3D.$scen.ncf

# Output
setenv MGNERS $OUTDIR/MGNERS.$GDNAM3D.$scen.${SDATE}.nostress.ncf

########################################################################
## Run MEGAN
setenv LOGFILE $LOGDIR/run.$PROG.$GDNAM3D.$scen.$SDATE.nostress.log
if ( -e $LOGFILE ) rm $LOGFILE
if ( $RUN_MEGAN == 'Y' ) then
   rm -f $MGNERS
   $EXE | tee $LOGDIR/log.run.$PROG.$GDNAM3D.$scen.$SDATE.nostress.txt
   #$EXE 
endif

@ JD++
end  # End while JD

end # scen
end # dom

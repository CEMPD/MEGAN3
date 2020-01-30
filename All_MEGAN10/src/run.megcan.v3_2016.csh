#! /bin/csh -f
########################################################################
#source setcase.csh
## Directory setups
setenv PROMPTFLAG N

# Program directory
setenv PROG   megcan
setenv EXE    ./$PROG

# Input map data directory
setenv INPDIR ../inputs

# MCIP input directory
#setenv METDIR $MGNINP/MGNMET

# Output directory
setenv OUTDIR ../outputs

# Log directory
setenv LOGDIR ./LOGS
mkdir -p $LOGDIR
########################################################################

foreach dom ( 12US1 )
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
setenv GDNAM3D  ${dom} 
setenf LOGFILE  $LOGDIR/run.$PROG.$GDNAM3D.$SDATE.log

# CANTYP
setenv CANTYP $INPDIR/CT3_$GDNAM3D.ncf

# LAIS46
setenv LAIS46 $INPDIR/LAI3_$GDNAM3D.ncf

# MGNMET
setenv MGNMET $INPDIR/MET.MEGAN.$GDNAM3D.${SDATE}.ncf

# Output
setenv CANMET $OUTDIR/CANMET.$GDNAM3D.${SDATE}.ncf

########################################################################
## Run MEGAN
if ( $RUN_MEGAN == 'Y' ) then
   rm -f $CANMET
   #$EXE | tee $LOGDIR/run.$PROG.$GDNAM3D.$SDATE.log
   $EXE 
endif

@ JD++
end  # End while JD

end  # dom

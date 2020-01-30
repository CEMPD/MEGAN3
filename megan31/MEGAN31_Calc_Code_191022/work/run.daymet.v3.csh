#! /bin/csh -f
########################################################################
## Common setups
source setcase.csh

## Directory setups
setenv PRJ USA 
setenv PROMPTFLAG N

# Program directory
setenv PROG daymet
setenv EXEDIR /disk8/MEGAN3/source_code/MEGAN3/src/DAYMET
setenv EXE    $EXEDIR/$PROG

# Input MCIP met directory
setenv METDIR $MGNINP/MGNMET

# Output directory
setenv OUTDIR $MGNINT

# Log directory
setenv LOGDIR $MGNLOG/daymet
########################################################################

# Grid definition
#foreach dom ( 36 12 )
foreach dom ( 36 )
setenv GRIDDESC $cwd/GRIDDESC
setenv GDNAM3D tceq_${dom}km 

setenv EPISDATE 2013145	# 2013/05/25        #Episode start date
setenv STIME 0		    #start time
setenv RLENG 240000	    # time step of meteorology files
setenv NDAYS 54	            # number of meteorology files

# loop over episode
set i = 1
while ( $i <= $NDAYS )
set j = `printf "%03i" $i`
@ JDATE = $EPISDATE + $i - 1

# MGNMET
setenv MGNMET$j $METDIR/MET.MEGAN.$GDNAM3D.rad45.${JDATE}.ncf
@ i++
end 

# Output
setenv DailyMET $OUTDIR/DAYMET.$GDNAM3D.ncf

########################################################################
## Run MEGAN
rm -f $DailyMET
if ( ! -e $LOGDIR ) mkdir -p $LOGDIR
$EXE | tee $LOGDIR/log.run.$PROG.$GDNAM3D.txt

end # dom

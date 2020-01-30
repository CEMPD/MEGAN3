#! /bin/csh -f
########################################################################
source setcase.csh
## Directory setups
setenv PRJ USA 
setenv PROMPTFLAG N

# Program directory
setenv PROG   megvea
setenv EXEDIR /disk8/MEGAN3_1/source_code/MEGAN3_1/src/MEGVEA
setenv EXE    $EXEDIR/$PROG

# Input map data directory
setenv INPDIR $MGNINP/MAP

# Output directory
setenv OUTDIR $MGNINT

# Log directory
setenv LOGDIR $MGNLOG/megvea
mkdir -p $LOGDIR

########################################################################

setenv Layers 5	         # canopy vertical layers, default is 5

# User's options to select specific emission activity factors to be applied
setenv GAMAQ_YN   Y      # [Y/N]: Y applies air quality stress; default is N
		         #        if set to Y, user needs to set AQFILE below
setenv GAMCO2_YN  Y      # [Y/N]: Y applies emission response to CO2; default is N
setenv GAMHW_YN   Y      # [Y/N]: Y applies emission response to high wind storm; default is N
setenv GAMHT_YN   Y      # [Y/N]: Y applies emission response to high temperature; default is N
setenv GAMLT_YN   Y      # [Y/N]: Y applies emission response to low temperature; default is N
setenv GAMSM_YN   Y      # [Y/N]: Y applies emission response to soil moisture; default is N
setenv GAMBD_YN   Y      # [Y/N]: Y applies bidirectional exchange LAI response; default is N

########################################################################

foreach dom ( 36 12)
foreach scen ( J4 )
set JD = 2013145
while ($JD <= 2013198)
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
setenv GDNAM3D tceq_${dom}km 

# number of previous days used for Max/Min temeprature/wind speed
if ( $JD == 2013145 ) then
setenv N_MaxT 1		# Number of past days for maximum temperature to be used
                        # [neglected if GAMHT_YN is set to N]
setenv N_MinT 1		# Number of past days for minimum temeprature to be used
                        # [neglected if GAMLT_YN is set to N]
setenv N_MaxWS 1        # number of past days for maximum wind speed to be used
                        # [neglected if GAMHW_YN is set to N]
else
setenv N_MaxT 2         # Number of past days for maximum temperature to be used
                        # [neglected if GAMHT_YN is set to N]
setenv N_MinT 2         # Number of past days for minimum temeprature to be used
                        # [neglected if GAMLT_YN is set to N]
setenv N_MaxWS 2        # Number of past days for maximum wind speed to be used
                        # [neglected if GAMHW_YN is set to N]
endif

# LAIS46
setenv LAIS46 $INPDIR/LAI3_$GDNAM3D.ncf

# CANMET
setenv CANMET $MGNINT/CANMET.$GDNAM3D.${SDATE}.ncf

# DailyMET
setenv DailyMET $MGNINT/DAYMET.$GDNAM3D.ncf

# MEGSEA output
setenv SMFILE $MGNINT/MGNSEA.$GDNAM3D.${SDATE}.ncf

# AQFILE (required if GAMAQ_YN is set to Y)
setenv AQFILE $INPDIR/W126_$GDNAM3D.ncf

# LDFILE
setenv LDFILE $INPDIR/LDF_$GDNAM3D.2019b.$scen.ncf

# Output
setenv MGNERS $MGNINT/MGNERS.$GDNAM3D.$scen.${SDATE}.ncf

########################################################################
## Run MEGAN
if ( $RUN_MEGAN == 'Y' ) then
   rm -f $MGNERS
   $EXE | tee $LOGDIR/log.run.$PROG.$GDNAM3D.$scen.$SDATE.txt
endif

@ JD++
end  # End while JD

end # scen
end # dom

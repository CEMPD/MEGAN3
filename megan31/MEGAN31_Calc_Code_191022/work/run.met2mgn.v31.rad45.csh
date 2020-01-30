#!/bin/csh
#
# MET2MGN v3
# --
#
#
# TPAR2IOAPI v2.03a 
# --added 26-category landuse capability for mm5camx (number of landuse categories defined by NLU) 
# --added capability for LATLON and UTM projections
# --added capability for MCIP v3.3 input (2m temperatures)
# --bug in PAR processing subroutine fixed where first few hours in GMT produced zero PAR
# --added code to fill missing par data (if valid data exists for the hours surrounding it)
#
# TPAR2IOAPI v2.0
# --added capability for MM5 or MCIP input
# 
#
#        RGRND/PAR options:
#           setenv MM5RAD  Y   Solar radiation obtained from MM5
#           OR 
#           setenv MCIPRAD Y   Solar radiation obtained from MCIP
#                  --MEGAN will internally calculate PAR for each of these options and user needs to  
#                    specify `setenv PAR_INPUT N' in the MEGAN runfile
#           OR
#           setenv SATPAR Y (satellite-derived PAR from UMD GCIP/SRB files)
#                  --user needs to specify `setenv PAR_INPUT Y' in the MEGAN runfile
#
#        TEMP options:
#           setenv CAMXTEMP Y         2m temperature, calculated from mm5camx output files
#           OR
#           setenv MM5MET  Y         2m temperature, calculated from MM5 output files
#                                     Note: 2m temperature is calculated since the P-X/ACM PBL
#                                     MM5 configuration (most commonly used LSM/PBL scheme for AQ 
#                                     modeling purposes) does not produce 2m temperatures.
#           OR
#           setenv MCIPMET Y         temperature obtained from MCIP
#              -setenv TMCIP  TEMP2   2m temperature, use for MCIP v3.3 or newer
#              -setenv TMCIP  TEMP1P5 1.5m temperature, use for MCIP v3.2 or older
#
#        TZONE   time zone for input mm5CAMx files 
#        NLAY    number of layers contained in input mm5CAMx files 
#        NLU     number of landuse categories contained in CAMx landuse file 
#

############################################################

# Setting up episode   IMPORTANT: you cannot run a period a few days at a time.  You must run the entire period at once.  If additional days are required, extend the period and run the entire thing again (this is due to the PFILE issue, see: ../Input/MGNMET/QA_steps/README.txt)
############################################################
foreach dom (36 12)
set sJDATE = (2013145)
set eJDATE = (2013198)
foreach ii  (1)
set STJD = $sJDATE[$ii]
set EDJD = $eJDATE[$ii]

setenv EPISODE_SDATE $sJDATE[$ii]
setenv EPISODE_STIME 000000    

#set for grid
setenv GRIDDESC GRIDDESC
setenv GDNAM3D tceq_${dom}km


# Setting up directories and common environment variable
############################################################
source ./setcase.csh

setenv PROG met2mgn
setenv EXE /disk8/MEGAN3_1/source_code/MEGAN3_1/src/MET2MGN_rad45/$PROG.rad45
#setenv EXE /disk8/MEGAN3_1/source_code/MEGAN3_1/src/MET2MGN/$PROG


set logdir = logdir/$PROG
if ( ! -e $logdir) mkdir -p $logdir

set INPPATH     = /acme5/disk43/aqrp_bvoc/camx/inputs/mcip_v4.2/aqrp_bvoc_${dom}km
set OUTPATH     = $MGNINP/MGNMET
if (! -e $OUTPATH) mkdir $OUTPATH

setenv PFILE $OUTPATH/PFILE
rm -fv $PFILE

# Looping
############################################################
set JDATE = $STJD
while ($JDATE <= $EDJD)

if ($JDATE == 2007366) set JDATE = 2008001
if ($JDATE == 2008367) set JDATE = 2009001
if ($JDATE == 2011366) set JDATE = 2012001
if ($JDATE == 2012367) set JDATE = 2013001
if ($JDATE == 2013366) set JDATE = 2014001
@ jdy  = $JDATE - 2000000
set Y4 = `yj2ymd $JDATE | awk '{print $1}'`
set Y2 = `echo $Y4 | cut -c 3-4`
set MM = `yj2ymd $JDATE | awk '{print $2}'`
set DD = `yj2ymd $JDATE | awk '{print $3}'`

@ JDATEm1 = $JDATE - 1
if ($JDATEm1 == 2007000) set JDATEm1 = 2006365
if ($JDATEm1 == 2008000) set JDATEm1 = 2007365
if ($JDATEm1 == 2009000) set JDATEm1 = 2008366
if ($JDATEm1 == 2012000) set JDATEm1 = 2011365
if ($JDATEm1 == 2013000) set JDATEm1 = 2012366
@ jdym1  = $JDATEm1 - 2000000
set Y4m1 = `yj2ymd $JDATEm1 | awk '{print $1}'`
set Y2m1 = `echo $Y4m1 | cut -c 3-4`
set MMm1 = `yj2ymd $JDATEm1 | awk '{print $2}'`
set DDm1 = `yj2ymd $JDATEm1 | awk '{print $3}'`

#set start/end dates
setenv STDATE ${jdy}00
setenv ENDATE ${jdy}24

#TEMP/PAR input choices
#
#set if using MM5 output files
setenv MM5MET N
setenv MM5RAD N
#setenv numMM5 2
#setenv MM5file1 /pete/pete5/fcorner/met/links/MMOUT_DOMAIN1_G$Y4$MM$DD
#setenv MM5file2 /pete/pete5/fcorner/met/links/MMOUT_DOMAIN1_G$Y4$MM$DD

#set if using UMD satellite PAR data
set PARDIR = $MGNINP/PAR
setenv SATPAR N
set satpar1 = "$PARDIR/$Y2m1${MMm1}par.h"
set satpar2 = "$PARDIR/$Y2${MM}par.h"

if ($satpar1 == $satpar2) then
  setenv numSATPAR 1
  setenv SATPARFILE1 $satpar2
else
  setenv numSATPAR 2
  setenv SATPARFILE1 $satpar1
  setenv SATPARFILE2 $satpar2
endif

#set if using MCIP output files
setenv MCIPMET Y
setenv TMCIP  TEMP2          #MCIP v3.3 or newer
#setenv TMCIP  TEMP1P5       #MCIP v3.2 or older
setenv SOICRO_YN  Y      # [Y/N] set Y for MCIP4.5+ to read soil moisture/soil
                         # temperature from SOI_CRO file

setenv MCIPRAD Y

if ($JDATE == $EPISODE_SDATE) then
  setenv METCRO2Dfile1 $INPPATH/METCRO2D_aqrp_bvoc_${dom}km.$Y4$MM$DD
  if ( $SOICRO_YN == Y ) then
   setenv SOICROfile1   $INPPATH/SOI_CRO_aqrp_bvoc_${dom}km.$Y4$MM$DD
  endif
else
  setenv METCRO2Dfile1 $INPPATH/METCRO2D_aqrp_bvoc_${dom}km.$Y4m1$MMm1$DDm1
  setenv METCRO2Dfile2 $INPPATH/METCRO2D_aqrp_bvoc_${dom}km.$Y4$MM$DD
  if ( $SOICRO_YN == Y ) then
   setenv SOICROfile1   $INPPATH/SOI_CRO_aqrp_bvoc_${dom}km.$Y4m1$MMm1$DDm1
   setenv SOICROfile2   $INPPATH/SOI_CRO_aqrp_bvoc_${dom}km.$Y4$MM$DD
  endif

endif
setenv METCRO3Dfile  $INPPATH/METCRO3D_aqrp_bvoc_${dom}km.$Y4$MM$DD
setenv METDOT3Dfile  $INPPATH/METDOT3D_aqrp_bvoc_${dom}km.$Y4$MM$DD

setenv OUTFILE $OUTPATH/MET.MEGAN.$GDNAM3D.rad45.$JDATE.ncf
rm -rf $OUTFILE

$EXE |tee $logdir/log.$PROG.$GDNAM3D.rad45.$JDATE.txt 

@ JDATE++
end  # End while JDATE
end  # End foreach ii
end  # End foreach dom

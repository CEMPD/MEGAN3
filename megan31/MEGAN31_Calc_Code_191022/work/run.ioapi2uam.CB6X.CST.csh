#! /bin/csh -f

######################################################################
#
#   IOAPI2UAM converts CMAQ 1-D emissions files (I/O API) to CAMx
#   low-level emissions files (UAM-IV). Since UAM-IV format limits
#   length of species names up tp 10 characters, last 6 characters
#   (normally blanks) of CMAQ species names are truncated. Emission
#   rate is converted from mol/s (or g/s) to mol/hr (or g/hr).It also
#   shifts time-zone from GMT to user-selected local time.
#
#   INPUT ENVIRONMENTAL VARIABLES:
#
#      INFILE1      - Logical name for input file 1 (current day)
#      INFILE2      - Logical name for input file 2 (next day;
#                      required only if additional data is needed
#                      due to time zone shifting; map projection
#                      consistency won't be checked)
#      OUTFILE      - Logical name for output file
#      TZONE        - Output time-zone (8 for PST, etc.)
#      SDATE        - Output start date (YYJJJ)
#      STIME        - Output start time (HHMMSS) in TZONE
#      RLENG        - Output run length (HHMMSS; 240000 for a CAMx
#                      daily emissions input)
#
######################################################################
## Directory setups

# Program directory
setenv EXEDIR /disk43/BOEM_GOMR/camx/inputs/MEGANv2.10/bin

# Input directory
setenv INPDIR ../Output

# Output directory
setenv OUTDIR ../CAMx_ready_cst
if ( ! -e $OUTDIR ) mkdir -p $OUTDIR

# Log directory
set LOGDIR = logdir/ioapi2uam_cst
mkdir -p $LOGDIR


setenv MECHANISM    CB6X

foreach dom (36 12)
foreach scen ( J4 )
setenv GDNAM3D tceq_${dom}km


set STJD = (2013145)
set EDJD = (2013198)

foreach i (1)

set jd = $STJD[$i]
while ( $jd <= $EDJD[$i])
  if ($jd == 2005366) set jd = 2006001

  setenv RLENG   240000

  @ nd = $jd + 1
  if ($nd == 2005366) set nd = 2006001

  if ($jd != $EDJD[$i]) then
    setenv INFILE1 $INPDIR/MEGANv31.$GDNAM3D.$scen.$MECHANISM.${jd}.BDSNP.ncf
    setenv INFILE2 $INPDIR/MEGANv31.$GDNAM3D.$scen.$MECHANISM.${nd}.BDSNP.ncf
    setenv OUTFILE $OUTDIR/MEGANv31.$GDNAM3D.2019b_$scen.rad45.$MECHANISM.BDSNP.${jd}.cst.camx
    setenv RLENG   240000
    setenv TZONE   6
    setenv SDATE   $jd
    setenv STIME   0
  else
    setenv INFILE1 $INPDIR/MEGANv31.$GDNAM3D.$scen.$MECHANISM.${jd}.BDSNP.ncf
    setenv OUTFILE $OUTDIR/MEGANv31.$GDNAM3D.2019b_$scen.rad45.$MECHANISM.BDSNP.${jd}.cst.camx
    setenv RLENG   80000
    setenv TZONE   6
    setenv SDATE   $jd
    setenv STIME   0
  endif

  rm -f $OUTFILE
  $EXEDIR/ioapi2uam | tee $LOGDIR/log.run.ioapi2uam.$MECHANISM.$GDNAM3D.$scen.$jd.txt

  @ jd++

end  # End while
end  # End foreach i

end # scen
end  # End foreach dom


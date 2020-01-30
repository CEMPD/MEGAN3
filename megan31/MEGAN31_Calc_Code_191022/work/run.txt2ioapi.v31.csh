#! /bin/csh -f
########################################################################
## Common setups
source setcase.csh

foreach dom ( 36km 12km )
foreach scen ( J4 )

setenv PROMPTFLAG N
setenv PROG   txt2ioapi
setenv EXEDIR /disk8/MEGAN3_1/source_code/MEGAN3_1/src/TXT2IOAPI
setenv EXEC   $EXEDIR/$PROG
setenv GRIDDESC $cwd/GRIDDESC
setenv GDNAM3D tceq_$dom

set dd = `echo $dom | cut -c 1-2`
## File setups
#############  Inputs #####################

# Emission factor file
setenv EFSTXTF $MGNINP/MAP/grid_EF.TCEQ$dd.2019b.$scen.csv

# Canopy type file
setenv CTTXTF $MGNINP/MAP/CT3.${GDNAM3D}.csv

# Leaf area index file
setenv LAITXTF $MGNINP/MAP/LAI3.${GDNAM3D}.csv

# W126 input file (if trurning on air pollution stress)
setenv W126TXTF $MGNINP/MAP/grid_W126.$GDNAM3D.csv

# Light dependent factor file 
setenv LDFTXTF $MGNINP/MAP/grid_LDF.TCEQ$dd.2019b.$scen.csv

# Files below are needed if BDSNP soil NO algorith is used
# Climate arid file
setenv ARIDTXTF $MGNINP/MAP/grid_arid.$GDNAM3D.csv
# Climate nonarid file
setenv NONARIDTXTF $MGNINP/MAP/grid_non_arid.$GDNAM3D.csv
# Biome type file
setenv LANDTYPETXTF $MGNINP/MAP/grid_LANDTYPE.$GDNAM3D.csv
# Fertilizer input file (in kg/m3)
setenv FERTTXTF $MGNINP/MAP/grid_FERT.$GDNAM3D.csv
# Nitrogen depostiion file (total nitrogen deposition in ng/m2/s)
setenv NDEPTXTF $MGNINP/MAP/grid_NITROGEN.$GDNAM3D.csv

################# Outputs #######################
setenv EFMAPS  $MGNINP/MAP/EFMAPS31.2019b.${GDNAM3D}.$scen.ncf
setenv CANTYP  $MGNINP/MAP/CT3_${GDNAM3D}.ncf
setenv LAIS46  $MGNINP/MAP/LAI3_${GDNAM3D}.ncf
setenv W126FILE $MGNINP/MAP/W126_${GDNAM3D}.ncf
setenv LDFILE $MGNINP/MAP/LDF_${GDNAM3D}.2019b.$scen.ncf
setenv ARIDFILE $MGNINP/MAP/ARID_${GDNAM3D}.ncf
setenv NONARIDFILE $MGNINP/MAP/NONARID_${GDNAM3D}.ncf
setenv LANDTYPEFILE $MGNINP/MAP/LANDTYPE_${GDNAM3D}.ncf
setenv FERTRESFILE $MGNINP/MAP/FERT_${GDNAM3D}.ncf
setenv NDEPFILE $MGNINP/MAP/NDEP_${GDNAM3D}.ncf

## Run control
setenv RUN_EFS F       # [T|F]
setenv RUN_LAI F       # [T|F]
setenv RUN_CANTYP F    # [T|F]
setenv RUN_W126 F      # [T|F]
setenv RUN_LDF F       # [T|F]
setenv RUN_ARID F      # [T|F] 
setenv RUN_NONARID F   # [T|F]
setenv RUN_FERT T      # [T|F] fertilizer data
setenv RUN_NITROGEN F  # [T|F] total nitrogen deposition data (ng/m2/s)
setenv RUN_LANDTYPE F  # [T|F] biome type
########################################################################


## Run TXT2IOAPI
#rm -f $LAIS46
#rm -f $EFMAPS
#rm -f $CANTYP
#rm -f $W126FILE
#rm -f $LDFILE
#rm -f $ARIDFILE
#rm -f $NONARIDFILE
#rm -f $LANDTYPEFILE
#rm -f $NDEPFILE
rm -f $FERTRESFILE
if ( ! -e $MGNLOG/$PROG ) mkdir -p $MGNLOG/$PROG
if ( $RUN_EFS == T || $RUN_LDF == T ) then
$EXEC | tee $MGNLOG/$PROG/log.run.$PROG.${GDNAM3D}.$scen.txt
else
$EXEC | tee $MGNLOG/$PROG/log.run.$PROG.${GDNAM3D}.txt
endif

end
end

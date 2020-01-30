#
setenv MGNHOME $cwd:h
setenv MGNSRC $MGNHOME/src
setenv MGNLIB $MGNHOME/lib
setenv MGNEXE $MGNHOME/bin
setenv MGNRUN $MGNHOME/work
setenv MGNINP $MGNHOME/Input
setenv MGNOUT $MGNHOME/Output
setenv MGNINT $MGNHOME/Output/INT
setenv MGNLOG $MGNHOME/work/logdir


if ( ! -e $MGNINP ) then
   mkdir -p $MGNINP/MAP
   mkdir -p $MGNINP/MGNMET
   mkdir -p $MGNINP/PAR
endif
if ( ! -e $MGNINT ) mkdir -p $MGNINT
if ( ! -e $MGNLOG ) mkdir -p $MGNLOG

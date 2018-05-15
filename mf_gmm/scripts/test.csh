#!/bin/tcsh -f
echo "hi"
set x = 5
echo $x
foreach file ( `cat ../etc`)
set ctlfile = ../etc/$file.ctl 
set outfile = ../etc/trainclass.$file.ctl

#/bin/rm $outfile ; touch $outfile
#foreach x (`cat $ctlfile`)
end
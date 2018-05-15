#!/bin/csh -f

if ($#argv != 2) then
    echo "USAGE $0 phone eta"
    exit 0
endif

set phone = $1
set eta = $2

set id = /projects/speech-raj/expts/TIMIT/bushderby
set pgm = /homes/bhiksha/work/estimate_mixgau/classify_mixgau/tst

$pgm $id/ctldir/$phone.ctl $id/modelctl ./lkhds/$phone.lkhds $eta >! ./log/$phone.$eta.log

exit 0

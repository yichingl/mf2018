#!/bin/tcsh -f

## Select one of the FMM training methods
#  cmn: cepstral mean normalization
#set pgm = ../src/trainclasses/linux/trainclass.cmn.delta
set pgm = ../src/trainclasses/linux/trainclass.nocmn.nodelta
#set pgm = ../src/trainclasses/linux/trainclass.nocmn.delta

## Train a different GMM for each of the trainset-list provided
foreach file (../etc/trainclass.*)
set name = $file:t:r:e
echo $name
$pgm $file 64 ../trained/$name.nocmn.nodelta.64 2
## Format: $program-name $input-filename [num of Gaussians to use]  $output-filename [not num iter - must be multiple of number of vals]

end

echo done...

exit 0

#!/bin/sh

CMDNAME=`basename $0`
if [ $# -ne 1 ]; then
  echo "USAGE: ${CMDNAME} ctlfile"
  exit 1
fi 

CTLFILE=$1

 INDIR=./
OUTDIR=./pnc
WAVDIR=./expand

DOWNSAMPLER=/net/catalpa/usr1/kyama/project/algorithm/Downsampling/dn31
WAV2PNCC=/net/catalpa/usr1/kyama/project/algorithm/PNCC_C/wav2pncc

MSCRIPT=`mktemp tmp.XXXXXX`

 NFILTS=40
FFTSIZE=1024
  LFREQ=100
  HFREQ=6800
  SRATE=16000
 FBFILE=GTFB_${NFILTS}_${FFTSIZE}_${LFREQ}_${HFREQ}_${SRATE}.bin

if [ ! -f ${FBFILE} ]; then
  echo "ComputeGammatoneFilterBankCoeff(${NFILTS}, ${FFTSIZE}, ${LFREQ}, ${HFREQ}, ${SRATE}, '${FBFILE}');" >${MSCRIPT}

  matlab -nodisplay -nodesktop -nojvm -nosplash <${MSCRIPT} >/dev/null 2>&1
fi

: <<'#__COMMENT_OUT__'
for file in `cat ${CTLFILE}`
do

  INFILE=${INDIR}/${file}.renamed
  WAVFILE=${WAVDIR}/${file}.sph
  OUTFILE=${OUTDIR}/${file}.${LFREQ}-${HFREQ}_${NFILTS}filts.pnc

  if [ ! -f ${WAVFILE} ]; then
    sox -t wav ${INFILE} -t sph -r ${SRATE} ${WAVFILE}
  fi

  ${WAV2PNCC} -i ${WAVFILE} -o ${OUTFILE} -f ${FBFILE}

done
#__COMMENT_OUT__

# rm ${MSCRIPT}

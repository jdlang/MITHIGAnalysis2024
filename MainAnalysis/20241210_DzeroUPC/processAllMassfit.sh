#!/bin/bash
# To use:
# bash processAllSamples.sh
# optional flag: -d <config_dir> , uses different config directory from default
# optional flag: -c , runs clean.sh before processing configs

MASSFIT_DIR="pt2-5_massfitSettings"
MASSFIT_VERSION="_useGammaNForNgammaForFitFunc"
MASSFIT_LIST=(
  "fullAnalysis$MASSFIT_VERSION.json"
  "systDalpha$MASSFIT_VERSION.json"
  "systDchi2cl$MASSFIT_VERSION.json"
  "systDsvpv$MASSFIT_VERSION.json"
  "systDtrkPt$MASSFIT_VERSION.json"
  "systFitComb$MASSFIT_VERSION.json"
  "systFitPkBg$MASSFIT_VERSION.json"
  "systFitSigAlpha$MASSFIT_VERSION.json"
  "systFitSigMean$MASSFIT_VERSION.json"
  "systMassWindow$MASSFIT_VERSION.json"
  "systRapGapLoose$MASSFIT_VERSION.json"
  "systRapGapTight$MASSFIT_VERSION.json"
)
DO_CLEAN=0

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -d)
      MASSFIT_DIR="$2"
      shift 2
      ;;
    -c|--clean)
      DO_CLEAN=1
      shift
      ;;
  esac
done
if [[ "$DO_CLEAN" -eq "1" ]]; then
  source clean.sh
  wait
else
  make
  wait
fi

# Process configs
echo ""
echo "Config directory: $MASSFIT_DIR"
echo ""
for MASSFIT_JSON in ${MASSFIT_LIST[@]}; do
  bash massfit.sh $MASSFIT_DIR/$MASSFIT_JSON
  wait
done
wait

bash plot.sh "pt2-5_plotSettings/fullAnalysis.json"
root -l -b -q plotCompareDataMCmass.cpp
root -l -b -q plotMassfitSignalStudy.cpp

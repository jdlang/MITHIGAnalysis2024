#!/bin/bash
# To use:
# bash processAllSamples.sh
# optional flag: -d <config_dir> , uses different config directory from default
# optional flag: -c , runs clean.sh before processing configs

SAMPLE_DIR="pt2-5_sampleSettings"
SAMPLE_VERSION="_skimV4_mcReweighting"
SAMPLE_LIST=(
  "fullAnalysis$SAMPLE_VERSION.json"
  "systDalpha$SAMPLE_VERSION.json"
  "systDchi2cl$SAMPLE_VERSION.json"
  "systDsvpv$SAMPLE_VERSION.json"
  "systDtrkPt$SAMPLE_VERSION.json"
  "systRapGapLoose$SAMPLE_VERSION.json"
  "systRapGapTight$SAMPLE_VERSION.json"
)
MASSFIT_DIR="pt2-5_fitSettings"
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
DO_SAMPLES=1
DO_MASSFIT=1
DO_PLOTS=1

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -s)
      DO_SAMPLES=1
      DO_MASSFIT=0
      DO_PLOTS=0
      shift
      ;;
    -m)
      DO_SAMPLES=0
      DO_MASSFIT=1
      DO_PLOTS=0
      shift
      ;;
    -p)
      DO_SAMPLES=0
      DO_MASSFIT=0
      DO_PLOTS=1
      shift
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

# Process sample configs
if [[ "$DO_SAMPLES" -eq "1" ]]; then
  echo ""
  echo "Config directory: $SAMPLE_DIR"
  echo ""
  for SAMPLE_JSON in ${SAMPLE_LIST[@]}; do
      echo "Processing: $SAMPLE_DIR/$SAMPLE_JSON"
      bash makeMicroTree.sh $SAMPLE_DIR/$SAMPLE_JSON
      wait
  done
  wait
else
  echo "Skipping sample processing."
fi

# Process massfit configs
if [[ "$DO_MASSFIT" -eq "1" ]]; then
  echo ""
  echo "MassFit config directory: $MASSFIT_DIR"
  echo ""
  for MASSFIT_JSON in ${MASSFIT_LIST[@]}; do
    echo "Processing: $MASSFIT_DIR/$MASSFIT_JSON"
    bash massfit.sh $MASSFIT_DIR/$MASSFIT_JSON
    wait
  done
  wait
else
  echo "Skipping massfit processing."
fi

if [[ "$DO_PLOTS" -eq "1" ]]; then
  bash plot.sh "pt2-5_plotSettings/fullAnalysis.json"
  root -l -b -q plotCompareDataMCmass.cpp
  root -l -b -q plotMassfitSignalStudy.cpp
else
  echo "Skipping plotting."
fi

echo "All done!"

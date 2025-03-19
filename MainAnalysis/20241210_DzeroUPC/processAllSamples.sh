#!/bin/bash
# To use:
# bash processAllSamples.sh
# optional flag: -d <config_dir> , uses different config directory from default
# optional flag: -c , runs clean.sh before processing configs

SAMPLE_DIR="pt2-5_sampleSettings"
SAMPLE_VERSION="_skimV4"
SAMPLE_LIST=(
  "fullAnalysis$SAMPLE_VERSION.json"
  "systDalpha$SAMPLE_VERSION.json"
  "systDchi2cl$SAMPLE_VERSION.json"
  "systDsvpv$SAMPLE_VERSION.json"
  "systDtrkPt$SAMPLE_VERSION.json"
  "systRapGapLoose$SAMPLE_VERSION.json"
  "systRapGapTight$SAMPLE_VERSION.json"
)
DO_CLEAN=0

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -d)
      SAMPLE_DIR="$2"
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
echo "Config directory: $SAMPLE_DIR"
echo ""
for SAMPLE_JSON in ${SAMPLELIST[@]}; do
    bash makeMicroTree.sh $SAMPLE_DIR/$SAMPLE_JSON
    wait
done
wait

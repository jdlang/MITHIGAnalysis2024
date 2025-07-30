#!/bin/bash

FILEPATH=${1}
COUNTER=${2}
OUTPUT=${3}
ISDATA=${4}

mkdir -p "${OUTPUT}/temp_inputs/"
FILE="${OUTPUT}/temp_inputs/job_${COUNTER}.root"
rm $FILE &> /dev/null
cp $FILEPATH $FILE
wait

echo "Processing $FILE"
./Execute --Input "$FILE" \
   --Output "${OUTPUT}/output_${COUNTER}.root" \
   --Year 2025 \
   --TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
   --DoGenLevel false \
   --IsData true \
   --IsPP true \
   --IsBackground false \
   --DoAlternateTrackSelection true \
   --AlternateTrackSelection 0 \
   --CheckZ true
wait

sleep 0.2
rm $FILE
wait

#!/bin/bash
DATE=$(date +%Y%m%d)

# Source of forest files (use xrootd path)
SOURCE_SERVER="root://eoscms.cern.ch/"
SOURCE_DIR="/store/group/phys_heavyions/jdlang/run3_2023Data_Jan2024ReReco/"
# Output of skimmed files (use xrootd path)
OUTPUT_SERVER="root://eoscms.cern.ch/"
OUTPUT_DIR="/store/group/phys_heavyions/$USER/run3_2023Data_Jan2024ReReco_Skims_$DATE"

# Job settings (memory and storage are in GB)
FILES_PER_JOB=500
JOB_MEMORY=5
JOB_STORAGE=20
CMSSW_VERSION="CMSSW_13_2_13"

# Local directory for condor configs
CONFIG_DIR="condorSkimConfigs_${DATE}"
MASTER_FILE_LIST="${CONFIG_DIR}/forestFilesForSkim.txt"

# Skim Config Settings
# (the arguments given to ./Execute in RunParallel.sh)
CONFIG_Year=2023
CONFIG_TriggerRejection=2
CONFIG_EventRejection=true
CONFIG_ZDCGapRejection=true
CONFIG_DRejection=or
CONFIG_ZDCMinus=1000
CONFIG_ZDCPlus=1100
CONFIG_IsData=true
CONFIG_PFTree=particleFlowAnalyser/pftree

# Include VOMS proxy in process
REFRESH_PROXY=0

# For testing/debugging (set to 0 to run all):
MAX_JOBS=0

$ProjectBase/SampleGeneration/20241203_CondorSkimUtils/InitCondorSkim.sh $SOURCE_SERVER $SOURCE_DIR $OUTPUT_SERVER $OUTPUT_DIR $FILES_PER_JOB $JOB_MEMORY $JOB_STORAGE $CMSSW_VERSION $CONFIG_DIR $MASTER_FILE_LIST $CONFIG_Year $CONFIG_TriggerRejection $CONFIG_EventRejection $CONFIG_ZDCGapRejection $CONFIG_DRejection $CONFIG_ZDCMinus $CONFIG_ZDCPlus $CONFIG_IsData $CONFIG_PFTree $REFRESH_PROXY $MAX_JOBS

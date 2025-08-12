#!/bin/bash

PD_LIST=(
  0
)

for PD_NUM in ${PD_LIST[@]}; do
  ls /eos/cms/store/group/phys_heavyions/jdlang/Run3_OxygenRAA/PromptForest_ZDiMu/IonPhysics${PD_NUM}/crab_NeNe_IonPhysics*_ZDiMu/25*/000*/*.root > filelist_2025_NeNe_PromptReco_ZDiMu_IonPhysics$PD_NUM.txt
  wait
  
  sed -i "s|/eos/cms||" filelist_2025_NeNe_PromptReco_ZDiMu_IonPhysics$PD_NUM.txt
  wait

  bash RunParallelData.sh $PD_NUM > log_RunParallelData_NeNe_IonPhysics${PD_NUM}.txt
  wait
done

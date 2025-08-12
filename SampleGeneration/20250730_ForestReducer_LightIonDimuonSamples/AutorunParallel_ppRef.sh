#!/bin/bash

PD_LIST=(
  3
)

for PD_NUM in ${PD_LIST[@]}; do
  ls /eos/cms/store/group/phys_heavyions/xirong/Run3_OxygenRAA/PromptForest_ZDiMu/PPRefSingleMuon${PD_NUM}/crab_pp_PPRefSingleMuon${PD_NUM}_ZDiMu/25*/000*/*.root > filelist_2024_ppRef_PromptReco_ZDiMu_SingleMuon$PD_NUM.txt
  wait
  
  sed -i "s|/eos/cms||" filelist_2024_ppRef_PromptReco_ZDiMu_SingleMuon$PD_NUM.txt
  wait

  bash RunParallelData_ppRef.sh $PD_NUM > log_RunParallelData_ppRef_SingleMuon${PD_NUM}.txt
  wait
done

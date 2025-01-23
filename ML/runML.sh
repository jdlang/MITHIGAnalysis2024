#!/bin/bash
python3 MLoptimization.py --input_file_mc /data00/UPCD0LowPtAnalysis_2023ZDCORData_2023reco/SkimsMC/20241216_v1_filelist20241216_Pthat2_ForceD0Decay100M_BeamA_v1/mergedfile.root --tree_name Tree --output_csv output_pt_p1p2_y_m1p1.cvs --output_root output_pt_p1p2_y_m1p1.root --output_model XGBoostMConly_pt_p1p2_y_m1p1.json --ptmin 1 --ptmax 2 --ymin -2 --ymax 2

#python3 MLoptimization.py --input_file_mc /data00/UPCD0LowPtAnalysis_2023ZDCORData_2023reco/SkimsMC/20241216_v1_filelist20241216_Pthat2_ForceD0Decay100M_BeamA_v1/mergedfile.root --tree_name Tree --output_csv output_pt_p1p2_y_m1p1.cvs --output_root output_pt_p1p2_y_m1p1.root --output_model XGBoostMConly_pt_p1p2_y_m1p1.json --ptmin 1 --ptmax 2 --ymin -2 --ymax 2

#python3 MLoptimization.py --input_file_mc /data00/UPCD0LowPtAnalysis_2023ZDCORData_2023reco/SkimsMC/20241216_v1_filelist20241216_Pthat2_ForceD0Decay100M_BeamA_v1/mergedfile.root --tree_name Tree --output_csv output_pt_p1p2_y_m1p1.cvs --output_root output_pt_p1p2_y_m1p1.root --output_model XGBoostMConly_pt_p1p2_y_m1p1.json --ptmin 1 --ptmax 2 --ymin -2 --ymax -1
#
#python3 MLoptimization.py --input_file_mc /data00/UPCD0LowPtAnalysis_2023ZDCORData_2023reco/SkimsMC/20241216_v1_filelist20241216_Pthat2_ForceD0Decay100M_BeamA_v1/mergedfile.root --tree_name Tree --output_csv output_pt_p1p2_y_m1p1.cvs --output_root output_pt_p1p2_y_m1p1.root --output_model XGBoostMConly_pt_p1p2_y_m1p1.json --ptmin 1 --ptmax 2 --ymin -1 --ymax 0
#
#python3 MLoptimization.py --input_file_mc /data00/UPCD0LowPtAnalysis_2023ZDCORData_2023reco/SkimsMC/20241216_v1_filelist20241216_Pthat2_ForceD0Decay100M_BeamA_v1/mergedfile.root --tree_name Tree --output_csv output_pt_p1p2_y_m1p1.cvs --output_root output_pt_p1p2_y_m1p1.root --output_model XGBoostMConly_pt_p1p2_y_m1p1.json --ptmin 1 --ptmax 2 --ymin 0 --ymax 1
#
#python3 MLoptimization.py --input_file_mc /data00/UPCD0LowPtAnalysis_2023ZDCORData_2023reco/SkimsMC/20241216_v1_filelist20241216_Pthat2_ForceD0Decay100M_BeamA_v1/mergedfile.root --tree_name Tree --output_csv output_pt_p1p2_y_m1p1.cvs --output_root output_pt_p1p2_y_m1p1.root --output_model XGBoostMConly_pt_p1p2_y_m1p1.json --ptmin 1 --ptmax 2 --ymin 1 --ymax 2

# ptmin: 1.0, ptmax: 2.0, ymin: -2.0, ymax: -1.0, ROC AUC: 0.7150
# ptmin: 1.0, ptmax: 2.0, ymin: -1.0, ymax:  0.0, ROC AUC: 0.8458
# ptmin: 1.0, ptmax: 2.0, ymin:  0.0, ymax:  1.0, ROC AUC: 0.8140
# ptmin: 1.0, ptmax: 2.0, ymin:  1.0, ymax:  2.0, ROC AUC: 0.8148

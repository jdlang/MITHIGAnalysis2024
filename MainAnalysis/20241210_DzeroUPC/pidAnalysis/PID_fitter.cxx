/* To compile, run this:
g++ PID_fitter.cpp -o ExecFitPID `root-config --cflags --glibs`
*/

#include <filesystem>
#include <iostream>
#include <cmath>

#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"

#include "PID_utils.h"

using namespace std;

/*
  MakePlot_logdedx_trkP()
  Acts as a standard dedx vs trkP 2D hist builder.
*/
TH2D* MakePlot_logdedx_trkP(
  TString plotName,
  TString plotTitle  = "dEdx vs. track p; track p (GeV); log(dEdx)"
  int   logdedx_nBins = logdedx_NBINS,
  float logdedx_min   = logdedx_MIN,
  float logdedx_max   = logdedx_MAX,
  int   trkP_nBins    = trkP_NBINS,
  float trkP_min      = trkP_MIN,
  float trkP_max      = trkP_MAX
) {
  TH2D* h2D_logdedx_trkP = new TH2D(
    plotName, plotTitle,
    logdedx_nBins, logdedx_min, logdedx_max,
    trkP_nBins, trkP_min, trkP_max
  );
  h2D_logdedx_trkP->Sumw2();
  return h2D_logdedx_trkP;
}

void FillDedxVsTrkP(
  TTree* tree,
  TString foutPath,
  vector<float> DptBins,
  vector<float> DyBins,
  bool isMC,
  bool absDy = true // Dy bins represent absolute values of y
) {
  // Get branches for tree
  int Dsize;
  vector<float>* Dpt          = new vector<float>();
  vector<float>* Dy           = new vector<float>();
  vector<float>* Dtrk1P       = new vector<float>();
  vector<float>* Dtrk1Pt      = new vector<float>();
  vector<float>* Dtrk1Eta     = new vector<float>();
  vector<float>* Dtrk1dedx    = new vector<float>();
  vector<int>*   Dgentk1pdgId = new vector<int>();
  tree->SetBranchStatus("*", 0);
  CheckAndSetBranch(tree, Dsize);
  CheckAndSetBranch(tree, Dpt);
  CheckAndSetBranch(tree, Dy);
  CheckAndSetBranch(tree, Dtrk1P);
  CheckAndSetBranch(tree, Dtrk1Pt);
  CheckAndSetBranch(tree, Dtrk1Eta);
  CheckAndSetBranch(tree, Dtrk1dedx);
  CheckAndSetBranch(tree, Dgentk1pdgId);
  
  // Construct hists for Dpt, Dy bins
  vector<vector<TH2D*>> hDedxVsTrkP;
  vector<vector<TH2D*>> hDedxVsTrkP_MC_pion;
  vector<vector<TH2D*>> hDedxVsTrkP_MC_kaon;
  vector<vector<TH2D*>> hDedxVsTrkP_MC_proton;
  
  for (int iDpt = 0; iDpt < (DptBins.size() - 1); iDpt++) {
    for (int iDy = 0; iDy < (DyBins.size() - 1); iDy++) {
      hDedxVsTrkP[iDpt][iDy] = MakePlot_logdedx_trkP(
        Form("hDedxVsTrkP_DptBin%d_DyBin%d", iDpt, iDy),
        Form("%.1f< Dpt < %.1f GeV, %.1f < Dy < %.1f; track p (GeV); log(dEdx)",
          DptBins[iDpt], DptBins[iDpt + 1], DyBins[iDy], DyBins[iDy + 1]));
      if (isMC) {
        hDedxVsTrkP_MC_pion[iDpt][iDy] = MakePlot_logdedx_trkP(
          Form("hDedxVsTrkP_MC_pion_DptBin%d_DyBin%d", iDpt, iDy),
          Form("%.1f< Dpt < %.1f GeV, %.1f < Dy < %.1f; track p (GeV); log(dEdx)",
            DptBins[iDpt], DptBins[iDpt + 1], DyBins[iDy], DyBins[iDy + 1]));
        hDedxVsTrkP_MC_kaon[iDpt][iDy] = MakePlot_logdedx_trkP(
          Form("hDedxVsTrkP_MC_kaon_DptBin%d_DyBin%d", iDpt, iDy),
          Form("%.1f< Dpt < %.1f GeV, %.1f < Dy < %.1f; track p (GeV); log(dEdx)",
            DptBins[iDpt], DptBins[iDpt + 1], DyBins[iDy], DyBins[iDy + 1]));
        hDedxVsTrkP_MC_proton[iDpt][iDy] = MakePlot_logdedx_trkP(
          Form("hDedxVsTrkP_MC_proton_DptBin%d_DyBin%d", iDpt, iDy),
          Form("%.1f< Dpt < %.1f GeV, %.1f < Dy < %.1f; track p (GeV); log(dEdx)",
            DptBins[iDpt], DptBins[iDpt + 1], DyBins[iDy], DyBins[iDy + 1]));
      }
    }// end loop over iDy
  }// end loop over iDpt
  
  for (int iEvt = 0; iEvt < tree->GetEntries(); iEvt++) {
    tree->GetEntry(iEvt);
    for (int iD = 0; iD < Dsize; iD++) {
      for (int iDpt = 0; iDpt < (DptBins.size() - 1); iDpt++) {
        for (int iDy = 0; iDy < (DyBins.size() - 1); iDy++) {
          float _Dy = Dy->at(iD);
          float _Dtrk1P = Dtrk1P->at(iD);
          float _logDtrk1dedx = TMath::Log(Dtrk1dedx->at(iD));
          if (absDy) _Dy = fabs(Dy->at(iD));
          if (Dpt->at(iD) < DptBins[iDpt]     ||
              Dpt->at(iD) > DptBins[iDpt + 1] ||
              _Dy < DyBins[iDy]               ||
              _Dy < DyBins[iDy + 1]
            ) continue;
          hDedxVsTrkP[iDpt][iDy]->Fill(_Dtrk1P, _logDtrk1dedx);
          if (isMC) {
            int _DtrkPdgId = Dgentk1pdgId->at(iD);
            if (_DtrkPdgId == pionID)
              hDedxVsTrkP_MC_pion[iDpt][iDy]->Fill(_Dtrk1P, _logDtrk1dedx);
            else if (_DtrkPdgId == kaonID)
              hDedxVsTrkP_MC_kaon[iDpt][iDy]->Fill(_Dtrk1P, _logDtrk1dedx);
            else if (_DtrkPdgId == protonID)
              hDedxVsTrkP_MC_proton[iDpt][iDy]->Fill(_Dtrk1P, _logDtrk1dedx);
          }
        }// end loop over iDy
      }// end loop over iDpt
    }// end loop over iD (Dsize)
  }// end loop over iEvt (nEntries)
  
  TFile* fout = TFile::Open(foutPath, "RECREATE");
  fout->cd();
  for (int iDpt = 0; iDpt < (DptBins.size() - 1); iDpt++) {
    for (int iDy = 0; iDy < (DyBins.size() - 1); iDy++) {
      hDedxVsTrkP[iDpt][iDy]->Write();
      if (isMC) {
        hDedxVsTrkP_MC_pion[iDpt][iDy]->Write();
        hDedxVsTrkP_MC_kaon[iDpt][iDy]->Write();
        hDedxVsTrkP_MC_proton[iDpt][iDy]->Write();
      }
    }// end loop over iDy
  }// end loop over iDpt
  fout->Close();
  fin->Close();
}

vector<TF1*> FitSingleParticleCurve(
  TH2D* hDedxVsTrkP_MC
) {
  int nTrackPBins = hDedxVsTrkP_MC->GetNbinsX();
  TString nameBase = hDedxVsTrkP_MC->GetName();
  TGraph* gDedxVsTrkP = new TGraphErrors(nTrackPBins);
  TGraph* gDedxVsTrkP_errHi = new TGraphErrors(nTrackPBins);
  TGraph* gDedxVsTrkP_errLo = new TGraphErrors(nTrackPBins);
  gDedxVsTrkP->SetName(nameBase + "_graph");
  for (int iBin = 1; iBin <= nTrackPBins; iBin++) {
    TH1D* hDedx = hDedxVsTrkP_MC->ProjectionY("hDedx", iBin, iBin);
    float binX = hDedxVsTrkP_MC->GetBinCenter(iBin);
    float binMean = hDedx->GetMean();
    float binErr = hDedx->GetStdDev();
    gDedxVsTrkP->SetPoint(iBin, binX, binMean);
    gDedxVsTrkP_errHi->SetPoint(iBin, binX, binMean + binErr);
    gDedxVsTrkP_errLo->SetPoint(iBin, binX, binMean - binErr);
    delete hDedx;
  }
  TF1* fDedxVsTrkP = new TF1(
    nameBase + "_fit", fSimpleBetheBloch, trkP_MIN, trkP_MAX, 4);
  fDedxVsTrkP->SetParameters(0.02, 0.1, 5., 1.);
  TF1* fDedxVsTrkP_errHi = new TF1(
    nameBase + "_fitErrHi", fSimpleBetheBloch, trkP_MIN, trkP_MAX, 4);
  fDedxVsTrkP_errHi->SetParameters(0.02, 0.1, 5., 1.);
  TF1* fDedxVsTrkP_errLo = new TF1(
    nameBase + "_fitErrLo", fSimpleBetheBloch, trkP_MIN, trkP_MAX, 4);
  fDedxVsTrkP_errLo->SetParameters(0.02, 0.1, 5., 1.);
  
  gDedxVsTrkP->Fit(fDedxVsTrkP);
  gDedxVsTrkP_errHi->Fit(fDedxVsTrkP_errHi);
  gDedxVsTrkP_errLo->Fit(fDedxVsTrkP_errLo);
  
  TH1D* hTrkP = hDedxVsTrkP_MC->ProjectionX("hTrkP");
  
  TF1* fTrkP = new TF1(nameBase + "_TrkP", "", trkP_MIN, trkP_MAX, 4);
  
  vector<TF1*> fDedxVsTrkP_array = {
    fDedxVsTrkP, fDedxVsTrkP_errHi, fDedxVsTrkP_errLo};
  return fDedxVsTrkP_array;
}

int main(
  int argc,
  char* argv[]
) {
  
  
  // Process tree
  
  // Plot dEdx per particle for MC -> set initial params for BetheBloch
  // Plot dEdx spread
  // Fit functions to dEdx spread
}

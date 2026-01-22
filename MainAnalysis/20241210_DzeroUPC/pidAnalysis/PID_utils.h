/* To compile, run this:
g++ fitPID.cpp -o ExecFitPID `root-config --cflags --glibs`
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

// Preprocessor constants
#define logdedx_MIN   0.
#define logdedx_MAX   10.
#define logdedx_NBINS 100
#define trkP_MIN      0.
#define trkP_MAX      5.
#define trkP_NBINS    100

#define pionId        211
#define kaonId        321
#define protonId      2212

#define DMAX          20000

// Preprocessor functions
#define CheckAndSetBranch(TREE, BRANCH)     \
if(TREE->GetBranch(#BRANCH))                \
{                                           \
  TREE->SetBranchStatus(#BRANCH, 1);        \
  TREE->SetBranchAddress(#BRANCH, &BRANCH); \
}

/*
  fSimpleBetheBloch()
  Acts as a function of <p> and returns value of <dedx>.
  Has 4 params.
*/
double fSimpleBetheBloch(
  double* p,
  double* par
) {
  double beta2 = pow(p[0] / par[2], 2);
  double dedx =
    -1 * par[0] * (1 / beta2) *
    (TMath::Log(2 / par[1] * beta2 / (1 - beta2)) - (2 * TMath::E()) )
    + par[3];
  return dedx;
}

double fTripleGauss(
  double* dedx,
  double* par
) {
  double gausPion   = par[0] * TMath::Exp(pow((dedx[0] - par[1])/(par[2]), 2));
  double gausKaon   = par[3] * TMath::Exp(pow((dedx[0] - par[4])/(par[5]), 2));
  double gausProton = par[6] * TMath::Exp(pow((dedx[0] - par[7])/(par[8]), 2));
  // Ensure reasonable ordering of means
  if (par[1] > par[4]) par[4] += fabs(par[1] - par[4]);
  if (par[4] > par[7]) par[7] += fabs(par[4] - par[7]);
  // Ensure reasonable ordering of peak size
  if (par[0] > par[3]) par[3] += fabs(par[0] - par[3]);
  if (par[3] > par[6]) par[6] += fabs(par[3] - par[6]);
}

void SkimDedxFromForestDfinder(
  TString finPath,
  TString foutPath,
  TString treeName
) {
  TFile* fin = TFile::Open(finPath, "READ");
  TTree* tree = (TTree*) fin->Get(treeName);
  tree->SetBranchStatus("*", 0);
  
  long nEntries = tree->GetEntries();
  int    Dsize;
  float* Dpt          = new float[DMAX];
  float* Dy           = new float[DMAX];
  float* Dtrk1P       = new float[DMAX];
  float* Dtrk1Pt      = new float[DMAX];
  float* Dtrk1Eta     = new float[DMAX];
  float* Dtrk1dedx    = new float[DMAX];
  int*   Dgentk1pdgId = new int[DMAX];
  float* Dtrk2P       = new float[DMAX];
  float* Dtrk2Pt      = new float[DMAX];
  float* Dtrk2Eta     = new float[DMAX];
  float* Dtrk2dedx    = new float[DMAX];
  int*   Dgentk2pdgId = new int[DMAX];
  CheckAndSetBranch(tree, Dsize);
  CheckAndSetBranch(tree, Dpt);
  CheckAndSetBranch(tree, Dy);
  CheckAndSetBranch(tree, Dtrk1P);
  CheckAndSetBranch(tree, Dtrk1Pt);
  CheckAndSetBranch(tree, Dtrk1Eta);
  CheckAndSetBranch(tree, Dtrk1dedx);
  CheckAndSetBranch(tree, Dgentk1pdgId);
  CheckAndSetBranch(tree, Dtrk2P);
  CheckAndSetBranch(tree, Dtrk2Pt);
  CheckAndSetBranch(tree, Dtrk2Eta
  CheckAndSetBranch(tree, Dtrk2dedx);
  CheckAndSetBranch(tree, Dgentk2pdgId);
  
  TFile* fout = TFile::Open(foutPath, "RECREATE");
  TTree* treeout = new TTree("tree", "Dfinder skim tree");
  vector<float>* vDpt          = new vector<float>();
  vector<float>* vDy           = new vector<float>();
  vector<float>* vDtrk1P       = new vector<float>();
  vector<float>* vDtrk1Pt      = new vector<float>();
  vector<float>* vDtrk1Eta     = new vector<float>();
  vector<float>* vDtrk1dedx    = new vector<float>();
  vector<int>*   vDgentk1pdgId = new vector<int>();
  vector<float>* vDtrk2P       = new vector<float>();
  vector<float>* vDtrk2Pt      = new vector<float>();
  vector<float>* vDtrk2Eta     = new vector<float>();
  vector<float>* vDtrk2dedx    = new vector<float>();
  vector<int>*   vDgentk2pdgId = new vector<int>();
  treeout->Branch("Dsize",        &Dsize, "Dsize/I");
  treeout->Branch("Dpt",          &vDpt);
  treeout->Branch("Dy",           &vDy);
  treeout->Branch("Dtrk1P",       &vDtrk1P);
  treeout->Branch("Dtrk1Pt",      &vDtrk1Pt);
  treeout->Branch("Dtrk1Eta",     &vDtrk1Eta);
  treeout->Branch("Dtrk1dedx",    &vDtrk1dedx);
  treeout->Branch("Dgentk1pdgId", &vDgentk1pdgId);
  treeout->Branch("Dtrk2P",       &vDtrk2P);
  treeout->Branch("Dtrk2Pt",      &vDtrk2Pt);
  treeout->Branch("Dtrk2Eta",     &vDtrk2Eta);
  treeout->Branch("Dtrk2dedx",    &vDtrk2dedx);
  treeout->Branch("Dgentk2pdgId", &vDgentk2pdgId);
  
  for (long int i = 0; i < nEntries; i++) {
    tree->GetEntry(i);
    for (int j = 0; j < Dsize; j++) {
      vDpt->push_back(vDpt->at(j));
      vDy->push_back(vDy->at(j));
      vDtrk1P->push_back(vDtrk1P->at(j));
      vDtrk1Pt->push_back(vDtrk1Pt->at(j));
      vDtrk1Eta->push_back(vDtrk1Eta->at(j));
      vDtrk1dedx->push_back(vDtrk1dedx->at(j));
      vDgentk1pdgId->push_back(vDgentk1pdgId->at(j));
      vDtrk2P->push_back(vDtrk2P->at(j));
      vDtrk2Pt->push_back(vDtrk2Pt->at(j));
      vDtrk2Eta->push_back(vDtrk2Eta->at(j));
      vDtrk2dedx->push_back(vDtrk2dedx->at(j));
      vDgentk2pdgId->push_back(vDgentk2pdgId->at(j));
    }
    treeout->Fill();
  }
  
  fout->cd();
  treeout->Write();
  fout->Close();
  fin->Close();
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

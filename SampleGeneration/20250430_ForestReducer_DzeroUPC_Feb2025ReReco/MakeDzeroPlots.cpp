#include <vector>
#include <filesystem>

#include "TFile.h"

using namespace std;

vector<int> validRuns2023 = {
  374804, 374810, 374828, 374833, 374834,
  374925, 374950, 374951, 374953, 374961,
  374970, 374997, 374998, 375001, 375002,
  375007, 375013, 375055, 375058, 375060,
  375061, 375064, 375110, 375145, 375164,
  375195, 375202, 375245, 375252, 375256,
  375259, 375300, 375317, 375391, 375413,
  375415, 375440, 375441, 375448, 375455,
  375463, 375465, 375483, 375491, 375507,
  375513, 375530, 375531, 375545, 375549,
  375658, 375659, 375665, 375666, 375695,
  375696, 375697, 375703, 375738, 375739,
  375740, 375744, 375746
};

bool checkIfValidRun(
  vector<int> validRuns,
  int run
) {
  bool pass = false;
  for (int i = 0; i < validRuns.size(); i++) {
    if (run == validRuns[i]) {
      pass = true;
      break;
    }
  }
  return pass;
}

void PlotComparisons (
  TH1F* hBase,
  TH1F* hJan24,
  TH1F* hFeb25,
  TString plotDir,
  bool doLogy = false,
  float hmin = 0.,
  float hmax = 0.
) {
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->SetMargin(0.15, 0.05, 0.10, 0.10);
  canvas->cd();
  
  if (hJan24 != nullptr) {
    hJan24->SetMarkerStyle(25);
    hJan24->SetMarkerSize(0.9);
    hJan24->SetMarkerColor(kAzure+9);
    hJan24->SetLineColor(kAzure+9);
    hJan24->SetLineWidth(2);
  }
  
  hFeb25->SetMarkerStyle(20);
  hFeb25->SetMarkerSize(1.0);
  hFeb25->SetMarkerColor(kViolet-1);
  hFeb25->SetLineColor(kViolet-1);
  hFeb25->SetLineWidth(2);
  
  if ( hmin == 0. && hmax == 0. ) {
    hmin = 0.;
    if (hJan24 != nullptr)
      hmax = TMath::Max(hJan24->GetMaximum(), hFeb25->GetMaximum());
    else hmax = hFeb25->GetMaximum();
  }
  if (!doLogy) {
    hBase->SetMinimum(hmin);
    hBase->SetMaximum(hmax*1.25);
  }
  else {
    hBase->SetMinimum(0.1);
    hBase->SetMaximum(hmax*10);
  }
  
  TLegend* legend = new TLegend(0.65, 0.75, 0.95, 0.9);
  legend->SetTextSize(0.03);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  if (hJan24 != nullptr) legend->AddEntry(hJan24, "Jan24 ReReco", "lp");
  legend->AddEntry(hFeb25, "Feb25 ReReco", "lp");
  
  hBase->Draw();
  if (hJan24 != nullptr) hJan24->Draw("same");
  hFeb25->Draw("same");
  legend->Draw();
  gStyle->SetOptStat(0);
  if (doLogy) gPad->SetLogy(1);
  
  canvas->Update();
  canvas->SaveAs(
    plotDir + "/ReRecoComparison_2023Data_Jan24_Feb25_" + hBase->GetName() + ".pdf"
  );
  
  if (hJan24 != nullptr) {
    canvas->Clear();
    canvas->cd();
    
    TH1F* hBaseNorm = (TH1F*) hBase->Clone(Form("%s_norm", hBase->GetName()));
    TH1F* hJan24Norm = (TH1F*) hJan24->Clone(Form("%s_norm", hJan24->GetName()));
    TH1F* hFeb25Norm = (TH1F*) hFeb25->Clone(Form("%s_norm", hFeb25->GetName()));
    
    hJan24Norm->Scale(1 / hJan24Norm->Integral());
    hFeb25Norm->Scale(1 / hFeb25Norm->Integral());
    
    if (!doLogy) {
      hBaseNorm->SetMinimum(0.);
      hBaseNorm->SetMaximum(
        1.25 * TMath::Max(hJan24Norm->GetMaximum(), hFeb25Norm->GetMaximum())
      );
    }
    else {
      hBaseNorm->SetMinimum(0.0001);
      hBaseNorm->SetMaximum(
        10 * TMath::Max(hJan24Norm->GetMaximum(), hFeb25Norm->GetMaximum())
      );
    }
    
    hBaseNorm->Draw();
    hJan24Norm->Draw("same");
    hFeb25Norm->Draw("same");
    legend->Draw();
    gStyle->SetOptStat(0);
    if (doLogy) gPad->SetLogy(1);
    
    canvas->Update();
    canvas->SaveAs(
      plotDir + "/ReRecoComparison_2023Data_Jan24_Feb25_" + hBase->GetName() + "_norm.pdf"
    );
  }
  delete canvas;
}

void MakeDzeroPlots (
  TString filepathJan24 = "/data00/jdlang/UPCD0LowPtAnalysis/SkimsData/20250508_Skim_2023Data_Jan2024ReReco_HIForward0.root",
  TString filepathFeb25 = "/data00/jdlang/UPCD0LowPtAnalysis/SkimsData/20250512_Skim_2023Data_Feb2025ReReco_HIForward0.root",
  TString plotDir = "skimComparisonPlots",
  float DptMin = 2.,
  float DptMax = 5.,
  float DyMin = -2.,
  float DyMax = 2.,
  bool selectgammaN = true,
  bool selectNgamma = true
) {
  cout << "Opening source files..." << endl;
  TFile* finJan24 = TFile::Open(filepathJan24, "READ");
  TFile* finFeb25 = TFile::Open(filepathFeb25, "READ");
  
  plotDir = plotDir + "_Dpt" + DptMin + "-" + DptMax + "_Dy" + DyMin + "-" + DyMax;
  if ( selectgammaN ) plotDir = plotDir + "_gammaN";
  if ( selectNgamma ) plotDir = plotDir + "_Ngamma";
  system(Form("mkdir -p %s", plotDir.Data()));
  TFile* fPlots = TFile::Open(
    plotDir+"/skimPlots.root", "RECREATE");
  
  TTree* treeJan24 = (TTree*) finJan24->Get("Tree");
  treeJan24->SetName("TreeJan24");
  TTree* treeFeb25 = (TTree*) finFeb25->Get("Tree");
  treeFeb25->SetName("TreeFeb25");
  
  cout << "Making histograms..." << endl;
  
  TH1F* hnVtx = new TH1F(
    "hnVtx", "nVtx; N_{Vertex}; Count",
    3, 1., 4.);
  TH1F* hnTrack = new TH1F(
    "hnTrack", "nTrack; N_{Tracks}; Count",
    10, 2., 12.);
  TH1F* hDsize = new TH1F(
    "hDsize", "Dsize; N_{D^{0} Candidates}; Count",
    10, 1., 11.);
  TH1F* hnVtxDSel = new TH1F(
    "hnVtxDSel", "nVtx (after D-selection); N_{Vertex}; Count",
    3, 1., 4.);
  TH1F* hnTrackDSel = new TH1F(
    "hnTrackDSel", "nTrack (after D-selection); N_{Tracks}; Count",
    10, 2., 12.);
  TH1F* hDsizeDSel = new TH1F(
    "hDsizeDSel", "Dsize (after D-selection); N_{D^{0} Candidates}; Count",
    10, 1., 11.);
  TH1F* hDsizeSelected = new TH1F(
    "hDsizeSelected", "Selected D's per event; N_{D^{0} Candidates}; Count",
    10, 1., 11.);
  TH1F* hDpt = new TH1F(
    "hDpt", "Dpt; D^{0} p_{T}; Count",
    24, 0., 12.);
  TH1F* hDy = new TH1F(
    "hDy", "Dy; D^{0} y; Count",
    16, -2., 2.);
  TH1F* hDmass = new TH1F(
    "hDmass", "Dmass; M_{K#pi}; Count",
    32, 1.46, 2.26);
  TH1F* hDalpha = new TH1F(
    "hDalpha", "Dalpha; D^{0} #alpha; Count",
    20, 0., 0.4);
  TH1F* hDchi2cl = new TH1F(
    "hDchi2cl", "Dchi2cl; D^{0} #chi^{2}; Count",
    20, 0., 1.0);
  TH1F* hDsvpvSig = new TH1F(
    "hDsvpvSig", "DsvpvSig; (Dsvpv / DsvpvErr); Count",
    20, 2., 12.);
  TH1F* hDtrk1Pt = new TH1F(
    "hDtrk1Pt", "Dtrk1Pt; D^{0}_{Track 1} p_{T}; Count",
    24, 0., 12.);
  TH1F* hDtrk2Pt = new TH1F(
    "hDtrk2Pt", "Dtrk2Pt; D^{0}_{Track 2} p_{T}; Count",
    24, 0., 12.);
  TH1F* hDtrk1Eta = new TH1F(
    "hDtrk1Eta", "Dtrk1Eta; D^{0}_{Track 1} #eta; Count",
    12, -3, 3.);
  TH1F* hDtrk2Eta = new TH1F(
    "hDtrk2Eta", "Dtrk2Eta; D^{0}_{Track 2} #eta; Count",
    12, -3, 3.);
  TH1F* hDtrk1dedx = new TH1F(
    "hDtrk1dedx", "Dtrk1dedx; D^{0}_{Track 1} Log(dE/dx); Count",
    80, 0., 2.);
  TH1F* hDtrk2dedx = new TH1F(
    "hDtrk2dedx", "Dtrk2Eta; D^{0}_{Track 2} Log(dE/dx); Count",
    80, 0., 2.);
  TH1F* hDtrk1MassHypo = new TH1F(
    "hDtrk1MassHypo", "Dtrk1MassHypo; D^{0}_{Track 1} Mass Hypothesis (GeV); Count",
    20, 0., 2.0);
  TH1F* hDtrk2MassHypo = new TH1F(
    "hDtrk2MassHypo", "Dtrk2MassHypo; D^{0}_{Track 2} Mass Hypothesis (GeV); Count",
    20, 0., 2.0);
    
  TH1F* hnVtx_Jan24 = (TH1F*) hnVtx->Clone("hnVtx_Jan24");
  TH1F* hnTrack_Jan24 = (TH1F*) hnTrack->Clone("hnTrack_Jan24");
  TH1F* hDsize_Jan24 = (TH1F*) hDsize->Clone("hDsize_Jan24");
  TH1F* hnVtxDSel_Jan24 = (TH1F*) hnVtxDSel->Clone("hnVtxDSel_Jan24");
  TH1F* hnTrackDSel_Jan24 = (TH1F*) hnTrackDSel->Clone("hnTrackDSel_Jan24");
  TH1F* hDsizeDSel_Jan24 = (TH1F*) hDsizeDSel->Clone("hDsizeDSel_Jan24");
  TH1F* hDsizeSelected_Jan24 = (TH1F*) hDsizeSelected->Clone("hDsizeSelected_Jan24");
  TH1F* hDpt_Jan24 = (TH1F*) hDpt->Clone("hDpt_Jan24");
  TH1F* hDy_Jan24 = (TH1F*) hDy->Clone("hDy_Jan24");
  TH1F* hDmass_Jan24 = (TH1F*) hDmass->Clone("hDmass_Jan24");
  TH1F* hDalpha_Jan24 = (TH1F*) hDalpha->Clone("hDalpha_Jan24");
  TH1F* hDchi2cl_Jan24 = (TH1F*) hDchi2cl->Clone("hDchi2cl_Jan24");
  TH1F* hDsvpvSig_Jan24 = (TH1F*) hDsvpvSig->Clone("hDsvpvSig_Jan24");
  TH1F* hDtrk1Pt_Jan24 = (TH1F*) hDtrk1Pt->Clone("hDtrk1Pt_Jan24");
  TH1F* hDtrk2Pt_Jan24 = (TH1F*) hDtrk2Pt->Clone("hDtrk2Pt_Jan24");
  hnVtx_Jan24->Sumw2();
  hnTrack_Jan24->Sumw2();
  hDsize_Jan24->Sumw2();
  hnVtxDSel_Jan24->Sumw2();
  hnTrackDSel_Jan24->Sumw2();
  hDsizeDSel_Jan24->Sumw2();
  hDsizeSelected_Jan24->Sumw2();
  hDpt_Jan24->Sumw2();
  hDy_Jan24->Sumw2();
  hDmass_Jan24->Sumw2();
  hDalpha_Jan24->Sumw2();
  hDchi2cl_Jan24->Sumw2();
  hDsvpvSig_Jan24->Sumw2();
  hDtrk1Pt_Jan24->Sumw2();
  hDtrk2Pt_Jan24->Sumw2();
  
  TH1F* hnVtx_Feb25 = (TH1F*) hnVtx->Clone("hnVtx_Feb25");
  TH1F* hnTrack_Feb25 = (TH1F*) hnTrack->Clone("hnTrack_Feb25");
  TH1F* hDsize_Feb25 = (TH1F*) hDsize->Clone("hDsize_Feb25");
  TH1F* hnVtxDSel_Feb25 = (TH1F*) hnVtxDSel->Clone("hnVtxDSel_Feb25");
  TH1F* hnTrackDSel_Feb25 = (TH1F*) hnTrackDSel->Clone("hnTrackDSel_Feb25");
  TH1F* hDsizeDSel_Feb25 = (TH1F*) hDsizeDSel->Clone("hDsizeDSel_Feb25");
  TH1F* hDsizeSelected_Feb25 = (TH1F*) hDsizeSelected->Clone("hDsizeSelected_Feb25");
  TH1F* hDpt_Feb25 = (TH1F*) hDpt->Clone("hDpt_Feb25");
  TH1F* hDy_Feb25 = (TH1F*) hDy->Clone("hDy_Feb25");
  TH1F* hDmass_Feb25 = (TH1F*) hDmass->Clone("hDmass_Feb25");
  TH1F* hDalpha_Feb25 = (TH1F*) hDalpha->Clone("hDalpha_Feb25");
  TH1F* hDchi2cl_Feb25 = (TH1F*) hDchi2cl->Clone("hDchi2cl_Feb25");
  TH1F* hDsvpvSig_Feb25 = (TH1F*) hDsvpvSig->Clone("hDsvpvSig_Feb25");
  TH1F* hDtrk1Pt_Feb25 = (TH1F*) hDtrk1Pt->Clone("hDtrk1Pt_Feb25");
  TH1F* hDtrk2Pt_Feb25 = (TH1F*) hDtrk2Pt->Clone("hDtrk2Pt_Feb25");
  TH1F* hDtrk1Eta_Feb25 = (TH1F*) hDtrk1Eta->Clone("hDtrk1Eta_Feb25");
  TH1F* hDtrk2Eta_Feb25 = (TH1F*) hDtrk2Eta->Clone("hDtrk2Eta_Feb25");
  TH1F* hDtrk1dedx_Feb25 = (TH1F*) hDtrk1dedx->Clone("hDtrk1dedx_Feb25");
  TH1F* hDtrk2dedx_Feb25 = (TH1F*) hDtrk2dedx->Clone("hDtrk2dedx_Feb25");
  TH1F* hDtrk1MassHypo_Feb25 = (TH1F*) hDtrk1MassHypo->Clone("hDtrk1MassHypo_Feb25");
  TH1F* hDtrk2MassHypo_Feb25 = (TH1F*) hDtrk2MassHypo->Clone("hDtrk2MassHypo_Feb25");
  hnVtx_Feb25->Sumw2();
  hnTrack_Feb25->Sumw2();
  hDsize_Feb25->Sumw2();
  hnVtxDSel_Feb25->Sumw2();
  hnTrackDSel_Feb25->Sumw2();
  hDsizeDSel_Feb25->Sumw2();
  hDsizeSelected_Feb25->Sumw2();
  hDpt_Feb25->Sumw2();
  hDy_Feb25->Sumw2();
  hDmass_Feb25->Sumw2();
  hDalpha_Feb25->Sumw2();
  hDchi2cl_Feb25->Sumw2();
  hDsvpvSig_Feb25->Sumw2();
  hDtrk1Pt_Feb25->Sumw2();
  hDtrk2Pt_Feb25->Sumw2();
  hDtrk1Eta_Feb25->Sumw2();
  hDtrk2Eta_Feb25->Sumw2();
  hDtrk1dedx_Feb25->Sumw2();
  hDtrk2dedx_Feb25->Sumw2();
  hDtrk1MassHypo_Feb25->Sumw2();
  hDtrk2MassHypo_Feb25->Sumw2();
  
  cout << "Setting Jan2024 branches..." << endl;
  int Run_24;
  int nVtx_24;
  int nTrack_24;
  int Dsize_24;
  bool selectedBkgFilter_24;
  bool selectedVtxFilter_24;
  bool isL1ZDCOr_24;
  bool isL1ZDCXORJet8_24;
  bool gapgammaN_24;
  bool gapNgamma_24;
  bool ZDCgammaN_24;
  bool ZDCNgamma_24;
  vector<bool>* DpassCut23PAS_24 = nullptr;
  vector<bool>* gammaN_24 = nullptr;
  vector<bool>* Ngamma_24 = nullptr;
  vector<float>* Dpt_24 = nullptr;
  vector<float>* Dy_24 = nullptr;
  vector<float>* Dmass_24 = nullptr;
  vector<float>* Dalpha_24 = nullptr;
  vector<float>* Dchi2cl_24 = nullptr;
  vector<float>* DsvpvDistance_24 = nullptr;
  vector<float>* DsvpvDisErr_24 = nullptr;
  vector<float>* Dtrk1Pt_24 = nullptr;
  vector<float>* Dtrk2Pt_24 = nullptr;
  treeJan24->SetBranchAddress("Run", &Run_24);
  treeJan24->SetBranchAddress("nVtx", &nVtx_24);
  treeJan24->SetBranchAddress("nTrackInAcceptanceHP", &nTrack_24);
  treeJan24->SetBranchAddress("Dsize", &Dsize_24);
  treeJan24->SetBranchAddress("selectedBkgFilter", &selectedBkgFilter_24);
  treeJan24->SetBranchAddress("selectedVtxFilter", &selectedVtxFilter_24);
  treeJan24->SetBranchAddress("isL1ZDCXORJet8", &isL1ZDCXORJet8_24);
  treeJan24->SetBranchAddress("isL1ZDCOr", &isL1ZDCOr_24);
  treeJan24->SetBranchAddress("gapgammaN", &gapgammaN_24);
  treeJan24->SetBranchAddress("gapNgamma", &gapNgamma_24);
  treeJan24->SetBranchAddress("ZDCgammaN", &ZDCgammaN_24);
  treeJan24->SetBranchAddress("ZDCNgamma", &ZDCNgamma_24);
  treeJan24->SetBranchAddress("DpassCut23PAS", &DpassCut23PAS_24);
  treeJan24->SetBranchAddress("gammaN", &gammaN_24);
  treeJan24->SetBranchAddress("Ngamma", &Ngamma_24);
  treeJan24->SetBranchAddress("Dpt", &Dpt_24);
  treeJan24->SetBranchAddress("Dy", &Dy_24);
  treeJan24->SetBranchAddress("Dmass", &Dmass_24);
  treeJan24->SetBranchAddress("Dalpha", &Dalpha_24);
  treeJan24->SetBranchAddress("Dchi2cl", &Dchi2cl_24);
  treeJan24->SetBranchAddress("DsvpvDistance", &DsvpvDistance_24);
  treeJan24->SetBranchAddress("DsvpvDisErr", &DsvpvDisErr_24);
  treeJan24->SetBranchAddress("Dtrk1Pt", &Dtrk1Pt_24);
  treeJan24->SetBranchAddress("Dtrk2Pt", &Dtrk2Pt_24);
  
  cout << "Filling Jan2024 histograms..." << endl;
  int nEntries_24 = treeJan24->GetEntries();
  int eventSel_24 = 0;
  int DSel_24 = 0;
  int DSel_Dsize3_24 = 0;
  int DSel_Dindex3_24 = 0;
  for (int iE = 0; iE < nEntries_24; iE++) {
    treeJan24->GetEntry(iE);
    if ( !isL1ZDCOr_24 ) continue;
    if ( !selectedBkgFilter_24 || !selectedVtxFilter_24 ) continue;
    if ( selectgammaN && selectNgamma ) {
      if (!(ZDCgammaN_24 && gapgammaN_24) && !(ZDCNgamma_24 && gapNgamma_24))
        continue;
    }
    else if ( selectgammaN && !(ZDCgammaN_24 && gapgammaN_24) ) continue;
    else if ( selectNgamma && !(ZDCNgamma_24 && gapNgamma_24) ) continue;
    if ( nVtx_24 >= 3 ) continue;
    if ( !checkIfValidRun(validRuns2023, Run_24) ) continue;
    eventSel_24++;
    hnVtx_Jan24->Fill(nVtx_24);
    hnTrack_Jan24->Fill(nTrack_24);
    hDsize_Jan24->Fill(Dsize_24);
    bool passDSel_24 = false;
    vector<bool> selectedDs_24(Dsize_24, false);
    for (int iD = 0; iD < Dsize_24; iD++) {
      if ( DpassCut23PAS_24->at(iD) == false ) continue;
      if ( Dpt_24->at(iD) < DptMin || Dpt_24->at(iD) > DptMax) continue;
      if ( Dy_24->at(iD) < DyMin || Dy_24->at(iD) > DyMax) continue;
      passDSel_24 = true;
      selectedDs_24[iD] = true;
    }
    if (!passDSel_24) continue;
    hnVtxDSel_Jan24->Fill(nVtx_24);
    hnTrackDSel_Jan24->Fill(nTrack_24);
    hDsizeDSel_Jan24->Fill(Dsize_24);
    int nDSelected_24 = 0;
    for (int iD = 0; iD < Dsize_24; iD++) {
      if (!selectedDs_24[iD]) continue;
      DSel_24++;
      nDSelected_24++;
      hDpt_Jan24->Fill(Dpt_24->at(iD));
      hDy_Jan24->Fill(Dy_24->at(iD));
      hDmass_Jan24->Fill(Dmass_24->at(iD));
      hDalpha_Jan24->Fill(Dalpha_24->at(iD));
      hDchi2cl_Jan24->Fill(Dchi2cl_24->at(iD));
      hDsvpvSig_Jan24->Fill(DsvpvDistance_24->at(iD)/DsvpvDisErr_24->at(iD));
      hDtrk1Pt_Jan24->Fill(Dtrk1Pt_24->at(iD));
      hDtrk2Pt_Jan24->Fill(Dtrk2Pt_24->at(iD));
      if ( Dsize_24 < 3 ) DSel_Dsize3_24++;
      if ( iD < 3 ) DSel_Dindex3_24++;
    }
    hDsizeSelected_Jan24->Fill(nDSelected_24);
  }
  
  cout << "Setting Feb2025 branches..." << endl;
  int Run_25;
  int nVtx_25;
  int nTrack_25;
  int Dsize_25;
  bool selectedBkgFilter_25;
  bool selectedVtxFilter_25;
  bool isL1ZDCOr_25;
  bool isL1ZDCXORJet8_25;
  bool gapgammaN_25;
  bool gapNgamma_25;
  bool ZDCgammaN_25;
  bool ZDCNgamma_25;
  vector<bool>* DpassCut23PAS_25 = nullptr;
  vector<bool>* gammaN_25 = nullptr;
  vector<bool>* Ngamma_25 = nullptr;
  vector<float>* Dpt_25 = nullptr;
  vector<float>* Dy_25 = nullptr;
  vector<float>* Dmass_25 = nullptr;
  vector<float>* Dalpha_25 = nullptr;
  vector<float>* Dchi2cl_25 = nullptr;
  vector<float>* DsvpvDistance_25 = nullptr;
  vector<float>* DsvpvDisErr_25 = nullptr;
  vector<float>* Dtrk1Pt_25 = nullptr;
  vector<float>* Dtrk2Pt_25 = nullptr;
  vector<float>* Dtrk1Eta_25 = nullptr;
  vector<float>* Dtrk2Eta_25 = nullptr;
  vector<float>* Dtrk1dedx_25 = nullptr;
  vector<float>* Dtrk2dedx_25 = nullptr;
  vector<float>* Dtrk1MassHypo_25 = nullptr;
  vector<float>* Dtrk2MassHypo_25 = nullptr;
  treeFeb25->SetBranchAddress("Run", &Run_25);
  treeFeb25->SetBranchAddress("nVtx", &nVtx_25);
  treeFeb25->SetBranchAddress("nTrackInAcceptanceHP", &nTrack_25);
  treeFeb25->SetBranchAddress("Dsize", &Dsize_25);
  treeFeb25->SetBranchAddress("selectedBkgFilter", &selectedBkgFilter_25);
  treeFeb25->SetBranchAddress("selectedVtxFilter", &selectedVtxFilter_25);
  treeFeb25->SetBranchAddress("isL1ZDCXORJet8", &isL1ZDCXORJet8_25);
  treeFeb25->SetBranchAddress("isL1ZDCOr", &isL1ZDCOr_25);
  treeFeb25->SetBranchAddress("gapgammaN", &gapgammaN_25);
  treeFeb25->SetBranchAddress("gapNgamma", &gapNgamma_25);
  treeFeb25->SetBranchAddress("ZDCgammaN", &ZDCgammaN_25);
  treeFeb25->SetBranchAddress("ZDCNgamma", &ZDCNgamma_25);
  treeFeb25->SetBranchAddress("DpassCut23PAS", &DpassCut23PAS_25);
  treeFeb25->SetBranchAddress("gammaN", &gammaN_25);
  treeFeb25->SetBranchAddress("Ngamma", &Ngamma_25);
  treeFeb25->SetBranchAddress("Dpt", &Dpt_25);
  treeFeb25->SetBranchAddress("Dy", &Dy_25);
  treeFeb25->SetBranchAddress("Dmass", &Dmass_25);
  treeFeb25->SetBranchAddress("Dalpha", &Dalpha_25);
  treeFeb25->SetBranchAddress("Dchi2cl", &Dchi2cl_25);
  treeFeb25->SetBranchAddress("DsvpvDistance", &DsvpvDistance_25);
  treeFeb25->SetBranchAddress("DsvpvDisErr", &DsvpvDisErr_25);
  treeFeb25->SetBranchAddress("Dtrk1Pt", &Dtrk1Pt_25);
  treeFeb25->SetBranchAddress("Dtrk2Pt", &Dtrk2Pt_25);
  treeFeb25->SetBranchAddress("Dtrk1Eta", &Dtrk1Eta_25);
  treeFeb25->SetBranchAddress("Dtrk2Eta", &Dtrk2Eta_25);
  treeFeb25->SetBranchAddress("Dtrk1dedx", &Dtrk1dedx_25);
  treeFeb25->SetBranchAddress("Dtrk2dedx", &Dtrk2dedx_25);
  treeFeb25->SetBranchAddress("Dtrk1MassHypo", &Dtrk1MassHypo_25);
  treeFeb25->SetBranchAddress("Dtrk2MassHypo", &Dtrk2MassHypo_25);
  
  cout << "Filling Feb2025 histograms..." << endl;
  int nEntries_25 = treeFeb25->GetEntries();
  int eventSel_25 = 0;
  int DSel_25 = 0;
  int DSel_Dsize3_25 = 0;
  int DSel_Dindex3_25 = 0;
  for (int iE = 0; iE < nEntries_25; iE++) {
    treeFeb25->GetEntry(iE);
    if ( !isL1ZDCOr_25) continue;
    if ( !selectedBkgFilter_25 || !selectedVtxFilter_25 ) continue;
    if ( selectgammaN && selectNgamma ) {
      if (!(ZDCgammaN_25 && gapgammaN_25) && !(ZDCNgamma_25 && gapNgamma_25))
        continue;
    }
    else if ( selectgammaN && !(ZDCgammaN_25 && gapgammaN_25) ) continue;
    else if ( selectNgamma && !(ZDCNgamma_25 && gapNgamma_25) ) continue;
    if ( nVtx_25 >= 3 ) continue;
    if ( !checkIfValidRun(validRuns2023, Run_25) ) continue;
    eventSel_25++;
    hnVtx_Feb25->Fill(nVtx_25);
    hnTrack_Feb25->Fill(nTrack_25);
    hDsize_Feb25->Fill(Dsize_25);
    bool passDSel_25 = false;
    vector<bool> selectedDs_25(Dsize_25, false);
    for (int iD = 0; iD < Dsize_25; iD++) {
      if ( DpassCut23PAS_25->at(iD) == false ) continue;
      if ( Dpt_25->at(iD) < DptMin || Dpt_25->at(iD) > DptMax) continue;
      if ( Dy_25->at(iD) < DyMin || Dy_25->at(iD) > DyMax) continue;
      passDSel_25 = true;
      selectedDs_25[iD] = true;
    }
    if (!passDSel_25) continue;
    hnVtxDSel_Feb25->Fill(nVtx_25);
    hnTrackDSel_Feb25->Fill(nTrack_25);
    hDsizeDSel_Feb25->Fill(Dsize_25);
    int nDSelected_25 = 0;
    for (int iD = 0; iD < Dsize_25; iD++) {
      if (!selectedDs_25[iD]) continue;
      DSel_25++;
      nDSelected_25++;
      hDpt_Feb25->Fill(Dpt_25->at(iD));
      hDy_Feb25->Fill(Dy_25->at(iD));
      hDmass_Feb25->Fill(Dmass_25->at(iD));
      hDalpha_Feb25->Fill(Dalpha_25->at(iD));
      hDchi2cl_Feb25->Fill(Dchi2cl_25->at(iD));
      hDsvpvSig_Feb25->Fill(DsvpvDistance_25->at(iD)/DsvpvDisErr_25->at(iD));
      hDtrk1Pt_Feb25->Fill(Dtrk1Pt_25->at(iD));
      hDtrk2Pt_Feb25->Fill(Dtrk2Pt_25->at(iD));
      hDtrk1Eta_Feb25->Fill(Dtrk1Eta_25->at(iD));
      hDtrk2Eta_Feb25->Fill(Dtrk2Eta_25->at(iD));
      hDtrk1dedx_Feb25->Fill(TMath::Log(Dtrk1dedx_25->at(iD)));
      hDtrk2dedx_Feb25->Fill(TMath::Log(Dtrk2dedx_25->at(iD)));
      hDtrk1MassHypo_Feb25->Fill(Dtrk1MassHypo_25->at(iD));
      hDtrk2MassHypo_Feb25->Fill(Dtrk2MassHypo_25->at(iD));
      if ( Dsize_25 < 3 ) DSel_Dsize3_25++;
      if ( iD < 3 ) DSel_Dindex3_25++;
    }
    hDsizeSelected_Feb25->Fill(nDSelected_25);
  }
  
  cout << "Jan24 Counts ----------" << endl;
  cout << "  Total Entries:   " << nEntries_24 << endl;
  cout << "  Event Selection: " << eventSel_24 <<
    Form(" -> %.5f", 1.0 * eventSel_24 / nEntries_24) <<
    " (evtSel/total)" << endl;
  cout << "  D Selection:     " << DSel_24 <<
    Form(" -> %.5f", 1.0 * DSel_24 / eventSel_24) << " (DSel/evtSel)" << endl;
  cout << "  Dindex < 3:      " << DSel_Dindex3_24 << endl;
  cout << "  Dsize < 3:       " << DSel_Dsize3_24 << endl;
  cout << "" << endl;
  
  cout << "Feb25 Counts ----------" << endl;
  cout << "  Total Entries:   " << nEntries_25 << endl;
  cout << "  Event Selection: " << eventSel_25 <<
    Form(" -> %.5f", 1.0 * eventSel_25 / nEntries_25) <<
    " (evtSel/total)" << endl;
  cout << "  D Selection:     " << DSel_25 <<
    Form(" -> %.5f", 1.0 * DSel_25 / eventSel_25) << " (DSel/evtSel)" << endl;
  cout << "  Dindex < 3:      " << DSel_Dindex3_25 << endl;
  cout << "  Dsize < 3:       " << DSel_Dsize3_25 << endl;
  cout << "" << endl;
  
  TCanvas* cTemp = new TCanvas();
  TString cutString = "DpassCut23PAS && Dpt > 2 && Dpt < 12 && selectedBkgFilter && selectedVtxFilter && ( (gammaN && gapgammaN && ZDCgammaN) || (Ngamma && gapNgamma && ZDCNgamma) )";
  
  cout << "Making plots..." << endl;
  PlotComparisons (hnVtx, hnVtx_Jan24, hnVtx_Feb25, plotDir, true);
  PlotComparisons (hnTrack, hnTrack_Jan24, hnTrack_Feb25, plotDir);
  PlotComparisons (hDsize, hDsize_Jan24, hDsize_Feb25, plotDir);
  PlotComparisons (hnVtxDSel, hnVtxDSel_Jan24, hnVtxDSel_Feb25, plotDir, true);
  PlotComparisons (hnTrackDSel, hnTrackDSel_Jan24, hnTrackDSel_Feb25, plotDir);
  PlotComparisons (hDsizeDSel, hDsizeDSel_Jan24, hDsizeDSel_Feb25, plotDir);
  PlotComparisons (hDsizeSelected, hDsizeSelected_Jan24, hDsizeSelected_Feb25, plotDir);
  PlotComparisons (hDpt, hDpt_Jan24, hDpt_Feb25, plotDir, true);
  PlotComparisons (hDy, hDy_Jan24, hDy_Feb25, plotDir);
  PlotComparisons (hDmass, hDmass_Jan24, hDmass_Feb25, plotDir);
  PlotComparisons (hDalpha, hDalpha_Jan24, hDalpha_Feb25, plotDir);
  PlotComparisons (hDchi2cl, hDchi2cl_Jan24, hDchi2cl_Feb25, plotDir);
  PlotComparisons (hDsvpvSig, hDsvpvSig_Jan24, hDsvpvSig_Feb25, plotDir);
  PlotComparisons (hDtrk1Pt, hDtrk1Pt_Jan24, hDtrk1Pt_Feb25, plotDir, true);
  PlotComparisons (hDtrk2Pt, hDtrk2Pt_Jan24, hDtrk2Pt_Feb25, plotDir, true);
  PlotComparisons (hDtrk1Eta, nullptr, hDtrk1Eta_Feb25, plotDir);
  PlotComparisons (hDtrk2Eta, nullptr, hDtrk2Eta_Feb25, plotDir);
  PlotComparisons (hDtrk1dedx, nullptr, hDtrk1dedx_Feb25, plotDir);
  PlotComparisons (hDtrk2dedx, nullptr, hDtrk2dedx_Feb25, plotDir);
  PlotComparisons (hDtrk1MassHypo, nullptr, hDtrk1MassHypo_Feb25, plotDir);
  PlotComparisons (hDtrk2MassHypo, nullptr, hDtrk2MassHypo_Feb25, plotDir);
  
  finJan24->Close();
  finFeb25->Close();
  
  cout << "Saving plots to root file..." << endl;
  fPlots->cd();
  hnVtx_Jan24->Write();
  hnVtx_Feb25->Write();
  hnVtxDSel_Jan24->Write();
  hnVtxDSel_Feb25->Write();
  hnTrack_Jan24->Write();
  hnTrack_Feb25->Write();
  hnTrackDSel_Jan24->Write();
  hnTrackDSel_Feb25->Write();
  hDsize_Jan24->Write();
  hDsize_Feb25->Write();
  hDsizeDSel_Jan24->Write();
  hDsizeDSel_Feb25->Write();
  hDsizeSelected_Jan24->Write();
  hDsizeSelected_Feb25->Write();
  hDpt_Jan24->Write();
  hDpt_Feb25->Write();
  hDy_Jan24->Write();
  hDy_Feb25->Write();
  hDmass_Jan24->Write();
  hDmass_Feb25->Write();
  hDalpha_Jan24->Write();
  hDalpha_Feb25->Write();
  hDchi2cl_Jan24->Write();
  hDchi2cl_Feb25->Write();
  hDsvpvSig_Jan24->Write();
  hDsvpvSig_Feb25->Write();
  hDtrk1Pt_Jan24->Write();
  hDtrk1Pt_Feb25->Write();
  hDtrk2Pt_Jan24->Write();
  hDtrk2Pt_Feb25->Write();
  hDtrk1Eta_Feb25->Write();
  hDtrk2Eta_Feb25->Write();
  hDtrk1dedx_Feb25->Write();
  hDtrk2dedx_Feb25->Write();
  hDtrk1MassHypo_Feb25->Write();
  hDtrk2MassHypo_Feb25->Write();
  fPlots->Close();
  
  cout << "Done!" << endl;
}

// run with:
// root -l -b -q plotSystematicsEval.cpp

#include <TTree.h>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TChain.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>

#include <RooAddPdf.h>
#include <RooAbsPdf.h>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooExponential.h>
#include <RooChebychev.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooWorkspace.h>
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooFormulaVar.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <regex>
#include <filesystem>

#include "include/structs.h"
#include "plotCrossSection.h"

using namespace std;
using namespace RooFit;

#define DMASS 1.86484
#define DMASSMIN 1.66
#define DMASSMAX 2.16
#define DMASSNBINS 40

const float DtrkPtCut = 0.5;
const float DsvpvCut = 2.5;
const float DalphaCut[4] = {0.2, 0.4, 0.4, 0.2};
const float Dchi2clCut = 0.1;

const bool RebuildHists = true;

const int nPtBins = 1;
double ptBins[2] = {2, 5};
const int nYBins = 4;
double yBins[5] = {-2, -1, 0, 1, 2};

void plot_yBin(
  TH1D* data_gN,
  TH1D* data_Ng,
  TH1D* MC_gN,
  TH1D* MC_Ng,
  string plotTitle,
  string plotXaxis,
  string plotYaxis,
  string fileName,
  float histMin = 0.,
  float histMax = 100.,
  float ratioMin = 0.0,
  float ratioMax = 2.0,
  bool useLogy = false
) {
  bool useMC = true;
  if (MC_gN == nullptr || MC_Ng == nullptr) useMC = false;
  
  cout << "Making hist template..." << endl;
  TH1D* histTemplate = new TH1D(
    "histTemplate",
    Form("%s;;%s", plotTitle.c_str(), plotYaxis.c_str()),
    data_gN->GetNbinsX(),
    data_gN->GetBinLowEdge(1),
    data_gN->GetBinLowEdge(data_gN->GetNbinsX() + 1)
  );
  histTemplate->SetMinimum(histMin);
  histTemplate->SetMaximum(histMax);
  
  cout << "Making ratio template..." << endl;
  TH1D* ratioTemplate = new TH1D(
    "ratioTemplate",
    Form(";%s;%s", plotXaxis.c_str(), plotYaxis.c_str()),
    data_gN->GetNbinsX(),
    data_gN->GetBinLowEdge(1),
    data_gN->GetBinLowEdge(data_gN->GetNbinsX() + 1)
  );
  ratioTemplate->SetMinimum(ratioMin);
  ratioTemplate->SetMaximum(ratioMax);
  
  cout << "Making data ratio hist..." << endl;
  TH1D* data_ratio = (TH1D*) data_gN->Clone("data_ratio");
  data_ratio->SetTitle(ratioTemplate->GetTitle());
  data_ratio->Divide(data_Ng);
  
  data_gN->SetLineWidth(2);
  data_gN->SetLineColor(kAzure+2);
  data_gN->SetMarkerColor(kAzure+2);
  data_gN->SetMarkerStyle(20);
  
  data_Ng->SetLineWidth(2);
  data_Ng->SetLineColor(kPink-8);
  data_Ng->SetMarkerColor(kPink-8);
  data_Ng->SetMarkerStyle(20);
  
  data_ratio->SetLineWidth(2);
  data_ratio->SetLineColor(kTeal+2);
  data_ratio->SetMarkerColor(kTeal+2);
  data_ratio->SetMarkerStyle(33);
  data_ratio->SetMarkerSize(1.4);
  
  TH1D* MC_ratio;
  if (useMC) {
    cout << "Making MC ratio hist..." << endl;
    MC_ratio = (TH1D*) MC_gN->Clone("MC_ratio");
    MC_ratio->SetTitle(ratioTemplate->GetTitle());
    MC_ratio->Divide(MC_Ng);
    
    MC_gN->SetLineWidth(2);
    MC_gN->SetLineColor(kAzure-9);
    MC_gN->SetMarkerColor(kAzure-9);
    MC_gN->SetMarkerStyle(24);
    
    MC_Ng->SetLineWidth(2);
    MC_Ng->SetLineColor(kPink+1);
    MC_Ng->SetMarkerColor(kPink+1);
    MC_Ng->SetMarkerStyle(24);
    
    MC_ratio->SetLineWidth(2);
    MC_ratio->SetLineColor(kTeal-9);
    MC_ratio->SetMarkerColor(kTeal-9);
    MC_ratio->SetMarkerStyle(27);
    MC_ratio->SetMarkerSize(1.4);
  }
  
  cout << "Making canvas..." << endl;
  TCanvas* canvas = new TCanvas("canvas", "", 600, 900);
  TPad* padTop = new TPad("padTop", "", 0.0, 0.3, 1.0, 1.0);
  TPad* padBot = new TPad("padBot", "", 0.0, 0.0, 1.0, 0.3);
  padTop->SetMargin(0.14, 0.06, 0.00, 0.16);
  padBot->SetMargin(0.14, 0.06, 0.30, 0.00);
  padTop->Draw();
  padBot->Draw();
  
  cout << "Making top plot..." << endl;
  padTop->cd();
  if (useLogy) gPad->SetLogy(1);
  histTemplate->Draw();
  histTemplate->SetLabelSize(0.025/0.7, "Y");
  histTemplate->SetTitleSize(0.03/0.7, "Y");
  histTemplate->SetTitleOffset(1.2, "Y");
  if (useMC) MC_gN->Draw("same");
  if (useMC) MC_Ng->Draw("same");
  data_gN->Draw("same");
  data_Ng->Draw("same");
  gStyle->SetOptStat(0);
  
  cout << "Making bottom plot..." << endl;
  padBot->cd();
  gPad->SetLogy(0);
  ratioTemplate->Draw();
  ratioTemplate->SetLabelSize(0.025/0.3, "XY");
  ratioTemplate->SetTitleSize(0.03/0.3, "XY");
  ratioTemplate->SetTitleOffset(0.5, "Y");
  ratioTemplate->SetTitleOffset(0.8, "X");
  TLine* unity = new TLine(
    ratioTemplate->GetBinLowEdge(1), 1.0,
    ratioTemplate->GetBinLowEdge(ratioTemplate->GetNbinsX() + 1), 1.0);
  unity->SetLineColor(kGray);
  unity->SetLineWidth(1);
  unity->SetLineStyle(9);
  unity->Draw();
  if (useMC) MC_ratio->Draw("same");
  data_ratio->Draw("same");
  gStyle->SetOptStat(0);
  
  cout << "Adding legend..." << endl;
  canvas->cd();
  gPad->SetLogy(0);
  canvas->Update();
  TLegend* legend;
  if (useMC) legend = new TLegend(0.25, 0.74, 0.75, 0.88);
  else legend = new TLegend(0.25, 0.80, 0.75, 0.88);
  legend->SetTextSize(0.022/0.7);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(data_gN, "Data #gammaN", "lp");
  if (useMC) legend->AddEntry(MC_gN, "MC #gammaN (BeamA)", "lp");
  legend->AddEntry(data_Ng, "Data N#gamma", "lp");
  if (useMC) legend->AddEntry(MC_Ng, "MC N#gamma (BeamB)", "lp");
  legend->AddEntry(data_ratio, "Ratio of data (#gammaN / N#gamma)", "lp");
  if (useMC) legend->AddEntry(MC_ratio, "Ratio of MC (#gammaN / N#gamma)", "lp");
  legend->Draw();
  
  cout << "Saving..." << endl;
  canvas->SaveAs(fileName.c_str());
  
  cout << "Done with plot!" << endl;
  delete canvas;
  delete histTemplate;
  delete ratioTemplate;
  delete data_ratio;
  delete MC_ratio;
  delete unity;
}

void plotEvalDatasets(
) {
  TFile* fData;
  TFile* fMCA;
  TFile* fMCB;
  TFile* histFile;
  TTree* tData;
  TTree* tMCA;
  TTree* tMCB;
  
  vector<TH1D*> zVtx_gN(8, nullptr);
  vector<TH1D*> zVtx_Ng(8, nullptr);
  vector<TH1D*> nVtx_gN(8, nullptr);
  vector<TH1D*> nVtx_Ng(8, nullptr);
  vector<TH1D*> HFEmaxP_gN(8, nullptr);
  vector<TH1D*> HFEmaxP_Ng(8, nullptr);
  vector<TH1D*> HFEmaxM_gN(8, nullptr);
  vector<TH1D*> HFEmaxM_Ng(8, nullptr);
  vector<TH1D*> ZDCsumP_gN(8, nullptr);
  vector<TH1D*> ZDCsumP_Ng(8, nullptr);
  vector<TH1D*> ZDCsumM_gN(8, nullptr);
  vector<TH1D*> ZDCsumM_Ng(8, nullptr);
  
  if (filesystem::exists("plot/DataEval/hists.root")) {
    histFile = TFile::Open("plot/DataEval/hists.root", "READ");
  }
  else {
    fData = TFile::Open(
      "/data00/UPCD0LowPtAnalysis_2023ZDCORData_2023reco/SkimsData/20250312_ForestDfinderData23Skim_v4.root",
      "READ"
    );
    fMCA = TFile::Open(
      "/data00/UPCD0LowPtAnalysis_2023ZDCORData_2023reco/SkimsMC/20250324_v4_Pthat0_ForceD0Decay_BeamA/mergedfile.root",
      "READ"
    );
    fMCB = TFile::Open(
      "/data00/UPCD0LowPtAnalysis_2023ZDCORData_2023reco/SkimsMC/20250324_v4_Pthat0_ForceD0Decay_BeamB/mergedfile.root",
      "READ"
    );
    tData = (TTree*) fData->Get("Tree");
    tMCA  = (TTree*) fMCA->Get("Tree");
    tMCB  = (TTree*) fMCB->Get("Tree");
    vector<TTree*> trees    = {  tData,   tMCA,   tMCB };
    
    for (int t = 0; t < trees.size(); t++) {
      int hIndex = 0;
      if (t > 0) hIndex = nYBins; // use back half of TH1 vectors for MC
      TTree* tree = trees[t];
      string treeName;
      if (t == 0) treeName = "Data";
      else if (t == 1) treeName = "MC-A";
      else if (t == 2) treeName = "MC-B";
      for (int iy = 0; iy < nYBins; iy++) {
        if (t == 0 || (t == 1)) {
          zVtx_gN[hIndex + iy] = new TH1D(
            Form("zVtx_%s_gN_pt2-5_y%.f-%.f",
            treeName.c_str(), yBins[iy], yBins[iy+1]),
            "zVtx_gN; z-vertex; Counts",
            25, -25, 25);
          nVtx_gN[hIndex + iy] = new TH1D(
            Form("nVtx_%s_gN_pt2-5_y%.f-%.f",
            treeName.c_str(), yBins[iy], yBins[iy+1]),
            "nVtx_gN; Number of Vertices; Counts",
            5, 0, 5);
          HFEmaxP_gN[hIndex + iy] = new TH1D(
            Form("HFEmaxP_%s_gN_pt2-5_y%.f-%.f",
            treeName.c_str(), yBins[iy], yBins[iy+1]),
            "HFEmaxP_gN; HF Energy; Counts",
            20, 0, 100);
          HFEmaxM_gN[hIndex + iy] = new TH1D(
            Form("HFEmaxM_%s_gN_pt2-5_y%.f-%.f",
            treeName.c_str(), yBins[iy], yBins[iy+1]),
            "HFEmaxM_gN; HF Energy; Counts",
            20, 0, 100);
          zVtx_gN[hIndex + iy]->Sumw2();
          nVtx_gN[hIndex + iy]->Sumw2();
          HFEmaxP_gN[hIndex + iy]->Sumw2();
          HFEmaxM_gN[hIndex + iy]->Sumw2();
        }
        if (t == 0 || (t == 2)) {
          zVtx_Ng[hIndex + iy] = new TH1D(
            Form("zVtx_%s_Ng_pt2-5_y%.f-%.f",
            treeName.c_str(), yBins[iy], yBins[iy+1]),
            "zVtx_Ng; z-vertex; Counts",
            25, -25, 25);
          nVtx_Ng[hIndex + iy] = new TH1D(
            Form("nVtx_%s_Ng_pt2-5_y%.f-%.f",
            treeName.c_str(), yBins[iy], yBins[iy+1]),
            "nVtx_Ng; Number of Vertices; Counts",
            5, 0, 5);
          HFEmaxP_Ng[hIndex + iy] = new TH1D(
            Form("HFEmaxP_%s_Ng_pt2-5_y%.f-%.f",
            treeName.c_str(), yBins[iy], yBins[iy+1]),
            "HFEmaxP_Ng; HF Energy; Counts",
            20, 0, 100);
          HFEmaxM_Ng[hIndex + iy] = new TH1D(
            Form("HFEmaxM_%s_Ng_pt2-5_y%.f-%.f",
            treeName.c_str(), yBins[iy], yBins[iy+1]),
            "HFEmaxM_Ng; HF Energy; Counts",
            20, 0, 100);
          zVtx_Ng[hIndex + iy]->Sumw2();
          nVtx_Ng[hIndex + iy]->Sumw2();
          HFEmaxP_Ng[hIndex + iy]->Sumw2();
          HFEmaxM_Ng[hIndex + iy]->Sumw2();
        }
        if (t == 0) {
          ZDCsumP_gN[hIndex + iy] = new TH1D(
            Form("ZDCsumP_%s_gN_pt2-5_y%.f-%.f",
            treeName.c_str(), yBins[iy], yBins[iy+1]),
            "ZDCsumP_gN; ZDC Sum; Counts",
            10, 1000, 13000);
          ZDCsumP_Ng[hIndex + iy] = new TH1D(
            Form("ZDCsumP_%s_Ng_pt2-5_y%.f-%.f",
            treeName.c_str(), yBins[iy], yBins[iy+1]),
            "ZDCsumP_Ng; ZDC Sum; Counts",
            10, 1000, 13000);
          ZDCsumM_gN[hIndex + iy] = new TH1D(
            Form("ZDCsumM_%s_gN_pt2-5_y%.f-%.f",
            treeName.c_str(), yBins[iy], yBins[iy+1]),
            "ZDCsumM_gN; ZDC Sum; Counts",
            10, 1000, 13000);
          ZDCsumM_Ng[hIndex + iy] = new TH1D(
            Form("ZDCsumM_%s_Ng_pt2-5_y%.f-%.f",
            treeName.c_str(), yBins[iy], yBins[iy+1]),
            "ZDCsumM_Ng; ZDC Sum; Counts",
            10, 1000, 13000);
          ZDCsumP_gN[hIndex + iy]->Sumw2();
          ZDCsumP_Ng[hIndex + iy]->Sumw2();
          ZDCsumM_gN[hIndex + iy]->Sumw2();
          ZDCsumM_Ng[hIndex + iy]->Sumw2();
        }
      }
      
      float VZ;
      int   nVtx;
      float HFEMaxPlus;
      float HFEMaxMinus;
      float ZDCsumPlus;
      float ZDCsumMinus;
      int   Dsize;
      bool  gapgammaN;
      bool  gapNgamma;
      vector<bool>*  gammaN = nullptr;
      vector<bool>*  Ngamma = nullptr;
      vector<float>* Dpt = nullptr;
      vector<float>* Dy = nullptr;
      vector<float>* Dmass = nullptr;
      vector<bool>*  DpassCut23PAS = nullptr;
      vector<float>* Dtrk1Pt = nullptr;
      vector<float>* Dtrk2Pt = nullptr;
      vector<float>* DsvpvDistance = nullptr;
      vector<float>* DsvpvDisErr = nullptr;
      vector<float>* Dchi2cl = nullptr;
      vector<float>* Dalpha = nullptr;
      
      tree->SetBranchAddress("gammaN", &gammaN);
      tree->SetBranchAddress("Ngamma", &Ngamma);
      tree->SetBranchAddress("gapgammaN", &gapgammaN);
      tree->SetBranchAddress("gapNgamma", &gapNgamma);
      tree->SetBranchAddress("VZ", &VZ);
      tree->SetBranchAddress("nVtx", &nVtx);
      tree->SetBranchAddress("HFEMaxPlus", &HFEMaxPlus);
      tree->SetBranchAddress("HFEMaxMinus", &HFEMaxMinus);
      tree->SetBranchAddress("ZDCsumPlus", &ZDCsumPlus);
      tree->SetBranchAddress("ZDCsumMinus", &ZDCsumMinus);
      tree->SetBranchAddress("Dsize", &Dsize);
      tree->SetBranchAddress("Dpt", &Dpt);
      tree->SetBranchAddress("Dy", &Dy);
      tree->SetBranchAddress("Dmass", &Dmass);
      tree->SetBranchAddress("DpassCut23PAS", &DpassCut23PAS);
      tree->SetBranchAddress("Dtrk1Pt", &Dtrk1Pt);
      tree->SetBranchAddress("Dtrk2Pt", &Dtrk2Pt);
      tree->SetBranchAddress("DsvpvDistance", &DsvpvDistance);
      tree->SetBranchAddress("DsvpvDisErr", &DsvpvDisErr);
      tree->SetBranchAddress("Dchi2cl", &Dchi2cl);
      tree->SetBranchAddress("Dalpha", &Dalpha);
      
      // Loop over events
      for (unsigned long int e = 0; e < tree->GetEntries(); e++) {
        if (e % 1000000 == 0) cout << "Filling from entry " << e << "..." << endl;
        tree->GetEntry(e);
        
        // Loop over Ds
        for (int iD = 0; iD < Dsize; iD++) {
          if (!DpassCut23PAS->at(iD)) continue;
          
          // Rap gap checks
          if (t == 0 && ((gammaN->at(iD) && !gapgammaN) ||
              (Ngamma->at(iD) && !gapNgamma))) continue;
          else if (t == 1 && !gapgammaN) continue;
          else if (t == 2 && !gapNgamma) continue;
          
          if (Dtrk1Pt->at(iD) < DtrkPtCut ||
              Dtrk2Pt->at(iD) < DtrkPtCut) continue;
          if ((DsvpvDistance->at(iD)/DsvpvDisErr->at(iD)) < DsvpvCut) continue;
          if (Dchi2cl->at(iD) < Dchi2clCut) continue;
          if (Dmass->at(iD) < DMASSMIN || Dmass->at(iD) > DMASSMAX) continue;
          if (Dpt->at(iD) < 2. || Dpt->at(iD) > 5.) continue;
          if (Dy->at(iD) < -2 || Dy->at(iD) > 2) continue;

          int yIndex = -1;
          if      (Dy->at(iD) >= -2 && Dy->at(iD) < -1) yIndex = 0;
          else if (Dy->at(iD) >= -1 && Dy->at(iD) <  0) yIndex = 1;
          else if (Dy->at(iD) >=  0 && Dy->at(iD) <  1) yIndex = 2;
          else if (Dy->at(iD) >=  1 && Dy->at(iD) <= 2) yIndex = 3;
          if (Dalpha->at(iD) > DalphaCut[yIndex]) continue;
          
          if ((t == 0 && gammaN->at(iD)) || t == 1) {
            zVtx_gN[hIndex + yIndex]->Fill(VZ);
            nVtx_gN[hIndex + yIndex]->Fill(nVtx);
            HFEmaxP_gN[hIndex + yIndex]->Fill(HFEMaxPlus);
            HFEmaxM_gN[hIndex + yIndex]->Fill(HFEMaxMinus);
            if (t == 0) ZDCsumP_gN[hIndex + yIndex]->Fill(ZDCsumPlus);
            if (t == 0) ZDCsumM_gN[hIndex + yIndex]->Fill(ZDCsumMinus);
          }
          else if ((t == 0  && Ngamma->at(iD)) || t == 2) {
            zVtx_Ng[hIndex + yIndex]->Fill(VZ);
            nVtx_Ng[hIndex + yIndex]->Fill(nVtx);
            HFEmaxP_Ng[hIndex + yIndex]->Fill(HFEMaxPlus);
            HFEmaxM_Ng[hIndex + yIndex]->Fill(HFEMaxMinus);
            if (t == 0) ZDCsumP_Ng[hIndex + yIndex]->Fill(ZDCsumPlus);
            if (t == 0) ZDCsumM_Ng[hIndex + yIndex]->Fill(ZDCsumMinus);
          }
        } // end loop over D's
      } // end loop over entries
    } // end loop over trees
    
    system("mkdir -p plot/DataEval");
    TFile* fOutput = new TFile("plot/DataEval/hists.root", "RECREATE");
    for (int i = 0; i < zVtx_gN.size(); i++) {
      if (zVtx_gN[i]->GetEntries() != 0) zVtx_gN[i]->Write();
      if (zVtx_Ng[i]->GetEntries() != 0) zVtx_Ng[i]->Write();
      if (nVtx_gN[i]->GetEntries() != 0) nVtx_gN[i]->Write();
      if (nVtx_Ng[i]->GetEntries() != 0) nVtx_Ng[i]->Write();
      if (HFEmaxP_gN[i]->GetEntries() != 0) HFEmaxP_gN[i]->Write();
      if (HFEmaxP_Ng[i]->GetEntries() != 0) HFEmaxP_Ng[i]->Write();
      if (HFEmaxM_gN[i]->GetEntries() != 0) HFEmaxM_gN[i]->Write();
      if (HFEmaxM_Ng[i]->GetEntries() != 0) HFEmaxM_Ng[i]->Write();
      if (i < 4) {
        if (ZDCsumP_gN[i]->GetEntries() != 0) ZDCsumP_gN[i]->Write();
        if (ZDCsumP_Ng[i]->GetEntries() != 0) ZDCsumP_Ng[i]->Write();
        if (ZDCsumM_gN[i]->GetEntries() != 0) ZDCsumM_gN[i]->Write();
        if (ZDCsumM_Ng[i]->GetEntries() != 0) ZDCsumM_Ng[i]->Write();
      }
    }
    fOutput->Close();
    fData->Close();
    fMCA->Close();
    fMCB->Close();
    
    histFile = TFile::Open("plot/DataEval/hists.root", "READ");
  }

  for (int iy = 0; iy < nYBins; iy++) {
    int iDatagN = iy;
    int iDataNg = nYBins - iy - 1;
    int iMCgN = iy;
    int iMCNg = nYBins - iy - 1;
    
    TH1D* zVtx_gN_Data = (TH1D*) histFile->Get(
      Form("zVtx_Data_gN_pt2-5_y%.f-%.f",
      yBins[iDatagN], yBins[iDatagN+1]));
    TH1D* zVtx_Ng_Data = (TH1D*) histFile->Get(
      Form("zVtx_Data_Ng_pt2-5_y%.f-%.f",
      yBins[iDataNg], yBins[iDataNg+1]));
    TH1D* zVtx_gN_MC = (TH1D*) histFile->Get(
      Form("zVtx_MC-A_gN_pt2-5_y%.f-%.f",
      yBins[iMCgN], yBins[iMCgN+1]));
    TH1D* zVtx_Ng_MC = (TH1D*) histFile->Get(
      Form("zVtx_MC-B_Ng_pt2-5_y%.f-%.f",
      yBins[iMCNg], yBins[iMCNg+1]));
    
    TH1D* nVtx_gN_Data = (TH1D*) histFile->Get(
      Form("nVtx_Data_gN_pt2-5_y%.f-%.f",
      yBins[iDatagN], yBins[iDatagN+1]));
    TH1D* nVtx_Ng_Data = (TH1D*) histFile->Get(
      Form("nVtx_Data_Ng_pt2-5_y%.f-%.f",
      yBins[iDataNg], yBins[iDataNg+1]));
    TH1D* nVtx_gN_MC = (TH1D*) histFile->Get(
      Form("nVtx_MC-A_gN_pt2-5_y%.f-%.f",
      yBins[iMCgN], yBins[iMCgN+1]));
    TH1D* nVtx_Ng_MC = (TH1D*) histFile->Get(
      Form("nVtx_MC-B_Ng_pt2-5_y%.f-%.f",
      yBins[iMCNg], yBins[iMCNg+1]));
    
    TH1D* HFEmaxP_gN_Data = (TH1D*) histFile->Get(
      Form("HFEmaxP_Data_gN_pt2-5_y%.f-%.f",
      yBins[iDatagN], yBins[iDatagN+1]));
    TH1D* HFEmaxP_Ng_Data = (TH1D*) histFile->Get(
      Form("HFEmaxP_Data_Ng_pt2-5_y%.f-%.f",
      yBins[iDataNg], yBins[iDataNg+1]));
    TH1D* HFEmaxP_gN_MC = (TH1D*) histFile->Get(
      Form("HFEmaxP_MC-A_gN_pt2-5_y%.f-%.f",
      yBins[iMCgN], yBins[iMCgN+1]));
    TH1D* HFEmaxP_Ng_MC = (TH1D*) histFile->Get(
      Form("HFEmaxP_MC-B_Ng_pt2-5_y%.f-%.f",
      yBins[iMCNg], yBins[iMCNg+1]));
      
    TH1D* HFEmaxM_gN_Data = (TH1D*) histFile->Get(
      Form("HFEmaxM_Data_gN_pt2-5_y%.f-%.f",
      yBins[iDatagN], yBins[iDatagN+1]));
    TH1D* HFEmaxM_Ng_Data = (TH1D*) histFile->Get(
      Form("HFEmaxM_Data_Ng_pt2-5_y%.f-%.f",
      yBins[iDataNg], yBins[iDataNg+1]));
    TH1D* HFEmaxM_gN_MC = (TH1D*) histFile->Get(
      Form("HFEmaxM_MC-A_gN_pt2-5_y%.f-%.f",
      yBins[iMCgN], yBins[iMCgN+1]));
    TH1D* HFEmaxM_Ng_MC = (TH1D*) histFile->Get(
      Form("HFEmaxM_MC-B_Ng_pt2-5_y%.f-%.f",
      yBins[iMCNg], yBins[iMCNg+1]));
    
    TH1D* ZDCsumP_gN_Data = (TH1D*) histFile->Get(
      Form("ZDCsumP_Data_gN_pt2-5_y%.f-%.f",
      yBins[iDatagN], yBins[iDatagN+1]));
    TH1D* ZDCsumP_Ng_Data = (TH1D*) histFile->Get(
      Form("ZDCsumP_Data_Ng_pt2-5_y%.f-%.f",
      yBins[iDataNg], yBins[iDataNg+1]));
      
    TH1D* ZDCsumM_gN_Data = (TH1D*) histFile->Get(
      Form("ZDCsumM_Data_gN_pt2-5_y%.f-%.f",
      yBins[iDatagN], yBins[iDatagN+1]));
    TH1D* ZDCsumM_Ng_Data = (TH1D*) histFile->Get(
      Form("ZDCsumM_Data_Ng_pt2-5_y%.f-%.f",
      yBins[iDataNg], yBins[iDataNg+1]));
    
    string plotTitle = Form(
      "2 < Dp_{T} < 5 (GeV), %.f < Dy < %.f", yBins[iy], yBins[iy+1]);
    string fileName;
    cout << "Making plots for: " << plotTitle << endl;
    
    cout << "Plotting zVtx..." << endl;
    fileName = Form(
      "plot/DataEval/zVtx_pt2-5_y%.f-%.f.pdf", yBins[iy], yBins[iy+1]);
    plot_yBin(
      zVtx_gN_Data, // data_gN
      zVtx_Ng_Data, // data_Ng
      nullptr, // MC_gN
      nullptr, // MC_Ng
      plotTitle, // plotTitle
      "z-vertex", // plotXaxis
      "Counts", // plotYaxis
      fileName, // fileName
      0., 250,
      0., 4.
    );
    fileName = "";
    
    cout << "Plotting nVtx..." << endl;
    fileName = Form(
      "plot/DataEval/nVtx_pt2-5_y%.f-%.f.pdf", yBins[iy], yBins[iy+1]);
    plot_yBin(
      nVtx_gN_Data, // data_gN
      nVtx_Ng_Data, // data_Ng
      nullptr, // MC_gN
      nullptr, // MC_Ng
      plotTitle, // plotTitle
      "Number of Vertices", // plotXaxis
      "Counts", // plotYaxis
      fileName, // fileName
      0.1, 10000,
      0., 4.,
      true
    );
    fileName = "";
    
    cout << "Plotting HFEMaxPlus..." << endl;
    fileName = Form(
      "plot/DataEval/HFEMax_gNPlus_NgMinus_pt2-5_y%.f-%.f.pdf", yBins[iy], yBins[iy+1]);
    plot_yBin(
      HFEmaxP_gN_Data, // data_gN
      HFEmaxM_Ng_Data, // data_Ng
      nullptr, // MC_gN
      nullptr, // MC_Ng
      plotTitle, // plotTitle
      "HF-Plus (HF-Minus) Energy", // plotXaxis
      "Counts", // plotYaxis
      fileName, // fileName
      0., 250,
      0., 4.
    );
    fileName = "";
    
    cout << "Plotting HFEMaxMinus..." << endl;
    fileName = Form(
      "plot/DataEval/HFEMax_gNMinus_NgPlus_pt2-5_y%.f-%.f.pdf", yBins[iy], yBins[iy+1]);
    plot_yBin(
      HFEmaxM_gN_Data, // data_gN
      HFEmaxP_Ng_Data, // data_Ng
      nullptr, // MC_gN
      nullptr, // MC_Ng
      plotTitle, // plotTitle
      "HF-Minus (HF-Plus) Energy", // plotXaxis
      "Counts", // plotYaxis
      fileName, // fileName
      0., 250,
      0., 4.
    );
    fileName = "";
    
    cout << "Plotting ZDCsumPlus..." << endl;
    fileName = Form(
      "plot/DataEval/ZDCsum_gNPlus_NgMinus_pt2-5_y%.f-%.f.pdf", yBins[iy], yBins[iy+1]);
    plot_yBin(
      ZDCsumP_gN_Data, // data_gN
      ZDCsumM_Ng_Data, // data_Ng
      nullptr, // MC_gN
      nullptr, // MC_Ng
      plotTitle, // plotTitle
      "ZDC-Plus (ZDC-Minus) Sum", // plotXaxis
      "Counts", // plotYaxis
      fileName, // fileName
      0., 20,
      0., 4.
    );
    fileName = "";
    
    cout << "Plotting ZDCsumMinus..." << endl;
    fileName = Form(
      "plot/DataEval/ZDCsum_gNMinus_NgPlus_pt2-5_y%.f-%.f.pdf", yBins[iy], yBins[iy+1]);
    plot_yBin(
      ZDCsumM_gN_Data, // data_gN
      ZDCsumP_Ng_Data, // data_Ng
      nullptr, // MC_gN
      nullptr, // MC_Ng
      plotTitle, // plotTitle
      "ZDC-Minus (ZDC-Plus) Sum", // plotXaxis
      "Counts", // plotYaxis
      fileName, // fileName
      0., 20,
      0., 4.
    );
    fileName = "";
    
    zVtx_gN_Data->Scale(1/zVtx_gN_Data->Integral());
    zVtx_Ng_Data->Scale(1/zVtx_Ng_Data->Integral());
    zVtx_gN_MC->Scale(1/zVtx_gN_MC->Integral());
    zVtx_Ng_MC->Scale(1/zVtx_Ng_MC->Integral());
    
    nVtx_gN_Data->Scale(1/nVtx_gN_Data->Integral());
    nVtx_Ng_Data->Scale(1/nVtx_Ng_Data->Integral());
    nVtx_gN_MC->Scale(1/nVtx_gN_MC->Integral());
    nVtx_Ng_MC->Scale(1/nVtx_Ng_MC->Integral());
    
    HFEmaxP_gN_Data->Scale(1/HFEmaxP_gN_Data->Integral());
    HFEmaxP_Ng_Data->Scale(1/HFEmaxP_Ng_Data->Integral());
    HFEmaxM_gN_Data->Scale(1/HFEmaxM_gN_Data->Integral());
    HFEmaxM_Ng_Data->Scale(1/HFEmaxM_Ng_Data->Integral());
    
    HFEmaxP_gN_MC->Scale(1/HFEmaxP_gN_MC->Integral());
    HFEmaxP_Ng_MC->Scale(1/HFEmaxP_Ng_MC->Integral());
    HFEmaxM_gN_MC->Scale(1/HFEmaxM_gN_MC->Integral());
    HFEmaxM_Ng_MC->Scale(1/HFEmaxM_Ng_MC->Integral());
    
    ZDCsumP_gN_Data->Scale(1/ZDCsumP_gN_Data->Integral());
    ZDCsumP_Ng_Data->Scale(1/ZDCsumP_Ng_Data->Integral());
    ZDCsumM_gN_Data->Scale(1/ZDCsumM_gN_Data->Integral());
    ZDCsumM_Ng_Data->Scale(1/ZDCsumM_Ng_Data->Integral());
    
    cout << "Plotting zVtx..." << endl;
    fileName = Form(
      "plot/DataEval/zVtx_pt2-5_y%.f-%.f_norm.pdf", yBins[iy], yBins[iy+1]);
    plot_yBin(
      zVtx_gN_Data, // data_gN
      zVtx_Ng_Data, // data_Ng
      zVtx_gN_MC, // MC_gN
      zVtx_Ng_MC, // MC_Ng
      plotTitle, // plotTitle
      "z-vertex", // plotXaxis
      "Normalized Counts", // plotYaxis
      fileName, // fileName
      0., 0.3,
      0., 4.
    );
    fileName = "";
    
    cout << "Plotting nVtx..." << endl;
    fileName = Form(
      "plot/DataEval/nVtx_pt2-5_y%.f-%.f_norm.pdf", yBins[iy], yBins[iy+1]);
    plot_yBin(
      nVtx_gN_Data, // data_gN
      nVtx_Ng_Data, // data_Ng
      nVtx_gN_MC, // MC_gN
      nVtx_Ng_MC, // MC_Ng
      plotTitle, // plotTitle
      "Number of Vertices", // plotXaxis
      "Normalized Counts", // plotYaxis
      fileName, // fileName
      0.001, 10.0,
      0., 4.,
      true
    );
    fileName = "";
    
    cout << "Plotting HFEMaxPlus..." << endl;
    fileName = Form(
      "plot/DataEval/HFEMax_gNPlus_NgMinus_pt2-5_y%.f-%.f_norm.pdf", yBins[iy], yBins[iy+1]);
    plot_yBin(
      HFEmaxP_gN_Data, // data_gN
      HFEmaxM_Ng_Data, // data_Ng
      HFEmaxP_gN_MC, // MC_gN
      HFEmaxM_Ng_MC, // MC_Ng
      plotTitle, // plotTitle
      "HF-Plus (HF-Minus) Energy", // plotXaxis
      "Normalized Counts", // plotYaxis
      fileName, // fileName
      0., 0.05,
      0., 4.
    );
    fileName = "";
    
    cout << "Plotting HFEMaxMinus..." << endl;
    fileName = Form(
      "plot/DataEval/HFEMax_gNMinus_NgPlus_pt2-5_y%.f-%.f_norm.pdf", yBins[iy], yBins[iy+1]);
    plot_yBin(
      HFEmaxM_gN_Data, // data_gN
      HFEmaxP_Ng_Data, // data_Ng
      HFEmaxM_gN_MC, // MC_gN
      HFEmaxP_Ng_MC, // MC_Ng
      plotTitle, // plotTitle
      "HF-Minus (HF-Plus) Energy", // plotXaxis
      "Normalized Counts", // plotYaxis
      fileName, // fileName
      0., 0.35,
      0., 4.
    );
    fileName = "";
    
    cout << "Plotting ZDCsumPlus..." << endl;
    fileName = Form(
      "plot/DataEval/ZDCsum_gNPlus_NgMinus_pt2-5_y%.f-%.f_norm.pdf", yBins[iy], yBins[iy+1]);
    plot_yBin(
      ZDCsumP_gN_Data, // data_gN
      ZDCsumM_Ng_Data, // data_Ng
      nullptr, // MC_gN
      nullptr, // MC_Ng
      plotTitle, // plotTitle
      "ZDC-Plus (ZDC-Minus) Sum", // plotXaxis
      "Normalized Counts", // plotYaxis
      fileName, // fileName
      0., 0.4,
      0., 4.
    );
    fileName = "";
    
    cout << "Plotting ZDCsumMinus..." << endl;
    fileName = Form(
      "plot/DataEval/ZDCsum_gNMinus_NgPlus_pt2-5_y%.f-%.f_norm.pdf", yBins[iy], yBins[iy+1]);
    plot_yBin(
      ZDCsumM_gN_Data, // data_gN
      ZDCsumP_Ng_Data, // data_Ng
      nullptr, // MC_gN
      nullptr, // MC_Ng
      plotTitle, // plotTitle
      "ZDC-Minus (ZDC-Plus) Sum", // plotXaxis
      "Normalized Counts", // plotYaxis
      fileName, // fileName
      0., 0.4,
      0., 4.
    );
    fileName = "";
  }
}

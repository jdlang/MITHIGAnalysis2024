#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLegend.h>
#include <TNtuple.h>
#include <TTree.h>

#include <iostream>

using namespace std;
#include "CommandLine.h" // Yi's Commandline bundle
#include "Messenger.h"   // Yi's Messengers for reading data files
#include "ProgressBar.h" // Yi's fish progress bar
#include "helpMessage.h" // Print out help message
#include "parameter.h"   // The parameters used in the analysis
#include "utilities.h"   // Yen-Jie's random utility functions

#include "WeightHandler2D.h"
#include "WeightHandler1D.h"

#define DMASS 1.86484
#define DMASSMIN 1.66
#define DMASSMAX 2.26
#define DMASSNBINS 48

//============================================================//
// Function to check for configuration errors
//============================================================//
bool checkError(const Parameters &par) { return false; }

//======= trackSelection =====================================//
// Check if the track pass selection criteria
//============================================================//
bool dzeroSelection(DzeroUPCTreeMessenger *b, Parameters par, int j) { return true; }

//======= eventSelection =====================================//
// Check if the event pass eventSelection criteria
//============================================================//
bool eventSelection(DzeroUPCTreeMessenger *b, const Parameters &par) {
  if (par.IsData)
  {
    if (par.TriggerChoice == 1 && b->isL1ZDCOr == false)
      return false;
    if (par.TriggerChoice == 2 && b->isL1ZDCXORJet8 == false)
      return false;
  }

  if (b->selectedBkgFilter == false || b->selectedVtxFilter == false)
    return false;

  if (par.DoSystRapGap==-1)
  {
    // alternative (loose) rapidity gap selection
    if (par.IsGammaN && b->gammaN_EThreshSyst15() == false)
      return false;
    if (!par.IsGammaN && b->Ngamma_EThreshSyst15() == false)
      return false;
  }
  else if (par.DoSystRapGap==1)
  {
    // alternative (tight) rapidity gap selection
    if (par.IsGammaN && b->gammaN_EThreshSyst5p5() == false)
      return false;
    if (!par.IsGammaN && b->Ngamma_EThreshSyst5p5() == false)
      return false;
  } 
  else if (par.DoSystRapGap > 9) {
    // Custom rapidity gap threshold decision
    if (par.IsGammaN && b->gammaN_EThreshCustom(((float)par.DoSystRapGap)/10.) == false)
      return false;
    if (!par.IsGammaN && b->Ngamma_EThreshCustom(((float)par.DoSystRapGap)/10.) == false)
      return false;
  } 
  else
  {
    // nominal rapidity gap selection
    if (par.IsGammaN && (b->ZDCgammaN && b->gapgammaN) == false)
      return false;
    if (!par.IsGammaN && (b->ZDCNgamma && b->gapNgamma) == false)
      return false;
  }

  if (b->nVtx >= 3) return false;
  return true;
}

class DataAnalyzer {
public:
  TFile *inf, *outf;
  TH1D *hDmass;
  TH1D *hDtrk1Pt;
  TH1D *hDtrk2Pt;
  TH1D *hDsvpvSig;
  TH1D *hDalpha;
  TH1D *hDchi2cl;
  DzeroUPCTreeMessenger *MDzeroUPC;
  TNtuple *nt;
  string title;
  TH2D *hHFEmaxPlus_vs_EvtMult;
  TH2D *hHFEmaxMinus_vs_EvtMult;
  TH1D *hDenEvtEff;
  TH1D *hNumEvtEff;
  TH1D *hRatioEvtEff;
  TH1D *hDenDEff;
  TH1D *hNumDEff;
  TH1D *hRatioDEff;

  DataAnalyzer(const char *filename, const char *outFilename, const char *mytitle = "")
      : inf(new TFile(filename)), MDzeroUPC(new DzeroUPCTreeMessenger(*inf, string("Tree"))), title(mytitle),
        outf(new TFile(outFilename, "recreate")) {
    outf->cd();
    nt = new TNtuple("nt", "D0 mass tree", "Dmass:Dgen");
  }

  ~DataAnalyzer() {
    deleteHistograms();
    delete nt;
    inf->Close();
    outf->Close();
    delete MDzeroUPC;
  }

  void analyze(Parameters &par) {
    outf->cd();
    hDmass = new TH1D(Form("hDmass%s", title.c_str()), "", DMASSNBINS, DMASSMIN, DMASSMAX);
    hDtrk1Pt = new TH1D(Form("hDtrk1Pt%s", title.c_str()), "", 50, 0.0, 5.0);
    hDtrk2Pt = new TH1D(Form("hDtrk2Pt%s", title.c_str()), "", 50, 0.0, 5.0);
    hDsvpvSig = new TH1D(Form("hDsvpvSig%s", title.c_str()), "", 50, 0.0, 5.0);
    hDalpha = new TH1D(Form("hDalpha%s", title.c_str()), "", 32, 0.0, 3.2);
    hDchi2cl = new TH1D(Form("hDchi2cl%s", title.c_str()), "", 20, 0.0, 1.0);
    hDenEvtEff = new TH1D(Form("hDenEvtEff%s", title.c_str()), "", 1, 0.5, 1.5);
    hNumEvtEff = new TH1D(Form("hNumEvtEff%s", title.c_str()), "", 1, 0.5, 1.5);
    hRatioEvtEff = (TH1D*) hNumEvtEff->Clone("hRatioEvtEff");
    hDenDEff = new TH1D(Form("hDenDEff%s", title.c_str()), "", 1, 0.5, 1.5);
    hNumDEff = new TH1D(Form("hNumDEff%s", title.c_str()), "", 1, 0.5, 1.5);
    hRatioDEff = (TH1D*) hNumDEff->Clone("hRatioDEff");
    hHFEmaxMinus_vs_EvtMult = nullptr;
    hHFEmaxPlus_vs_EvtMult = nullptr;

    bool doHFEmaxDistributions=(par.DoSystRapGap > 9);
    if (doHFEmaxDistributions) {
      hHFEmaxMinus_vs_EvtMult = new TH2D(Form("hHFEmaxMinus_vs_EvtMult%s", title.c_str()), "", 80, 0, 20, 200, 0, 1000);
      hHFEmaxPlus_vs_EvtMult = new TH2D(Form("hHFEmaxPlus_vs_EvtMult%s", title.c_str()), "", 80, 0, 20, 200, 0, 1000);

      hHFEmaxMinus_vs_EvtMult->Sumw2();
      hHFEmaxPlus_vs_EvtMult->Sumw2();
    }

    hDmass->Sumw2();
    hDtrk1Pt->Sumw2();
    hDtrk2Pt->Sumw2();
    hDsvpvSig->Sumw2();
    hDalpha->Sumw2();
    hDchi2cl->Sumw2();
    hDenEvtEff->Sumw2();
    hNumEvtEff->Sumw2();
    hDenDEff->Sumw2();
    hNumDEff->Sumw2();

    WeightHandler2D GptGyWH;
    if (par.DoGptGyReweighting)
    {
      GptGyWH.LoadFromFile(par.GptGyWeightFileName);
    }

    WeightHandler1D MultWH;
    if (par.DoMultReweighting)
    {
      MultWH.LoadFromFile(par.MultWeightFileName);
    }

    par.printParameters();
    unsigned long nEntry = MDzeroUPC->GetEntries() * par.scaleFactor;
    ProgressBar Bar(cout, nEntry);
    Bar.SetStyle(1);

    for (unsigned long i = 0; i < nEntry; i++) {
      MDzeroUPC->GetEntry(i);
      if (i % 1000 == 0) {
        Bar.Update(i);
        Bar.Print();
      }

      // Check if the event is a signal MC event, a.k.a., having at least one gen-level D candidate that falls into the (pt,y) bin of interest
      bool isSigMCEvt = false;
      double leadingGpt = -999.;
      double leadingGy  = -999.;
      if (!par.IsData){
        for (unsigned long j = 0; j < MDzeroUPC->Gpt->size(); j++) {
          if (MDzeroUPC->GisSignalCalc->at(j) == false)
            continue;
          if (MDzeroUPC->Gpt->at(j) > leadingGpt)
          {
            leadingGpt = MDzeroUPC->Gpt->at(j);
            leadingGy  = MDzeroUPC->Gy ->at(j);
          }

          if (MDzeroUPC->Gpt->at(j) < par.MinDzeroPT)
            continue;
          if (MDzeroUPC->Gpt->at(j) > par.MaxDzeroPT)
            continue;
          if (MDzeroUPC->Gy->at(j) < par.MinDzeroY)
            continue;
          if (MDzeroUPC->Gy->at(j) > par.MaxDzeroY)
            continue;
          isSigMCEvt = true;
        }
      }

      double GptGyWeight = 1.0;
      double MultWeight  = 1.0;
      if (!par.IsData)
      {
        if (par.DoGptGyReweighting) GptGyWeight = GptGyWH.GetWeight(leadingGpt, leadingGy, (!par.IsGammaN));
        if (par.DoMultReweighting ) MultWeight  = MultWH.GetWeight(MDzeroUPC->nTrackInAcceptanceHP);
      }

      if (!par.IsData && isSigMCEvt) hDenEvtEff->Fill(1, GptGyWeight*MultWeight);

      // Check if the event passes the selection criteria
      if (eventSelection(MDzeroUPC, par)) {
        if (!par.IsData && isSigMCEvt) hNumEvtEff->Fill(1, GptGyWeight*MultWeight);

        for (unsigned long j = 0; j < MDzeroUPC->Dalpha->size(); j++) {
          if (MDzeroUPC->Dpt->at(j) < par.MinDzeroPT)
            continue;
          if (MDzeroUPC->Dpt->at(j) > par.MaxDzeroPT)
            continue;
          if (MDzeroUPC->Dy->at(j) < par.MinDzeroY)
            continue;
          if (MDzeroUPC->Dy->at(j) > par.MaxDzeroY)
            continue;
          // if (MDzeroUPC->DpassCut23PAS->at(j) == false)
          if (par.DoSystD==0 && MDzeroUPC->DpassCut23PAS->at(j) == false) continue;
          // if (par.DoSystD==0 && MDzeroUPC->DpassCut23LowPt->at(j) == false) continue;
          if (par.DoSystD==1 && MDzeroUPC->DpassCut23PASSystDsvpvSig->at(j) == false) continue;
          if (par.DoSystD==2 && MDzeroUPC->DpassCut23PASSystDtrkPt->at(j) == false) continue;
          if (par.DoSystD==3 && MDzeroUPC->DpassCut23PASSystDalpha->at(j) == false) continue;
          if (par.DoSystD==4 && MDzeroUPC->DpassCut23PASSystDchi2cl->at(j) == false) continue;
          
          // explicit DtrkPt filter (should not trigger if DtrkPtCut is not set)
          if (MDzeroUPC->Dtrk1Pt->at(j) < par.DtrkPtCut || MDzeroUPC->Dtrk2Pt->at(j) < par.DtrkPtCut) continue;
          
          // explicit Dsvpv filter (should not trigger if DsvpvCut is not set)
          if ((MDzeroUPC->DsvpvDistance->at(j)/MDzeroUPC->DsvpvDisErr->at(j)) < par.DsvpvCut) continue;

          // explicit Dalpha filter (should not trigger if DalphaCut is not set)
          if (par.DalphaCut > 0. && MDzeroUPC->Dalpha->at(j) > par.DalphaCut) continue;
          
          // explicit Dchi2cl filter (should not trigger if DalphaCut is not set)
          if (MDzeroUPC->Dchi2cl->at(j) < par.Dchi2Cut) continue;

          hDmass->Fill((*MDzeroUPC->Dmass)[j]);
          hDtrk1Pt->Fill((*MDzeroUPC->Dtrk1Pt)[j]);
          hDtrk2Pt->Fill((*MDzeroUPC->Dtrk2Pt)[j]);
          hDsvpvSig->Fill((*MDzeroUPC->DsvpvDistance)[j]/(*MDzeroUPC->DsvpvDisErr)[j]);
          hDalpha->Fill((*MDzeroUPC->Dalpha)[j]);
          hDchi2cl->Fill((*MDzeroUPC->Dchi2cl)[j]);
          if (!par.IsData) {
            nt->Fill((*MDzeroUPC->Dmass)[j], (*MDzeroUPC->Dgen)[j]);
            if (MDzeroUPC->Dgen->at(j) == 23333) {
              hNumDEff->Fill(1, GptGyWeight*MultWeight);
            }
          } else
            nt->Fill((*MDzeroUPC->Dmass)[j], 0);

          // Fill HF E_max distributions for data
          if(doHFEmaxDistributions && par.IsData) {
            hHFEmaxMinus_vs_EvtMult->Fill(MDzeroUPC->HFEMaxMinus, MDzeroUPC->nTrackInAcceptanceHP);
            hHFEmaxPlus_vs_EvtMult->Fill(MDzeroUPC->HFEMaxPlus, MDzeroUPC->nTrackInAcceptanceHP);
          }
        } // end of reco-level Dzero loop

        if (!par.IsData && isSigMCEvt) {
          for (unsigned long j = 0; j < MDzeroUPC->Gpt->size(); j++) {
            if (MDzeroUPC->Gpt->at(j) < par.MinDzeroPT)
              continue;
            if (MDzeroUPC->Gpt->at(j) > par.MaxDzeroPT)
              continue;
            if (MDzeroUPC->Gy->at(j) < par.MinDzeroY)
              continue;
            if (MDzeroUPC->Gy->at(j) > par.MaxDzeroY)
              continue;
            if (MDzeroUPC->GisSignalCalc->at(j) == false)
              continue;
            hDenDEff->Fill(1, GptGyWeight*MultWeight);
            // Fill HF E_max distributions for MC
            if(doHFEmaxDistributions) {
              hHFEmaxMinus_vs_EvtMult->Fill(MDzeroUPC->HFEMaxMinus, MDzeroUPC->nTrackInAcceptanceHP);
              hHFEmaxPlus_vs_EvtMult->Fill(MDzeroUPC->HFEMaxPlus, MDzeroUPC->nTrackInAcceptanceHP);
            }
          } // end of gen-level Dzero loop
        }   // end of gen-level Dzero loop
      }   // end of event selection
    }     // end of event loop
  }       // end of analyze

  void writeHistograms(TFile *outf) {
    outf->cd();
    smartWrite(hDmass);
    smartWrite(hDtrk1Pt);
    smartWrite(hDtrk2Pt);
    smartWrite(hDsvpvSig);
    smartWrite(hDalpha);
    smartWrite(hDchi2cl);
    hRatioEvtEff->Divide(hNumEvtEff, hDenEvtEff, 1, 1, "B");
    hRatioDEff->Divide(hNumDEff, hDenDEff, 1, 1, "B");
    hDenEvtEff->Write();
    hNumEvtEff->Write();
    hRatioEvtEff->Write();
    hDenDEff->Write();
    hNumDEff->Write();
    hRatioDEff->Write();
    smartWrite(nt);
    smartWrite(hHFEmaxMinus_vs_EvtMult);
    smartWrite(hHFEmaxPlus_vs_EvtMult);
  }

private:
  void deleteHistograms() {
    delete hDmass;
    delete hDtrk1Pt;
    delete hDtrk2Pt;
    delete hDsvpvSig;
    delete hDalpha;
    delete hDchi2cl;
    delete hDenEvtEff;
    delete hNumEvtEff;
    delete hRatioEvtEff;
    delete hDenDEff;
    delete hNumDEff;
    delete hRatioDEff;
    if (hHFEmaxMinus_vs_EvtMult != nullptr) {
      delete hHFEmaxMinus_vs_EvtMult;
    }
    if (hHFEmaxPlus_vs_EvtMult != nullptr) {
      delete hHFEmaxPlus_vs_EvtMult;
    }
  }
};

//============================================================//
// Main analysis
//============================================================//
int main(int argc, char *argv[]) {
  if (printHelpMessage(argc, argv))
    return 0;
  CommandLine CL(argc, argv);
  float MinDzeroPT = CL.GetDouble("MinDzeroPT", 2);  // Minimum Dzero transverse momentum threshold for Dzero selection.
  float MaxDzeroPT = CL.GetDouble("MaxDzeroPT", 5);  // Maximum Dzero transverse momentum threshold for Dzero selection.
  float MinDzeroY = CL.GetDouble("MinDzeroY", -2);   // Minimum Dzero rapidity threshold for Dzero selection.
  float MaxDzeroY = CL.GetDouble("MaxDzeroY", +2);   // Maximum Dzero rapidity threshold for Dzero selection.
  bool IsGammaN = CL.GetBool("IsGammaN", true);      // GammaN analysis (or NGamma)
  int TriggerChoice = CL.GetInt("TriggerChoice", 2); // 0 = no trigger sel, 1 = isL1ZDCOr, 2 = isL1ZDCXORJet8
  float scaleFactor = CL.GetDouble("scaleFactor", 1); // Scale factor for the number of events to be processed.
  int DoSystRapGap = CL.GetInt("DoSystRapGap", 0);   // Systematic study: apply the alternative event selections
                                                     // 0 = nominal, 1 = tight, -1: loose
                                                     // 9 < DoSystRapGap: use custom HF energy threshold, the threshold value will be DoSystRapGap/10.
  int DoSystD = CL.GetInt("DoSystD", 0);             // Systematic study: apply the alternative D selections
                                                     // 0 = nominal, 1 = Dsvpv variation, 2: DtrkPt variation
                                                     // 3 = Dalpha variation, 4: Dchi2cl variation
  bool DoGptGyReweighting   = CL.GetBool  ("DoGptGyReweighting", false);
  string GptGyWeightFileName= CL.Get      ("GptGyWeightFileName", "../../WeightHandler/20250305_DzeroUPC_GptGyWeight/Weights/testWeight.root");
  bool DoMultReweighting   = CL.GetBool  ("DoMultReweighting", false);
  string MultWeightFileName= CL.Get      ("MultWeightFileName", "../../WeightHandler/20250312_DzeroUPC_multiplicityWeight/Weights/testWeight.root");
  
  float DtrkPtCut = CL.GetDouble("DtrkPtCut", -1.0);  // Track pt cut. If >=0, overwrites DoSystD.
  float DsvpvCut = CL.GetDouble("DsvpvCut", -1.0);    // Track decay length cut. If >=0, overwrites DoSystD.
  float DalphaCut = CL.GetDouble("DalphaCut", -1.0);  // Alpha cut. If > 0, overwrites DoSystD.
  float Dchi2Cut = CL.GetDouble("Dchi2Cut", -1.0);    // Dchi2cl cut. If > 0, overwrites DoSystD.
  if (DsvpvCut >= 0.) DoSystD = 1;
  if (DtrkPtCut >= 0.) DoSystD = 2;
  if (DalphaCut > 0.) DoSystD = 3;
  if (Dchi2Cut >= 0.) DoSystD = 4;
  
  bool IsData = CL.GetBool("IsData", 0);              // Data or MC
  
  Parameters par(MinDzeroPT, MaxDzeroPT, MinDzeroY, MaxDzeroY, IsGammaN, TriggerChoice, IsData, scaleFactor,
                 DoSystRapGap, DoSystD,
                 DoGptGyReweighting, GptGyWeightFileName,
                 DoMultReweighting, MultWeightFileName,
                 DtrkPtCut, DsvpvCut, DalphaCut, Dchi2Cut);
  par.input = CL.Get("Input", "mergedSample.root"); // Input file
  par.output = CL.Get("Output", "output.root");     // Output file
  par.nThread = CL.GetInt("nThread", 1);            // The number of threads to be used for parallel processing.
  par.nChunk =
      CL.GetInt("nChunk", 1); // Specifies which chunk (segment) of the data to process, used in parallel processing.
  if (checkError(par))
    return -1;
  std::cout << "Parameters are set" << std::endl;
  // Analyze Data
  DataAnalyzer analyzer(par.input.c_str(), par.output.c_str(), "");
  analyzer.analyze(par);
  analyzer.writeHistograms(analyzer.outf);
  saveParametersToHistograms(par, analyzer.outf);
  cout << "done!" << analyzer.outf->GetName() << endl;

  return 0;
}

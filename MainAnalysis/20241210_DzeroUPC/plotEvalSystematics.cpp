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

#include "include/structs.h"
#include "plotCrossSection.h"

using namespace std;
using namespace RooFit;

#define DMASS 1.86484
#define DMASSMIN 1.66
#define DMASSMAX 2.26
#define DMASSNBINS 48

// Color settings
int kNominal = kBlack;
// MassFits
int kFitMean = kPink+2;
int kFitAlpha = kOrange-8;
int kFitComb = kSpring+2;
int kFitPkbg = kTeal-8;
int kMassWin = kAzure-8;
// Cuts
int kDalpha = kOrange+2;
int kDchi2cl = kSpring-8;
int kDsvpv = kTeal+2;
int kDtrkPt = kAzure+2;
int kRapGapLoose = kViolet+2;
int kRapGapTight = kPink-8;

void styleHist_MC(TH1D* mcHist)
{
  mcHist->SetMarkerColor(kRed);
  mcHist->SetMarkerStyle(24);
  mcHist->SetMarkerSize(1.);
  mcHist->SetLineColor(kRed);
  mcHist->SetLineWidth(2);
}

void styleHist_Data(TH1D* dataHist)
{
  dataHist->SetMarkerColor(kBlack);
  dataHist->SetMarkerStyle(20);
  dataHist->SetMarkerSize(1.);
  dataHist->SetLineColor(kBlack);
  dataHist->SetLineWidth(2);
}

void massfit_rawYield_vs_y(
  vector<string> directories,
  string plotTitle,
  string plotLabel,
  int nYBins,
  double* yBins,
  double ptMin,
  double ptMax,
  int isGammaN
) {
  vector<string> input_nominal;
  vector<string> input_fitMean;
  vector<string> input_fitAlpha;
  vector<string> input_fitComb;
  vector<string> input_fitPkbg;
  vector<string> input_massWin;
  
  for (const string& directory : directories) {
    input_nominal.push_back(Form("%s/MassFit/correctedYields.md", directory.c_str()));
    input_fitMean.push_back(Form("%s/MassFit_systFitSigMean/correctedYields.md", directory.c_str()));
    input_fitAlpha.push_back(Form("%s/MassFit_systFitSigAlpha/correctedYields.md", directory.c_str()));
    input_fitComb.push_back(Form("%s/MassFit_systComb/correctedYields.md", directory.c_str()));
    input_fitPkbg.push_back(Form("%s/MassFit_systPkBg/correctedYields.md", directory.c_str()));
    input_massWin.push_back(Form("%s/MassFit_systMassWindow040/correctedYields.md", directory.c_str()));
  }
  
  vector<Point> points_nominal  = getPointArr(ptMin, ptMax, isGammaN, input_nominal);
  vector<Point> points_fitMean  = getPointArr(ptMin, ptMax, isGammaN, input_fitMean);
  vector<Point> points_fitAlpha = getPointArr(ptMin, ptMax, isGammaN, input_fitAlpha);
  vector<Point> points_fitComb  = getPointArr(ptMin, ptMax, isGammaN, input_fitComb);
  vector<Point> points_fitPkbg  = getPointArr(ptMin, ptMax, isGammaN, input_fitPkbg);
  vector<Point> points_massWin  = getPointArr(ptMin, ptMax, isGammaN, input_massWin);
  
  vector<double> rawYieldVal_nominal  = getDoubleArr(points_nominal, [](Point& p) -> double { return p.rawYield;} );
  vector<double> rawYieldErr_nominal  = getDoubleArr(points_nominal, [](Point& p) -> double { return p.rawYieldError;} );
  vector<double> rawYieldVal_fitMean  = getDoubleArr(points_fitMean, [](Point& p) -> double { return p.rawYield;} );
  vector<double> rawYieldErr_fitMean  = getDoubleArr(points_fitMean, [](Point& p) -> double { return p.rawYieldError;} );
  vector<double> rawYieldVal_fitAlpha = getDoubleArr(points_fitAlpha, [](Point& p) -> double { return p.rawYield;} );
  vector<double> rawYieldErr_fitAlpha = getDoubleArr(points_fitAlpha, [](Point& p) -> double { return p.rawYieldError;} );
  vector<double> rawYieldVal_fitComb  = getDoubleArr(points_fitComb, [](Point& p) -> double { return p.rawYield;} );
  vector<double> rawYieldErr_fitComb  = getDoubleArr(points_fitComb, [](Point& p) -> double { return p.rawYieldError;} );
  vector<double> rawYieldVal_fitPkbg  = getDoubleArr(points_fitPkbg, [](Point& p) -> double { return p.rawYield;} );
  vector<double> rawYieldErr_fitPkbg  = getDoubleArr(points_fitPkbg, [](Point& p) -> double { return p.rawYieldError;} );
  vector<double> rawYieldVal_massWin  = getDoubleArr(points_massWin, [](Point& p) -> double { return p.rawYield;} );
  vector<double> rawYieldErr_massWin  = getDoubleArr(points_massWin, [](Point& p) -> double { return p.rawYieldError;} );

  TH1D* rawYieldTemplate  = new TH1D(
    "rawYieldTemplate", (plotTitle + "; y; Raw Yield").c_str(), nYBins, yBins);
  TGraphAsymmErrors* rawYield_nominal  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* rawYield_fitMean  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* rawYield_fitAlpha = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* rawYield_fitComb  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* rawYield_fitPkbg  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* rawYield_massWin  = new TGraphAsymmErrors(nYBins);
  
  TH1D* ratioTemplate  = new TH1D(
    "ratioTemplate", "; y; syst/nominal", nYBins, yBins);
  TGraphAsymmErrors* ratio_fitMean  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_fitAlpha = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_fitComb  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_fitPkbg  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_massWin  = new TGraphAsymmErrors(nYBins);
    
  double offset = 0.02;
  double offset_nominal  = 5 * offset;
  double offset_fitMean  = 3 * offset;
  double offset_fitAlpha = 1 * offset;
  double offset_fitComb  = -1 * offset;
  double offset_fitPkbg  = -3 * offset;
  double offset_massWin  = -5 * offset;
  for (int yBin = 0; yBin < nYBins; yBin++) {
    double yBinCenter = 0.5 * (yBins[yBin] + yBins[yBin+1]);
    // Set yield plots
    rawYield_nominal->SetPoint(yBin,
      yBinCenter - offset_nominal, rawYieldVal_nominal[yBin]);
    rawYield_nominal->SetPointError(yBin,
      0.5 - offset_nominal, 0.5 + offset_nominal,
      rawYieldErr_nominal[yBin], rawYieldErr_nominal[yBin]
    );
    rawYield_fitMean->SetPoint(yBin,
      yBinCenter - offset_fitMean, rawYieldVal_fitMean[yBin]);
    rawYield_fitMean->SetPointError(yBin,
      0.5 - offset_fitMean, 0.5 + offset_fitMean,
      rawYieldErr_fitMean[yBin], rawYieldErr_fitMean[yBin]
    );
    rawYield_fitAlpha->SetPoint(yBin,
      yBinCenter - offset_fitAlpha, rawYieldVal_fitAlpha[yBin]);
    rawYield_fitAlpha->SetPointError(yBin,
      0.5 - offset_fitAlpha, 0.5 + offset_fitAlpha,
      rawYieldErr_fitAlpha[yBin], rawYieldErr_fitAlpha[yBin]
    );
    rawYield_fitComb->SetPoint(yBin,
      yBinCenter - offset_fitComb, rawYieldVal_fitComb[yBin]);
    rawYield_fitComb->SetPointError(yBin,
      0.5 - offset_fitComb, 0.5 + offset_fitComb,
      rawYieldErr_fitComb[yBin], rawYieldErr_fitComb[yBin]
    );
    rawYield_fitPkbg->SetPoint(yBin,
      yBinCenter - offset_fitPkbg, rawYieldVal_fitPkbg[yBin]);
    rawYield_fitPkbg->SetPointError(yBin,
      0.5 - offset_fitPkbg, 0.5 + offset_fitPkbg,
      rawYieldErr_fitPkbg[yBin], rawYieldErr_fitPkbg[yBin]
    );
    rawYield_massWin->SetPoint(yBin,
      yBinCenter - offset_massWin, rawYieldVal_massWin[yBin]);
    rawYield_massWin->SetPointError(yBin,
      0.5 - offset_massWin, 0.5 + offset_massWin,
      rawYieldErr_massWin[yBin], rawYieldErr_massWin[yBin]
    );
    // Set ratio plots
    double ratioErr_fitMean = (
      TMath::Abs(rawYieldErr_fitMean[yBin] * rawYieldVal_nominal[yBin] -
      rawYieldVal_fitMean[yBin] * rawYieldErr_nominal[yBin]) /
      (rawYieldVal_nominal[yBin] * rawYieldVal_nominal[yBin])
    );
    ratio_fitMean->SetPoint(yBin,
      yBinCenter - offset_fitMean,
      rawYieldVal_fitMean[yBin]/rawYieldVal_nominal[yBin]);
    ratio_fitMean->SetPointError(yBin,
      0.5 - offset_fitMean, 0.5 + offset_fitMean,
      ratioErr_fitMean, ratioErr_fitMean
    );
    double ratioErr_fitAlpha = (
      TMath::Abs(rawYieldErr_fitAlpha[yBin] * rawYieldVal_nominal[yBin] -
      rawYieldVal_fitAlpha[yBin] * rawYieldErr_nominal[yBin]) /
      (rawYieldVal_nominal[yBin] * rawYieldVal_nominal[yBin])
    );
    ratio_fitAlpha->SetPoint(yBin,
      yBinCenter - offset_fitAlpha,
      rawYieldVal_fitAlpha[yBin]/rawYieldVal_nominal[yBin]);
    ratio_fitAlpha->SetPointError(yBin,
      0.5 - offset_fitAlpha, 0.5 + offset_fitAlpha,
      ratioErr_fitAlpha, ratioErr_fitAlpha
    );
    double ratioErr_fitComb = (
      TMath::Abs(rawYieldErr_fitComb[yBin] * rawYieldVal_nominal[yBin] -
      rawYieldVal_fitComb[yBin] * rawYieldErr_nominal[yBin]) /
      (rawYieldVal_nominal[yBin] * rawYieldVal_nominal[yBin])
    );
    ratio_fitComb->SetPoint(yBin,
      yBinCenter - offset_fitComb,
      rawYieldVal_fitComb[yBin]/rawYieldVal_nominal[yBin]);
    ratio_fitComb->SetPointError(yBin,
      0.5 - offset_fitComb, 0.5 + offset_fitComb,
      ratioErr_fitComb, ratioErr_fitComb
    );
    double ratioErr_fitPkbg = (
      TMath::Abs(rawYieldErr_fitPkbg[yBin] * rawYieldVal_nominal[yBin] -
      rawYieldVal_fitPkbg[yBin] * rawYieldErr_nominal[yBin]) /
      (rawYieldVal_nominal[yBin] * rawYieldVal_nominal[yBin])
    );
    ratio_fitPkbg->SetPoint(yBin,
      yBinCenter - offset_fitPkbg,
      rawYieldVal_fitPkbg[yBin]/rawYieldVal_nominal[yBin]);
    ratio_fitPkbg->SetPointError(yBin,
      0.5 - offset_fitPkbg, 0.5 + offset_fitPkbg,
      ratioErr_fitPkbg, ratioErr_fitPkbg
    );
    double ratioErr_massWin = (
      TMath::Abs(rawYieldErr_massWin[yBin] * rawYieldVal_nominal[yBin] -
      rawYieldVal_massWin[yBin] * rawYieldErr_nominal[yBin]) /
      (rawYieldVal_nominal[yBin] * rawYieldVal_nominal[yBin])
    );
    ratio_massWin->SetPoint(yBin,
      yBinCenter - offset_massWin,
      rawYieldVal_massWin[yBin]/rawYieldVal_nominal[yBin]);
    ratio_massWin->SetPointError(yBin,
      0.5 - offset_massWin, 0.5 + offset_massWin,
      ratioErr_massWin, ratioErr_massWin
    );
  }
  rawYieldTemplate->SetMinimum(0.);
  rawYieldTemplate->SetMaximum(160.);
  ratioTemplate->SetMinimum(0.7);
  ratioTemplate->SetMaximum(1.3);
  
  // Nominal
  rawYield_nominal->SetLineColor(kNominal);
  rawYield_nominal->SetLineWidth(2);
  // SystFitMean
  rawYield_fitMean->SetLineColor(kFitMean);
  rawYield_fitMean->SetLineWidth(2);
  ratio_fitMean->SetLineColor(kFitMean);
  ratio_fitMean->SetLineWidth(2);
  // SystFitAlpha
  rawYield_fitAlpha->SetLineColor(kFitAlpha);
  rawYield_fitAlpha->SetLineWidth(2);
  ratio_fitAlpha->SetLineColor(kFitAlpha);
  ratio_fitAlpha->SetLineWidth(2);
  // SystFitComb
  rawYield_fitComb->SetLineColor(kFitComb);
  rawYield_fitComb->SetLineWidth(2);
  ratio_fitComb->SetLineColor(kFitComb);
  ratio_fitComb->SetLineWidth(2);
  // SystFitPkbg
  rawYield_fitPkbg->SetLineColor(kFitPkbg);
  rawYield_fitPkbg->SetLineWidth(2);
  ratio_fitPkbg->SetLineColor(kFitPkbg);
  ratio_fitPkbg->SetLineWidth(2);
  // SystMassWin
  rawYield_massWin->SetLineColor(kMassWin);
  rawYield_massWin->SetLineWidth(2);
  ratio_massWin->SetLineColor(kMassWin);
  ratio_massWin->SetLineWidth(2);
  
  TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
  TPad* padTop = new TPad("padTop", "", 0.0, 0.3, 1.0, 1.0);
  TPad* padBot = new TPad("padBot", "", 0.0, 0.0, 1.0, 0.3);
  padTop->SetMargin(0.18, 0.02, 0.00, 0.16);
  padBot->SetMargin(0.18, 0.02, 0.30, 0.00);
  padTop->Draw();
  padBot->Draw();
  
  padTop->cd();
  rawYieldTemplate->Draw();
  rawYieldTemplate->SetLabelSize(0.03/0.7, "Y");
  rawYieldTemplate->SetTitleSize(0.035/0.7, "Y");
  rawYield_nominal->Draw("same p");
  rawYield_fitMean->Draw("same p");
  rawYield_fitAlpha->Draw("same p");
  rawYield_fitComb->Draw("same p");
  rawYield_fitPkbg->Draw("same p");
  rawYield_massWin->Draw("same p");
  gStyle->SetOptStat(0);
  
  padBot->cd();
  ratioTemplate->Draw();
  ratioTemplate->SetLabelSize(0.03/0.3, "XY");
  ratioTemplate->SetTitleSize(0.035/0.3, "XY");
  ratioTemplate->SetTitleOffset(0.4, "Y");
  TLine* unity = new TLine(yBins[0], 1.0, yBins[nYBins], 1.0);
  unity->SetLineColor(kGray);
  unity->SetLineWidth(1);
  unity->SetLineStyle(9);
  unity->Draw();
  ratio_fitMean->Draw("same p");
  ratio_fitAlpha->Draw("same p");
  ratio_fitComb->Draw("same p");
  ratio_fitPkbg->Draw("same p");
  ratio_massWin->Draw("same p");
  gStyle->SetOptStat(0);
  
  canvas->cd();
  canvas->Update();
  float legShift = 0.;
  if (isGammaN) legShift = 0.4;
  TLegend* legend = new TLegend(0.2 + legShift, 0.68, 0.6 + legShift, 0.88);
  legend->SetTextSize(0.015/0.7);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(rawYield_nominal,  "Nominal Raw Yield", "l");
  legend->AddEntry(rawYield_fitMean,  "Fit Syst: Signal mean", "l");
  legend->AddEntry(rawYield_fitAlpha, "Fit Syst: Signal alpha", "l");
  legend->AddEntry(rawYield_fitComb,  "Fit Syst: Comb. Background", "l");
  legend->AddEntry(rawYield_fitPkbg,  "Fit Syst: KK + #pi#pi Peaks", "l");
  legend->AddEntry(rawYield_massWin,  "Fit Syst: Mass Window", "l");
  legend->Draw();
  
  canvas->SaveAs(Form("plot/SystEval_MassFit/%s_rawYield_vs_y.pdf", plotLabel.c_str()));
  
  delete canvas;
  delete rawYieldTemplate;
  delete rawYield_nominal;
  delete rawYield_fitMean;
  delete rawYield_fitAlpha;
  delete rawYield_fitComb;
  delete rawYield_fitPkbg;
  delete rawYield_massWin;
  delete ratioTemplate;
  delete ratio_fitMean;
  delete ratio_fitAlpha;
  delete ratio_fitComb;
  delete ratio_fitPkbg;
  delete ratio_massWin;
}

void massfit_lambda_vs_y(
  vector<string> directories,
  string plotTitle,
  string plotLabel,
  int nYBins,
  double* yBins,
  double ptMin,
  double ptMax,
  int isGammaN
) {
  vector<string> input_nominal;
  vector<string> input_fitMean;
  vector<string> input_fitAlpha;
  // NO fitComb since it does not use exponential
  vector<string> input_fitPkbg;
  vector<string> input_massWin;
  
  for (const string& directory : directories) {
    input_nominal.push_back(Form("%s/MassFit/combdata.dat", directory.c_str()));
    input_fitMean.push_back(Form("%s/MassFit_systFitSigMean/combdata.dat", directory.c_str()));
    input_fitAlpha.push_back(Form("%s/MassFit_systFitSigAlpha/combdata.dat", directory.c_str()));
    input_fitPkbg.push_back(Form("%s/MassFit_systPkBg/combdata.dat", directory.c_str()));
    input_massWin.push_back(Form("%s/MassFit_systMassWindow040/combdata.dat", directory.c_str()));
  }
    
  vector<double> lambdaVal_nominal;
  vector<double> lambdaErr_nominal;
  vector<double> lambdaVal_fitMean;
  vector<double> lambdaErr_fitMean;
  vector<double> lambdaVal_fitAlpha;
  vector<double> lambdaErr_fitAlpha;
  vector<double> lambdaVal_fitPkbg;
  vector<double> lambdaErr_fitPkbg;
  vector<double> lambdaVal_massWin;
  vector<double> lambdaErr_massWin;
  
  for (int i = 0; i < directories.size(); i++) {
    CombinatoricsBkgParams comb_nominal = CombinatoricsBkgParams(input_nominal[i]);
    lambdaVal_nominal.push_back(comb_nominal.lambda.getVal());
    lambdaErr_nominal.push_back(comb_nominal.lambda.getError());
    CombinatoricsBkgParams comb_fitMean = CombinatoricsBkgParams(input_fitMean[i]);
    lambdaVal_fitMean.push_back(comb_fitMean.lambda.getVal());
    lambdaErr_fitMean.push_back(comb_fitMean.lambda.getError());
    CombinatoricsBkgParams comb_fitAlpha = CombinatoricsBkgParams(input_fitAlpha[i]);
    lambdaVal_fitAlpha.push_back(comb_fitAlpha.lambda.getVal());
    lambdaErr_fitAlpha.push_back(comb_fitAlpha.lambda.getError());
    CombinatoricsBkgParams comb_fitPkbg = CombinatoricsBkgParams(input_fitPkbg[i]);
    lambdaVal_fitPkbg.push_back(comb_fitPkbg.lambda.getVal());
    lambdaErr_fitPkbg.push_back(comb_fitPkbg.lambda.getError());
    CombinatoricsBkgParams comb_massWin = CombinatoricsBkgParams(input_massWin[i]);
    lambdaVal_massWin.push_back(comb_massWin.lambda.getVal());
    lambdaErr_massWin.push_back(comb_massWin.lambda.getError());
  }

  TH1D* lambdaTemplate  = new TH1D(
    "lambdaTemplate", (plotTitle + "; y; Raw Yield").c_str(), nYBins, yBins);
  TGraphAsymmErrors* lambda_nominal  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* lambda_fitMean  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* lambda_fitAlpha = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* lambda_fitPkbg  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* lambda_massWin  = new TGraphAsymmErrors(nYBins);
  
  TH1D* ratioTemplate  = new TH1D(
    "ratioTemplate", "; y; syst/nominal", nYBins, yBins);
  TGraphAsymmErrors* ratio_fitMean  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_fitAlpha = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_fitPkbg  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_massWin  = new TGraphAsymmErrors(nYBins);
    
  double offset = 0.02;
  double offset_nominal  = 5 * offset;
  double offset_fitMean  = 3 * offset;
  double offset_fitAlpha = 1 * offset;
  double offset_fitPkbg  = -3 * offset;
  double offset_massWin  = -5 * offset;
  for (int yBin = 0; yBin < nYBins; yBin++) {
    double yBinCenter = 0.5 * (yBins[yBin] + yBins[yBin+1]);
    // Set yield plots
    lambda_nominal->SetPoint(yBin,
      yBinCenter - offset_nominal, lambdaVal_nominal[yBin]);
    lambda_nominal->SetPointError(yBin,
      0.5 - offset_nominal, 0.5 + offset_nominal,
      lambdaErr_nominal[yBin], lambdaErr_nominal[yBin]
    );
    lambda_fitMean->SetPoint(yBin,
      yBinCenter - offset_fitMean, lambdaVal_fitMean[yBin]);
    lambda_fitMean->SetPointError(yBin,
      0.5 - offset_fitMean, 0.5 + offset_fitMean,
      lambdaErr_fitMean[yBin], lambdaErr_fitMean[yBin]
    );
    lambda_fitAlpha->SetPoint(yBin,
      yBinCenter - offset_fitAlpha, lambdaVal_fitAlpha[yBin]);
    lambda_fitAlpha->SetPointError(yBin,
      0.5 - offset_fitAlpha, 0.5 + offset_fitAlpha,
      lambdaErr_fitAlpha[yBin], lambdaErr_fitAlpha[yBin]
    );
    lambda_fitPkbg->SetPoint(yBin,
      yBinCenter - offset_fitPkbg, lambdaVal_fitPkbg[yBin]);
    lambda_fitPkbg->SetPointError(yBin,
      0.5 - offset_fitPkbg, 0.5 + offset_fitPkbg,
      lambdaErr_fitPkbg[yBin], lambdaErr_fitPkbg[yBin]
    );
    lambda_massWin->SetPoint(yBin,
      yBinCenter - offset_massWin, lambdaVal_massWin[yBin]);
    lambda_massWin->SetPointError(yBin,
      0.5 - offset_massWin, 0.5 + offset_massWin,
      lambdaErr_massWin[yBin], lambdaErr_massWin[yBin]
    );
    // Set ratio plots
    double ratioErr_fitMean = (
      TMath::Abs(lambdaErr_fitMean[yBin] * lambdaVal_nominal[yBin] -
      lambdaVal_fitMean[yBin] * lambdaErr_nominal[yBin]) /
      (lambdaVal_nominal[yBin] * lambdaVal_nominal[yBin])
    );
    ratio_fitMean->SetPoint(yBin,
      yBinCenter - offset_fitMean,
      lambdaVal_fitMean[yBin]/lambdaVal_nominal[yBin]);
    ratio_fitMean->SetPointError(yBin,
      0.5 - offset_fitMean, 0.5 + offset_fitMean,
      ratioErr_fitMean, ratioErr_fitMean
    );
    double ratioErr_fitAlpha = (
      TMath::Abs(lambdaErr_fitAlpha[yBin] * lambdaVal_nominal[yBin] -
      lambdaVal_fitAlpha[yBin] * lambdaErr_nominal[yBin]) /
      (lambdaVal_nominal[yBin] * lambdaVal_nominal[yBin])
    );
    ratio_fitAlpha->SetPoint(yBin,
      yBinCenter - offset_fitAlpha,
      lambdaVal_fitAlpha[yBin]/lambdaVal_nominal[yBin]);
    ratio_fitAlpha->SetPointError(yBin,
      0.5 - offset_fitAlpha, 0.5 + offset_fitAlpha,
      ratioErr_fitAlpha, ratioErr_fitAlpha
    );
    double ratioErr_fitPkbg = (
      TMath::Abs(lambdaErr_fitPkbg[yBin] * lambdaVal_nominal[yBin] -
      lambdaVal_fitPkbg[yBin] * lambdaErr_nominal[yBin]) /
      (lambdaVal_nominal[yBin] * lambdaVal_nominal[yBin])
    );
    ratio_fitPkbg->SetPoint(yBin,
      yBinCenter - offset_fitPkbg,
      lambdaVal_fitPkbg[yBin]/lambdaVal_nominal[yBin]);
    ratio_fitPkbg->SetPointError(yBin,
      0.5 - offset_fitPkbg, 0.5 + offset_fitPkbg,
      ratioErr_fitPkbg, ratioErr_fitPkbg
    );
    double ratioErr_massWin = (
      TMath::Abs(lambdaErr_massWin[yBin] * lambdaVal_nominal[yBin] -
      lambdaVal_massWin[yBin] * lambdaErr_nominal[yBin]) /
      (lambdaVal_nominal[yBin] * lambdaVal_nominal[yBin])
    );
    ratio_massWin->SetPoint(yBin,
      yBinCenter - offset_massWin,
      lambdaVal_massWin[yBin]/lambdaVal_nominal[yBin]);
    ratio_massWin->SetPointError(yBin,
      0.5 - offset_massWin, 0.5 + offset_massWin,
      ratioErr_massWin, ratioErr_massWin
    );
  }
  lambdaTemplate->SetMinimum(-10.);
  lambdaTemplate->SetMaximum(0.0);
  ratioTemplate->SetMinimum(0.7);
  ratioTemplate->SetMaximum(1.3);
  
  // Nominal
  lambda_nominal->SetLineColor(kNominal);
  lambda_nominal->SetLineWidth(2);
  // SystFitMean
  lambda_fitMean->SetLineColor(kFitMean);
  lambda_fitMean->SetLineWidth(2);
  ratio_fitMean->SetLineColor(kFitMean);
  ratio_fitMean->SetLineWidth(2);
  // SystFitAlpha
  lambda_fitAlpha->SetLineColor(kFitAlpha);
  lambda_fitAlpha->SetLineWidth(2);
  ratio_fitAlpha->SetLineColor(kFitAlpha);
  ratio_fitAlpha->SetLineWidth(2);
  // SystFitPkbg
  lambda_fitPkbg->SetLineColor(kFitPkbg);
  lambda_fitPkbg->SetLineWidth(2);
  ratio_fitPkbg->SetLineColor(kFitPkbg);
  ratio_fitPkbg->SetLineWidth(2);
  // SystMassWin
  lambda_massWin->SetLineColor(kMassWin);
  lambda_massWin->SetLineWidth(2);
  ratio_massWin->SetLineColor(kMassWin);
  ratio_massWin->SetLineWidth(2);
  
  TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
  TPad* padTop = new TPad("padTop", "", 0.0, 0.3, 1.0, 1.0);
  TPad* padBot = new TPad("padBot", "", 0.0, 0.0, 1.0, 0.3);
  padTop->SetMargin(0.18, 0.02, 0.00, 0.16);
  padBot->SetMargin(0.18, 0.02, 0.30, 0.00);
  padTop->Draw();
  padBot->Draw();
  
  padTop->cd();
  lambdaTemplate->Draw();
  lambdaTemplate->SetLabelSize(0.03/0.7, "Y");
  lambdaTemplate->SetTitleSize(0.035/0.7, "Y");
  lambda_nominal->Draw("same p");
  lambda_fitMean->Draw("same p");
  lambda_fitAlpha->Draw("same p");
  lambda_fitPkbg->Draw("same p");
  lambda_massWin->Draw("same p");
  gStyle->SetOptStat(0);
  
  padBot->cd();
  ratioTemplate->Draw();
  ratioTemplate->SetLabelSize(0.03/0.3, "XY");
  ratioTemplate->SetTitleSize(0.035/0.3, "XY");
  ratioTemplate->SetTitleOffset(0.4, "Y");
  TLine* unity = new TLine(yBins[0], 1.0, yBins[nYBins], 1.0);
  unity->SetLineColor(kGray);
  unity->SetLineWidth(1);
  unity->SetLineStyle(9);
  unity->Draw();
  ratio_fitMean->Draw("same p");
  ratio_fitAlpha->Draw("same p");
  ratio_fitPkbg->Draw("same p");
  ratio_massWin->Draw("same p");
  gStyle->SetOptStat(0);
  
  canvas->cd();
  canvas->Update();
  float legShift = 0.;
  if (isGammaN) legShift = 0.4;
  TLegend* legend = new TLegend(0.2 + legShift, 0.72, 0.6 + legShift, 0.88);
  legend->SetTextSize(0.015/0.7);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(lambda_nominal,  "Nominal Exponential #lambda", "l");
  legend->AddEntry(lambda_fitMean,  "Fit Syst: Signal mean", "l");
  legend->AddEntry(lambda_fitAlpha, "Fit Syst: Signal alpha", "l");
  legend->AddEntry(lambda_fitPkbg,  "Fit Syst: KK + #pi#pi Peaks", "l");
  legend->AddEntry(lambda_massWin,  "Fit Syst: Mass Window", "l");
  legend->Draw();
  
  canvas->SaveAs(Form("plot/SystEval_MassFit/%s_lambda_vs_y.pdf", plotLabel.c_str()));
  
  delete canvas;
  delete lambdaTemplate;
  delete lambda_nominal;
  delete lambda_fitMean;
  delete lambda_fitAlpha;
  delete lambda_fitPkbg;
  delete lambda_massWin;
  delete ratioTemplate;
  delete ratio_fitMean;
  delete ratio_fitAlpha;
  delete ratio_fitPkbg;
  delete ratio_massWin;
}

void cuts_rawYield_vs_y(
  vector<string> nominalDirs,
  vector<string> systDalphaDirs,
  vector<string> systDchi2clDirs,
  vector<string> systDsvpvDirs,
  vector<string> systDtrkPtDirs,
  vector<string> systRapGapLooseDirs,
  vector<string> systRapGapTightDirs,
  string plotTitle,
  string plotLabel,
  int nYBins,
  double* yBins,
  double ptMin,
  double ptMax,
  int isGammaN
) {
  vector<string> input_nominal;
  vector<string> input_Dalpha;
  vector<string> input_Dchi2cl;
  vector<string> input_Dsvpv;
  vector<string> input_DtrkPt;
  vector<string> input_RapGapLoose;
  vector<string> input_RapGapTight;
  
  for (int i = 0; i < nominalDirs.size(); i++) {
    input_nominal.push_back(Form("%s/MassFit/correctedYields.md", nominalDirs[i].c_str()));
    input_Dalpha.push_back(Form("%s/MassFit/correctedYields.md", systDalphaDirs[i].c_str()));
    input_Dchi2cl.push_back(Form("%s/MassFit/correctedYields.md", systDchi2clDirs[i].c_str()));
    input_Dsvpv.push_back(Form("%s/MassFit/correctedYields.md", systDsvpvDirs[i].c_str()));
    input_DtrkPt.push_back(Form("%s/MassFit/correctedYields.md", systDtrkPtDirs[i].c_str()));
    input_RapGapLoose.push_back(Form("%s/MassFit/correctedYields.md", systRapGapLooseDirs[i].c_str()));
    input_RapGapTight.push_back(Form("%s/MassFit/correctedYields.md", systRapGapTightDirs[i].c_str()));
  }
  
  vector<Point> points_nominal  = getPointArr(ptMin, ptMax, isGammaN, input_nominal);
  vector<Point> points_Dalpha  = getPointArr(ptMin, ptMax, isGammaN, input_Dalpha);
  vector<Point> points_Dchi2cl = getPointArr(ptMin, ptMax, isGammaN, input_Dchi2cl);
  vector<Point> points_Dsvpv  = getPointArr(ptMin, ptMax, isGammaN, input_Dsvpv);
  vector<Point> points_DtrkPt  = getPointArr(ptMin, ptMax, isGammaN, input_DtrkPt);
  vector<Point> points_RapGapLoose  = getPointArr(ptMin, ptMax, isGammaN, input_RapGapLoose);
  vector<Point> points_RapGapTight  = getPointArr(ptMin, ptMax, isGammaN, input_RapGapTight);
  
  vector<double> rawYieldVal_nominal  = getDoubleArr(points_nominal, [](Point& p) -> double { return p.rawYield;});
  vector<double> rawYieldErr_nominal  = getDoubleArr(points_nominal, [](Point& p) -> double { return p.rawYieldError;});
  vector<double> rawYieldVal_Dalpha  = getDoubleArr(points_Dalpha, [](Point& p) -> double { return p.rawYield;});
  vector<double> rawYieldErr_Dalpha  = getDoubleArr(points_Dalpha, [](Point& p) -> double { return p.rawYieldError;});
  vector<double> rawYieldVal_Dchi2cl = getDoubleArr(points_Dchi2cl, [](Point& p) -> double { return p.rawYield;});
  vector<double> rawYieldErr_Dchi2cl = getDoubleArr(points_Dchi2cl, [](Point& p) -> double { return p.rawYieldError;});
  vector<double> rawYieldVal_Dsvpv  = getDoubleArr(points_Dsvpv, [](Point& p) -> double { return p.rawYield;});
  vector<double> rawYieldErr_Dsvpv  = getDoubleArr(points_Dsvpv, [](Point& p) -> double { return p.rawYieldError;});
  vector<double> rawYieldVal_DtrkPt  = getDoubleArr(points_DtrkPt, [](Point& p) -> double { return p.rawYield;});
  vector<double> rawYieldErr_DtrkPt  = getDoubleArr(points_DtrkPt, [](Point& p) -> double { return p.rawYieldError;});
  vector<double> rawYieldVal_RapGapLoose = getDoubleArr(points_RapGapLoose, [](Point& p) -> double { return p.rawYield;});
  vector<double> rawYieldErr_RapGapLoose = getDoubleArr(points_RapGapLoose, [](Point& p) -> double { return p.rawYieldError;});
  vector<double> rawYieldVal_RapGapTight = getDoubleArr(points_RapGapLoose, [](Point& p) -> double { return p.rawYield;});
  vector<double> rawYieldErr_RapGapTight = getDoubleArr(points_RapGapLoose, [](Point& p) -> double { return p.rawYieldError;});

  TH1D* rawYieldTemplate  = new TH1D(
    "rawYieldTemplate", (plotTitle + "; y; Raw Yield").c_str(), nYBins, yBins);
  TGraphAsymmErrors* rawYield_nominal  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* rawYield_Dalpha  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* rawYield_Dchi2cl = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* rawYield_Dsvpv  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* rawYield_DtrkPt  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* rawYield_RapGapLoose  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* rawYield_RapGapTight  = new TGraphAsymmErrors(nYBins);

  TH1D* ratioTemplate  = new TH1D(
    "ratioTemplate", "; y; syst/nominal", nYBins, yBins);
  TGraphAsymmErrors* ratio_Dalpha  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_Dchi2cl = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_Dsvpv  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_DtrkPt  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_RapGapLoose  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_RapGapTight  = new TGraphAsymmErrors(nYBins);
    
  double offset = 0.02;
  double offset_nominal  = 6 * offset;
  double offset_Dalpha  = 4 * offset;
  double offset_Dchi2cl = 2 * offset;
  double offset_Dsvpv  = 0 * offset;
  double offset_DtrkPt  = -3 * offset;
  double offset_RapGapLoose  = -4 * offset;
  double offset_RapGapTight  = -6 * offset;
  for (int yBin = 0; yBin < nYBins; yBin++) {
    double yBinCenter = 0.5 * (yBins[yBin] + yBins[yBin+1]);
    // Set yield plots
    rawYield_nominal->SetPoint(yBin,
      yBinCenter - offset_nominal, rawYieldVal_nominal[yBin]);
    rawYield_nominal->SetPointError(yBin,
      0.5 - offset_nominal, 0.5 + offset_nominal,
      rawYieldErr_nominal[yBin], rawYieldErr_nominal[yBin]
    );
    rawYield_Dalpha->SetPoint(yBin,
      yBinCenter - offset_Dalpha, rawYieldVal_Dalpha[yBin]);
    rawYield_Dalpha->SetPointError(yBin,
      0.5 - offset_Dalpha, 0.5 + offset_Dalpha,
      rawYieldErr_Dalpha[yBin], rawYieldErr_Dalpha[yBin]
    );
    rawYield_Dchi2cl->SetPoint(yBin,
      yBinCenter - offset_Dchi2cl, rawYieldVal_Dchi2cl[yBin]);
    rawYield_Dchi2cl->SetPointError(yBin,
      0.5 - offset_Dchi2cl, 0.5 + offset_Dchi2cl,
      rawYieldErr_Dchi2cl[yBin], rawYieldErr_Dchi2cl[yBin]
    );
    rawYield_Dsvpv->SetPoint(yBin,
      yBinCenter - offset_Dsvpv, rawYieldVal_Dsvpv[yBin]);
    rawYield_Dsvpv->SetPointError(yBin,
      0.5 - offset_Dsvpv, 0.5 + offset_Dsvpv,
      rawYieldErr_Dsvpv[yBin], rawYieldErr_Dsvpv[yBin]
    );
    rawYield_DtrkPt->SetPoint(yBin,
      yBinCenter - offset_DtrkPt, rawYieldVal_DtrkPt[yBin]);
    rawYield_DtrkPt->SetPointError(yBin,
      0.5 - offset_DtrkPt, 0.5 + offset_DtrkPt,
      rawYieldErr_DtrkPt[yBin], rawYieldErr_DtrkPt[yBin]
    );
    rawYield_RapGapLoose->SetPoint(yBin,
      yBinCenter - offset_RapGapLoose, rawYieldVal_RapGapLoose[yBin]);
    rawYield_RapGapLoose->SetPointError(yBin,
      0.5 - offset_RapGapLoose, 0.5 + offset_RapGapLoose,
      rawYieldErr_RapGapLoose[yBin], rawYieldErr_RapGapLoose[yBin]
    );
    rawYield_RapGapTight->SetPoint(yBin,
      yBinCenter - offset_RapGapTight, rawYieldVal_RapGapTight[yBin]);
    rawYield_RapGapTight->SetPointError(yBin,
      0.5 - offset_RapGapTight, 0.5 + offset_RapGapTight,
      rawYieldErr_RapGapTight[yBin], rawYieldErr_RapGapTight[yBin]
    );
    // Set ratio plots
    double ratioErr_Dalpha = (
      TMath::Abs(rawYieldErr_Dalpha[yBin] * rawYieldVal_nominal[yBin] -
      rawYieldVal_Dalpha[yBin] * rawYieldErr_nominal[yBin]) /
      (rawYieldVal_nominal[yBin] * rawYieldVal_nominal[yBin])
    );
    ratio_Dalpha->SetPoint(yBin,
      yBinCenter - offset_Dalpha,
      rawYieldVal_Dalpha[yBin]/rawYieldVal_nominal[yBin]);
    ratio_Dalpha->SetPointError(yBin,
      0.5 - offset_Dalpha, 0.5 + offset_Dalpha,
      ratioErr_Dalpha, ratioErr_Dalpha
    );
    double ratioErr_Dchi2cl = (
      TMath::Abs(rawYieldErr_Dchi2cl[yBin] * rawYieldVal_nominal[yBin] -
      rawYieldVal_Dchi2cl[yBin] * rawYieldErr_nominal[yBin]) /
      (rawYieldVal_nominal[yBin] * rawYieldVal_nominal[yBin])
    );
    ratio_Dchi2cl->SetPoint(yBin,
      yBinCenter - offset_Dchi2cl,
      rawYieldVal_Dchi2cl[yBin]/rawYieldVal_nominal[yBin]);
    ratio_Dchi2cl->SetPointError(yBin,
      0.5 - offset_Dchi2cl, 0.5 + offset_Dchi2cl,
      ratioErr_Dchi2cl, ratioErr_Dchi2cl
    );
    double ratioErr_Dsvpv = (
      TMath::Abs(rawYieldErr_Dsvpv[yBin] * rawYieldVal_nominal[yBin] -
      rawYieldVal_Dsvpv[yBin] * rawYieldErr_nominal[yBin]) /
      (rawYieldVal_nominal[yBin] * rawYieldVal_nominal[yBin])
    );
    ratio_Dsvpv->SetPoint(yBin,
      yBinCenter - offset_Dsvpv,
      rawYieldVal_Dsvpv[yBin]/rawYieldVal_nominal[yBin]);
    ratio_Dsvpv->SetPointError(yBin,
      0.5 - offset_Dsvpv, 0.5 + offset_Dsvpv,
      ratioErr_Dsvpv, ratioErr_Dsvpv
    );
    double ratioErr_DtrkPt = (
      TMath::Abs(rawYieldErr_DtrkPt[yBin] * rawYieldVal_nominal[yBin] -
      rawYieldVal_DtrkPt[yBin] * rawYieldErr_nominal[yBin]) /
      (rawYieldVal_nominal[yBin] * rawYieldVal_nominal[yBin])
    );
    ratio_DtrkPt->SetPoint(yBin,
      yBinCenter - offset_DtrkPt,
      rawYieldVal_DtrkPt[yBin]/rawYieldVal_nominal[yBin]);
    ratio_DtrkPt->SetPointError(yBin,
      0.5 - offset_DtrkPt, 0.5 + offset_DtrkPt,
      ratioErr_DtrkPt, ratioErr_DtrkPt
    );
    double ratioErr_RapGapLoose = (
      TMath::Abs(rawYieldErr_RapGapLoose[yBin] * rawYieldVal_nominal[yBin] -
      rawYieldVal_RapGapLoose[yBin] * rawYieldErr_nominal[yBin]) /
      (rawYieldVal_nominal[yBin] * rawYieldVal_nominal[yBin])
    );
    ratio_RapGapLoose->SetPoint(yBin,
      yBinCenter - offset_RapGapLoose,
      rawYieldVal_RapGapLoose[yBin]/rawYieldVal_nominal[yBin]);
    ratio_RapGapLoose->SetPointError(yBin,
      0.5 - offset_RapGapLoose, 0.5 + offset_RapGapLoose,
      ratioErr_RapGapLoose, ratioErr_RapGapLoose
    );
    double ratioErr_RapGapTight = (
      TMath::Abs(rawYieldErr_RapGapTight[yBin] * rawYieldVal_nominal[yBin] -
      rawYieldVal_RapGapTight[yBin] * rawYieldErr_nominal[yBin]) /
      (rawYieldVal_nominal[yBin] * rawYieldVal_nominal[yBin])
    );
    ratio_RapGapTight->SetPoint(yBin,
      yBinCenter - offset_RapGapTight,
      rawYieldVal_RapGapTight[yBin]/rawYieldVal_nominal[yBin]);
    ratio_RapGapTight->SetPointError(yBin,
      0.5 - offset_RapGapTight, 0.5 + offset_RapGapTight,
      ratioErr_RapGapTight, ratioErr_RapGapTight
    );
  }
  rawYieldTemplate->SetMinimum(0.);
  rawYieldTemplate->SetMaximum(200.);
  ratioTemplate->SetMinimum(0.9);
  ratioTemplate->SetMaximum(1.4);
  
  // Nominal
  rawYield_nominal->SetLineColor(kNominal);
  rawYield_nominal->SetLineWidth(2);
  // SystFitDalpha
  rawYield_Dalpha->SetLineColor(kDalpha);
  rawYield_Dalpha->SetLineWidth(2);
  ratio_Dalpha->SetLineColor(kDalpha);
  ratio_Dalpha->SetLineWidth(2);
  // SystDchi2cl
  rawYield_Dchi2cl->SetLineColor(kDchi2cl);
  rawYield_Dchi2cl->SetLineWidth(2);
  ratio_Dchi2cl->SetLineColor(kDchi2cl);
  ratio_Dchi2cl->SetLineWidth(2);
  // SystDsvpv
  rawYield_Dsvpv->SetLineColor(kDsvpv);
  rawYield_Dsvpv->SetLineWidth(2);
  ratio_Dsvpv->SetLineColor(kDsvpv);
  ratio_Dsvpv->SetLineWidth(2);
  // SystDtrkPt
  rawYield_DtrkPt->SetLineColor(kDtrkPt);
  rawYield_DtrkPt->SetLineWidth(2);
  ratio_DtrkPt->SetLineColor(kDtrkPt);
  ratio_DtrkPt->SetLineWidth(2);
  // SystRapGapLoose
  rawYield_RapGapLoose->SetLineColor(kRapGapLoose);
  rawYield_RapGapLoose->SetLineWidth(2);
  ratio_RapGapLoose->SetLineColor(kRapGapLoose);
  ratio_RapGapLoose->SetLineWidth(2);
  // SystRapGapTight
  rawYield_RapGapTight->SetLineColor(kRapGapTight);
  rawYield_RapGapTight->SetLineWidth(2);
  ratio_RapGapTight->SetLineColor(kRapGapTight);
  ratio_RapGapTight->SetLineWidth(2);
  
  TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
  TPad* padTop = new TPad("padTop", "", 0.0, 0.3, 1.0, 1.0);
  TPad* padBot = new TPad("padBot", "", 0.0, 0.0, 1.0, 0.3);
  padTop->SetMargin(0.18, 0.02, 0.00, 0.16);
  padBot->SetMargin(0.18, 0.02, 0.30, 0.00);
  padTop->Draw();
  padBot->Draw();
  
  padTop->cd();
  rawYieldTemplate->Draw();
  rawYieldTemplate->SetLabelSize(0.03/0.7, "Y");
  rawYieldTemplate->SetTitleSize(0.035/0.7, "Y");
  rawYield_nominal->Draw("same p");
  rawYield_Dalpha->Draw("same p");
  rawYield_Dchi2cl->Draw("same p");
  rawYield_Dsvpv->Draw("same p");
  rawYield_DtrkPt->Draw("same p");
  rawYield_RapGapLoose->Draw("same p");
  rawYield_RapGapTight->Draw("same p");
  gStyle->SetOptStat(0);
  
  padBot->cd();
  ratioTemplate->Draw();
  ratioTemplate->SetLabelSize(0.03/0.3, "XY");
  ratioTemplate->SetTitleSize(0.035/0.3, "XY");
  ratioTemplate->SetTitleOffset(0.4, "Y");
  TLine* unity = new TLine(yBins[0], 1.0, yBins[nYBins], 1.0);
  unity->SetLineColor(kGray);
  unity->SetLineWidth(1);
  unity->SetLineStyle(9);
  unity->Draw();
  ratio_Dalpha->Draw("same p");
  ratio_Dchi2cl->Draw("same p");
  ratio_Dsvpv->Draw("same p");
  ratio_DtrkPt->Draw("same p");
  ratio_RapGapLoose->Draw("same p");
  ratio_RapGapTight->Draw("same p");
  gStyle->SetOptStat(0);
  
  canvas->cd();
  canvas->Update();
  float legShift = 0.;
  if (isGammaN) legShift = 0.4;
  TLegend* legend = new TLegend(0.2 + legShift, 0.65, 0.6 + legShift, 0.88);
  legend->SetTextSize(0.015/0.7);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(rawYield_nominal,  "Nominal Raw Yield", "l");
  legend->AddEntry(rawYield_Dalpha,  "Cut Syst: Dalpha", "l");
  legend->AddEntry(rawYield_Dchi2cl, "Cut Syst: Dchi2cl", "l");
  legend->AddEntry(rawYield_Dsvpv,  "Cut Syst: Dsvpv", "l");
  legend->AddEntry(rawYield_DtrkPt,  "Cut Syst: DtrkPt", "l");
  legend->AddEntry(rawYield_RapGapLoose,  "Cut Syst: RapGapLoose", "l");
  legend->AddEntry(rawYield_RapGapTight,  "Cut Syst: RapGapTight", "l");
  legend->Draw();
  
  canvas->SaveAs(Form("plot/SystEval_Cuts/%s_rawYield_vs_y.pdf", plotLabel.c_str()));
  
  delete canvas;
  delete rawYieldTemplate;
  delete rawYield_nominal;
  delete rawYield_Dalpha;
  delete rawYield_Dchi2cl;
  delete rawYield_Dsvpv;
  delete rawYield_DtrkPt;
  delete rawYield_RapGapLoose;
  delete rawYield_RapGapTight;
  delete ratioTemplate;
  delete ratio_Dalpha;
  delete ratio_Dchi2cl;
  delete ratio_Dsvpv;
  delete ratio_DtrkPt;
  delete ratio_RapGapLoose;
  delete ratio_RapGapTight;
}

void cuts_effEvt_vs_y(
  vector<string> nominalDirs,
  vector<string> systDalphaDirs,
  vector<string> systDchi2clDirs,
  vector<string> systDsvpvDirs,
  vector<string> systDtrkPtDirs,
  vector<string> systRapGapLooseDirs,
  vector<string> systRapGapTightDirs,
  string plotTitle,
  string plotLabel,
  int nYBins,
  double* yBins,
  double ptMin,
  double ptMax,
  int isGammaN
) {
  vector<string> input_nominal;
  vector<string> input_Dalpha;
  vector<string> input_Dchi2cl;
  vector<string> input_Dsvpv;
  vector<string> input_DtrkPt;
  vector<string> input_RapGapLoose;
  vector<string> input_RapGapTight;
  
  for (int i = 0; i < nominalDirs.size(); i++) {
    input_nominal.push_back(Form("%s/MassFit/correctedYields.md", nominalDirs[i].c_str()));
    input_Dalpha.push_back(Form("%s/MassFit/correctedYields.md", systDalphaDirs[i].c_str()));
    input_Dchi2cl.push_back(Form("%s/MassFit/correctedYields.md", systDchi2clDirs[i].c_str()));
    input_Dsvpv.push_back(Form("%s/MassFit/correctedYields.md", systDsvpvDirs[i].c_str()));
    input_DtrkPt.push_back(Form("%s/MassFit/correctedYields.md", systDtrkPtDirs[i].c_str()));
    input_RapGapLoose.push_back(Form("%s/MassFit/correctedYields.md", systRapGapLooseDirs[i].c_str()));
    input_RapGapTight.push_back(Form("%s/MassFit/correctedYields.md", systRapGapTightDirs[i].c_str()));
  }
  
  vector<Point> points_nominal  = getPointArr(ptMin, ptMax, isGammaN, input_nominal);
  vector<Point> points_Dalpha  = getPointArr(ptMin, ptMax, isGammaN, input_Dalpha);
  vector<Point> points_Dchi2cl = getPointArr(ptMin, ptMax, isGammaN, input_Dchi2cl);
  vector<Point> points_Dsvpv  = getPointArr(ptMin, ptMax, isGammaN, input_Dsvpv);
  vector<Point> points_DtrkPt  = getPointArr(ptMin, ptMax, isGammaN, input_DtrkPt);
  vector<Point> points_RapGapLoose  = getPointArr(ptMin, ptMax, isGammaN, input_RapGapLoose);
  vector<Point> points_RapGapTight  = getPointArr(ptMin, ptMax, isGammaN, input_RapGapTight);
  
  vector<double> effEvtVal_nominal  = getDoubleArr(points_nominal, [](Point& p) -> double { return p.effEvt;});
  vector<double> effEvtErr_nominal  = getDoubleArr(points_nominal, [](Point& p) -> double { return p.effEvtError;});
  vector<double> effEvtVal_Dalpha  = getDoubleArr(points_Dalpha, [](Point& p) -> double { return p.effEvt;});
  vector<double> effEvtErr_Dalpha  = getDoubleArr(points_Dalpha, [](Point& p) -> double { return p.effEvtError;});
  vector<double> effEvtVal_Dchi2cl = getDoubleArr(points_Dchi2cl, [](Point& p) -> double { return p.effEvt;});
  vector<double> effEvtErr_Dchi2cl = getDoubleArr(points_Dchi2cl, [](Point& p) -> double { return p.effEvtError;});
  vector<double> effEvtVal_Dsvpv  = getDoubleArr(points_Dsvpv, [](Point& p) -> double { return p.effEvt;});
  vector<double> effEvtErr_Dsvpv  = getDoubleArr(points_Dsvpv, [](Point& p) -> double { return p.effEvtError;});
  vector<double> effEvtVal_DtrkPt  = getDoubleArr(points_DtrkPt, [](Point& p) -> double { return p.effEvt;});
  vector<double> effEvtErr_DtrkPt  = getDoubleArr(points_DtrkPt, [](Point& p) -> double { return p.effEvtError;});
  vector<double> effEvtVal_RapGapLoose = getDoubleArr(points_RapGapLoose, [](Point& p) -> double { return p.effEvt;});
  vector<double> effEvtErr_RapGapLoose = getDoubleArr(points_RapGapLoose, [](Point& p) -> double { return p.effEvtError;});
  vector<double> effEvtVal_RapGapTight = getDoubleArr(points_RapGapLoose, [](Point& p) -> double { return p.effEvt;});
  vector<double> effEvtErr_RapGapTight = getDoubleArr(points_RapGapLoose, [](Point& p) -> double { return p.effEvtError;});

  TH1D* effEvtTemplate  = new TH1D(
    "effEvtTemplate", (plotTitle + "; y; Raw Yield").c_str(), nYBins, yBins);
  TGraphAsymmErrors* effEvt_nominal  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* effEvt_Dalpha  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* effEvt_Dchi2cl = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* effEvt_Dsvpv  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* effEvt_DtrkPt  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* effEvt_RapGapLoose  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* effEvt_RapGapTight  = new TGraphAsymmErrors(nYBins);

  TH1D* ratioTemplate  = new TH1D(
    "ratioTemplate", "; y; syst/nominal", nYBins, yBins);
  TGraphAsymmErrors* ratio_Dalpha  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_Dchi2cl = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_Dsvpv  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_DtrkPt  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_RapGapLoose  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_RapGapTight  = new TGraphAsymmErrors(nYBins);
    
  double offset = 0.02;
  double offset_nominal  = 6 * offset;
  double offset_Dalpha  = 4 * offset;
  double offset_Dchi2cl = 2 * offset;
  double offset_Dsvpv  = 0 * offset;
  double offset_DtrkPt  = -3 * offset;
  double offset_RapGapLoose  = -4 * offset;
  double offset_RapGapTight  = -6 * offset;
  for (int yBin = 0; yBin < nYBins; yBin++) {
    double yBinCenter = 0.5 * (yBins[yBin] + yBins[yBin+1]);
    // Set yield plots
    effEvt_nominal->SetPoint(yBin,
      yBinCenter - offset_nominal, effEvtVal_nominal[yBin]);
    effEvt_nominal->SetPointError(yBin,
      0.5 - offset_nominal, 0.5 + offset_nominal,
      effEvtErr_nominal[yBin], effEvtErr_nominal[yBin]
    );
    effEvt_Dalpha->SetPoint(yBin,
      yBinCenter - offset_Dalpha, effEvtVal_Dalpha[yBin]);
    effEvt_Dalpha->SetPointError(yBin,
      0.5 - offset_Dalpha, 0.5 + offset_Dalpha,
      effEvtErr_Dalpha[yBin], effEvtErr_Dalpha[yBin]
    );
    effEvt_Dchi2cl->SetPoint(yBin,
      yBinCenter - offset_Dchi2cl, effEvtVal_Dchi2cl[yBin]);
    effEvt_Dchi2cl->SetPointError(yBin,
      0.5 - offset_Dchi2cl, 0.5 + offset_Dchi2cl,
      effEvtErr_Dchi2cl[yBin], effEvtErr_Dchi2cl[yBin]
    );
    effEvt_Dsvpv->SetPoint(yBin,
      yBinCenter - offset_Dsvpv, effEvtVal_Dsvpv[yBin]);
    effEvt_Dsvpv->SetPointError(yBin,
      0.5 - offset_Dsvpv, 0.5 + offset_Dsvpv,
      effEvtErr_Dsvpv[yBin], effEvtErr_Dsvpv[yBin]
    );
    effEvt_DtrkPt->SetPoint(yBin,
      yBinCenter - offset_DtrkPt, effEvtVal_DtrkPt[yBin]);
    effEvt_DtrkPt->SetPointError(yBin,
      0.5 - offset_DtrkPt, 0.5 + offset_DtrkPt,
      effEvtErr_DtrkPt[yBin], effEvtErr_DtrkPt[yBin]
    );
    effEvt_RapGapLoose->SetPoint(yBin,
      yBinCenter - offset_RapGapLoose, effEvtVal_RapGapLoose[yBin]);
    effEvt_RapGapLoose->SetPointError(yBin,
      0.5 - offset_RapGapLoose, 0.5 + offset_RapGapLoose,
      effEvtErr_RapGapLoose[yBin], effEvtErr_RapGapLoose[yBin]
    );
    effEvt_RapGapTight->SetPoint(yBin,
      yBinCenter - offset_RapGapTight, effEvtVal_RapGapTight[yBin]);
    effEvt_RapGapTight->SetPointError(yBin,
      0.5 - offset_RapGapTight, 0.5 + offset_RapGapTight,
      effEvtErr_RapGapTight[yBin], effEvtErr_RapGapTight[yBin]
    );
    // Set ratio plots
    double ratioErr_Dalpha = (
      TMath::Abs(effEvtErr_Dalpha[yBin] * effEvtVal_nominal[yBin] -
      effEvtVal_Dalpha[yBin] * effEvtErr_nominal[yBin]) /
      (effEvtVal_nominal[yBin] * effEvtVal_nominal[yBin])
    );
    ratio_Dalpha->SetPoint(yBin,
      yBinCenter - offset_Dalpha,
      effEvtVal_Dalpha[yBin]/effEvtVal_nominal[yBin]);
    ratio_Dalpha->SetPointError(yBin,
      0.5 - offset_Dalpha, 0.5 + offset_Dalpha,
      effEvtErr_Dalpha[yBin], effEvtErr_Dalpha[yBin]
    );
    double ratioErr_Dchi2cl = (
      TMath::Abs(effEvtErr_Dchi2cl[yBin] * effEvtVal_nominal[yBin] -
      effEvtVal_Dchi2cl[yBin] * effEvtErr_nominal[yBin]) /
      (effEvtVal_nominal[yBin] * effEvtVal_nominal[yBin])
    );
    ratio_Dchi2cl->SetPoint(yBin,
      yBinCenter - offset_Dchi2cl,
      effEvtVal_Dchi2cl[yBin]/effEvtVal_nominal[yBin]);
    ratio_Dchi2cl->SetPointError(yBin,
      0.5 - offset_Dchi2cl, 0.5 + offset_Dchi2cl,
      effEvtErr_Dchi2cl[yBin], effEvtErr_Dchi2cl[yBin]
    );
    double ratioErr_Dsvpv = (
      TMath::Abs(effEvtErr_Dsvpv[yBin] * effEvtVal_nominal[yBin] -
      effEvtVal_Dsvpv[yBin] * effEvtErr_nominal[yBin]) /
      (effEvtVal_nominal[yBin] * effEvtVal_nominal[yBin])
    );
    ratio_Dsvpv->SetPoint(yBin,
      yBinCenter - offset_Dsvpv,
      effEvtVal_Dsvpv[yBin]/effEvtVal_nominal[yBin]);
    ratio_Dsvpv->SetPointError(yBin,
      0.5 - offset_Dsvpv, 0.5 + offset_Dsvpv,
      effEvtErr_Dsvpv[yBin], effEvtErr_Dsvpv[yBin]
    );
    double ratioErr_DtrkPt = (
      TMath::Abs(effEvtErr_DtrkPt[yBin] * effEvtVal_nominal[yBin] -
      effEvtVal_DtrkPt[yBin] * effEvtErr_nominal[yBin]) /
      (effEvtVal_nominal[yBin] * effEvtVal_nominal[yBin])
    );
    ratio_DtrkPt->SetPoint(yBin,
      yBinCenter - offset_DtrkPt,
      effEvtVal_DtrkPt[yBin]/effEvtVal_nominal[yBin]);
    ratio_DtrkPt->SetPointError(yBin,
      0.5 - offset_DtrkPt, 0.5 + offset_DtrkPt,
      effEvtErr_DtrkPt[yBin], effEvtErr_DtrkPt[yBin]
    );
    double ratioErr_RapGapLoose = (
      TMath::Abs(effEvtErr_RapGapLoose[yBin] * effEvtVal_nominal[yBin] -
      effEvtVal_RapGapLoose[yBin] * effEvtErr_nominal[yBin]) /
      (effEvtVal_nominal[yBin] * effEvtVal_nominal[yBin])
    );
    ratio_RapGapLoose->SetPoint(yBin,
      yBinCenter - offset_RapGapLoose,
      effEvtVal_RapGapLoose[yBin]/effEvtVal_nominal[yBin]);
    ratio_RapGapLoose->SetPointError(yBin,
      0.5 - offset_RapGapLoose, 0.5 + offset_RapGapLoose,
      ratioErr_RapGapLoose, ratioErr_RapGapLoose
    );
    double ratioErr_RapGapTight = (
      TMath::Abs(effEvtErr_RapGapTight[yBin] * effEvtVal_nominal[yBin] -
      effEvtVal_RapGapTight[yBin] * effEvtErr_nominal[yBin]) /
      (effEvtVal_nominal[yBin] * effEvtVal_nominal[yBin])
    );
    ratio_RapGapTight->SetPoint(yBin,
      yBinCenter - offset_RapGapTight,
      effEvtVal_RapGapTight[yBin]/effEvtVal_nominal[yBin]);
    ratio_RapGapTight->SetPointError(yBin,
      0.5 - offset_RapGapTight, 0.5 + offset_RapGapTight,
      ratioErr_RapGapTight, ratioErr_RapGapTight
    );
  }
  effEvtTemplate->SetMinimum(0.95);
  effEvtTemplate->SetMaximum(1.05);
  ratioTemplate->SetMinimum(0.95);
  ratioTemplate->SetMaximum(1.05);
  
  // Nominal
  effEvt_nominal->SetLineColor(kNominal);
  effEvt_nominal->SetLineWidth(2);
  // SystFitDalpha
  effEvt_Dalpha->SetLineColor(kDalpha);
  effEvt_Dalpha->SetLineWidth(2);
  ratio_Dalpha->SetLineColor(kDalpha);
  ratio_Dalpha->SetLineWidth(2);
  // SystDchi2cl
  effEvt_Dchi2cl->SetLineColor(kDchi2cl);
  effEvt_Dchi2cl->SetLineWidth(2);
  ratio_Dchi2cl->SetLineColor(kDchi2cl);
  ratio_Dchi2cl->SetLineWidth(2);
  // SystDsvpv
  effEvt_Dsvpv->SetLineColor(kDsvpv);
  effEvt_Dsvpv->SetLineWidth(2);
  ratio_Dsvpv->SetLineColor(kDsvpv);
  ratio_Dsvpv->SetLineWidth(2);
  // SystDtrkPt
  effEvt_DtrkPt->SetLineColor(kDtrkPt);
  effEvt_DtrkPt->SetLineWidth(2);
  ratio_DtrkPt->SetLineColor(kDtrkPt);
  ratio_DtrkPt->SetLineWidth(2);
  // SystRapGapLoose
  effEvt_RapGapLoose->SetLineColor(kRapGapLoose);
  effEvt_RapGapLoose->SetLineWidth(2);
  ratio_RapGapLoose->SetLineColor(kRapGapLoose);
  ratio_RapGapLoose->SetLineWidth(2);
  // SystRapGapTight
  effEvt_RapGapTight->SetLineColor(kRapGapTight);
  effEvt_RapGapTight->SetLineWidth(2);
  ratio_RapGapTight->SetLineColor(kRapGapTight);
  ratio_RapGapTight->SetLineWidth(2);
  
  TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
  TPad* padTop = new TPad("padTop", "", 0.0, 0.3, 1.0, 1.0);
  TPad* padBot = new TPad("padBot", "", 0.0, 0.0, 1.0, 0.3);
  padTop->SetMargin(0.18, 0.02, 0.00, 0.16);
  padBot->SetMargin(0.18, 0.02, 0.30, 0.00);
  padTop->Draw();
  padBot->Draw();
  
  padTop->cd();
  effEvtTemplate->Draw();
  effEvtTemplate->SetLabelSize(0.03/0.7, "Y");
  effEvtTemplate->SetTitleSize(0.035/0.7, "Y");
  effEvt_nominal->Draw("same p");
  effEvt_Dalpha->Draw("same p");
  effEvt_Dchi2cl->Draw("same p");
  effEvt_Dsvpv->Draw("same p");
  effEvt_DtrkPt->Draw("same p");
  effEvt_RapGapLoose->Draw("same p");
  effEvt_RapGapTight->Draw("same p");
  gStyle->SetOptStat(0);
  
  padBot->cd();
  ratioTemplate->Draw();
  ratioTemplate->SetLabelSize(0.03/0.3, "XY");
  ratioTemplate->SetTitleSize(0.035/0.3, "XY");
  ratioTemplate->SetTitleOffset(0.4, "Y");
  TLine* unity = new TLine(yBins[0], 1.0, yBins[nYBins], 1.0);
  unity->SetLineColor(kGray);
  unity->SetLineWidth(1);
  unity->SetLineStyle(9);
  unity->Draw();
  ratio_Dalpha->Draw("same p");
  ratio_Dchi2cl->Draw("same p");
  ratio_Dsvpv->Draw("same p");
  ratio_DtrkPt->Draw("same p");
  ratio_RapGapLoose->Draw("same p");
  ratio_RapGapTight->Draw("same p");
  gStyle->SetOptStat(0);
  
  canvas->cd();
  canvas->Update();
  float legShift = 0.;
  if (isGammaN) legShift = 0.4;
  TLegend* legend = new TLegend(0.2 + legShift, 0.65, 0.6 + legShift, 0.88);
  legend->SetTextSize(0.015/0.7);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(effEvt_nominal,  "Nominal Event Efficiency", "l");
  legend->AddEntry(effEvt_Dalpha,  "Cut Syst: Dalpha", "l");
  legend->AddEntry(effEvt_Dchi2cl, "Cut Syst: Dchi2cl", "l");
  legend->AddEntry(effEvt_Dsvpv,  "Cut Syst: Dsvpv", "l");
  legend->AddEntry(effEvt_DtrkPt,  "Cut Syst: DtrkPt", "l");
  legend->AddEntry(effEvt_RapGapLoose,  "Cut Syst: RapGapLoose", "l");
  legend->AddEntry(effEvt_RapGapTight,  "Cut Syst: RapGapTight", "l");
  legend->Draw();
  
  canvas->SaveAs(Form("plot/SystEval_Cuts/%s_effEvt_vs_y.pdf", plotLabel.c_str()));
  
  delete canvas;
  delete effEvtTemplate;
  delete effEvt_nominal;
  delete effEvt_Dalpha;
  delete effEvt_Dchi2cl;
  delete effEvt_Dsvpv;
  delete effEvt_DtrkPt;
  delete effEvt_RapGapLoose;
  delete effEvt_RapGapTight;
  delete ratioTemplate;
  delete ratio_Dalpha;
  delete ratio_Dchi2cl;
  delete ratio_Dsvpv;
  delete ratio_DtrkPt;
  delete ratio_RapGapLoose;
  delete ratio_RapGapTight;
}

void cuts_effD_vs_y(
  vector<string> nominalDirs,
  vector<string> systDalphaDirs,
  vector<string> systDchi2clDirs,
  vector<string> systDsvpvDirs,
  vector<string> systDtrkPtDirs,
  vector<string> systRapGapLooseDirs,
  vector<string> systRapGapTightDirs,
  string plotTitle,
  string plotLabel,
  int nYBins,
  double* yBins,
  double ptMin,
  double ptMax,
  int isGammaN
) {
  vector<string> input_nominal;
  vector<string> input_Dalpha;
  vector<string> input_Dchi2cl;
  vector<string> input_Dsvpv;
  vector<string> input_DtrkPt;
  vector<string> input_RapGapLoose;
  vector<string> input_RapGapTight;
  
  for (int i = 0; i < nominalDirs.size(); i++) {
    input_nominal.push_back(Form("%s/MassFit/correctedYields.md", nominalDirs[i].c_str()));
    input_Dalpha.push_back(Form("%s/MassFit/correctedYields.md", systDalphaDirs[i].c_str()));
    input_Dchi2cl.push_back(Form("%s/MassFit/correctedYields.md", systDchi2clDirs[i].c_str()));
    input_Dsvpv.push_back(Form("%s/MassFit/correctedYields.md", systDsvpvDirs[i].c_str()));
    input_DtrkPt.push_back(Form("%s/MassFit/correctedYields.md", systDtrkPtDirs[i].c_str()));
    input_RapGapLoose.push_back(Form("%s/MassFit/correctedYields.md", systRapGapLooseDirs[i].c_str()));
    input_RapGapTight.push_back(Form("%s/MassFit/correctedYields.md", systRapGapTightDirs[i].c_str()));
  }
  
  vector<Point> points_nominal  = getPointArr(ptMin, ptMax, isGammaN, input_nominal);
  vector<Point> points_Dalpha  = getPointArr(ptMin, ptMax, isGammaN, input_Dalpha);
  vector<Point> points_Dchi2cl = getPointArr(ptMin, ptMax, isGammaN, input_Dchi2cl);
  vector<Point> points_Dsvpv  = getPointArr(ptMin, ptMax, isGammaN, input_Dsvpv);
  vector<Point> points_DtrkPt  = getPointArr(ptMin, ptMax, isGammaN, input_DtrkPt);
  vector<Point> points_RapGapLoose  = getPointArr(ptMin, ptMax, isGammaN, input_RapGapLoose);
  vector<Point> points_RapGapTight  = getPointArr(ptMin, ptMax, isGammaN, input_RapGapTight);
  
  vector<double> effDVal_nominal  = getDoubleArr(points_nominal, [](Point& p) -> double { return p.effD;});
  vector<double> effDErr_nominal  = getDoubleArr(points_nominal, [](Point& p) -> double { return p.effDError;});
  vector<double> effDVal_Dalpha  = getDoubleArr(points_Dalpha, [](Point& p) -> double { return p.effD;});
  vector<double> effDErr_Dalpha  = getDoubleArr(points_Dalpha, [](Point& p) -> double { return p.effDError;});
  vector<double> effDVal_Dchi2cl = getDoubleArr(points_Dchi2cl, [](Point& p) -> double { return p.effD;});
  vector<double> effDErr_Dchi2cl = getDoubleArr(points_Dchi2cl, [](Point& p) -> double { return p.effDError;});
  vector<double> effDVal_Dsvpv  = getDoubleArr(points_Dsvpv, [](Point& p) -> double { return p.effD;});
  vector<double> effDErr_Dsvpv  = getDoubleArr(points_Dsvpv, [](Point& p) -> double { return p.effDError;});
  vector<double> effDVal_DtrkPt  = getDoubleArr(points_DtrkPt, [](Point& p) -> double { return p.effD;});
  vector<double> effDErr_DtrkPt  = getDoubleArr(points_DtrkPt, [](Point& p) -> double { return p.effDError;});
  vector<double> effDVal_RapGapLoose = getDoubleArr(points_RapGapLoose, [](Point& p) -> double { return p.effD;});
  vector<double> effDErr_RapGapLoose = getDoubleArr(points_RapGapLoose, [](Point& p) -> double { return p.effDError;});
  vector<double> effDVal_RapGapTight = getDoubleArr(points_RapGapLoose, [](Point& p) -> double { return p.effD;});
  vector<double> effDErr_RapGapTight = getDoubleArr(points_RapGapLoose, [](Point& p) -> double { return p.effDError;});

  TH1D* effDTemplate  = new TH1D(
    "effDTemplate", (plotTitle + "; y; Raw Yield").c_str(), nYBins, yBins);
  TGraphAsymmErrors* effD_nominal  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* effD_Dalpha  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* effD_Dchi2cl = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* effD_Dsvpv  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* effD_DtrkPt  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* effD_RapGapLoose  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* effD_RapGapTight  = new TGraphAsymmErrors(nYBins);

  TH1D* ratioTemplate  = new TH1D(
    "ratioTemplate", "; y; syst/nominal", nYBins, yBins);
  TGraphAsymmErrors* ratio_Dalpha  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_Dchi2cl = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_Dsvpv  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_DtrkPt  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_RapGapLoose  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_RapGapTight  = new TGraphAsymmErrors(nYBins);
    
  double offset = 0.02;
  double offset_nominal  = 6 * offset;
  double offset_Dalpha  = 4 * offset;
  double offset_Dchi2cl = 2 * offset;
  double offset_Dsvpv  = 0 * offset;
  double offset_DtrkPt  = -3 * offset;
  double offset_RapGapLoose  = -4 * offset;
  double offset_RapGapTight  = -6 * offset;
  for (int yBin = 0; yBin < nYBins; yBin++) {
    double yBinCenter = 0.5 * (yBins[yBin] + yBins[yBin+1]);
    // Set yield plots
    effD_nominal->SetPoint(yBin,
      yBinCenter - offset_nominal, effDVal_nominal[yBin]);
    effD_nominal->SetPointError(yBin,
      0.5 - offset_nominal, 0.5 + offset_nominal,
      effDErr_nominal[yBin], effDErr_nominal[yBin]
    );
    effD_Dalpha->SetPoint(yBin,
      yBinCenter - offset_Dalpha, effDVal_Dalpha[yBin]);
    effD_Dalpha->SetPointError(yBin,
      0.5 - offset_Dalpha, 0.5 + offset_Dalpha,
      effDErr_Dalpha[yBin], effDErr_Dalpha[yBin]
    );
    effD_Dchi2cl->SetPoint(yBin,
      yBinCenter - offset_Dchi2cl, effDVal_Dchi2cl[yBin]);
    effD_Dchi2cl->SetPointError(yBin,
      0.5 - offset_Dchi2cl, 0.5 + offset_Dchi2cl,
      effDErr_Dchi2cl[yBin], effDErr_Dchi2cl[yBin]
    );
    effD_Dsvpv->SetPoint(yBin,
      yBinCenter - offset_Dsvpv, effDVal_Dsvpv[yBin]);
    effD_Dsvpv->SetPointError(yBin,
      0.5 - offset_Dsvpv, 0.5 + offset_Dsvpv,
      effDErr_Dsvpv[yBin], effDErr_Dsvpv[yBin]
    );
    effD_DtrkPt->SetPoint(yBin,
      yBinCenter - offset_DtrkPt, effDVal_DtrkPt[yBin]);
    effD_DtrkPt->SetPointError(yBin,
      0.5 - offset_DtrkPt, 0.5 + offset_DtrkPt,
      effDErr_DtrkPt[yBin], effDErr_DtrkPt[yBin]
    );
    effD_RapGapLoose->SetPoint(yBin,
      yBinCenter - offset_RapGapLoose, effDVal_RapGapLoose[yBin]);
    effD_RapGapLoose->SetPointError(yBin,
      0.5 - offset_RapGapLoose, 0.5 + offset_RapGapLoose,
      effDErr_RapGapLoose[yBin], effDErr_RapGapLoose[yBin]
    );
    effD_RapGapTight->SetPoint(yBin,
      yBinCenter - offset_RapGapTight, effDVal_RapGapTight[yBin]);
    effD_RapGapTight->SetPointError(yBin,
      0.5 - offset_RapGapTight, 0.5 + offset_RapGapTight,
      effDErr_RapGapTight[yBin], effDErr_RapGapTight[yBin]
    );
    // Set ratio plots
    double ratioErr_Dalpha = (
      TMath::Abs(effDErr_Dalpha[yBin] * effDVal_nominal[yBin] -
      effDVal_Dalpha[yBin] * effDErr_nominal[yBin]) /
      (effDVal_nominal[yBin] * effDVal_nominal[yBin])
    );
    ratio_Dalpha->SetPoint(yBin,
      yBinCenter - offset_Dalpha,
      effDVal_Dalpha[yBin]/effDVal_nominal[yBin]);
    ratio_Dalpha->SetPointError(yBin,
      0.5 - offset_Dalpha, 0.5 + offset_Dalpha,
      ratioErr_Dalpha, ratioErr_Dalpha
    );
    double ratioErr_Dchi2cl = (
      TMath::Abs(effDErr_Dchi2cl[yBin] * effDVal_nominal[yBin] -
      effDVal_Dchi2cl[yBin] * effDErr_nominal[yBin]) /
      (effDVal_nominal[yBin] * effDVal_nominal[yBin])
    );
    ratio_Dchi2cl->SetPoint(yBin,
      yBinCenter - offset_Dchi2cl,
      effDVal_Dchi2cl[yBin]/effDVal_nominal[yBin]);
    ratio_Dchi2cl->SetPointError(yBin,
      0.5 - offset_Dchi2cl, 0.5 + offset_Dchi2cl,
      ratioErr_Dchi2cl, ratioErr_Dchi2cl
    );
    double ratioErr_Dsvpv = (
      TMath::Abs(effDErr_Dsvpv[yBin] * effDVal_nominal[yBin] -
      effDVal_Dsvpv[yBin] * effDErr_nominal[yBin]) /
      (effDVal_nominal[yBin] * effDVal_nominal[yBin])
    );
    ratio_Dsvpv->SetPoint(yBin,
      yBinCenter - offset_Dsvpv,
      effDVal_Dsvpv[yBin]/effDVal_nominal[yBin]);
    ratio_Dsvpv->SetPointError(yBin,
      0.5 - offset_Dsvpv, 0.5 + offset_Dsvpv,
      ratioErr_Dsvpv, ratioErr_Dsvpv
    );
    double ratioErr_DtrkPt = (
      TMath::Abs(effDErr_DtrkPt[yBin] * effDVal_nominal[yBin] -
      effDVal_DtrkPt[yBin] * effDErr_nominal[yBin]) /
      (effDVal_nominal[yBin] * effDVal_nominal[yBin])
    );
    ratio_DtrkPt->SetPoint(yBin,
      yBinCenter - offset_DtrkPt,
      effDVal_DtrkPt[yBin]/effDVal_nominal[yBin]);
    ratio_DtrkPt->SetPointError(yBin,
      0.5 - offset_DtrkPt, 0.5 + offset_DtrkPt,
      ratioErr_DtrkPt, ratioErr_DtrkPt
    );
    double ratioErr_RapGapLoose = (
      TMath::Abs(effDErr_RapGapLoose[yBin] * effDVal_nominal[yBin] -
      effDVal_RapGapLoose[yBin] * effDErr_nominal[yBin]) /
      (effDVal_nominal[yBin] * effDVal_nominal[yBin])
    );
    ratio_RapGapLoose->SetPoint(yBin,
      yBinCenter - offset_RapGapLoose,
      effDVal_RapGapLoose[yBin]/effDVal_nominal[yBin]);
    ratio_RapGapLoose->SetPointError(yBin,
      0.5 - offset_RapGapLoose, 0.5 + offset_RapGapLoose,
      ratioErr_RapGapLoose, ratioErr_RapGapLoose
    );
    double ratioErr_RapGapTight = (
      TMath::Abs(effDErr_RapGapTight[yBin] * effDVal_nominal[yBin] -
      effDVal_RapGapTight[yBin] * effDErr_nominal[yBin]) /
      (effDVal_nominal[yBin] * effDVal_nominal[yBin])
    );
    ratio_RapGapTight->SetPoint(yBin,
      yBinCenter - offset_RapGapTight,
      effDVal_RapGapTight[yBin]/effDVal_nominal[yBin]);
    ratio_RapGapTight->SetPointError(yBin,
      0.5 - offset_RapGapTight, 0.5 + offset_RapGapTight,
      ratioErr_RapGapTight, ratioErr_RapGapTight
    );
  }
  effDTemplate->SetMinimum(0.);
  effDTemplate->SetMaximum(0.2);
  ratioTemplate->SetMinimum(0.9);
  ratioTemplate->SetMaximum(1.5);
  
  // Nominal
  effD_nominal->SetLineColor(kNominal);
  effD_nominal->SetLineWidth(2);
  // SystFitDalpha
  effD_Dalpha->SetLineColor(kDalpha);
  effD_Dalpha->SetLineWidth(2);
  ratio_Dalpha->SetLineColor(kDalpha);
  ratio_Dalpha->SetLineWidth(2);
  // SystDchi2cl
  effD_Dchi2cl->SetLineColor(kDchi2cl);
  effD_Dchi2cl->SetLineWidth(2);
  ratio_Dchi2cl->SetLineColor(kDchi2cl);
  ratio_Dchi2cl->SetLineWidth(2);
  // SystDsvpv
  effD_Dsvpv->SetLineColor(kDsvpv);
  effD_Dsvpv->SetLineWidth(2);
  ratio_Dsvpv->SetLineColor(kDsvpv);
  ratio_Dsvpv->SetLineWidth(2);
  // SystDtrkPt
  effD_DtrkPt->SetLineColor(kDtrkPt);
  effD_DtrkPt->SetLineWidth(2);
  ratio_DtrkPt->SetLineColor(kDtrkPt);
  ratio_DtrkPt->SetLineWidth(2);
  // SystRapGapLoose
  effD_RapGapLoose->SetLineColor(kRapGapLoose);
  effD_RapGapLoose->SetLineWidth(2);
  ratio_RapGapLoose->SetLineColor(kRapGapLoose);
  ratio_RapGapLoose->SetLineWidth(2);
  // SystRapGapTight
  effD_RapGapTight->SetLineColor(kRapGapTight);
  effD_RapGapTight->SetLineWidth(2);
  ratio_RapGapTight->SetLineColor(kRapGapTight);
  ratio_RapGapTight->SetLineWidth(2);
  
  TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
  TPad* padTop = new TPad("padTop", "", 0.0, 0.3, 1.0, 1.0);
  TPad* padBot = new TPad("padBot", "", 0.0, 0.0, 1.0, 0.3);
  padTop->SetMargin(0.18, 0.02, 0.00, 0.16);
  padBot->SetMargin(0.18, 0.02, 0.30, 0.00);
  padTop->Draw();
  padBot->Draw();
  
  padTop->cd();
  effDTemplate->Draw();
  effDTemplate->SetLabelSize(0.03/0.7, "Y");
  effDTemplate->SetTitleSize(0.035/0.7, "Y");
  effD_nominal->Draw("same p");
  effD_Dalpha->Draw("same p");
  effD_Dchi2cl->Draw("same p");
  effD_Dsvpv->Draw("same p");
  effD_DtrkPt->Draw("same p");
  effD_RapGapLoose->Draw("same p");
  effD_RapGapTight->Draw("same p");
  gStyle->SetOptStat(0);
  
  padBot->cd();
  ratioTemplate->Draw();
  ratioTemplate->SetLabelSize(0.03/0.3, "XY");
  ratioTemplate->SetTitleSize(0.035/0.3, "XY");
  ratioTemplate->SetTitleOffset(0.4, "Y");
  TLine* unity = new TLine(yBins[0], 1.0, yBins[nYBins], 1.0);
  unity->SetLineColor(kGray);
  unity->SetLineWidth(1);
  unity->SetLineStyle(9);
  unity->Draw();
  ratio_Dalpha->Draw("same p");
  ratio_Dchi2cl->Draw("same p");
  ratio_Dsvpv->Draw("same p");
  ratio_DtrkPt->Draw("same p");
  ratio_RapGapLoose->Draw("same p");
  ratio_RapGapTight->Draw("same p");
  gStyle->SetOptStat(0);
  
  canvas->cd();
  canvas->Update();
  float legShift = 0.;
  if (isGammaN) legShift = 0.4;
  TLegend* legend = new TLegend(0.2 + legShift, 0.65, 0.6 + legShift, 0.88);
  legend->SetTextSize(0.015/0.7);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(effD_nominal,  "Nominal D^{0} Efficiency", "l");
  legend->AddEntry(effD_Dalpha,  "Cut Syst: Dalpha", "l");
  legend->AddEntry(effD_Dchi2cl, "Cut Syst: Dchi2cl", "l");
  legend->AddEntry(effD_Dsvpv,  "Cut Syst: Dsvpv", "l");
  legend->AddEntry(effD_DtrkPt,  "Cut Syst: DtrkPt", "l");
  legend->AddEntry(effD_RapGapLoose,  "Cut Syst: RapGapLoose", "l");
  legend->AddEntry(effD_RapGapTight,  "Cut Syst: RapGapTight", "l");
  legend->Draw();
  
  canvas->SaveAs(Form("plot/SystEval_Cuts/%s_effD_vs_y.pdf", plotLabel.c_str()));
  
  delete canvas;
  delete effDTemplate;
  delete effD_nominal;
  delete effD_Dalpha;
  delete effD_Dchi2cl;
  delete effD_Dsvpv;
  delete effD_DtrkPt;
  delete effD_RapGapLoose;
  delete effD_RapGapTight;
  delete ratioTemplate;
  delete ratio_Dalpha;
  delete ratio_Dchi2cl;
  delete ratio_Dsvpv;
  delete ratio_DtrkPt;
  delete ratio_RapGapLoose;
  delete ratio_RapGapTight;
}

void cuts_lambda_vs_y(
  vector<string> nominalDirs,
  vector<string> systDalphaDirs,
  vector<string> systDchi2clDirs,
  vector<string> systDsvpvDirs,
  vector<string> systDtrkPtDirs,
  vector<string> systRapGapLooseDirs,
  vector<string> systRapGapTightDirs,
  string plotTitle,
  string plotLabel,
  int nYBins,
  double* yBins,
  double ptMin,
  double ptMax,
  int isGammaN
) {
  vector<string> input_nominal;
  vector<string> input_Dalpha;
  vector<string> input_Dchi2cl;
  vector<string> input_Dsvpv;
  vector<string> input_DtrkPt;
  vector<string> input_RapGapLoose;
  vector<string> input_RapGapTight;
  
  for (int i = 0; i < nominalDirs.size(); i++) {
    input_nominal.push_back(Form("%s/MassFit/combdata.dat", nominalDirs[i].c_str()));
    input_Dalpha.push_back(Form("%s/MassFit/combdata.dat", systDalphaDirs[i].c_str()));
    input_Dchi2cl.push_back(Form("%s/MassFit/combdata.dat", systDchi2clDirs[i].c_str()));
    input_Dsvpv.push_back(Form("%s/MassFit/combdata.dat", systDsvpvDirs[i].c_str()));
    input_DtrkPt.push_back(Form("%s/MassFit/combdata.dat", systDtrkPtDirs[i].c_str()));
    input_RapGapLoose.push_back(Form("%s/MassFit/combdata.dat", systRapGapLooseDirs[i].c_str()));
    input_RapGapTight.push_back(Form("%s/MassFit/combdata.dat", systRapGapTightDirs[i].c_str()));
  }
    
  vector<double> lambdaVal_nominal;
  vector<double> lambdaErr_nominal;
  vector<double> lambdaVal_Dalpha;
  vector<double> lambdaErr_Dalpha;
  vector<double> lambdaVal_Dchi2cl;
  vector<double> lambdaErr_Dchi2cl;
  vector<double> lambdaVal_Dsvpv;
  vector<double> lambdaErr_Dsvpv;
  vector<double> lambdaVal_DtrkPt;
  vector<double> lambdaErr_DtrkPt;
  vector<double> lambdaVal_RapGapLoose;
  vector<double> lambdaErr_RapGapLoose;
  vector<double> lambdaVal_RapGapTight;
  vector<double> lambdaErr_RapGapTight;
  
  for (int i = 0; i < input_nominal.size(); i++) {
    CombinatoricsBkgParams comb_nominal = CombinatoricsBkgParams(input_nominal[i]);
    lambdaVal_nominal.push_back(comb_nominal.lambda.getVal());
    lambdaErr_nominal.push_back(comb_nominal.lambda.getError());
    CombinatoricsBkgParams comb_Dalpha = CombinatoricsBkgParams(input_Dalpha[i]);
    lambdaVal_Dalpha.push_back(comb_Dalpha.lambda.getVal());
    lambdaErr_Dalpha.push_back(comb_Dalpha.lambda.getError());
    CombinatoricsBkgParams comb_Dchi2cl = CombinatoricsBkgParams(input_Dchi2cl[i]);
    lambdaVal_Dchi2cl.push_back(comb_Dchi2cl.lambda.getVal());
    lambdaErr_Dchi2cl.push_back(comb_Dchi2cl.lambda.getError());
    CombinatoricsBkgParams comb_Dsvpv = CombinatoricsBkgParams(input_Dsvpv[i]);
    lambdaVal_Dsvpv.push_back(comb_Dsvpv.lambda.getVal());
    lambdaErr_Dsvpv.push_back(comb_Dsvpv.lambda.getError());
    CombinatoricsBkgParams comb_DtrkPt = CombinatoricsBkgParams(input_DtrkPt[i]);
    lambdaVal_DtrkPt.push_back(comb_DtrkPt.lambda.getVal());
    lambdaErr_DtrkPt.push_back(comb_DtrkPt.lambda.getError());
    CombinatoricsBkgParams comb_RapGapLoose = CombinatoricsBkgParams(input_RapGapLoose[i]);
    lambdaVal_RapGapLoose.push_back(comb_RapGapLoose.lambda.getVal());
    lambdaErr_RapGapLoose.push_back(comb_RapGapLoose.lambda.getError());
    CombinatoricsBkgParams comb_RapGapTight = CombinatoricsBkgParams(input_RapGapTight[i]);
    lambdaVal_RapGapTight.push_back(comb_RapGapTight.lambda.getVal());
    lambdaErr_RapGapTight.push_back(comb_RapGapTight.lambda.getError());
  }

  TH1D* lambdaTemplate  = new TH1D(
    "lambdaTemplate", (plotTitle + "; y; Raw Yield").c_str(), nYBins, yBins);
  TGraphAsymmErrors* lambda_nominal  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* lambda_Dalpha  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* lambda_Dchi2cl = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* lambda_Dsvpv  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* lambda_DtrkPt  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* lambda_RapGapLoose  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* lambda_RapGapTight  = new TGraphAsymmErrors(nYBins);
  
  TH1D* ratioTemplate  = new TH1D(
    "ratioTemplate", "; y; syst/nominal", nYBins, yBins);
  TGraphAsymmErrors* ratio_Dalpha  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_Dchi2cl = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_Dsvpv  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_DtrkPt  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_RapGapLoose  = new TGraphAsymmErrors(nYBins);
  TGraphAsymmErrors* ratio_RapGapTight  = new TGraphAsymmErrors(nYBins);
    
  double offset = 0.02;
  double offset_nominal     = 6 * offset;
  double offset_Dalpha      = 4 * offset;
  double offset_Dchi2cl     = 2 * offset;
  double offset_Dsvpv       = 0 * offset;
  double offset_DtrkPt      = -2 * offset;
  double offset_RapGapLoose = -4 * offset;
  double offset_RapGapTight = -6 * offset;
  for (int yBin = 0; yBin < nYBins; yBin++) {
    double yBinCenter = 0.5 * (yBins[yBin] + yBins[yBin+1]);
    // Set yield plots
    lambda_nominal->SetPoint(yBin,
      yBinCenter - offset_nominal, lambdaVal_nominal[yBin]);
    lambda_nominal->SetPointError(yBin,
      0.5 - offset_nominal, 0.5 + offset_nominal,
      lambdaErr_nominal[yBin], lambdaErr_nominal[yBin]
    );
    lambda_Dalpha->SetPoint(yBin,
      yBinCenter - offset_Dalpha, lambdaVal_Dalpha[yBin]);
    lambda_Dalpha->SetPointError(yBin,
      0.5 - offset_Dalpha, 0.5 + offset_Dalpha,
      lambdaErr_Dalpha[yBin], lambdaErr_Dalpha[yBin]
    );
    lambda_Dchi2cl->SetPoint(yBin,
      yBinCenter - offset_Dchi2cl, lambdaVal_Dchi2cl[yBin]);
    lambda_Dchi2cl->SetPointError(yBin,
      0.5 - offset_Dchi2cl, 0.5 + offset_Dchi2cl,
      lambdaErr_Dchi2cl[yBin], lambdaErr_Dchi2cl[yBin]
    );
    lambda_Dsvpv->SetPoint(yBin,
      yBinCenter - offset_Dsvpv, lambdaVal_Dsvpv[yBin]);
    lambda_Dsvpv->SetPointError(yBin,
      0.5 - offset_Dsvpv, 0.5 + offset_Dsvpv,
      lambdaErr_Dsvpv[yBin], lambdaErr_Dsvpv[yBin]
    );
    lambda_DtrkPt->SetPoint(yBin,
      yBinCenter - offset_DtrkPt, lambdaVal_DtrkPt[yBin]);
    lambda_DtrkPt->SetPointError(yBin,
      0.5 - offset_DtrkPt, 0.5 + offset_DtrkPt,
      lambdaErr_DtrkPt[yBin], lambdaErr_DtrkPt[yBin]
    );
    lambda_RapGapLoose->SetPoint(yBin,
      yBinCenter - offset_RapGapLoose, lambdaVal_RapGapLoose[yBin]);
    lambda_RapGapLoose->SetPointError(yBin,
      0.5 - offset_RapGapLoose, 0.5 + offset_RapGapLoose,
      lambdaErr_RapGapLoose[yBin], lambdaErr_RapGapLoose[yBin]
    );
    lambda_RapGapTight->SetPoint(yBin,
      yBinCenter - offset_RapGapTight, lambdaVal_RapGapTight[yBin]);
    lambda_RapGapTight->SetPointError(yBin,
      0.5 - offset_RapGapTight, 0.5 + offset_RapGapTight,
      lambdaErr_RapGapTight[yBin], lambdaErr_RapGapTight[yBin]
    );
    // Set ratio plots
    double ratioErr_Dalpha = (
      TMath::Abs(lambdaErr_Dalpha[yBin] * lambdaVal_nominal[yBin] -
      lambdaVal_Dalpha[yBin] * lambdaErr_nominal[yBin]) /
      (lambdaVal_nominal[yBin] * lambdaVal_nominal[yBin])
    );
    ratio_Dalpha->SetPoint(yBin,
      yBinCenter - offset_Dalpha,
      lambdaVal_Dalpha[yBin]/lambdaVal_nominal[yBin]);
    ratio_Dalpha->SetPointError(yBin,
      0.5 - offset_Dalpha, 0.5 + offset_Dalpha,
      ratioErr_Dalpha, ratioErr_Dalpha
    );
    double ratioErr_Dchi2cl = (
      TMath::Abs(lambdaErr_Dchi2cl[yBin] * lambdaVal_nominal[yBin] -
      lambdaVal_Dchi2cl[yBin] * lambdaErr_nominal[yBin]) /
      (lambdaVal_nominal[yBin] * lambdaVal_nominal[yBin])
    );
    ratio_Dchi2cl->SetPoint(yBin,
      yBinCenter - offset_Dchi2cl,
      lambdaVal_Dchi2cl[yBin]/lambdaVal_nominal[yBin]);
    ratio_Dchi2cl->SetPointError(yBin,
      0.5 - offset_Dchi2cl, 0.5 + offset_Dchi2cl,
      ratioErr_Dchi2cl, ratioErr_Dchi2cl
    );
    double ratioErr_Dsvpv = (
      TMath::Abs(lambdaErr_Dsvpv[yBin] * lambdaVal_nominal[yBin] -
      lambdaVal_Dsvpv[yBin] * lambdaErr_nominal[yBin]) /
      (lambdaVal_nominal[yBin] * lambdaVal_nominal[yBin])
    );
    ratio_Dsvpv->SetPoint(yBin,
      yBinCenter - offset_Dsvpv,
      lambdaVal_Dsvpv[yBin]/lambdaVal_nominal[yBin]);
    ratio_Dsvpv->SetPointError(yBin,
      0.5 - offset_Dsvpv, 0.5 + offset_Dsvpv,
      ratioErr_Dsvpv, ratioErr_Dsvpv
    );
    double ratioErr_DtrkPt = (
      TMath::Abs(lambdaErr_DtrkPt[yBin] * lambdaVal_nominal[yBin] -
      lambdaVal_DtrkPt[yBin] * lambdaErr_nominal[yBin]) /
      (lambdaVal_nominal[yBin] * lambdaVal_nominal[yBin])
    );
    ratio_DtrkPt->SetPoint(yBin,
      yBinCenter - offset_DtrkPt,
      lambdaVal_DtrkPt[yBin]/lambdaVal_nominal[yBin]);
    ratio_DtrkPt->SetPointError(yBin,
      0.5 - offset_DtrkPt, 0.5 + offset_DtrkPt,
      ratioErr_DtrkPt, ratioErr_DtrkPt
    );
    double ratioErr_RapGapLoose = (
      TMath::Abs(lambdaErr_RapGapLoose[yBin] * lambdaVal_nominal[yBin] -
      lambdaVal_RapGapLoose[yBin] * lambdaErr_nominal[yBin]) /
      (lambdaVal_nominal[yBin] * lambdaVal_nominal[yBin])
    );
    ratio_RapGapLoose->SetPoint(yBin,
      yBinCenter - offset_RapGapLoose,
      lambdaVal_RapGapLoose[yBin]/lambdaVal_nominal[yBin]);
    ratio_RapGapLoose->SetPointError(yBin,
      0.5 - offset_RapGapLoose, 0.5 + offset_RapGapLoose,
      ratioErr_RapGapLoose, ratioErr_RapGapLoose
    );
    double ratioErr_RapGapTight = (
      TMath::Abs(lambdaErr_RapGapTight[yBin] * lambdaVal_nominal[yBin] -
      lambdaVal_RapGapTight[yBin] * lambdaErr_nominal[yBin]) /
      (lambdaVal_nominal[yBin] * lambdaVal_nominal[yBin])
    );
    ratio_RapGapTight->SetPoint(yBin,
      yBinCenter - offset_RapGapTight,
      lambdaVal_RapGapTight[yBin]/lambdaVal_nominal[yBin]);
    ratio_RapGapTight->SetPointError(yBin,
      0.5 - offset_RapGapTight, 0.5 + offset_RapGapTight,
      ratioErr_RapGapTight, ratioErr_RapGapTight
    );
  }
  lambdaTemplate->SetMinimum(-10.);
  lambdaTemplate->SetMaximum(0.0);
  ratioTemplate->SetMinimum(0.7);
  ratioTemplate->SetMaximum(1.2);
  
  // Nominal
  lambda_nominal->SetLineColor(kNominal);
  lambda_nominal->SetLineWidth(2);
  // SystFitMean
  lambda_Dalpha->SetLineColor(kDalpha);
  lambda_Dalpha->SetLineWidth(2);
  ratio_Dalpha->SetLineColor(kDalpha);
  ratio_Dalpha->SetLineWidth(2);
  // SystFitAlpha
  lambda_Dchi2cl->SetLineColor(kDchi2cl);
  lambda_Dchi2cl->SetLineWidth(2);
  ratio_Dchi2cl->SetLineColor(kDchi2cl);
  ratio_Dchi2cl->SetLineWidth(2);
  // SystFitPkbg
  lambda_Dsvpv->SetLineColor(kDsvpv);
  lambda_Dsvpv->SetLineWidth(2);
  ratio_Dsvpv->SetLineColor(kDsvpv);
  ratio_Dsvpv->SetLineWidth(2);
  // SystMassWin
  lambda_DtrkPt->SetLineColor(kDtrkPt);
  lambda_DtrkPt->SetLineWidth(2);
  ratio_DtrkPt->SetLineColor(kDtrkPt);
  ratio_DtrkPt->SetLineWidth(2);
  // SystRapGapLoose
  lambda_RapGapLoose->SetLineColor(kRapGapLoose);
  lambda_RapGapLoose->SetLineWidth(2);
  ratio_RapGapLoose->SetLineColor(kRapGapLoose);
  ratio_RapGapLoose->SetLineWidth(2);
  // SystRapGapTight
  lambda_RapGapTight->SetLineColor(kRapGapTight);
  lambda_RapGapTight->SetLineWidth(2);
  ratio_RapGapTight->SetLineColor(kRapGapTight);
  ratio_RapGapTight->SetLineWidth(2);
  
  TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
  TPad* padTop = new TPad("padTop", "", 0.0, 0.3, 1.0, 1.0);
  TPad* padBot = new TPad("padBot", "", 0.0, 0.0, 1.0, 0.3);
  padTop->SetMargin(0.18, 0.02, 0.00, 0.16);
  padBot->SetMargin(0.18, 0.02, 0.30, 0.00);
  padTop->Draw();
  padBot->Draw();
  
  padTop->cd();
  lambdaTemplate->Draw();
  lambdaTemplate->SetLabelSize(0.03/0.7, "Y");
  lambdaTemplate->SetTitleSize(0.035/0.7, "Y");
  lambda_nominal->Draw("same p");
  lambda_Dalpha->Draw("same p");
  lambda_Dchi2cl->Draw("same p");
  lambda_Dsvpv->Draw("same p");
  lambda_DtrkPt->Draw("same p");
  lambda_RapGapLoose->Draw("same p");
  lambda_RapGapTight->Draw("same p");
  gStyle->SetOptStat(0);
  
  padBot->cd();
  ratioTemplate->Draw();
  ratioTemplate->SetLabelSize(0.03/0.3, "XY");
  ratioTemplate->SetTitleSize(0.035/0.3, "XY");
  ratioTemplate->SetTitleOffset(0.4, "Y");
  TLine* unity = new TLine(yBins[0], 1.0, yBins[nYBins], 1.0);
  unity->SetLineColor(kGray);
  unity->SetLineWidth(1);
  unity->SetLineStyle(9);
  unity->Draw();
  ratio_Dalpha->Draw("same p");
  ratio_Dchi2cl->Draw("same p");
  ratio_Dsvpv->Draw("same p");
  ratio_DtrkPt->Draw("same p");
  ratio_RapGapLoose->Draw("same p");
  ratio_RapGapTight->Draw("same p");
  gStyle->SetOptStat(0);
  
  canvas->cd();
  canvas->Update();
  float legShift = 0.;
  if (!isGammaN) legShift = 0.4;
  TLegend* legend = new TLegend(0.2 + legShift, 0.32, 0.6 + legShift, 0.55);
  legend->SetTextSize(0.015/0.7);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(lambda_nominal,  "Nominal Exponential #lambda", "l");
  legend->AddEntry(lambda_Dalpha,  "Cut Syst: Dalpha", "l");
  legend->AddEntry(lambda_Dchi2cl, "Cut Syst: Dchi2cl", "l");
  legend->AddEntry(lambda_Dsvpv,  "Cut Syst: Dsvpv", "l");
  legend->AddEntry(lambda_DtrkPt,  "Cut Syst: DtrkPt", "l");
  legend->AddEntry(lambda_RapGapLoose,  "Cut Syst: RapGapLoose", "l");
  legend->AddEntry(lambda_RapGapTight,  "Cut Syst: RapGapTight", "l");
  legend->Draw();
  
  canvas->SaveAs(Form("plot/SystEval_Cuts/%s_lambda_vs_y.pdf", plotLabel.c_str()));
  
  delete canvas;
  delete lambdaTemplate;
  delete lambda_nominal;
  delete lambda_Dalpha;
  delete lambda_Dchi2cl;
  delete lambda_Dsvpv;
  delete lambda_DtrkPt;
  delete lambda_RapGapLoose;
  delete lambda_RapGapTight;
  delete ratioTemplate;
  delete ratio_Dalpha;
  delete ratio_Dchi2cl;
  delete ratio_Dsvpv;
  delete ratio_DtrkPt;
  delete ratio_RapGapLoose;
  delete ratio_RapGapTight;
}

void cuts_rawYield_vs_m(
  vector<string> nominalDirs,
  vector<string> systDalphaDirs,
  vector<string> systDchi2clDirs,
  vector<string> systDsvpvDirs,
  vector<string> systDtrkPtDirs,
  vector<string> systRapGapLooseDirs,
  vector<string> systRapGapTightDirs,
  string plotTitle,
  string plotLabel,
  int nYBins,
  double* yBins,
  double ptMin,
  double ptMax,
  int isGammaN
) {
  const int nMassBins = 6;
  double massWinUp[nMassBins+1] = {
    0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425
  };
  vector<string> massWinLabel = {
    "MassWin_166-201", "MassWin_166-206",
    "MassWin_166-211", "MassWin_166-216",
    "MassWin_166-221", "MassWin_166-226"
  };
  
  for (int i = 0; i < nYBins; i++) {
    vector<string> input_nominal;
    vector<string> input_Dalpha;
    vector<string> input_Dchi2cl;
    vector<string> input_Dsvpv;
    vector<string> input_DtrkPt;
    vector<string> input_RapGapLoose;
    vector<string> input_RapGapTight;
    for (int j = 0; j < nMassBins; j++) {
      input_nominal.push_back(Form("%s/MassFit_%s/correctedYields.md", nominalDirs[i].c_str(), massWinLabel[j].c_str()));
      input_Dalpha.push_back(Form("%s/MassFit_%s/correctedYields.md", systDalphaDirs[i].c_str(), massWinLabel[j].c_str()));
      input_Dchi2cl.push_back(Form("%s/MassFit_%s/correctedYields.md", systDchi2clDirs[i].c_str(), massWinLabel[j].c_str()));
      input_Dsvpv.push_back(Form("%s/MassFit_%s/correctedYields.md", systDsvpvDirs[i].c_str(), massWinLabel[j].c_str()));
      input_DtrkPt.push_back(Form("%s/MassFit_%s/correctedYields.md", systDtrkPtDirs[i].c_str(), massWinLabel[j].c_str()));
      input_RapGapLoose.push_back(Form("%s/MassFit_%s/correctedYields.md", systRapGapLooseDirs[i].c_str(), massWinLabel[j].c_str()));
      input_RapGapTight.push_back(Form("%s/MassFit_%s/correctedYields.md", systRapGapTightDirs[i].c_str(), massWinLabel[j].c_str()));
    }
    vector<Point> points_nominal = getPointArr(ptMin, ptMax, isGammaN, input_nominal);
    vector<Point> points_Dalpha = getPointArr(ptMin, ptMax, isGammaN, input_Dalpha);
    vector<Point> points_Dchi2cl = getPointArr(ptMin, ptMax, isGammaN, input_Dchi2cl);
    vector<Point> points_Dsvpv = getPointArr(ptMin, ptMax, isGammaN, input_Dsvpv);
    vector<Point> points_DtrkPt = getPointArr(ptMin, ptMax, isGammaN, input_DtrkPt);
    vector<Point> points_RapGapLoose = getPointArr(ptMin, ptMax, isGammaN, input_RapGapLoose);
    vector<Point> points_RapGapTight = getPointArr(ptMin, ptMax, isGammaN, input_RapGapTight);
    
    vector<double> rawYieldVal_nominal  = getDoubleArr(points_nominal, [](Point& p) -> double { return p.rawYield;});
    vector<double> rawYieldErr_nominal  = getDoubleArr(points_nominal, [](Point& p) -> double { return p.rawYieldError;});
    vector<double> rawYieldVal_Dalpha  = getDoubleArr(points_Dalpha, [](Point& p) -> double { return p.rawYield;});
    vector<double> rawYieldErr_Dalpha  = getDoubleArr(points_Dalpha, [](Point& p) -> double { return p.rawYieldError;});
    vector<double> rawYieldVal_Dchi2cl = getDoubleArr(points_Dchi2cl, [](Point& p) -> double { return p.rawYield;});
    vector<double> rawYieldErr_Dchi2cl = getDoubleArr(points_Dchi2cl, [](Point& p) -> double { return p.rawYieldError;});
    vector<double> rawYieldVal_Dsvpv  = getDoubleArr(points_Dsvpv, [](Point& p) -> double { return p.rawYield;});
    vector<double> rawYieldErr_Dsvpv  = getDoubleArr(points_Dsvpv, [](Point& p) -> double { return p.rawYieldError;});
    vector<double> rawYieldVal_DtrkPt  = getDoubleArr(points_DtrkPt, [](Point& p) -> double { return p.rawYield;});
    vector<double> rawYieldErr_DtrkPt  = getDoubleArr(points_DtrkPt, [](Point& p) -> double { return p.rawYieldError;});
    vector<double> rawYieldVal_RapGapLoose = getDoubleArr(points_RapGapLoose, [](Point& p) -> double { return p.rawYield;});
    vector<double> rawYieldErr_RapGapLoose = getDoubleArr(points_RapGapLoose, [](Point& p) -> double { return p.rawYieldError;});
    vector<double> rawYieldVal_RapGapTight = getDoubleArr(points_RapGapLoose, [](Point& p) -> double { return p.rawYield;});
    vector<double> rawYieldErr_RapGapTight = getDoubleArr(points_RapGapLoose, [](Point& p) -> double { return p.rawYieldError;});
    
    plotTitle = plotTitle + Form(", %.0f < Dy < %.0f", yBins[i], yBins[i+1]);
    TH1D* rawYieldTemplate  = new TH1D(
      "rawYieldTemplate", (plotTitle + "; Upper Mass Window (GeV); Raw Yield").c_str(), nMassBins, massWinUp);
    TGraphAsymmErrors* rawYield_nominal = new TGraphAsymmErrors(nMassBins);
    TGraphAsymmErrors* rawYield_Dalpha = new TGraphAsymmErrors(nMassBins);
    TGraphAsymmErrors* rawYield_Dchi2cl = new TGraphAsymmErrors(nMassBins);
    TGraphAsymmErrors* rawYield_Dsvpv = new TGraphAsymmErrors(nMassBins);
    TGraphAsymmErrors* rawYield_DtrkPt = new TGraphAsymmErrors(nMassBins);
    TGraphAsymmErrors* rawYield_RapGapLoose = new TGraphAsymmErrors(nMassBins);
    TGraphAsymmErrors* rawYield_RapGapTight = new TGraphAsymmErrors(nMassBins);

    TH1D* ratioTemplate  = new TH1D(
      "ratioTemplate", "; Upper Mass Window (GeV); syst/nominal", nMassBins, massWinUp);
    TGraphAsymmErrors* ratio_Dalpha = new TGraphAsymmErrors(nMassBins);
    TGraphAsymmErrors* ratio_Dchi2cl = new TGraphAsymmErrors(nMassBins);
    TGraphAsymmErrors* ratio_Dsvpv = new TGraphAsymmErrors(nMassBins);
    TGraphAsymmErrors* ratio_DtrkPt = new TGraphAsymmErrors(nMassBins);
    TGraphAsymmErrors* ratio_RapGapLoose = new TGraphAsymmErrors(nMassBins);
    TGraphAsymmErrors* ratio_RapGapTight = new TGraphAsymmErrors(nMassBins);
    
    double offset = 0.001;
    double offset_nominal  = 6 * offset;
    double offset_Dalpha  = 4 * offset;
    double offset_Dchi2cl = 2 * offset;
    double offset_Dsvpv  = 0 * offset;
    double offset_DtrkPt  = -3 * offset;
    double offset_RapGapLoose  = -4 * offset;
    double offset_RapGapTight  = -6 * offset;
    for (int mBin = 0; mBin < nMassBins; mBin++) {
      double mBinCenter = 0.5 * (massWinUp[mBin] + massWinUp[mBin+1]);
      double mBinWidth = 0.5 * (massWinUp[mBin+1] - massWinUp[mBin]);
      // Set yield plots
      rawYield_nominal->SetPoint(mBin,
        mBinCenter - offset_nominal, rawYieldVal_nominal[mBin]);
      rawYield_nominal->SetPointError(mBin,
        mBinWidth - offset_nominal, mBinWidth + offset_nominal,
        rawYieldErr_nominal[mBin], rawYieldErr_nominal[mBin]
      );
      rawYield_Dalpha->SetPoint(mBin,
        mBinCenter - offset_Dalpha, rawYieldVal_Dalpha[mBin]);
      rawYield_Dalpha->SetPointError(mBin,
        mBinWidth - offset_Dalpha, mBinWidth + offset_Dalpha,
        rawYieldErr_Dalpha[mBin], rawYieldErr_Dalpha[mBin]
      );
      rawYield_Dchi2cl->SetPoint(mBin,
        mBinCenter - offset_Dchi2cl, rawYieldVal_Dchi2cl[mBin]);
      rawYield_Dchi2cl->SetPointError(mBin,
        mBinWidth - offset_Dchi2cl, mBinWidth + offset_Dchi2cl,
        rawYieldErr_Dchi2cl[mBin], rawYieldErr_Dchi2cl[mBin]
      );
      rawYield_Dsvpv->SetPoint(mBin,
        mBinCenter - offset_Dsvpv, rawYieldVal_Dsvpv[mBin]);
      rawYield_Dsvpv->SetPointError(mBin,
        mBinWidth - offset_Dsvpv, mBinWidth + offset_Dsvpv,
        rawYieldErr_Dsvpv[mBin], rawYieldErr_Dsvpv[mBin]
      );
      rawYield_DtrkPt->SetPoint(mBin,
        mBinCenter - offset_DtrkPt, rawYieldVal_DtrkPt[mBin]);
      rawYield_DtrkPt->SetPointError(mBin,
        mBinWidth - offset_DtrkPt, mBinWidth + offset_DtrkPt,
        rawYieldErr_DtrkPt[mBin], rawYieldErr_DtrkPt[mBin]
      );
      rawYield_RapGapLoose->SetPoint(mBin,
        mBinCenter - offset_RapGapLoose, rawYieldVal_RapGapLoose[mBin]);
      rawYield_RapGapLoose->SetPointError(mBin,
        mBinWidth - offset_RapGapLoose, mBinWidth + offset_RapGapLoose,
        rawYieldErr_RapGapLoose[mBin], rawYieldErr_RapGapLoose[mBin]
      );
      rawYield_RapGapTight->SetPoint(mBin,
        mBinCenter - offset_RapGapTight, rawYieldVal_RapGapTight[mBin]);
      rawYield_RapGapTight->SetPointError(mBin,
        mBinWidth - offset_RapGapTight, mBinWidth + offset_RapGapTight,
        rawYieldErr_RapGapTight[mBin], rawYieldErr_RapGapTight[mBin]
      );
      // Set ratio plots
      double ratioErr_Dalpha = (
        TMath::Abs(rawYieldErr_Dalpha[mBin] * rawYieldVal_nominal[mBin] -
        rawYieldVal_Dalpha[mBin] * rawYieldErr_nominal[mBin]) /
        (rawYieldVal_nominal[mBin] * rawYieldVal_nominal[mBin])
      );
      ratio_Dalpha->SetPoint(mBin,
        mBinCenter - offset_Dalpha,
        rawYieldVal_Dalpha[mBin]/rawYieldVal_nominal[mBin]);
      ratio_Dalpha->SetPointError(mBin,
        mBinWidth - offset_Dalpha, mBinWidth + offset_Dalpha,
        ratioErr_Dalpha, ratioErr_Dalpha
      );
      double ratioErr_Dchi2cl = (
        TMath::Abs(rawYieldErr_Dchi2cl[mBin] * rawYieldVal_nominal[mBin] -
        rawYieldVal_Dchi2cl[mBin] * rawYieldErr_nominal[mBin]) /
        (rawYieldVal_nominal[mBin] * rawYieldVal_nominal[mBin])
      );
      ratio_Dchi2cl->SetPoint(mBin,
        mBinCenter - offset_Dchi2cl,
        rawYieldVal_Dchi2cl[mBin]/rawYieldVal_nominal[mBin]);
      ratio_Dchi2cl->SetPointError(mBin,
        mBinWidth - offset_Dchi2cl, mBinWidth + offset_Dchi2cl,
        ratioErr_Dchi2cl, ratioErr_Dchi2cl
      );
      double ratioErr_Dsvpv = (
        TMath::Abs(rawYieldErr_Dsvpv[mBin] * rawYieldVal_nominal[mBin] -
        rawYieldVal_Dsvpv[mBin] * rawYieldErr_nominal[mBin]) /
        (rawYieldVal_nominal[mBin] * rawYieldVal_nominal[mBin])
      );
      ratio_Dsvpv->SetPoint(mBin,
        mBinCenter - offset_Dsvpv,
        rawYieldVal_Dsvpv[mBin]/rawYieldVal_nominal[mBin]);
      ratio_Dsvpv->SetPointError(mBin,
        mBinWidth - offset_Dsvpv, mBinWidth + offset_Dsvpv,
        ratioErr_Dsvpv, ratioErr_Dsvpv
      );
      double ratioErr_DtrkPt = (
        TMath::Abs(rawYieldErr_DtrkPt[mBin] * rawYieldVal_nominal[mBin] -
        rawYieldVal_DtrkPt[mBin] * rawYieldErr_nominal[mBin]) /
        (rawYieldVal_nominal[mBin] * rawYieldVal_nominal[mBin])
      );
      ratio_DtrkPt->SetPoint(mBin,
        mBinCenter - offset_DtrkPt,
        rawYieldVal_DtrkPt[mBin]/rawYieldVal_nominal[mBin]);
      ratio_DtrkPt->SetPointError(mBin,
        mBinWidth - offset_DtrkPt, mBinWidth + offset_DtrkPt,
        ratioErr_DtrkPt, ratioErr_DtrkPt
      );
      double ratioErr_RapGapLoose = (
        TMath::Abs(rawYieldErr_RapGapLoose[mBin] * rawYieldVal_nominal[mBin] -
        rawYieldVal_RapGapLoose[mBin] * rawYieldErr_nominal[mBin]) /
        (rawYieldVal_nominal[mBin] * rawYieldVal_nominal[mBin])
      );
      ratio_RapGapLoose->SetPoint(mBin,
        mBinCenter - offset_RapGapLoose,
        rawYieldVal_RapGapLoose[mBin]/rawYieldVal_nominal[mBin]);
      ratio_RapGapLoose->SetPointError(mBin,
        mBinWidth - offset_RapGapLoose, mBinWidth + offset_RapGapLoose,
        ratioErr_RapGapLoose, ratioErr_RapGapLoose
      );
      double ratioErr_RapGapTight = (
        TMath::Abs(rawYieldErr_RapGapTight[mBin] * rawYieldVal_nominal[mBin] -
        rawYieldVal_RapGapTight[mBin] * rawYieldErr_nominal[mBin]) /
        (rawYieldVal_nominal[mBin] * rawYieldVal_nominal[mBin])
      );
      ratio_RapGapTight->SetPoint(mBin,
        mBinCenter - offset_RapGapTight,
        rawYieldVal_RapGapTight[mBin]/rawYieldVal_nominal[mBin]);
      ratio_RapGapTight->SetPointError(mBin,
        mBinWidth - offset_RapGapTight, mBinWidth + offset_RapGapTight,
        ratioErr_RapGapTight, ratioErr_RapGapTight
      );
    }
    rawYieldTemplate->SetMinimum(0.);
    rawYieldTemplate->SetMaximum(240.);
    ratioTemplate->SetMinimum(0.8);
    ratioTemplate->SetMaximum(1.6);
  
    // Nominal
    rawYield_nominal->SetLineColor(kNominal);
    rawYield_nominal->SetLineWidth(2);
    // SystFitDalpha
    rawYield_Dalpha->SetLineColor(kDalpha);
    rawYield_Dalpha->SetLineWidth(2);
    ratio_Dalpha->SetLineColor(kDalpha);
    ratio_Dalpha->SetLineWidth(2);
    // SystDchi2cl
    rawYield_Dchi2cl->SetLineColor(kDchi2cl);
    rawYield_Dchi2cl->SetLineWidth(2);
    ratio_Dchi2cl->SetLineColor(kDchi2cl);
    ratio_Dchi2cl->SetLineWidth(2);
    // SystDsvpv
    rawYield_Dsvpv->SetLineColor(kDsvpv);
    rawYield_Dsvpv->SetLineWidth(2);
    ratio_Dsvpv->SetLineColor(kDsvpv);
    ratio_Dsvpv->SetLineWidth(2);
    // SystDtrkPt
    rawYield_DtrkPt->SetLineColor(kDtrkPt);
    rawYield_DtrkPt->SetLineWidth(2);
    ratio_DtrkPt->SetLineColor(kDtrkPt);
    ratio_DtrkPt->SetLineWidth(2);
    // SystRapGapLoose
    rawYield_RapGapLoose->SetLineColor(kRapGapLoose);
    rawYield_RapGapLoose->SetLineWidth(2);
    ratio_RapGapLoose->SetLineColor(kRapGapLoose);
    ratio_RapGapLoose->SetLineWidth(2);
    // SystRapGapTight
    rawYield_RapGapTight->SetLineColor(kRapGapTight);
    rawYield_RapGapTight->SetLineWidth(2);
    ratio_RapGapTight->SetLineColor(kRapGapTight);
    ratio_RapGapTight->SetLineWidth(2);
    
    TCanvas* canvas = new TCanvas("canvas", "", 600, 900);
    TPad* padTop = new TPad("padTop", "", 0.0, 0.3, 1.0, 1.0);
    TPad* padBot = new TPad("padBot", "", 0.0, 0.0, 1.0, 0.3);
    padTop->SetMargin(0.15, 0.05, 0.00, 0.16);
    padBot->SetMargin(0.15, 0.05, 0.30, 0.00);
    padTop->Draw();
    padBot->Draw();
    
    padTop->cd();
    rawYieldTemplate->Draw();
    rawYieldTemplate->SetLabelSize(0.03/0.7, "Y");
    rawYieldTemplate->SetTitleSize(0.035/0.7, "Y");
    rawYieldTemplate->SetTitleOffset(1.6, "Y");
    rawYield_nominal->Draw("same p");
    rawYield_Dalpha->Draw("same p");
    rawYield_Dchi2cl->Draw("same p");
    rawYield_Dsvpv->Draw("same p");
    rawYield_DtrkPt->Draw("same p");
    rawYield_RapGapLoose->Draw("same p");
    rawYield_RapGapTight->Draw("same p");
    gStyle->SetOptStat(0);
    
    padBot->cd();
    ratioTemplate->Draw();
    ratioTemplate->SetLabelSize(0.03/0.3, "XY");
    ratioTemplate->SetTitleSize(0.035/0.3, "XY");
    ratioTemplate->SetTitleOffset(0.65, "Y");
    TLine* unity = new TLine(massWinUp[0], 1.0, massWinUp[nMassBins], 1.0);
    unity->SetLineColor(kGray);
    unity->SetLineWidth(1);
    unity->SetLineStyle(9);
    unity->Draw();
    ratio_Dalpha->Draw("same p");
    ratio_Dchi2cl->Draw("same p");
    ratio_Dsvpv->Draw("same p");
    ratio_DtrkPt->Draw("same p");
    ratio_RapGapLoose->Draw("same p");
    ratio_RapGapTight->Draw("same p");
    gStyle->SetOptStat(0);
    
    canvas->cd();
    canvas->Update();
    float legShift = 0.;
    TLegend* legend = new TLegend(0.2, 0.73, 0.6, 0.88);
    legend->SetTextSize(0.02/0.7);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->AddEntry(rawYield_nominal,  "Nominal Raw Yield", "l");
    legend->AddEntry(rawYield_Dalpha,  "Cut Syst: Dalpha", "l");
    legend->AddEntry(rawYield_Dchi2cl, "Cut Syst: Dchi2cl", "l");
    legend->AddEntry(rawYield_Dsvpv,  "Cut Syst: Dsvpv", "l");
    legend->AddEntry(rawYield_DtrkPt,  "Cut Syst: DtrkPt", "l");
    legend->AddEntry(rawYield_RapGapLoose,  "Cut Syst: RapGapLoose", "l");
    legend->AddEntry(rawYield_RapGapTight,  "Cut Syst: RapGapTight", "l");
    legend->Draw();
    
    canvas->SaveAs(Form("plot/SystEval_Cuts/%s_rawYield_vs_m_y%.0f-%.0f.pdf",
      plotLabel.c_str(), yBins[i], yBins[i+1]));
    
    delete canvas;
    delete rawYieldTemplate;
    delete rawYield_nominal;
    delete rawYield_Dalpha;
    delete rawYield_Dchi2cl;
    delete rawYield_Dsvpv;
    delete rawYield_DtrkPt;
    delete rawYield_RapGapLoose;
    delete rawYield_RapGapTight;
    delete ratioTemplate;
    delete ratio_Dalpha;
    delete ratio_Dchi2cl;
    delete ratio_Dsvpv;
    delete ratio_DtrkPt;
    delete ratio_RapGapLoose;
    delete ratio_RapGapTight;
  }
}

void massfit_rawYield_vs_m(
  vector<string> nominalDirs,
  string plotTitle,
  string plotLabel,
  int nYBins,
  double* yBins,
  double ptMin,
  double ptMax,
  int isGammaN
) {
  const int nMassBins = 6;
  double massWinUp[nMassBins+1] = {
    0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425
  };
  vector<string> massWinLabel = {
    "MassWin_166-201", "MassWin_166-206",
    "MassWin_166-211", "MassWin_166-216",
    "MassWin_166-221", "MassWin_166-226"
  };
  
  for (int i = 0; i < nYBins; i++) {
    vector<string> input_nominal;
    vector<string> input_fitMean;
    vector<string> input_fitAlpha;
    vector<string> input_fitPkbg;
    
    for (int j = 0; j < nMassBins; j++) {
      input_nominal.push_back(Form("%s/MassFit_%s/correctedYields.md", nominalDirs[i].c_str(), massWinLabel[j].c_str()));
      input_fitMean.push_back(Form("%s/MassFit_systFitSigMean_%s/correctedYields.md", nominalDirs[i].c_str(), massWinLabel[j].c_str()));
      input_fitAlpha.push_back(Form("%s/MassFit_systFitSigAlpha_%s/correctedYields.md", nominalDirs[i].c_str(), massWinLabel[j].c_str()));
      input_fitPkbg.push_back(Form("%s/MassFit_systPkBg_%s/correctedYields.md", nominalDirs[i].c_str(), massWinLabel[j].c_str()));
    }
    vector<Point> points_nominal = getPointArr(ptMin, ptMax, isGammaN, input_nominal);
    vector<Point> points_fitMean = getPointArr(ptMin, ptMax, isGammaN, input_fitMean);
    vector<Point> points_fitAlpha = getPointArr(ptMin, ptMax, isGammaN, input_fitAlpha);
    vector<Point> points_fitPkbg = getPointArr(ptMin, ptMax, isGammaN, input_fitPkbg);
    
    vector<double> rawYieldVal_nominal  = getDoubleArr(points_nominal, [](Point& p) -> double { return p.rawYield;});
    vector<double> rawYieldErr_nominal  = getDoubleArr(points_nominal, [](Point& p) -> double { return p.rawYieldError;});
    vector<double> rawYieldVal_fitMean  = getDoubleArr(points_fitMean, [](Point& p) -> double { return p.rawYield;});
    vector<double> rawYieldErr_fitMean  = getDoubleArr(points_fitMean, [](Point& p) -> double { return p.rawYieldError;});
    vector<double> rawYieldVal_fitAlpha = getDoubleArr(points_fitAlpha, [](Point& p) -> double { return p.rawYield;});
    vector<double> rawYieldErr_fitAlpha = getDoubleArr(points_fitAlpha, [](Point& p) -> double { return p.rawYieldError;});
    vector<double> rawYieldVal_fitPkbg  = getDoubleArr(points_fitPkbg, [](Point& p) -> double { return p.rawYield;});
    vector<double> rawYieldErr_fitPkbg  = getDoubleArr(points_fitPkbg, [](Point& p) -> double { return p.rawYieldError;});
    
    plotTitle = plotTitle + Form(", %.0f < Dy < %.0f", yBins[i], yBins[i+1]);
    TH1D* rawYieldTemplate  = new TH1D(
      "rawYieldTemplate", (plotTitle + "; Upper Mass Window (GeV); Raw Yield").c_str(), nMassBins, massWinUp);
    TGraphAsymmErrors* rawYield_nominal = new TGraphAsymmErrors(nMassBins);
    TGraphAsymmErrors* rawYield_fitMean = new TGraphAsymmErrors(nMassBins);
    TGraphAsymmErrors* rawYield_fitAlpha = new TGraphAsymmErrors(nMassBins);
    TGraphAsymmErrors* rawYield_fitPkbg = new TGraphAsymmErrors(nMassBins);

    TH1D* ratioTemplate  = new TH1D(
      "ratioTemplate", "; Upper Mass Window (GeV); syst/nominal", nMassBins, massWinUp);
    TGraphAsymmErrors* ratio_fitMean = new TGraphAsymmErrors(nMassBins);
    TGraphAsymmErrors* ratio_fitAlpha = new TGraphAsymmErrors(nMassBins);
    TGraphAsymmErrors* ratio_fitPkbg = new TGraphAsymmErrors(nMassBins);
    
    double offset = 0.001;
    double offset_nominal  = 6 * offset;
    double offset_fitMean  = 4 * offset;
    double offset_fitAlpha = 2 * offset;
    double offset_fitPkbg  = 0 * offset;
    for (int mBin = 0; mBin < nMassBins; mBin++) {
      double mBinCenter = 0.5 * (massWinUp[mBin] + massWinUp[mBin+1]);
      double mBinWidth = 0.5 * (massWinUp[mBin+1] - massWinUp[mBin]);
      // Set yield plots
      rawYield_nominal->SetPoint(mBin,
        mBinCenter - offset_nominal, rawYieldVal_nominal[mBin]);
      rawYield_nominal->SetPointError(mBin,
        mBinWidth - offset_nominal, mBinWidth + offset_nominal,
        rawYieldErr_nominal[mBin], rawYieldErr_nominal[mBin]
      );
      rawYield_fitMean->SetPoint(mBin,
        mBinCenter - offset_fitMean, rawYieldVal_fitMean[mBin]);
      rawYield_fitMean->SetPointError(mBin,
        mBinWidth - offset_fitMean, mBinWidth + offset_fitMean,
        rawYieldErr_fitMean[mBin], rawYieldErr_fitMean[mBin]
      );
      rawYield_fitAlpha->SetPoint(mBin,
        mBinCenter - offset_fitAlpha, rawYieldVal_fitAlpha[mBin]);
      rawYield_fitAlpha->SetPointError(mBin,
        mBinWidth - offset_fitAlpha, mBinWidth + offset_fitAlpha,
        rawYieldErr_fitAlpha[mBin], rawYieldErr_fitAlpha[mBin]
      );
      rawYield_fitPkbg->SetPoint(mBin,
        mBinCenter - offset_fitPkbg, rawYieldVal_fitPkbg[mBin]);
      rawYield_fitPkbg->SetPointError(mBin,
        mBinWidth - offset_fitPkbg, mBinWidth + offset_fitPkbg,
        rawYieldErr_fitPkbg[mBin], rawYieldErr_fitPkbg[mBin]
      );
      // Set ratio plots
      double ratioErr_fitMean = (
        TMath::Abs(rawYieldErr_fitMean[mBin] * rawYieldVal_nominal[mBin] -
        rawYieldVal_fitMean[mBin] * rawYieldErr_nominal[mBin]) /
        (rawYieldVal_nominal[mBin] * rawYieldVal_nominal[mBin])
      );
      ratio_fitMean->SetPoint(mBin,
        mBinCenter - offset_fitMean,
        rawYieldVal_fitMean[mBin]/rawYieldVal_nominal[mBin]);
      ratio_fitMean->SetPointError(mBin,
        mBinWidth - offset_fitMean, mBinWidth + offset_fitMean,
        ratioErr_fitMean, ratioErr_fitMean
      );
      double ratioErr_fitAlpha = (
        TMath::Abs(rawYieldErr_fitAlpha[mBin] * rawYieldVal_nominal[mBin] -
        rawYieldVal_fitAlpha[mBin] * rawYieldErr_nominal[mBin]) /
        (rawYieldVal_nominal[mBin] * rawYieldVal_nominal[mBin])
      );
      ratio_fitAlpha->SetPoint(mBin,
        mBinCenter - offset_fitAlpha,
        rawYieldVal_fitAlpha[mBin]/rawYieldVal_nominal[mBin]);
      ratio_fitAlpha->SetPointError(mBin,
        mBinWidth - offset_fitAlpha, mBinWidth + offset_fitAlpha,
        ratioErr_fitAlpha, ratioErr_fitAlpha
      );
      double ratioErr_fitPkbg = (
        TMath::Abs(rawYieldErr_fitPkbg[mBin] * rawYieldVal_nominal[mBin] -
        rawYieldVal_fitPkbg[mBin] * rawYieldErr_nominal[mBin]) /
        (rawYieldVal_nominal[mBin] * rawYieldVal_nominal[mBin])
      );
      ratio_fitPkbg->SetPoint(mBin,
        mBinCenter - offset_fitPkbg,
        rawYieldVal_fitPkbg[mBin]/rawYieldVal_nominal[mBin]);
      ratio_fitPkbg->SetPointError(mBin,
        mBinWidth - offset_fitPkbg, mBinWidth + offset_fitPkbg,
        ratioErr_fitPkbg, ratioErr_fitPkbg
      );
    }
    rawYieldTemplate->SetMinimum(0.);
    rawYieldTemplate->SetMaximum(200.);
    ratioTemplate->SetMinimum(0.7);
    ratioTemplate->SetMaximum(1.1);
  
    // Nominal
    rawYield_nominal->SetLineColor(kNominal);
    rawYield_nominal->SetLineWidth(2);
    // SystFitMean
    rawYield_fitMean->SetLineColor(kFitMean);
    rawYield_fitMean->SetLineWidth(2);
    ratio_fitMean->SetLineColor(kFitMean);
    ratio_fitMean->SetLineWidth(2);
    // SystFitAlpha
    rawYield_fitAlpha->SetLineColor(kFitAlpha);
    rawYield_fitAlpha->SetLineWidth(2);
    ratio_fitAlpha->SetLineColor(kFitAlpha);
    ratio_fitAlpha->SetLineWidth(2);
    // SystFitPkbg
    rawYield_fitPkbg->SetLineColor(kFitPkbg);
    rawYield_fitPkbg->SetLineWidth(2);
    ratio_fitPkbg->SetLineColor(kFitPkbg);
    ratio_fitPkbg->SetLineWidth(2);
    
    TCanvas* canvas = new TCanvas("canvas", "", 600, 900);
    TPad* padTop = new TPad("padTop", "", 0.0, 0.3, 1.0, 1.0);
    TPad* padBot = new TPad("padBot", "", 0.0, 0.0, 1.0, 0.3);
    padTop->SetMargin(0.15, 0.05, 0.00, 0.16);
    padBot->SetMargin(0.15, 0.05, 0.30, 0.00);
    padTop->Draw();
    padBot->Draw();
    
    padTop->cd();
    rawYieldTemplate->Draw();
    rawYieldTemplate->SetLabelSize(0.03/0.7, "Y");
    rawYieldTemplate->SetTitleSize(0.035/0.7, "Y");
    rawYieldTemplate->SetTitleOffset(1.6, "Y");
    rawYield_nominal->Draw("same p");
    rawYield_fitMean->Draw("same p");
    rawYield_fitAlpha->Draw("same p");
    rawYield_fitPkbg->Draw("same p");
    gStyle->SetOptStat(0);
    
    padBot->cd();
    ratioTemplate->Draw();
    ratioTemplate->SetLabelSize(0.03/0.3, "XY");
    ratioTemplate->SetTitleSize(0.035/0.3, "XY");
    ratioTemplate->SetTitleOffset(0.65, "Y");
    TLine* unity = new TLine(massWinUp[0], 1.0, massWinUp[nMassBins], 1.0);
    unity->SetLineColor(kGray);
    unity->SetLineWidth(1);
    unity->SetLineStyle(9);
    unity->Draw();
    ratio_fitMean->Draw("same p");
    ratio_fitAlpha->Draw("same p");
    ratio_fitPkbg->Draw("same p");
    gStyle->SetOptStat(0);
    
    canvas->cd();
    canvas->Update();
    float legShift = 0.;
    TLegend* legend = new TLegend(0.2, 0.73, 0.6, 0.88);
    legend->SetTextSize(0.02/0.7);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->AddEntry(rawYield_nominal,  "Nominal Raw Yield", "l");
    legend->AddEntry(rawYield_fitMean,  "Fit Syst: Signal mean", "l");
    legend->AddEntry(rawYield_fitAlpha, "Fit Syst: Signal alpha", "l");
    legend->AddEntry(rawYield_fitPkbg,  "Fit Syst: KK + #pi#pi Peaks", "l");
    legend->Draw();
    
    canvas->SaveAs(Form("plot/SystEval_MassFit/%s_rawYield_vs_m_y%.0f-%.0f.pdf",
      plotLabel.c_str(), yBins[i], yBins[i+1]));
    
    delete canvas;
    delete rawYieldTemplate;
    delete rawYield_nominal;
    delete rawYield_fitMean;
    delete rawYield_fitAlpha;
    delete rawYield_fitPkbg;
    delete ratioTemplate;
    delete ratio_fitMean;
    delete ratio_fitAlpha;
    delete ratio_fitPkbg;
  }
}

void save_corrYieldPlot(
  TH1D* corrYieldTemplate,
  TH1D* corrYield_nominal,
  TH1D* ratioTemplate,
  TH1D* ratio_nominal,
  string filename,
  float ratioMin = 0.8, float ratioMax = 1.2,
  float yieldMin = 0.0, float yieldMax = 3.0
) {
  ratioTemplate->SetMinimum(ratioMin);
  ratioTemplate->SetMaximum(ratioMax);
  corrYieldTemplate->SetMinimum(yieldMin);
  corrYieldTemplate->SetMaximum(yieldMax);

  // Styling
  corrYield_nominal->SetMarkerColor(kNominal);
  corrYield_nominal->SetLineColor(kNominal);
  corrYield_nominal->SetLineWidth(2);
  ratio_nominal->SetMarkerColor(kNominal);
  ratio_nominal->SetLineColor(kNominal);
  ratio_nominal->SetLineWidth(2);
  
  TCanvas* canvas = new TCanvas("canvas", "", 600, 900);
  TPad* padTop = new TPad("padTop", "", 0.0, 0.3, 1.0, 1.0);
  TPad* padBot = new TPad("padBot", "", 0.0, 0.0, 1.0, 0.3);
  padTop->SetMargin(0.10, 0.10, 0.00, 0.16);
  padBot->SetMargin(0.10, 0.10, 0.30, 0.00);
  padTop->Draw();
  padBot->Draw();
  
  padTop->cd();
  corrYieldTemplate->Draw();
  corrYieldTemplate->SetLabelSize(0.025/0.7, "Y");
  corrYieldTemplate->SetTitleSize(0.03/0.7, "Y");
  corrYieldTemplate->SetTitleOffset(1.6, "Y");
  corrYield_nominal->Draw("same");
  gStyle->SetOptStat(0);
  
  padBot->cd();
  ratioTemplate->Draw();
  ratioTemplate->SetLabelSize(0.025/0.3, "XY");
  ratioTemplate->SetTitleSize(0.03/0.3, "XY");
  ratioTemplate->SetTitleOffset(0.65, "Y");
  TLine* unity = new TLine(
    ratioTemplate->GetBinLowEdge(1), 1.0,
    ratioTemplate->GetBinLowEdge(ratioTemplate->GetNbinsX()), 1.0);
  unity->SetLineColor(kGray);
  unity->SetLineWidth(1);
  unity->SetLineStyle(9);
  unity->Draw();
  ratio_nominal->Draw("same");
  gStyle->SetOptStat(0);
  
  canvas->cd();
  canvas->Update();
  TLegend* legend = new TLegend(0.4, 0.84, 0.88, 0.88);
  legend->SetTextSize(0.025/0.7);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(corrYield_nominal,  "Nominal Corr. Yield", "l");
  legend->Draw();
  
  canvas->SaveAs(filename.c_str());
  
  delete canvas;
}

void cuts_corrYield_vs_DtrkPtCut(
  vector<string> nominalDirs,
  vector<string> systSubdir,
  string plotTitle,
  string plotLabel,
  int nYBins,
  double* yBins,
  double ptMin,
  double ptMax,
  int isGammaN
) {
  const int nDtrkPtBins = 6;
  double DtrkPtBins[nDtrkPtBins+1] = {
    0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25
  };
  vector<string> DtrkPtLabel = {
    "_070", "_080", "_090", "_100", "_110", "_120"
  };
  
  for (int i = 0; i < nYBins; i++) {
    vector<string> input_nominal;
    for (int j = 0; j < nDtrkPtBins; j++) {
      input_nominal.push_back(Form("systDtrkPt%s/%s/MassFit/correctedYields.md", DtrkPtLabel[j].c_str(), systSubdir[i].c_str()));
    }
    vector<Point> points_nominal = getPointArr(ptMin, ptMax, isGammaN, input_nominal);
    
    vector<double> corrYieldVal_nominal  = getDoubleArr(points_nominal, [](Point& p) -> double { return p.correctedYield;});
    vector<double> corrYieldErr_nominal  = getDoubleArr(points_nominal, [](Point& p) -> double { return p.correctedYieldError;});
          
    TH1D* corrYieldTemplate  = new TH1D(
      "corrYieldTemplate",
      Form("%s, %.0f < Dy < %.0f; DtrkPt Cut (GeV); Corrected Yield",
        plotTitle.c_str(), yBins[i], yBins[i+1]),
      nDtrkPtBins, DtrkPtBins
    );
    TH1D* corrYield_nominal = new TH1D(
      "corrYield_nominal",
      corrYieldTemplate->GetTitle(),
      nDtrkPtBins, DtrkPtBins
    );
    TH1D* corrYield_ref = new TH1D(
      "corrYield_ref", plotTitle.c_str(), nDtrkPtBins, DtrkPtBins);
    
    for (int bin = 1; bin <= nDtrkPtBins; bin++) {
      corrYield_nominal->SetBinContent(bin, corrYieldVal_nominal[bin-1]);
      corrYield_nominal->SetBinError(bin, corrYieldErr_nominal[bin-1]);
      
      corrYield_ref->SetBinContent(bin, corrYieldVal_nominal[2]);
      corrYield_ref->SetBinError(bin, corrYieldErr_nominal[2]);
    }
    corrYield_nominal->Sumw2();
    corrYield_ref->Sumw2();
    TH1D* ratioTemplate  = new TH1D(
      "ratioTemplate", "; DtrkPt Cut (GeV); Scan / [DtrkPt > 1.0]",
      nDtrkPtBins, DtrkPtBins
    );
    TH1D* ratio_nominal = (TH1D*) corrYield_nominal->Clone("ratio_nominal");
    ratio_nominal->SetTitle(ratioTemplate->GetTitle());
    ratio_nominal->Sumw2();
    ratio_nominal->Divide(corrYield_ref);
    
    save_corrYieldPlot(
      corrYieldTemplate,
      corrYield_nominal,
      ratioTemplate,
      ratio_nominal,
      Form("plot/SystEval_CutScans/%s_corrYield_vs_DtrkPt_y%.0f-%.0f.pdf",
        plotLabel.c_str(), yBins[i], yBins[i+1]),
      0.4, // ratioMin
      1.6  // ratioMax
    );
    
    delete corrYieldTemplate;
    delete corrYield_nominal;
    delete corrYield_ref;
    delete ratioTemplate;
    delete ratio_nominal;
  }
}

void cuts_corrYield_vs_DsvpvCut(
  vector<string> nominalDirs,
  vector<string> systSubdir,
  string plotTitle,
  string plotLabel,
  int nYBins,
  double* yBins,
  double ptMin,
  double ptMax,
  int isGammaN
) {
  const int nDsvpvBins = 5;
  double DsvpvBins[nDsvpvBins+1] = {
    1.875, 2.125, 2.375, 2.625, 2.875, 3.125
  };
  vector<string> DsvpvLabel = {
    "_200", "_225", "_250", "_275", "_300"
  };
  
  for (int i = 0; i < nYBins; i++) {
    vector<string> input_nominal;
    for (int j = 0; j < nDsvpvBins; j++) {
      input_nominal.push_back(Form("systDsvpv%s/%s/MassFit/correctedYields.md", DsvpvLabel[j].c_str(), systSubdir[i].c_str()));
    }
    vector<Point> points_nominal = getPointArr(ptMin, ptMax, isGammaN, input_nominal);
    
    vector<double> corrYieldVal_nominal  = getDoubleArr(points_nominal, [](Point& p) -> double { return p.correctedYield;});
    vector<double> corrYieldErr_nominal  = getDoubleArr(points_nominal, [](Point& p) -> double { return p.correctedYieldError;});
    
    TH1D* corrYieldTemplate  = new TH1D(
      "corrYieldTemplate",
      Form(
        "%s, %.0f < Dy < %.0f; Dsvpv Significance Cut; Corrected Yield",
        plotTitle.c_str(), yBins[i], yBins[i+1]),
      nDsvpvBins, DsvpvBins
    );
    TH1D* corrYield_nominal = new TH1D(
      "corrYield_nominal",
      corrYieldTemplate->GetTitle(),
      nDsvpvBins, DsvpvBins
    );
    TH1D* corrYield_ref = new TH1D(
      "corrYield_ref", plotTitle.c_str(), nDsvpvBins, DsvpvBins);
    
    for (int bin = 1; bin <= nDsvpvBins; bin++) {
      corrYield_nominal->SetBinContent(bin, corrYieldVal_nominal[bin-1]);
      corrYield_nominal->SetBinError(bin, corrYieldErr_nominal[bin-1]);
      
      corrYield_ref->SetBinContent(bin, corrYieldVal_nominal[2]);
      corrYield_ref->SetBinError(bin, corrYieldErr_nominal[2]);
    }
    corrYield_nominal->Sumw2();
    corrYield_ref->Sumw2();
    TH1D* ratioTemplate  = new TH1D(
      "ratioTemplate", "; Dsvpv Significance Cut; Scan / [Dsvpv > 2.5]",
      nDsvpvBins, DsvpvBins
    );
    TH1D* ratio_nominal = (TH1D*) corrYield_nominal->Clone("ratio_nominal");
    ratio_nominal->SetTitle(ratioTemplate->GetTitle());
    ratio_nominal->Sumw2();
    ratio_nominal->Divide(corrYield_ref);
    
    save_corrYieldPlot(
      corrYieldTemplate,
      corrYield_nominal,
      ratioTemplate,
      ratio_nominal,
      Form("plot/SystEval_CutScans/%s_corrYield_vs_Dsvpv_y%.0f-%.0f.pdf",
        plotLabel.c_str(), yBins[i], yBins[i+1]),
      0.4, // ratioMin
      1.6  // ratioMax
    );
    
    delete corrYieldTemplate;
    delete corrYield_nominal;
    delete corrYield_ref;
    delete ratioTemplate;
    delete ratio_nominal;
  }
}

void cuts_corrYield_vs_DalphaCut(
  vector<string> nominalDirs,
  vector<string> systSubdir,
  string plotTitle,
  string plotLabel,
  int nYBins,
  double* yBins,
  double ptMin,
  double ptMax,
  int isGammaN
) {
  const int nDalphaBins = 7;
  double DalphaBins[nDalphaBins+1] = {
    0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85
  };
  vector<string> DalphaLabel = {
    "_020", "_030", "_040", "_050", "_060", "_070", "_080"
  };
  
  for (int i = 0; i < nYBins; i++) {
    vector<string> input_nominal;
    for (int j = 0; j < nDalphaBins; j++) {
      input_nominal.push_back(Form("systDalpha%s/%s/MassFit/correctedYields.md", DalphaLabel[j].c_str(), systSubdir[i].c_str()));
    }
    vector<Point> points_nominal = getPointArr(ptMin, ptMax, isGammaN, input_nominal);
    
    vector<double> corrYieldVal_nominal  = getDoubleArr(points_nominal, [](Point& p) -> double { return p.correctedYield;});
    vector<double> corrYieldErr_nominal  = getDoubleArr(points_nominal, [](Point& p) -> double { return p.correctedYieldError;});
    
    TH1D* corrYieldTemplate  = new TH1D(
      "corrYieldTemplate",
      Form(
        "%s, %.0f < Dy < %.0f; Dalpha Significance Cut; Corrected Yield",
        plotTitle.c_str(), yBins[i], yBins[i+1]),
      nDalphaBins, DalphaBins
    );
    TH1D* corrYield_nominal = new TH1D(
      "corrYield_nominal",
      corrYieldTemplate->GetTitle(),
      nDalphaBins, DalphaBins
    );
    TH1D* corrYield_ref = new TH1D(
      "corrYield_ref", plotTitle.c_str(), nDalphaBins, DalphaBins);
    
    for (int bin = 1; bin <= nDalphaBins; bin++) {
      corrYield_nominal->SetBinContent(bin, corrYieldVal_nominal[bin-1]);
      corrYield_nominal->SetBinError(bin, corrYieldErr_nominal[bin-1]);
      
      if (i == 0 || i == 3) { // ends use < 0.2
        corrYield_ref->SetBinContent(bin, corrYieldVal_nominal[0]);
        corrYield_ref->SetBinError(bin, corrYieldErr_nominal[0]);
      }
      else { // barrel uses < 0.4
        corrYield_ref->SetBinContent(bin, corrYieldVal_nominal[2]);
        corrYield_ref->SetBinError(bin, corrYieldErr_nominal[2]);
      }
    }
    corrYield_nominal->Sumw2();
    corrYield_ref->Sumw2();
    TH1D* ratioTemplate  = new TH1D(
      "ratioTemplate", "",
      nDalphaBins, DalphaBins
    );
    if (i == 0 || i == 3) ratioTemplate->SetTitle(
      "; Dalpha Significance Cut; Scan / [Dalpha < 0.2]");
    else ratioTemplate->SetTitle(
      "; Dalpha Significance Cut; Scan / [Dalpha < 0.4]");
    
    TH1D* ratio_nominal = (TH1D*) corrYield_nominal->Clone("ratio_nominal");
    ratio_nominal->SetTitle(ratioTemplate->GetTitle());
    ratio_nominal->Sumw2();
    ratio_nominal->Divide(corrYield_ref);
    
    save_corrYieldPlot(
      corrYieldTemplate,
      corrYield_nominal,
      ratioTemplate,
      ratio_nominal,
      Form("plot/SystEval_CutScans/%s_corrYield_vs_Dalpha_y%.0f-%.0f.pdf",
        plotLabel.c_str(), yBins[i], yBins[i+1]),
      0.4, // ratioMin
      1.6  // ratioMax
    );
    
    delete corrYieldTemplate;
    delete corrYield_nominal;
    delete corrYield_ref;
    delete ratioTemplate;
    delete ratio_nominal;
  }
}

  const int nDBins = 6;
  double DBins[nDBins+1] = {
    0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325
  };
  vector<string> Dchi2Label = {
    "_005", "_010", "_015", "_020", "_025", "_030"
  };

void cuts_corrYield_rawYield_vs_DCut(
  vector<string> systSubdir_gN,
  vector<string> systSubdir_Ng,
  string plotTitle,
  string plotTag,
  string DCutDir,
  string DCutXaxis,
  int nDBins,
  double* DBins,
  string* DLabels,
  int nYBins,
  double* yBins,
  int* refBins,
  double ratioMin = 0.2,
  double ratioMax = 1.8,
  double corrYieldMin = 0.0,
  double corrYieldMax = 3.0,
  double rawYieldMin = 0.0,
  double rawYieldMax = 200.0,
  double ptMin = 2.0,
  double ptMax = 5.0
) {
  
  for (int i = 0; i < nYBins; i++) {
    vector<string> input_gN;
    vector<string> input_Ng;
    for (int j = 0; j < nDBins; j++) {
      input_gN.push_back(Form("%s%s/%s/MassFit/correctedYields.md",
        DCutDir.c_str(), DLabels[j].c_str(), systSubdir_gN[i].c_str()));
      input_Ng.push_back(Form("%s%s/%s/MassFit/correctedYields.md",
        DCutDir.c_str(), DLabels[j].c_str(), systSubdir_Ng[nYBins-i-1].c_str()));
    }
    vector<Point> points_gN = getPointArr(ptMin, ptMax, 1, input_gN);
    vector<Point> points_Ng = getPointArr(ptMin, ptMax, 0, input_Ng);
    
    vector<double> corrYieldVal_gN = getDoubleArr(
      points_gN, [](Point& p) -> double { return p.correctedYield;});
    vector<double> corrYieldErr_gN = getDoubleArr(
      points_gN, [](Point& p) -> double { return p.correctedYieldError;});
    vector<double> corrYieldVal_Ng = getDoubleArr(
      points_Ng, [](Point& p) -> double { return p.correctedYield;});
    vector<double> corrYieldErr_Ng = getDoubleArr(
      points_Ng, [](Point& p) -> double { return p.correctedYieldError;});
    
    vector<double> rawYieldVal_gN  = getDoubleArr(
      points_gN, [](Point& p) -> double { return p.rawYield;});
    vector<double> rawYieldErr_gN  = getDoubleArr(
      points_gN, [](Point& p) -> double { return p.rawYieldError;});
    vector<double> rawYieldVal_Ng  = getDoubleArr(
      points_Ng, [](Point& p) -> double { return p.rawYield;});
    vector<double> rawYieldErr_Ng  = getDoubleArr(
      points_Ng, [](Point& p) -> double { return p.rawYieldError;});
    
    vector<double> effEvtVal_gN  = getDoubleArr(
      points_gN, [](Point& p) -> double { return p.effEvt;});
    vector<double> effEvtErr_gN  = getDoubleArr(
      points_gN, [](Point& p) -> double { return p.effEvtError;});
    vector<double> effEvtVal_Ng  = getDoubleArr(
      points_Ng, [](Point& p) -> double { return p.effEvt;});
    vector<double> effEvtErr_Ng  = getDoubleArr(
      points_Ng, [](Point& p) -> double { return p.effEvtError;});
      
    vector<double> effDVal_gN  = getDoubleArr(
      points_gN, [](Point& p) -> double { return p.effD;});
    vector<double> effDErr_gN  = getDoubleArr(
      points_gN, [](Point& p) -> double { return p.effDError;});
    vector<double> effDVal_Ng  = getDoubleArr(
      points_Ng, [](Point& p) -> double { return p.effD;});
    vector<double> effDErr_Ng  = getDoubleArr(
      points_Ng, [](Point& p) -> double { return p.effDError;});
    
    TH1D* corrYieldTemplate  = new TH1D(
      "corrYieldTemplate",
      Form(
        "%s, %.0f < Dy < %.0f; %s; Corrected Yield",
        plotTitle.c_str(), yBins[i], yBins[i+1], DCutXaxis.c_str()),
      nDBins, DBins
    );
    TH1D* corrYield_gN = new TH1D(
      "corrYield_gN", corrYieldTemplate->GetTitle(), nDBins, DBins);
    TH1D* corrYield_Ng = new TH1D(
      "corrYield_Ng", corrYieldTemplate->GetTitle(), nDBins, DBins);
    TH1D* corrYield_ref_gN = new TH1D(
      "corrYield_ref_gN", corrYieldTemplate->GetTitle(), nDBins, DBins);
    TH1D* corrYield_ref_Ng = new TH1D(
      "corrYield_ref_Ng", corrYieldTemplate->GetTitle(), nDBins, DBins);
    
    TH1D* rawYieldTemplate  = new TH1D(
      "rawYieldTemplate",
      Form(
        "%s, %.0f < Dy < %.0f; %s; Raw Yield",
        plotTitle.c_str(), yBins[i], yBins[i+1], DCutXaxis.c_str()),
      nDBins, DBins
    );
    TH1D* rawYield_gN = new TH1D(
      "rawYield_gN", rawYieldTemplate->GetTitle(), nDBins, DBins);
    TH1D* rawYield_Ng = new TH1D(
      "rawYield_Ng", rawYieldTemplate->GetTitle(), nDBins, DBins);
    TH1D* rawYield_ref_gN = new TH1D(
      "rawYield_ref_gN", rawYieldTemplate->GetTitle(), nDBins, DBins);
    TH1D* rawYield_ref_Ng = new TH1D(
      "rawYield_ref_Ng", rawYieldTemplate->GetTitle(), nDBins, DBins);
    
    TH1D* effEvtRatio = new TH1D(
      "effEvtRatio",
      Form(
        "%s, %.0f < Dy < %.0f; %s; #varepsilon^{#gammaN}_{Event} / #varepsilon^{N#gamma}_{Event}",
        plotTitle.c_str(), yBins[i], yBins[i+1], DCutXaxis.c_str()),
      nDBins, DBins
    );
    TH1D* effEvt_gN = new TH1D(
      "effEvt_gN", effEvtRatio->GetTitle(), nDBins, DBins);
    TH1D* effEvt_Ng = new TH1D(
      "effEvt_Ng", effEvtRatio->GetTitle(), nDBins, DBins);
    
    TH1D* effDRatio = new TH1D(
      "effDRatio",
      Form(
        "%s, %.0f < Dy < %.0f; %s; #varepsilon^{#gammaN}_{D^{0}} / #varepsilon^{N#gamma}_{D^{0}}",
        plotTitle.c_str(), yBins[i], yBins[i+1], DCutXaxis.c_str()),
      nDBins, DBins
    );
    TH1D* effD_gN = new TH1D(
      "effD_gN", effDRatio->GetTitle(), nDBins, DBins);
    TH1D* effD_Ng = new TH1D(
      "effD_Ng", effDRatio->GetTitle(), nDBins, DBins);
    
    for (int bin = 1; bin <= nDBins; bin++) {
      corrYield_gN->SetBinContent(bin, corrYieldVal_gN[bin-1]);
      corrYield_gN->SetBinError(bin, corrYieldErr_gN[bin-1]);
      corrYield_Ng->SetBinContent(bin, corrYieldVal_Ng[bin-1]);
      corrYield_Ng->SetBinError(bin, corrYieldErr_Ng[bin-1]);
      rawYield_gN->SetBinContent(bin, rawYieldVal_gN[bin-1]);
      rawYield_gN->SetBinError(bin, rawYieldErr_gN[bin-1]);
      rawYield_Ng->SetBinContent(bin, rawYieldVal_Ng[bin-1]);
      rawYield_Ng->SetBinError(bin, rawYieldErr_Ng[bin-1]);
      
      corrYield_ref_gN->SetBinContent(bin, corrYieldVal_gN[refBins[i]]);
      corrYield_ref_gN->SetBinError(bin, corrYieldErr_gN[refBins[i]]);
      corrYield_ref_Ng->SetBinContent(bin, corrYieldVal_Ng[refBins[i]]);
      corrYield_ref_Ng->SetBinError(bin, corrYieldErr_Ng[refBins[i]]);
      rawYield_ref_gN->SetBinContent(bin, rawYieldVal_gN[refBins[i]]);
      rawYield_ref_gN->SetBinError(bin, rawYieldErr_gN[refBins[i]]);
      rawYield_ref_Ng->SetBinContent(bin, rawYieldVal_Ng[refBins[i]]);
      rawYield_ref_Ng->SetBinError(bin, rawYieldErr_Ng[refBins[i]]);
      
      effEvt_gN->SetBinContent(bin, effEvtVal_gN[bin-1]);
      effEvt_gN->SetBinError(bin, effEvtErr_gN[bin-1]);
      effEvt_Ng->SetBinContent(bin, effEvtVal_Ng[bin-1]);
      effEvt_Ng->SetBinError(bin, effEvtErr_Ng[bin-1]);
      effD_gN->SetBinContent(bin, effDVal_gN[bin-1]);
      effD_gN->SetBinError(bin, effDErr_gN[bin-1]);
      effD_Ng->SetBinContent(bin, effDVal_Ng[bin-1]);
      effD_Ng->SetBinError(bin, effDErr_Ng[bin-1]);
    }
    corrYield_gN->Sumw2();
    corrYield_Ng->Sumw2();
    rawYield_gN->Sumw2();
    rawYield_Ng->Sumw2();
    corrYield_ref_gN->Sumw2();
    corrYield_ref_Ng->Sumw2();
    rawYield_ref_gN->Sumw2();
    rawYield_ref_Ng->Sumw2();
    effEvt_gN->Sumw2();
    effEvt_Ng->Sumw2();
    effD_gN->Sumw2();
    effD_Ng->Sumw2();
    
    effEvtRatio->Divide(effEvt_gN, effEvt_Ng);
    effDRatio->Divide(effD_gN, effD_Ng);
    
    TH1D* ratioTemplate  = new TH1D(
      "ratioTemplate",
      "",
      nDBins, DBins);
    ratioTemplate->SetTitle(Form(
      "; %s; Scan / [Yield at %.1f]",
      DCutXaxis.c_str(), 0.5*(DBins[refBins[i]] + DBins[refBins[i]+1])));
    
    TH1D* ratio_corrYield_gN = (TH1D*) corrYield_gN->Clone("ratio_corrYield_gN");
    ratio_corrYield_gN->SetTitle(ratioTemplate->GetTitle());
    ratio_corrYield_gN->Sumw2();
    ratio_corrYield_gN->Divide(corrYield_ref_gN);
    TH1D* ratio_corrYield_Ng = (TH1D*) corrYield_Ng->Clone("ratio_corrYield_Ng");
    ratio_corrYield_Ng->SetTitle(ratioTemplate->GetTitle());
    ratio_corrYield_Ng->Sumw2();
    ratio_corrYield_Ng->Divide(corrYield_ref_Ng);
    
    TH1D* ratio_rawYield_gN = (TH1D*) rawYield_gN->Clone("ratio_rawYield_gN");
    ratio_rawYield_gN->SetTitle(ratioTemplate->GetTitle());
    ratio_rawYield_gN->Sumw2();
    ratio_rawYield_gN->Divide(rawYield_ref_gN);
    TH1D* ratio_rawYield_Ng = (TH1D*) rawYield_Ng->Clone("ratio_rawYield_Ng");
    ratio_rawYield_Ng->SetTitle(ratioTemplate->GetTitle());
    ratio_rawYield_Ng->Sumw2();
    ratio_rawYield_Ng->Divide(rawYield_ref_Ng);
    
    TH1D* ratio_corrYields = (TH1D*) corrYield_gN->Clone("ratio_corrYields");
    ratio_corrYields->Divide(corrYield_Ng);

    // Styling gammaN
    corrYield_gN->SetMarkerColor(kAzure+2);
    corrYield_gN->SetLineColor(kAzure+2);
    corrYield_gN->SetLineWidth(2);
    corrYield_gN->SetMarkerStyle(20);
    ratio_corrYield_gN->SetMarkerColor(kAzure+2);
    ratio_corrYield_gN->SetLineColor(kAzure+2);
    ratio_corrYield_gN->SetLineWidth(2);
    rawYield_gN->SetMarkerColor(kAzure-9);
    rawYield_gN->SetLineColor(kAzure-9);
    rawYield_gN->SetLineWidth(1);
    rawYield_gN->SetMarkerStyle(24);
    ratio_rawYield_gN->SetMarkerColor(kAzure-9);
    ratio_rawYield_gN->SetLineColor(kAzure-9);
    ratio_rawYield_gN->SetLineWidth(1);
    
    // Styling Ngamma
    corrYield_Ng->SetMarkerColor(kPink-8);
    corrYield_Ng->SetLineColor(kPink-8);
    corrYield_Ng->SetLineWidth(2);
    corrYield_Ng->SetMarkerStyle(20);
    ratio_corrYield_Ng->SetMarkerColor(kPink-8);
    ratio_corrYield_Ng->SetLineColor(kPink-8);
    ratio_corrYield_Ng->SetLineWidth(2);
    rawYield_Ng->SetMarkerColor(kPink+1);
    rawYield_Ng->SetLineColor(kPink+1);
    rawYield_Ng->SetLineWidth(1);
    rawYield_Ng->SetMarkerStyle(24);
    ratio_rawYield_Ng->SetMarkerColor(kPink+1);
    ratio_rawYield_Ng->SetLineColor(kPink+1);
    ratio_rawYield_Ng->SetLineWidth(1);
    
    ratio_corrYields->SetMarkerColor(kTeal+2);
    ratio_corrYields->SetLineColor(kTeal+2);
    ratio_corrYields->SetLineWidth(2);
    ratio_corrYields->SetMarkerStyle(33);
    ratio_corrYields->SetMarkerSize(1.4);
    
    // Scaling plots
    corrYieldTemplate->SetMinimum(corrYieldMin);
    corrYieldTemplate->SetMaximum(corrYieldMax);
    ratioTemplate->SetMinimum(ratioMin);
    ratioTemplate->SetMaximum(ratioMax);
    rawYield_Ng->Scale(
      (corrYieldMax - corrYieldMin) / (rawYieldMax - rawYieldMin));
    rawYield_gN->Scale(
      (corrYieldMax - corrYieldMin) / (rawYieldMax - rawYieldMin));
    
    TCanvas* canvas = new TCanvas("canvas", "", 600, 900);
    TPad* padTop = new TPad("padTop", "", 0.0, 0.3, 1.0, 1.0);
    TPad* padBot = new TPad("padBot", "", 0.0, 0.0, 1.0, 0.3);
    padTop->SetMargin(0.14, 0.10, 0.00, 0.16);
    padBot->SetMargin(0.14, 0.10, 0.30, 0.00);
    padTop->Draw();
    padBot->Draw();
    
    padTop->cd();
    corrYieldTemplate->Draw();
    corrYieldTemplate->SetLabelSize(0.025/0.7, "Y");
    corrYieldTemplate->SetTitleSize(0.03/0.7, "Y");
    corrYieldTemplate->SetTitleOffset(1.2, "Y");
    TGaxis *Raxis = new TGaxis(
      corrYield_gN->GetBinLowEdge(corrYield_gN->GetNbinsX() + 1),
      corrYieldMin,
      corrYield_gN->GetBinLowEdge(corrYield_gN->GetNbinsX() + 1),
      corrYieldMax,
      rawYieldMin,
      rawYieldMax,
      510, "+L"
    );
    Raxis->SetLineColor(kGray+1);
    Raxis->SetLabelFont(corrYieldTemplate->GetYaxis()->GetLabelFont());
    Raxis->SetLabelSize(corrYieldTemplate->GetYaxis()->GetLabelSize());
    Raxis->SetLabelColor(kGray+1);
    Raxis->SetTitle("Raw Yield");
    Raxis->SetTitleFont(corrYieldTemplate->GetYaxis()->GetTitleFont());
    Raxis->SetTitleSize(0.03/0.7);
    Raxis->SetTitleOffset(1.2);
    Raxis->SetTitleColor(kGray+1);
    Raxis->Draw();
    rawYield_Ng->Draw("same");
    rawYield_gN->Draw("same");
    corrYield_Ng->Draw("same");
    corrYield_gN->Draw("same");
    gStyle->SetOptStat(0);
    
    padBot->cd();
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
    ratio_rawYield_Ng->Draw("same");
    ratio_rawYield_gN->Draw("same");
    ratio_corrYield_Ng->Draw("same");
    ratio_corrYield_gN->Draw("same");
    ratio_corrYields->Draw("same");
    gStyle->SetOptStat(0);
    
    canvas->cd();
    canvas->Update();
    TLegend* legend = new TLegend(0.25, 0.74, 0.75, 0.88);
    legend->SetTextSize(0.022/0.7);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->AddEntry(corrYield_gN, "Corrected Yield, #gammaN", "lp");
    legend->AddEntry(corrYield_Ng, "Corrected Yield, N#gamma", "lp");
    legend->AddEntry(rawYield_gN, "Raw Yield, #gammaN", "lp");
    legend->AddEntry(rawYield_Ng, "Raw Yield, N#gamma", "lp");
    legend->AddEntry(ratio_corrYields, "Ratio of Corr. Yields (#gammaN / N#gamma)", "lp");
    legend->Draw();
    
    canvas->SaveAs(Form("plot/SystEval_CutScans/%s_y%.0f-%.0f.pdf",
        plotTag.c_str(), yBins[i], yBins[i+1]));
  
    effEvtRatio->SetMarkerColor(kOrange+2);
    effEvtRatio->SetLineColor(kOrange+2);
    effEvtRatio->SetLineWidth(2);
    effDRatio->SetMarkerColor(kViolet+2);
    effDRatio->SetLineColor(kViolet+2);
    effDRatio->SetLineWidth(2);
  
    effEvtRatio->SetMinimum(0.99);
    effEvtRatio->SetMaximum(1.01);
    effDRatio->SetMinimum(0.9);
    effDRatio->SetMaximum(1.1);
  
    TCanvas* canvas2 = new TCanvas("canvas2", "", 600, 900);
    TPad* padTop2 = new TPad("padTop2", "", 0.0, 0.5, 1.0, 1.0);
    TPad* padBot2 = new TPad("padBot2", "", 0.0, 0.0, 1.0, 0.5);
    padTop2->SetMargin(0.12, 0.12, 0.00, 0.20);
    padBot2->SetMargin(0.12, 0.12, 0.20, 0.00);
    padTop2->Draw();
    padBot2->Draw();
    
    padTop2->cd();
    effEvtRatio->Draw();
    effEvtRatio->SetLabelSize(0.025/0.5, "Y");
    effEvtRatio->SetTitleSize(0.03/0.5, "Y");
    effEvtRatio->SetTitleOffset(1.0, "Y");
    TLine* unityTop = new TLine(
      effEvtRatio->GetBinLowEdge(1), 1.0,
      effEvtRatio->GetBinLowEdge(effEvtRatio->GetNbinsX()+1), 1.0);
    unityTop->SetLineColor(kGray);
    unityTop->SetLineWidth(1);
    unityTop->SetLineStyle(9);
    unityTop->Draw();
    effEvtRatio->Draw("same");
    gStyle->SetOptStat(0);
  
    padBot2->cd();
    effDRatio->Draw();
    effDRatio->SetLabelSize(0.025/0.5, "XY");
    effDRatio->SetTitleSize(0.03/0.5, "XY");
    effDRatio->SetTitleOffset(1.0, "Y");
    TLine* unityBot = new TLine(
      effDRatio->GetBinLowEdge(1), 1.0,
      effDRatio->GetBinLowEdge(effDRatio->GetNbinsX()+1), 1.0);
    unityBot->SetLineColor(kGray);
    unityBot->SetLineWidth(1);
    unityBot->SetLineStyle(9);
    unityBot->Draw();
    effDRatio->Draw("same");
    gStyle->SetOptStat(0);
  
    canvas2->cd();
    canvas2->Update();
    TLegend* legend2 = new TLegend(0.25, 0.82, 0.75, 0.88);
    legend2->SetTextSize(0.022/0.5);
    legend2->SetFillStyle(0);
    legend2->SetBorderSize(0);
    legend2->AddEntry(effEvtRatio, "Ratio of Event Eff. (#gammaN / N#gamma)", "l");
    legend2->AddEntry(effDRatio, "Ratio of D^{0} Eff. (#gammaN / N#gamma)", "l");
    legend2->Draw();
  
    canvas2->SaveAs(Form("plot/SystEval_CutScans/%s_effRatios_y%.0f-%.0f.pdf",
        plotTag.c_str(), yBins[i], yBins[i+1]));
    
    delete canvas;
    delete corrYieldTemplate;
    delete rawYieldTemplate;
    delete ratioTemplate;
    delete corrYield_gN;
    delete corrYield_Ng;
    delete corrYield_ref_gN;
    delete corrYield_ref_Ng;
    delete ratio_corrYield_gN;
    delete ratio_corrYield_Ng;
    delete ratio_corrYields;
    delete rawYield_gN;
    delete rawYield_Ng;
    delete rawYield_ref_gN;
    delete rawYield_ref_Ng;
    delete ratio_rawYield_gN;
    delete ratio_rawYield_Ng;
  }
}

void plotSystematicsEval()
{
  int isGammaN[2] = {1, 0};
  const int nPtBins = 1;
  double ptBins[nPtBins+1] = {2, 5};
  const int nYBins = 4;
  double yBins[nYBins+1]  = {-2, -1, 0, 1, 2};
  
  system("mkdir -p plot/SystEval_MassFit/");
  system("mkdir -p plot/SystEval_Cuts/");
  system("mkdir -p plot/SystEval_CutScans/");
  
  // Iterate through gammaN options
  vector<string> fullAnalysis_gN;
  vector<string> fullAnalysis_Ng;
  vector<string> systDalpha_gN;
  vector<string> systDalpha_Ng;
  vector<string> systDchi2cl_gN;
  vector<string> systDchi2cl_Ng;
  vector<string> systDsvpv_gN;
  vector<string> systDsvpv_Ng;
  vector<string> systDtrkPt_gN;
  vector<string> systDtrkPt_Ng;
  vector<string> systSubdir_gN;
  vector<string> systSubdir_Ng;
  for (const int& gammaN : isGammaN) {
    // Iterate through pt bins
    for (int ipt = 0; ipt < nPtBins; ipt++) {
      double ptBinLo = ptBins[ipt];
      double ptBinHi = ptBins[ipt+1];
      // Iterate through y bins to make list of source directories
      vector<string> fullAnalysis;
      vector<string> systDalpha;
      vector<string> systDchi2cl;
      vector<string> systDsvpv;
      vector<string> systDtrkPt;
      vector<string> systRapGapLoose;
      vector<string> systRapGapTight;
      vector<string> systSubdir;
      for (int iy = 0; iy < nYBins; iy++) {
        double yBinLo = yBins[iy];
        double yBinHi = yBins[iy+1];
        string subdir = Form(
          "pt%.0f-%.0f_y%.0f-%.0f_IsGammaN%d",
          ptBinLo, ptBinHi, yBinLo, yBinHi, gammaN
        );
        fullAnalysis.push_back("fullAnalysis/" + subdir + "/");
        systDalpha.push_back("systDalpha/" + subdir + "/");
        systDchi2cl.push_back("systDchi2cl/" + subdir + "/");
        systDsvpv.push_back("systDsvpv/" + subdir + "/");
        systDtrkPt.push_back("systDtrkPt/" + subdir + "/");
        systRapGapLoose.push_back("systRapGapLoose/" + subdir + "/");
        systRapGapTight.push_back("systRapGapTight/" + subdir + "/");
        systSubdir.push_back(subdir);
        if (gammaN) {
          fullAnalysis_gN.push_back("fullAnalysis/" + subdir + "/");
          systDalpha_gN.push_back("systDalpha/" + subdir + "/");
          systDchi2cl_gN.push_back("systDchi2cl/" + subdir + "/");
          systDsvpv_gN.push_back("systDsvpv/" + subdir + "/");
          systDtrkPt_gN.push_back("systDtrkPt/" + subdir + "/");
          systSubdir_gN.push_back(subdir);
        }
        else {
          fullAnalysis_Ng.push_back("fullAnalysis/" + subdir + "/");
          systDalpha_Ng.push_back("systDalpha/" + subdir + "/");
          systDchi2cl_Ng.push_back("systDchi2cl/" + subdir + "/");
          systDsvpv_Ng.push_back("systDsvpv/" + subdir + "/");
          systDtrkPt_Ng.push_back("systDtrkPt/" + subdir + "/");
          systSubdir_Ng.push_back(subdir);
        }
      }
      string plotLabel = Form("pt%.0f-%.0f_IsGammaN%d", ptBinLo, ptBinHi, gammaN);
      string plotTitle = Form("%.1f #leq D_{p_{T}} < %.1f (GeV), %s",
        ptBins[ipt], ptBins[ipt+1], (gammaN == 1 ? "#gammaN" : "N#gamma")
      );
//      massfit_rawYield_vs_y(fullAnalysis, plotTitle, plotLabel, nYBins, yBins, ptBinLo, ptBinHi, gammaN);
//      massfit_lambda_vs_y(fullAnalysis, plotTitle, plotLabel, nYBins, yBins, ptBinLo, ptBinHi, gammaN);
//      massfit_rawYield_vs_m(fullAnalysis, plotTitle, plotLabel, nYBins, yBins, ptBinLo, ptBinHi, gammaN);
//      cuts_rawYield_vs_y(
//        fullAnalysis, systDalpha, systDchi2cl, systDsvpv, systDtrkPt, systRapGapLoose, systRapGapTight,
//        plotTitle, plotLabel, nYBins, yBins, ptBinLo, ptBinHi, gammaN);
//      cuts_effEvt_vs_y(
//        fullAnalysis, systDalpha, systDchi2cl, systDsvpv, systDtrkPt, systRapGapLoose, systRapGapTight,
//        plotTitle, plotLabel, nYBins, yBins, ptBinLo, ptBinHi, gammaN);
//      cuts_effD_vs_y(
//        fullAnalysis, systDalpha, systDchi2cl, systDsvpv, systDtrkPt, systRapGapLoose, systRapGapTight,
//        plotTitle, plotLabel, nYBins, yBins, ptBinLo, ptBinHi, gammaN);
//      cuts_lambda_vs_y(
//        fullAnalysis, systDalpha, systDchi2cl, systDsvpv, systDtrkPt, systRapGapLoose, systRapGapTight,
//        plotTitle, plotLabel, nYBins, yBins, ptBinLo, ptBinHi, gammaN);
//      cuts_rawYield_vs_m(
//        fullAnalysis, systDalpha, systDchi2cl, systDsvpv, systDtrkPt, systRapGapLoose, systRapGapTight,
//        plotTitle, plotLabel, nYBins, yBins, ptBinLo, ptBinHi, gammaN);
//  cuts_corrYield_vs_DtrkPtCut(
//    fullAnalysis_gN, fullAnalysis_Ng,
//    systSubdir_gN, systSubdir_Ng,
//    plotTitle, plotLabel, nYBins, yBins, ptBinLo, ptBinHi
//  );
//  cuts_corrYield_vs_DsvpvCut(
//    fullAnalysis_gN, fullAnalysis_Ng,
//    systSubdir_gN, systSubdir_Ng,
//    plotTitle, plotLabel, nYBins, yBins, ptBinLo, ptBinHi
//  );
//  cuts_corrYield_vs_DalphaCut(
//    fullAnalysis_gN, fullAnalysis_Ng,
//    systSubdir_gN, systSubdir_Ng,
//    plotTitle, plotLabel, nYBins, yBins, ptBinLo, ptBinHi
//  );
//  cuts_corrYield_vs_Dchi2Cut(
//    fullAnalysis_gN, fullAnalysis_Ng,
//    systSubdir_gN, systSubdir_Ng,
//    plotTitle, plotLabel, nYBins, yBins, ptBinLo, ptBinHi
//  );
    }
  }
  const int nDtrkPtBins = 6;
  double DtrkPtBins[nDtrkPtBins+1] = {
    0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25};
  string DtrkPtLabels[nDtrkPtBins+1] = {
    "_070", "_080", "_090", "_100", "_110", "_120"};
  int DtrkPtRefBins[nYBins] = {3, 3, 3, 3};
  cuts_corrYield_rawYield_vs_DCut(
    systSubdir_gN,
    systSubdir_Ng,
    "2.0 #leq D_{p_{T}} < 5.0 (GeV)", // plot title
    "pt2-5_yields_vs_DtrkPt", // plot tag (for filename)
    "systDtrkPt", // DCut dir
    "DtrkPt Cut (GeV)", // DCut X-axis label
    nDtrkPtBins, // nDCutBins
    DtrkPtBins, // DCutBins
    DtrkPtLabels, // DCutLabels
    nYBins,
    yBins,
    DtrkPtRefBins
  );
  
  const int nDsvpvBins = 5;
  double DsvpvBins[nDsvpvBins+1] = {
    1.875, 2.125, 2.375, 2.635, 2.875, 3.125};
  string DsvpvLabels[nDsvpvBins+1] = {
    "_200", "_225", "_250", "_275", "_300"};
  int DsvpvRefBins[nYBins] = {2, 2, 2, 2};
  cuts_corrYield_rawYield_vs_DCut(
    systSubdir_gN,
    systSubdir_Ng,
    "2.0 #leq D_{p_{T}} < 5.0 (GeV)", // plot title
    "pt2-5_yields_vs_DsvpvSig", // plot tag (for filename)
    "systDsvpv", // DCut dir
    "Dsvpv Significance Cut", // DCut X-axis label
    nDsvpvBins, // nDCutBins
    DsvpvBins, // DCutBins
    DsvpvLabels, // DCutLabels
    nYBins,
    yBins,
    DsvpvRefBins
  );
  
  const int nDalphaBins = 5;
  double DalphaBins[nDalphaBins+1] = {
    0.15, 0.25, 0.35, 0.45, 0.55, 0.65};
  string DalphaLabels[nDalphaBins+1] = {
    "_020", "_030", "_040", "_050", "_060"};
  int DalphaRefBins[nYBins] = {1, 4, 4, 1};
  cuts_corrYield_rawYield_vs_DCut(
    systSubdir_gN,
    systSubdir_Ng,
    "2.0 #leq D_{p_{T}} < 5.0 (GeV)", // plot title
    "pt2-5_yields_vs_Dalpha", // plot tag (for filename)
    "systDalpha", // DCut dir
    "Dalpha Cut (GeV)", // DCut X-axis label
    nDalphaBins, // nDCutBins
    DalphaBins, // DCutBins
    DalphaLabels, // DCutLabels
    nYBins,
    yBins,
    DalphaRefBins
  );
  
  const int nDchi2clBins = 6;
  double Dchi2clBins[nDchi2clBins+1] = {
    0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325};
  string Dchi2clLabels[nDchi2clBins+1] = {
    "_005", "_010", "_015", "_020", "_025", "_030"};
  int Dchi2clRefBins[nYBins] = {1, 1, 1, 1};
  cuts_corrYield_rawYield_vs_DCut(
    systSubdir_gN,
    systSubdir_Ng,
    "2.0 #leq D_{p_{T}} < 5.0 (GeV)", // plot title
    "pt2-5_yields_vs_Dchi2cl", // plot tag (for filename)
    "systDchi2cl", // DCut dir
    "Dchi2cl Cut", // DCut X-axis label
    nDchi2clBins, // nDCutBins
    Dchi2clBins, // DCutBins
    Dchi2clLabels, // DCutLabels
    nYBins,
    yBins,
    Dchi2clRefBins
  );
}

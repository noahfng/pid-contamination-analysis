#include <algorithm>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TString.h"

#include "AddTrees.h"

void nSigma_Plot_dx() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    const Char_t* baseDir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";
    TChain chain("twotauchain");
    AddTrees(chain, baseDir);

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);
    const Char_t* subs[5] = {"El","Mu","Pi","Ka","Pr"};
    for (Int_t i = 0; i < 5; ++i) {
        chain.SetBranchStatus(Form("fTrkTPCnSigma%s", subs[i]), 1);
    }

    Float_t inner[2];
    Float_t tpcNS[5][2];
    chain.SetBranchAddress("fTrkTPCinnerParam", inner);
    for (Int_t i = 0; i < 5; ++i) {
        chain.SetBranchAddress(Form("fTrkTPCnSigma%s", subs[i]), tpcNS[i]);
    }

    const Int_t   nBins  = 200;
    const Double_t pMin  = 0.55, pMax = 0.70;
    const Double_t xMin  = -15.0, xMax = 15.0;
    const Int_t   nParts = 5;
    const TString names[nParts]  = {"e","#mu","#pi","K","p"};
    const Int_t   colors[nParts] = {kBlue, kRed, kGreen+2, kOrange+7, kViolet};


    TH1F* hRes[nParts];
    TH1F* hDeriv[nParts];
    for (Int_t h = 0; h < nParts; ++h) {

        hRes[h] = new TH1F(
            Form("res_%s", names[h].Data()),
            Form("n#sigma_{%s} for %.2f<p<%.2f; n#sigma; Counts",
                 names[h].Data(), pMin, pMax),
            nBins, xMin, xMax
        );
        hDeriv[h] = (TH1F*)hRes[h]->Clone(Form("deriv_%s", names[h].Data()));
        hDeriv[h]->SetTitle(
            Form("d/dx n#sigma_{%s} for %.2f<p<%.2f; n#sigma; dN/dx",
                 names[h].Data(), pMin, pMax)
        );
        hDeriv[h]->SetLineColor(colors[h]);
        hDeriv[h]->SetMarkerColor(colors[h]);
        hDeriv[h]->SetMarkerStyle(kFullCircle);
        hDeriv[h]->SetMarkerSize(0.8);
        hDeriv[h]->GetYaxis()->SetRangeUser(-10, 10);
    }

    Long64_t nEntries = std::min(chain.GetEntries(), static_cast<Long64_t>(1e6));
    for (Long64_t i = 0; i < nEntries; ++i) {
        chain.GetEntry(i);
        for (Int_t tr = 0; tr < 2; ++tr) {
            Float_t pG = inner[tr];
            if (pG < pMin || pG > pMax) continue;
            for (Int_t h = 0; h < nParts; ++h) {
                Float_t val = tpcNS[h][tr];
                if (!TMath::IsNaN(val)) {
                    hRes[h]->Fill(val);
                }
            }
        }
    }

    for (Int_t h = 0; h < nParts; ++h) {
        Int_t nb = hRes[h]->GetNbinsX();
        for (Int_t ib = 1; ib <= nb; ++ib) {
            if (ib == 1 || ib == nb) {
                hDeriv[h]->SetBinContent(ib, 0.);
            } else {
                Double_t xLo = hRes[h]->GetBinCenter(ib - 1);
                Double_t xHi = hRes[h]->GetBinCenter(ib + 1);
                Double_t yLo = hRes[h]->GetBinContent(ib - 1);
                Double_t yHi = hRes[h]->GetBinContent(ib + 1);
                hDeriv[h]->SetBinContent(ib, (yHi - yLo)/(xHi - xLo));
            }
        }
    }

    TCanvas* c = new TCanvas("c","d/dx n#sigma", 800, 600);
    c->SetLeftMargin(0.15);
    c->SetGrid();
    c->Print("nSigmaTPC_derivative.pdf[");
    for (Int_t h = 0; h < nParts; ++h) {
        c->Clear();
        hDeriv[h]->Draw("E1"); 
        c->Print("nSigmaTPC_derivative.pdf");
    }
    c->Print("nSigmaTPC_derivative.pdf]");
}

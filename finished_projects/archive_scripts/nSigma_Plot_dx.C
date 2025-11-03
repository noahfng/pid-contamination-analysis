#include <algorithm>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TString.h"

#include <AddTrees.h>
#include <helpers.h>

void nSigma_Plot_dx() {
    auto help = new helper();
    const Int_t nParts = helper::nParts;
    const Int_t NtrkMax = help->NtrkMax;
    const Int_t   nBins  = 200;
    const Double_t pMin  = 0.55, pMax = 0.70;
    const Double_t xMin  = -15.0, xMax = 15.0;
    const Double_t nEntriesLimit = 1e7;
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    TChain chain("twotauchain");
    AddTrees(chain, help->base_dir);

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);
    for (Int_t i = 0; i < nParts; ++i) {
        chain.SetBranchStatus(Form("fTrkTPCnSigma%s", help->pNames[i]), 1);
    }
    std::vector<Float_t> inner(NtrkMax);
    std::vector<std::vector<Float_t>> tpcNS(nParts, std::vector<Float_t>(NtrkMax));
    chain.SetBranchAddress("fTrkTPCinnerParam", inner.data());
    for (Int_t i = 0; i < nParts; ++i) {
        chain.SetBranchAddress(Form("fTrkTPCnSigma%s", help->pNames[i]), tpcNS[i].data());
    }


    TH1D* hRes[nParts];
    TH1D* hDeriv[nParts];
    for (Int_t h = 0; h < nParts; ++h) {

        hRes[h] = new TH1D(
            Form("res_%s", help->pCodes[h]),
            Form("n#sigma_{%s} for %.2f<p<%.2f; n#sigma; Counts",
                 help->pCodes[h], pMin, pMax),
            nBins, xMin, xMax
        );
        hDeriv[h] = (TH1D*)hRes[h]->Clone(Form("deriv_%s", help->pCodes[h]));
        hDeriv[h]->SetTitle(
            Form("d/dx n#sigma_{%s} for %.2f<p<%.2f; n#sigma; dN/dx",
                 help->pCodes[h], pMin, pMax)
        );
        hDeriv[h]->SetLineColor(help->colors[h]);
        hDeriv[h]->SetMarkerColor(help->colors[h]);
        hDeriv[h]->SetMarkerStyle(kFullCircle);
        hDeriv[h]->SetMarkerSize(0.8);
        hDeriv[h]->GetYaxis()->SetRangeUser(-10, 10);
    }

    Long64_t nEntries = std::min(chain.GetEntries(), static_cast<Long64_t>(nEntriesLimit));
    for (Long64_t i = 0; i < nEntries; ++i) {
        chain.GetEntry(i);
        for (Int_t tr = 0; tr < NtrkMax; ++tr) {
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

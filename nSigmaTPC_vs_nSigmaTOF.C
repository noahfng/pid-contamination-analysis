#include <TChain.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>

#include "AddTrees.h"

void nSigmaTPC_vs_nSigmaTOF()
{   gStyle->SetPalette(kRainBow);
    const Char_t* baseDir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";
    const Int_t    nBins = 200;
    const Double_t xMin  = -20, xMax = 20;
    const Double_t yMin  = -20, yMax = 20;
    const Double_t nEntriesLimit = 1e7;
    TChain chain("twotauchain");
    AddTrees(chain, baseDir);

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);
    chain.SetBranchStatus("fTrkTOFexpMom",     1);
    const Char_t* subs[5] = {"El","Mu","Pi","Ka","Pr"};
    for (Int_t i = 0; i < 5; ++i) {
        chain.SetBranchStatus(Form("fTrkTPCnSigma%s", subs[i]), 1);
        chain.SetBranchStatus(Form("fTrkTOFnSigma%s", subs[i]), 1);
    }

    Float_t inner[2], tofExpMom[2];
    Float_t tpcNS[5][2], tofNS[5][2];
    chain.SetBranchAddress("fTrkTPCinnerParam", inner);
    chain.SetBranchAddress("fTrkTOFexpMom",      tofExpMom);
    for (Int_t i = 0; i < 5; ++i) {
        chain.SetBranchAddress(Form("fTrkTPCnSigma%s", subs[i]), tpcNS[i]);
        chain.SetBranchAddress(Form("fTrkTOFnSigma%s", subs[i]), tofNS[i]);
    }

    TH2F* h2[5];
    for (Int_t i = 0; i < 5; ++i) {
        h2[i] = new TH2F(Form("h2_%s", subs[i]),
                         Form("n#sigma_{TPC} vs n#sigma_{TOF} (%s);n#sigma_{TPC};n#sigma_{TOF}", subs[i]),
                         nBins, xMin, xMax,
                         nBins, yMin, yMax);
        h2[i]->Sumw2();
    }

    Long64_t nEntries = std::min(chain.GetEntries(), static_cast<Long64_t>(nEntriesLimit));
    for (Long64_t ev = 0; ev < nEntries; ++ev) {
        chain.GetEntry(ev);
        for (Int_t t = 0; t < 2; ++t) {
            if (tofExpMom[t] < 0) continue;
            for (Int_t i = 0; i < 5; ++i) {
                Double_t ntpc = tpcNS[i][t];
                Double_t ntof = tofNS[i][t];
                if (!TMath::IsNaN(ntpc) && !TMath::IsNaN(ntof)) {
                    h2[i]->Fill(ntpc, ntof);
                }
            }
        }
    }

    TCanvas* c = new TCanvas("c2d","nSigma TPC vs TOF",800,700);
    gStyle->SetOptStat(0);
    for (Int_t i = 0; i < 5; ++i) {
        c->Clear();
        h2[i]->Draw("COLZ");
        c->SetLogz();
        c->SetRightMargin(0.15);
        c->SaveAs(Form("nSigma2D_%s.pdf", subs[i]));
    }
    delete c;
}
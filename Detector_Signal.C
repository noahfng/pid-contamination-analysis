#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"

#include <AddTrees.h>
#include <get_expected_signal.h>

void Detector_Signal() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1);

    const Bool_t plotTPC       = true;
    const Bool_t plotTOF       = true;
    const Bool_t applyNorm     = false;
    const Double_t pStart = 0.1, pEnd = 1.5, step = 0.1;
    const Int_t    nSteps = Int_t(std::floor((pEnd - pStart) / step + 0.5));
    const Int_t   nBins = 200;
    const Double_t xMin  =   0, xMax = 250;
    const Double_t nEntriesMax = 1e7;

    const char* baseDir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";
    TChain chain("twotauchain");
    AddTrees(chain, baseDir);

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkTPCsignal", 1);
    chain.SetBranchStatus("fTrkTOFsignal", 1);
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);

    Float_t TPCsignal[2], TOFsignal[2], inner[2];
    chain.SetBranchAddress("fTrkTPCsignal", TPCsignal);
    chain.SetBranchAddress("fTrkTOFsignal", TOFsignal);
    chain.SetBranchAddress("fTrkTPCinnerParam", inner);

    Long64_t nEntries = std::min(chain.GetEntries(), static_cast<Long64_t>(nEntriesMax));

    std::vector<TH1F*> hTPCs, hTOFs;
    if (plotTPC) hTPCs.reserve(nSteps);
    if (plotTOF) hTOFs.reserve(nSteps);


    for (Int_t i = 0; i < nSteps; ++i) {
        Double_t low  = pStart + i*step;
        Double_t high = low + step;
        if (plotTPC) {
            hTPCs.push_back(
              new TH1F(Form("hTPC_%.1f_%.1f", low, high),
                       Form("TPC signal; TPC Signal ; Entries"),
                       nBins, xMin, xMax)
            );
            hTPCs.back()->SetMarkerColor(kBlack);
            hTPCs.back()->SetLineColor(kBlack);
            hTPCs.back()->SetMarkerSize(0.75);
            hTPCs.back()->SetMarkerStyle(kFullCircle);
        }
        if (plotTOF) {
            hTOFs.push_back(
              new TH1F(Form("hTOF_%.1f_%.1f", low, high),
                       Form("TOF signal; TOF Signal ; Entries"),
                       nBins, xMin, xMax)
            );
            hTOFs.back()->SetMarkerColor(kBlack);
            hTOFs.back()->SetLineColor(kBlack);
            hTOFs.back()->SetMarkerSize(0.75);
            hTOFs.back()->SetMarkerStyle(kFullCircle);
        }
    }

    for (Long64_t i = 0; i < nEntries; ++i) {
        chain.GetEntry(i);
        for (int trk = 0; trk < 2; ++trk) {
            Double_t p = inner[trk];
            Int_t    ibin = Int_t((p - pStart)/step);
            if (ibin < 0 || ibin >= nSteps) continue;  
            if (plotTPC) hTPCs[ibin]->Fill(TPCsignal[trk]);
            if (plotTOF) hTOFs[ibin]->Fill(TOFsignal[trk]);
        }
    }

    if (applyNorm) {
        for (auto* h : hTPCs) if (h->GetEntries()>0) h->Scale(1./h->GetEntries());
        for (auto* h : hTOFs) if (h->GetEntries()>0) h->Scale(1./h->GetEntries());
    }

    TCanvas* cTPC = nullptr;
    TCanvas* cTOF = nullptr;
    const char* fTPC = "TPCSignalDistribution.pdf";
    const char* fTOF = "TOFSignalDistribution.pdf";

    if (plotTPC) {
        cTPC = new TCanvas("cTPC","TPC slices",800,600);
        cTPC->SetLogy();
        cTPC->Print(Form("%s[", fTPC)); 
    }
    if (plotTOF) {
        cTOF = new TCanvas("cTOF","TOF slices",800,600);
        cTOF->SetLogy();
        cTOF->Print(Form("%s[", fTOF));  
    }

    for (Int_t i = 0; i < nSteps; ++i) {
        Double_t low  = pStart + i*step;
        Double_t high = low + step;

        if (plotTPC) {
            cTPC->Clear();
            hTPCs[i]->SetTitle(Form("TPC: %.2f < p < %.2f GeV/c", low, high));
            hTPCs[i]->Draw("E1");
            cTPC->Print(fTPC);          
        }
        if (plotTOF) {
            cTOF->Clear();
            hTOFs[i]->SetTitle(Form("TOF: %.2f < p < %.2f GeV/c", low, high));
            hTOFs[i]->Draw("E1");
            cTOF->Print(fTOF);          
        }
    }

    if (plotTPC) cTPC->Print(Form("%s]", fTPC));
    if (plotTOF) cTOF->Print(Form("%s]", fTOF));
}
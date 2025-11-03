#include "TChain.h"  
#include "TH2D.h"      
#include "TGraph.h"    
#include "TCanvas.h"    
#include "TLegend.h"   
#include "TMath.h"      
#include "TStyle.h"     
#include "TString.h"

#include <AddTrees.h>   // project-specific: file discovery/chain fill
#include <helpers.h>    // project-specific: masses, charges, colors, Bethe–Bloch

// Plot TPC dE/dx vs momentum and overlay Bethe–Bloch curves.
void dEdx_vs_p() {
    auto help = new helper();
    gStyle->SetPalette(kRainBow);

    // basic config    
    const Double_t nEntriesLimit = 1e7;
    const Int_t nPoints = 500;
    const Bool_t KaExclusion = false; // TOF-based Kaon veto
    const Bool_t PrExclusion = true; // TOF-based Proton veto
    const Bool_t tofFilter = false; // require TOF info
    const Double_t pMin = 0.3, pMax = 5.0; // GeV/c range for axes/curves
    
    const Int_t nParts = helper::nParts;
    const Int_t NtrkMax = help->NtrkMax;
    
    // input data chain
    TChain chain("twotauchain");
    AddTrees(chain, help->base_dir);

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);
    chain.SetBranchStatus("fTrkTPCsignal",     1);
    chain.SetBranchStatus("fTrkTOFnSigmaKa",   1);
    chain.SetBranchStatus("fTrkTOFnSigmaPr",   1);
    chain.SetBranchStatus("fTrkTOFexpMom",     1);

    std::vector<Float_t> inner(NtrkMax);
    std::vector<Float_t> signal(NtrkMax);
    std::vector<Float_t> TOFnSigmaKa(NtrkMax);
    std::vector<Float_t> TOFnSigmaPr(NtrkMax);
    std::vector<Float_t> expMom(NtrkMax);

    chain.SetBranchAddress("fTrkTPCinnerParam", inner.data());
    chain.SetBranchAddress("fTrkTPCsignal",     signal.data());
    chain.SetBranchAddress("fTrkTOFnSigmaKa",   TOFnSigmaKa.data());
    chain.SetBranchAddress("fTrkTOFnSigmaPr",   TOFnSigmaPr.data());
    chain.SetBranchAddress("fTrkTOFexpMom",     expMom.data());

    // histogram
    const Double_t step = (pMax - pMin) / nPoints;
    Long64_t nEntries = std::min(chain.GetEntries(), static_cast<Long64_t>(nEntriesLimit));
    TH2D *hist = new TH2D("dedx_vs_p1",
                          "TPC dE/dx vs p;p [GeV/c];dE/dx [arb.u.]",
                          250, pMin, pMax,
                          100,   0, 120);
    
    // event/track loop
    for (Long64_t i = 0; i < nEntries; ++i) {
        chain.GetEntry(i);
        for (Int_t t = 0; t < NtrkMax; ++t) {
            if (expMom[t] < 0 && tofFilter)
                continue;
            if (inner[t] <= 0 || signal[t] <= 0)
                continue;
            if (KaExclusion && !TMath::IsNaN(TOFnSigmaKa[t]) && TMath::Abs(TOFnSigmaKa[t]) < 3.0) 
                continue;
            if (PrExclusion && !TMath::IsNaN(TOFnSigmaPr[t]) && TMath::Abs(TOFnSigmaPr[t]) < 3.0) 
                continue;
            
            hist->Fill(inner[t], signal[t]);
        }
    }

    // draw heatmap
    TCanvas* c = new TCanvas("c","dE/dx vs p (tracks)",800,600); 
    c->SetLogz(); 
    c->SetLogx();
    c->SetGrid();
    hist->SetStats(false);
    hist->Draw("COLZ");

    // overlay Bethe–Bloch bands
    TLegend *leg = new TLegend(0, 0.10, 0.15, 0.30);
    leg->SetBorderSize(0); 
    leg->SetFillStyle(0);

    for (Int_t i = 0; i < nParts; ++i) {
        TGraph *g = new TGraph();
        g->SetLineColor(help->colors[i]);
        g->SetLineWidth(2);
        Int_t idx=0;
        for (Int_t j = 0; j <= nPoints; ++j) {
            Double_t pG = pMin + j*step; // GeV/c
            Double_t pM = pG * 1000.; // helper expects MeV/c
            Double_t d  = help->getTPCSignal(pM, help->pMasses[i], help->pCharges[i]);
            if (d<0 || TMath::IsNaN(d)) continue;
            g->SetPoint(idx++, pG, d);
        }
        g->Draw("L SAME");
        leg->AddEntry(g, help->pNames[i], "l");
    }
    leg->Draw();

    c->Print("Bethe_Bloch.pdf");
}

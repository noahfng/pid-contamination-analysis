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
    gStyle->SetNumberContours(256);

    // basic config    
    const Double_t nEntriesLimit = 1e10;
    const Int_t nPoints = 500;
    const Bool_t KaExclusion = false; // TOF-based Kaon veto
    const Bool_t PrExclusion = false; // TOF-based Proton veto
    const Bool_t requireTOF = (KaExclusion || PrExclusion); // automatically require TOF info if any TOF-based veto is used
    const Double_t pMin = 0.1, pMax = 10.0; // GeV/c range for axes/curves
    
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
    TH2D *hist = new TH2D("dedx_vs_p1",";#it{p} (GeV/#it{c});#LTd#it{E}/d#it{x}#GT (a.u.)",1000, pMin, pMax, 1000, 20, 120);
    hist->GetXaxis()->SetTitleOffset(1.3);
    
    // event/track loop
    for (Long64_t i = 0; i < nEntries; ++i) {
        chain.GetEntry(i);
        for (Int_t t = 0; t < NtrkMax; ++t) {
            if (expMom[t] < 0 && requireTOF)
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
    c->SetRightMargin(0.12);
    c->SetLogz(); 
    c->SetLogx();
    c->SetGrid(0,0);      // no grid
    c->SetTickx();
    c->SetTicky();
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.14); 
    c->SetTopMargin(0.06);

    hist->SetStats(false);
    hist->Draw("COLZ");

    // overlay Bethe–Bloch bands
    TLegend *leg = new TLegend(0.82, 0.74, 0.86, 0.94);
    leg->SetBorderSize(1); 
    leg->SetMargin(0.45);
    leg->SetFillColorAlpha(kWhite, 0.8);

    for (Int_t i = 0; i < nParts; ++i) {
        TGraph *g = new TGraph();
        Int_t idx=0;
        for (Int_t j = 0; j <= nPoints; ++j) {
            Double_t pG = pMin + j*step; // GeV/c
            Double_t pM = pG * 1000.; // helper expects MeV/c
            Double_t d  = help->getTPCSignal(pM, help->pMasses[i], help->pCharges[i]);
            if (d<0 || TMath::IsNaN(d)) continue;
            g->SetPoint(idx++, pG, d);
        }

        g->SetLineColor(help->colors[i]);   
        g->SetLineWidth(2);
        g->Draw("L SAME");

        leg->AddEntry(g, help->pCodes[i], "l");
    }
    leg->Draw();

    c->Print("Bethe_Bloch.pdf");
}

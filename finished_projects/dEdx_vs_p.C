#include "TChain.h"  
#include "TH2F.h"      
#include "TGraph.h"    
#include "TCanvas.h"    
#include "TLegend.h"   
#include "TMath.h"      
#include "TStyle.h"     
#include "TString.h"

#include <AddTrees.h> 
#include <get_expected_signal.h>

void dEdx_vs_p() {
    gStyle->SetPalette(kRainBow);
    const Char_t* base_dir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";
    TChain chain("twotauchain");
    AddTrees(chain, base_dir);

    chain.SetBranchStatus("*", 0);
    Float_t inner[2];
    Float_t signal[2];
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);
    chain.SetBranchStatus("fTrkTPCsignal",     1);
    chain.SetBranchAddress("fTrkTPCinnerParam", inner);
    chain.SetBranchAddress("fTrkTPCsignal",     signal);

    const Double_t nEntriesLimit = 1e6;
    const Int_t nPoints = 500;
    const Double_t pMin = 0.3, pMax = 5.0;
    const Double_t step = (pMax - pMin) / nPoints;
    Long64_t nEntries = std::min(chain.GetEntries(), static_cast<Long64_t>(nEntriesLimit));
    TH2F *hist = new TH2F("dedx_vs_p1",
                          "TPC dE/dx vs p;p [GeV/c];dE/dx [arb.u.]",
                          250, pMin, pMax,
                          100,   0, 120);
    for (Int_t i=0; i<nEntries; ++i) {
        chain.GetEntry(i);
        if (inner[0]>0 && signal[0]>0) hist->Fill(inner[0], signal[0]);
        if (inner[1]>0 && signal[1]>0) hist->Fill(inner[1], signal[1]);
    }

    const Int_t nParts = 5;
    const TString names[nParts]    = {"e", "#mu", "#pi", "K", "p"};
    const Double_t masses[nParts]  = {
        0.51099895,   // e
        105.6583755,  // μ
        139.57039,    // π
        493.677,      // K
        938.27208816  // p
    };
    const Double_t charges[nParts] = {1, 1, 1, 1, 1};
    const Int_t colors[nParts] = {kBlue, kGreen+2, kOrange+7, kMagenta+2, kCyan+1};

    TCanvas* c = new TCanvas("c","dE/dx vs p (tracks)",800,600); 
    c->SetLogz(); 
    c->SetLogx();
    c->SetGrid();
    hist->SetStats(false);
    hist->Draw("COLZ");

    TLegend *leg = new TLegend(0, 0.10, 0.15, 0.30);
    leg->SetBorderSize(0); 
    leg->SetFillStyle(0);

    for (Int_t i=0; i<nParts; ++i) {
        TGraph *g = new TGraph();
        g->SetLineColor(colors[i]);
        g->SetLineWidth(2);
        Int_t idx=0;
        for (Int_t j=0; j<=nPoints; ++j) {
            Double_t pG = pMin + j*step;
            Double_t pM = pG * 1000.;
            Double_t d  = get_expected_signal(pM, masses[i], charges[i]);
            if (d<0 || TMath::IsNaN(d)) continue;
            g->SetPoint(idx++, pG, d);
        }
        g->Draw("L SAME");
        leg->AddEntry(g, names[i], "l");
    }
    leg->Draw();

    c->Print("Bethe_Bloch.pdf");
}

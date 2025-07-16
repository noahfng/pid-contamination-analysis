#include <algorithm>
#include <vector>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TF1.h"  
#include "TSpectrum.h"
#include "TSystem.h"

#include <AddTrees.h>

void nSigma_Plot_manual(){
    gROOT->SetBatch(kFALSE);
    gStyle->SetOptStat(1);
    const char* baseDir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";
    TChain chain("twotauchain");
    AddTrees(chain, baseDir);

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);
    const char* subs[5] = {"El","Mu","Pi","Ka","Pr"};
    for (int i = 0; i < 5; ++i) {
        chain.SetBranchStatus(Form("fTrkTPCnSigma%s", subs[i]), 1);}
    Float_t inner[2];
    Float_t tpcNS[5][2];
    chain.SetBranchAddress("fTrkTPCinnerParam", inner);
    for (int i = 0; i < 5; ++i) {
        chain.SetBranchAddress(Form("fTrkTPCnSigma%s", subs[i]), tpcNS[i]);}

    const Int_t nBins = 200;
    const Double_t pMin = 0.45, pMax = 0.6; 
    const Double_t xMin = -15.0, xMax = 15.0;
    const Int_t nParts = 5;
    const TString names[nParts] = {"e","#mu","#pi","K","p"};
    const Int_t colors[nParts] = {kBlue, kGreen+2, kOrange+7, kMagenta+2, kCyan+1};

    TH1F *hRes[nParts];
    for (int h = 0; h < nParts; ++h) {
        hRes[h] = new TH1F(Form("n#sigma_{%s}",names[h].Data()), 
        Form("n#sigma_{%s} for %.2f < p < %.2f GeV/c; n#sigma_{%s}; Counts", names[h].Data(),pMin, pMax, names[h].Data()),nBins, xMin, xMax);
        hRes[h]->SetLineColor(colors[h]);
        hRes[h]->SetMarkerColor(colors[h]);}
    
    Long64_t nEntries = chain.GetEntries();
    nEntries = std::min(nEntries, static_cast<Long64_t>(1e6));
    for(Long64_t i = 0; i < nEntries; ++i){
        chain.GetEntry(i);
        for(int j = 0; j < 2; ++j){
            Float_t pG = inner[j]; 
            if(pG < pMin || pG > pMax)   continue;
            
            for (int h = 0; h < nParts; ++h) {
                if (!TMath::IsNaN(tpcNS[h][j])){
                    hRes[h]-> Fill(tpcNS[h][j]);
                }
            }
        }
    }

    TCanvas *c = new TCanvas("c", "Residuals", 800, 600);
    c->SetLeftMargin(0.15);
    c->SetGrid();
    c->Print("nSigmaTPC_plot.pdf[");
    for (int h = 0; h < nParts; ++h) {
        c->Clear();
        hRes[h]->SetMarkerStyle(kFullCircle);
        hRes[h]->SetMarkerColor(kBlack);
        hRes[h]->SetLineColor(kRed);
        hRes[h]->SetMarkerSize(0.75); 
        hRes[h]->Draw("E1");
        c->Update();
        gSystem->ProcessEvents();

        int N = 0;
        std::cout << "\nHistogram "     << h
                << " (" << names[h]   << "): how many Gaussians? ";
        std::cin  >> N;
        if (N <= 0) {                   
            c->Print("nSigmaTPC_plot.pdf");
            continue;
        }

        std::vector<double> means;
        means.reserve(N);
        std::cout << "  Enter " << N << " mean positions (space-separated): ";
        for (int i = 0; i < N; ++i) {
            double m;  std::cin >> m;
            means.push_back(m);
        }
        
            
        std::ostringstream form;
        for (int i = 0; i < N; ++i) {
            if (i) form << "+";
            form << "gaus(" << 3*i << ")";
        }
        TF1* fSum = new TF1(Form("sum_%d",h), form.str().c_str(), xMin, xMax);
        fSum->SetNpx(500);

        for (int i = 0; i < N; ++i) {
            double mu  = means[i];
            int    bin = hRes[h]->FindBin(mu);
            double amp = hRes[h]->GetBinContent(bin);

            fSum->SetParameter(3*i+0, amp);   // amplitude
            fSum->SetParameter(3*i+1, mu);    // mean (user value)
            fSum->SetParameter(3*i+2, 0.3);   // sigma (generic guess)
        }
        hRes[h]->Fit(fSum,"RQ0");

        hRes[h]->Fit(fSum,"RQ");
        double chi2 = fSum->GetChisquare();
        int NDF = fSum->GetNDF();

        hRes[h]->Draw("E1");
        
        TPaveText *pt = new TPaveText(0.0,0.9,0.2,1,"NDC");
        pt->AddText(Form("#chi^{2}/NDF = %.3f",chi2/NDF));
        pt->SetFillColorAlpha(0,0);       // transparent box
        pt->Draw("same");

        fSum->SetLineColor(kRed);
        fSum->SetLineWidth(2);
        fSum->Draw("same");                  

        for (int i = 0; i < N; ++i) {
            TF1* g = new TF1(Form("g_%d_%d",h,i),"gaus",xMin,xMax);
            for (int p = 0; p < 3; ++p)
                g->SetParameter(p, fSum->GetParameter(3*i+p));
            g->SetLineStyle(2);              
            g->SetLineColor(colors[i % nParts]);
            g->Draw("same");
        }

        c->Print("nSigmaTPC_plot.pdf");} 
    
    c->Print("nSigmaTPC_plot.pdf]");
}
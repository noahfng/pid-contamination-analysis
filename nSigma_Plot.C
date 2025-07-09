#include <TChain.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TMath.h>
#include <TF1.h>          
#include "finished_projects/AddTrees.h"

void nSigma_Plot(){
    gROOT->SetBatch(kTRUE);
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
    const Double_t pMin = 0.07, pMax = 0.55; 
    const Double_t xMin = -15.0, xMax = 15.0;
    const Int_t nParts = 5;
    const TString names[nParts] = {"e","#mu","#pi","K","p"};
    const Int_t colors[nParts] = {kBlue, kRed, kGreen+2, kOrange+7, kViolet};

    TH1F *hRes[nParts];
    for (int h = 0; h < nParts; ++h) {
        hRes[h] = new TH1F(Form("n#sigma_{%s}",names[h].Data()), 
        Form("Residuals for %.2f < p < %.2f GeV/c; n#sigma; Counts",pMin, pMax),nBins, xMin, xMax);
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
    c->Print("Residuals_of_Bethe_Bloch.pdf[");
    for (int h = 0; h < nParts; ++h) {
        c->Clear();
        hRes[h]->SetMarkerStyle(kFullCircle);
        hRes[h]->SetMarkerColor(kBlack);
        hRes[h]->SetLineColor(kRed);
        hRes[h]->SetMarkerSize(0.75); 
        hRes[h]->Draw("E1");
        c->Print("Residuals_of_Bethe_Bloch.pdf");} 
    
    c->Print("Residuals_of_Bethe_Bloch.pdf]");
}
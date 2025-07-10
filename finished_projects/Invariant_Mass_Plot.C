#include <TChain.h>          
#include <TCanvas.h>        
#include <TH1.h>                      
#include <TLegend.h>        
#include <TLine.h>         
#include <TMath.h>            
#include <TStyle.h>          
#include <TString.h>           
#include <TLorentzVector.h>    
#include <algorithm>        
#include "AddTrees.h"

void Invariant_Mass_Plot() {
    gROOT->SetBatch(kTRUE); 
    gStyle->SetOptStat(1);
    const char* baseDir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";
    const bool applyTPCnSigmaFilter = true;
    const Float_t nSigmaTPC = 3.0; 
    const bool applyTOFEventfilter = false; 
    const bool applyTOFnSigmaFilter = false; 
    const Float_t nSigmaTOF = 3.0;
  
    TChain chain("twotauchain");
    AddTrees(chain, baseDir);
    Long64_t nEntries = std::min(chain.GetEntries(), static_cast<Long64_t>(1e6));
    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkPx", 1);
    chain.SetBranchStatus("fTrkPy", 1);
    chain.SetBranchStatus("fTrkPz", 1);
    chain.SetBranchStatus("fTrkTOFexpMom", 1);
    const char* subs[5] = {"El","Mu","Pi","Ka","Pr"};
    for (int i = 0; i < 5; ++i) {
        chain.SetBranchStatus(Form("fTrkTPCnSigma%s", subs[i]), 1);
        chain.SetBranchStatus(Form("fTrkTOFnSigma%s", subs[i]), 1);}
    Float_t tpcNS[5][2], tofNS[5][2];
    Float_t px[2], py[2], pz[2];
    Float_t tofExpMom[2];
    chain.SetBranchAddress("fTrkTOFexpMom", tofExpMom);
    chain.SetBranchAddress("fTrkPx",px);
    chain.SetBranchAddress("fTrkPy",py);
    chain.SetBranchAddress("fTrkPz",pz);
    for (int i = 0; i < 5; ++i) {
        chain.SetBranchAddress(Form("fTrkTPCnSigma%s", subs[i]), tpcNS[i]);
        chain.SetBranchAddress(Form("fTrkTOFnSigma%s", subs[i]), tofNS[i]);}

    Double_t masses[5] = {
        0.00051099895,   // e
        0.1056583755,    // μ
        0.13957039,      // π
        0.493677,        // K
        0.93827208816    // p
    };
    const char* names[6] = {
    "e^{+}e^{-}",
    "#mu^{+}#mu^{-}",
    "#pi^{+}#pi^{-}",
    "K^{+}K^{-}",
    "p^{+}p^{-}",
    "K#pi"};
    Int_t colors[6] = {kBlue, kBlue, kBlue, kBlue, kBlue, kBlue};

    const Int_t   nPtBins = 100;
    const Float_t ptMax   = 5.0;
    TH1D* hM[6];
    for (int i = 0; i < 6; ++i) {
        hM[i] = new TH1D(Form("Invariant mass %s", names[i]),
                         Form("Invariant mass %s;M (GeV/#it{c}^{2});Entries", names[i]),
                         nPtBins, 0.0, ptMax);
        hM[i]->SetLineColor(colors[i]);
        hM[i]->SetLineWidth(2);
    }
    
    for (Long64_t i = 0; i < nEntries; ++i) {
        chain.GetEntry(i);
        if(applyTOFEventfilter && (tofExpMom[0] < 0.0 || tofExpMom[1] < 0.0)) continue; 

        for (int j = 0; j < 5; ++j) {
            if (applyTPCnSigmaFilter && (TMath::Abs(tpcNS[j][0]) > nSigmaTPC || TMath::Abs(tpcNS[j][1]) > nSigmaTPC))  continue;
            if (applyTOFnSigmaFilter && (TMath::Abs(tofNS[j][0]) > nSigmaTOF || TMath::Abs(tofNS[j][1]) > nSigmaTOF)) continue;

            Float_t p1 = TMath::Sqrt(px[0]*px[0] + py[0]*py[0] + pz[0]*pz[0]);
            Float_t p2 = TMath::Sqrt(px[1]*px[1] + py[1]*py[1] + pz[1]*pz[1]);
            Float_t e1 = TMath::Sqrt(p1*p1 + masses[j]*masses[j]);
            Float_t e2 = TMath::Sqrt(p2*p2 + masses[j]*masses[j]);

            Float_t Esum = e1 + e2;
            Float_t pxsum = px[0] + px[1];
            Float_t pysum = py[0] + py[1];
            Float_t pzsum = pz[0] + pz[1];
            Float_t M2 = Esum*Esum - (pxsum*pxsum + pysum*pysum + pzsum*pzsum);
            
            if (M2 > 0) {
                hM[j]->Fill(std::sqrt(M2));}
        }
        bool piK = true, Kpi = true;
        if (applyTPCnSigmaFilter) {
            piK &= (TMath::Abs(tpcNS[2][0]) < nSigmaTPC && TMath::Abs(tpcNS[3][1]) < nSigmaTPC);
            Kpi &= (TMath::Abs(tpcNS[3][0]) < nSigmaTPC && TMath::Abs(tpcNS[2][1]) < nSigmaTPC);}
        
        if (applyTOFnSigmaFilter) {
                piK &= (TMath::Abs(tofNS[2][0]) < nSigmaTOF && TMath::Abs(tofNS[3][1]) < nSigmaTOF);
                Kpi &= (TMath::Abs(tofNS[3][0]) < nSigmaTOF && TMath::Abs(tofNS[2][1]) < nSigmaTOF);}
        
        if ((applyTPCnSigmaFilter || applyTOFnSigmaFilter) && (piK == Kpi)) continue;

        int i1 = 2; // π
        int i2 = 3; // K
        if (Kpi) {
            i1 = 3; 
            i2 = 2;
        }
        
        TLorentzVector tl1;
        TLorentzVector tl2;

        double p1 = std::sqrt(px[0]*px[0] + py[0]*py[0] + pz[0]*pz[0]);
        double p2 = std::sqrt(px[1]*px[1] + py[1]*py[1] + pz[1]*pz[1]);

        double m1 = masses[i1];
        double m2 = masses[i2];

        double e1 = std::sqrt(p1*p1 + m1*m1);
        double e2 = std::sqrt(p2*p2 + m2*m2);

        tl1.SetPxPyPzE(px[0], py[0], pz[0], e1);
        tl2.SetPxPyPzE(px[1], py[1], pz[1], e2);

        TLorentzVector tlSum = tl1 + tl2;
        double ivm = tlSum.M();

        if(ivm <= 0) continue;
        hM[5]->Fill(ivm);
    }
    for (int ih = 0; ih < 6; ++ih) {
        Double_t integ = hM[ih]->GetEntries();
        if (integ > 0) hM[ih]->Scale(1.0/integ);
    }

    TCanvas* c = new TCanvas("c","Invariant Mass Pages", 800, 600);
    c->Print("InvariantMass.pdf[");      

    for (int i = 0; i < 6; ++i) {
        c->Clear();
        c->SetLogy();
        hM[i]->Draw("HIST");
        
        Double_t y1 = 0, y2 = hM[i]->GetMaximum()* 1.65;
        const Double_t mRho   = 0.763;  
        const Double_t mKstar = 0.890;  
        const Double_t mJpsi  = 3.0969; 
      
        TLine* l1= new TLine(mRho, y1, mRho, y2);
        l1->SetLineColor(kBlue);
        l1->SetLineStyle(2);
        l1->SetLineWidth(2);
        l1->Draw("SAME");
      
        TLine* l2= new TLine(mKstar, y1, mKstar, y2);
        l2->SetLineColor(kGreen+2);
        l2->SetLineStyle(2);
        l2->SetLineWidth(2);
        l2->Draw("SAME");
      
        TLine* l3 = new TLine(mJpsi, y1, mJpsi, y2);
        l3->SetLineColor(kRed);
        l3->SetLineStyle(2);
        l3->SetLineWidth(2);
        l3->Draw("SAME");

        TLegend* leg = new TLegend(0.8, 0.65, 0.88, 0.75);
        leg->AddEntry(l1, "#rho(770)",  "l");   
        leg->AddEntry(l2, "K*(892)", "l");
        leg->AddEntry(l3, "J/#psi",     "l");
        leg->SetBorderSize(0);
        leg->SetTextSize(0.03);
        leg->Draw("SAME");
        c->Print("InvariantMass.pdf");
    }

    c->Print("InvariantMass.pdf]");    
 
}


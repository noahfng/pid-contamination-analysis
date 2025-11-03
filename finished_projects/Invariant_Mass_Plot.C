#include <algorithm>

#include "TCanvas.h"
#include "TChain.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TROOT.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TStyle.h"

#include <AddTrees.h>   // project-specific: file discovery/chain fill
#include <helpers.h>    // project-specific: masses, charges, colors, Bethe–Bloch

void Invariant_Mass_Plot() {
    auto help = new helper();
    gROOT->SetBatch(kTRUE); 
    gStyle->SetOptStat(1);

    // basic config  
    const Bool_t applyTPCnSigmaFilter = true; // TPC PID gate
    const Float_t nSigmaTPC = 3.0; 
    const Bool_t applyTOFEventfilter = false; // require TOF info
    const Bool_t applyTOFnSigmaFilter = false; // TOF PID gate
    const Bool_t applyNorm = false; // normalize histograms by entries
    const Float_t nSigmaTOF = 3.0;

    const Double_t nEntriesLimit = 1e7;
    const Int_t   nPtBins = 100;
    const Float_t ptMax   = 5.0;

    const Int_t nParts = helper::nParts;
    const Int_t NtrkMax = help->NtrkMax;

    // input data chain
    TChain chain("twotauchain");
    AddTrees(chain, help->base_dir);

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkPx", 1);
    chain.SetBranchStatus("fTrkPy", 1);
    chain.SetBranchStatus("fTrkPz", 1);
    chain.SetBranchStatus("fTrkTOFexpMom", 1);
    for (Int_t i = 0; i < nParts; ++i) {
        chain.SetBranchStatus(Form("fTrkTPCnSigma%s", help->pNames[i]), 1);
        chain.SetBranchStatus(Form("fTrkTOFnSigma%s", help->pNames[i]), 1);
    }

    std::vector<Float_t> px(NtrkMax);
    std::vector<Float_t> py(NtrkMax);
    std::vector<Float_t> pz(NtrkMax);
    std::vector<Float_t> tofExpMom(NtrkMax);
    std::vector<std::vector<Float_t>> tpcNS(nParts, std::vector<Float_t>(NtrkMax));
    std::vector<std::vector<Float_t>> tofNS(nParts, std::vector<Float_t>(NtrkMax));

    chain.SetBranchAddress("fTrkTOFexpMom", tofExpMom.data());
    chain.SetBranchAddress("fTrkPx",px.data());
    chain.SetBranchAddress("fTrkPy",py.data());
    chain.SetBranchAddress("fTrkPz",pz.data());
    for (Int_t i = 0; i < nParts; ++i) {
        chain.SetBranchAddress(Form("fTrkTPCnSigma%s", help->pNames[i]), tpcNS[i].data());
        chain.SetBranchAddress(Form("fTrkTOFnSigma%s", help->pNames[i]), tofNS[i].data());
    }
    
    // channels: same-mass (e,e), (μ,μ), (π,π), (K,K), (p,p) and mixed (K,π)
    const Char_t* names[6] = {"e^{+}e^{-}", "#mu^{+}#mu^{-}", "#pi^{+}#pi^{-}", "K^{+}K^{-}", "p^{+}p^{-}", "K#pi"};
    Int_t colors[6] = {kBlue, kBlue, kBlue, kBlue, kBlue, kBlue};
    Long64_t nEntries = std::min(chain.GetEntries(), static_cast<Long64_t>(nEntriesLimit));

    // output histograms
    TH1D* hM[6];
    for (Int_t i = 0; i < 6; ++i) {
        hM[i] = new TH1D(Form("Invariant mass %s", names[i]),
                         Form("Invariant mass %s;M (GeV/#it{c}^{2});Entries", names[i]),
                         nPtBins, 0.0, ptMax);
        hM[i]->SetLineColor(colors[i]);
        hM[i]->SetLineWidth(2);
    }
    
    for (Long64_t i = 0; i < nEntries; ++i) {
        chain.GetEntry(i);

        // require valid TOF momentum for both tracks
        if(applyTOFEventfilter && (tofExpMom[0] < 0.0 || tofExpMom[1] < 0.0)) continue; 

        for (Int_t j = 0; j < 5; ++j) {
            // PID gates per track
            if (applyTPCnSigmaFilter && (TMath::Abs(tpcNS[j][0]) > nSigmaTPC || TMath::Abs(tpcNS[j][1]) > nSigmaTPC))  continue;
            if (applyTOFnSigmaFilter && (TMath::Abs(tofNS[j][0]) > nSigmaTOF || TMath::Abs(tofNS[j][1]) > nSigmaTOF)) continue;

            // build energies using track momenta (GeV) and species mass (convert MeV->GeV)
            Float_t p1 = TMath::Sqrt(px[0]*px[0] + py[0]*py[0] + pz[0]*pz[0]);
            Float_t p2 = TMath::Sqrt(px[1]*px[1] + py[1]*py[1] + pz[1]*pz[1]);
            Float_t e1 = TMath::Sqrt(p1*p1 + help->pMasses[j]*help->pMasses[j]/1e6);
            Float_t e2 = TMath::Sqrt(p2*p2 + help->pMasses[j]*help->pMasses[j]/1e6);

            // invariant mass M^2 = (E1+E2)^2 - |p1+p2|^2
            Float_t Esum = e1 + e2;
            Float_t pxsum = px[0] + px[1];
            Float_t pysum = py[0] + py[1];
            Float_t pzsum = pz[0] + pz[1];
            Float_t M2 = Esum*Esum - (pxsum*pxsum + pysum*pysum + pzsum*pzsum);
            
            if (M2 > 0) {
                hM[j]->Fill(std::sqrt(M2));}
        }

        // mixed Kπ hypothesis: decide assignment from PID (π-K vs K-π)
        Bool_t piK = true, Kpi = true;
        if (applyTPCnSigmaFilter) {
            piK &= (TMath::Abs(tpcNS[2][0]) < nSigmaTPC && TMath::Abs(tpcNS[3][1]) < nSigmaTPC);
            Kpi &= (TMath::Abs(tpcNS[3][0]) < nSigmaTPC && TMath::Abs(tpcNS[2][1]) < nSigmaTPC);
        }
        
        if (applyTOFnSigmaFilter) {
            piK &= (TMath::Abs(tofNS[2][0]) < nSigmaTOF && TMath::Abs(tofNS[3][1]) < nSigmaTOF);
            Kpi &= (TMath::Abs(tofNS[3][0]) < nSigmaTOF && TMath::Abs(tofNS[2][1]) < nSigmaTOF);
        }
        
        // if both or neither pass, skip (ambiguous assignment)
        if ((applyTPCnSigmaFilter || applyTOFnSigmaFilter) && (piK == Kpi)) continue;

        // pick mass assignment
        Int_t i1 = 2; // π
        Int_t i2 = 3; // K
        if (Kpi) {
            i1 = 3; 
            i2 = 2;
        }

        // four-vectors for mixed hypothesis (masses in GeV)
        TLorentzVector tl1, tl2;
        Double_t p1 = std::sqrt(px[0]*px[0] + py[0]*py[0] + pz[0]*pz[0]);
        Double_t p2 = std::sqrt(px[1]*px[1] + py[1]*py[1] + pz[1]*pz[1]);
        Double_t m1 = help->pMasses[i1]/1e3;
        Double_t m2 = help->pMasses[i2]/1e3;
        Double_t e1 = std::sqrt(p1*p1 + m1*m1);
        Double_t e2 = std::sqrt(p2*p2 + m2*m2);

        tl1.SetPxPyPzE(px[0], py[0], pz[0], e1);
        tl2.SetPxPyPzE(px[1], py[1], pz[1], e2);
        TLorentzVector tlSum = tl1 + tl2;
        Double_t ivm = tlSum.M();
        if(ivm > 0) hM[5]->Fill(ivm);
    }

    // optional per-hist normalization by entries
    if (applyNorm){
        for (Int_t ih = 0; ih < 6; ++ih) {
            Double_t integ = hM[ih]->GetEntries();
            if (integ > 0) hM[ih]->Scale(1.0/integ);
        }
    }
    
    // Plotting (PDF with one page per channel)
    TCanvas* c = new TCanvas("c","Invariant Mass Pages", 800, 600);
    c->Print("InvariantMass.pdf[");      

    // expected particle markers (GeV/c^2)
    const Double_t mRho   = 0.763;  
    const Double_t mKstar = 0.890;  
    const Double_t mJpsi  = 3.0969; 
    
    TLine* l1= new TLine(mRho, 0, mRho, 1);
    l1->SetLineColor(kBlue);
    l1->SetLineStyle(2);
    l1->SetLineWidth(2);
    
    TLine* l2= new TLine(mKstar, 0, mKstar, 1);
    l2->SetLineColor(kGreen+2);
    l2->SetLineStyle(2);
    l2->SetLineWidth(2);
    
    TLine* l3 = new TLine(mJpsi, 0, mJpsi, 1);
    l3->SetLineColor(kRed);
    l3->SetLineStyle(2);
    l3->SetLineWidth(2);

    TLegend* leg = new TLegend(0.8, 0.65, 0.88, 0.75);
    leg->AddEntry(l1, "#rho(770)", "l");   
    leg->AddEntry(l2, "K*(892)", "l");
    leg->AddEntry(l3, "J/#psi", "l");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);

    for (Int_t i = 0; i < 6; ++i) {
        c->Clear();
        c->SetLogy();
        Double_t y2 = hM[i]->GetMaximum()*1.65;  // extend marker lines to current ymax
        l1->SetY2(y2);
        l2->SetY2(y2);
        l3->SetY2(y2);
        hM[i]->Draw("HIST");
        l1->Draw("SAME");
        l2->Draw("SAME");
        l3->Draw("SAME");
        leg->Draw("SAME");
        c->Print("InvariantMass.pdf");
    }

    c->Print("InvariantMass.pdf]");    
 
}


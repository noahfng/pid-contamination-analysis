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
    gStyle->SetOptStat(0);

    // basic config
    const Double_t nEntriesLimit = 1e10;
    const Int_t   nPtBins = 300;
    const Float_t  ptMin   = 0.0;
    const Float_t ptMax   = 10.0;

    const Int_t nParts = helper::nParts;
    const Int_t NtrkMax = help->NtrkMax;
    
    const Bool_t KaExclusion = false;     // TOF-based Kaon veto
    const Bool_t PrExclusion = false;     // TOF-based Proton veto
    const Bool_t requireTOF  = (KaExclusion || PrExclusion);
    const Float_t nSigmaTOF = 3.0;        // TOF nσ cut for veto

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
                         nPtBins, ptMin, ptMax);
        hM[i]->SetMarkerStyle(kFullCircle);
        hM[i]->SetMarkerSize(0.5);
        hM[i]->SetMarkerColor(colors[i]);
        hM[i]->SetLineColor(colors[i]);
    }
    
    for (Long64_t i = 0; i < nEntries; ++i) {
        chain.GetEntry(i);

        if (requireTOF) {
            if (tofExpMom[0] < 0 || tofExpMom[1] < 0)
                continue;  // TOF missing → reject only when needed

            // Kaon veto (index 3 in helpers)
            if (KaExclusion) {
                if (!TMath::IsNaN(tofNS[3][0]) && TMath::Abs(tofNS[3][0]) < nSigmaTOF) continue;
                if (!TMath::IsNaN(tofNS[3][1]) && TMath::Abs(tofNS[3][1]) < nSigmaTOF) continue;
            }

            // Proton veto (index 4 in helpers) 
            if (PrExclusion) {
                if (!TMath::IsNaN(tofNS[4][0]) && TMath::Abs(tofNS[4][0]) < nSigmaTOF) continue;
                if (!TMath::IsNaN(tofNS[4][1]) && TMath::Abs(tofNS[4][1]) < nSigmaTOF) continue;
            }
        }

        for (Int_t j = 0; j < 5; ++j) {

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
        // mixed Kπ case: fill both assignments (track0=π,track1=K AND track0=K,track1=π)
        {
            const Int_t iPi = 2; // π
            const Int_t iK  = 3; // K

            auto fillKpi = [&](int idxPi, int idxK) {
                const double pPi = std::sqrt(px[idxPi]*px[idxPi] + py[idxPi]*py[idxPi] + pz[idxPi]*pz[idxPi]);
                const double pK  = std::sqrt(px[idxK]*px[idxK] + py[idxK]*py[idxK] + pz[idxK]*pz[idxK]);

                // masses are in MeV -> convert consistently via m^2 / 1e6 (GeV^2)
                const double ePi = std::sqrt(pPi*pPi + (help->pMasses[iPi]*help->pMasses[iPi])/1e6);
                const double eK  = std::sqrt(pK*pK + (help->pMasses[iK]*help->pMasses[iK])/1e6);

                TLorentzVector vPi, vK;
                vPi.SetPxPyPzE(px[idxPi], py[idxPi], pz[idxPi], ePi);
                vK.SetPxPyPzE(px[idxK], py[idxK], pz[idxK], eK);

                const double m = (vPi + vK).M();
                if (m > 0.0) hM[5]->Fill(m);
            };

            fillKpi(0, 1);
            fillKpi(1, 0);
        }
    }
    
    // Plotting (PDF with one page per channel)
    TCanvas* c = new TCanvas("c","Invariant Mass Pages", 800, 600);
    c->Print("InvariantMass.pdf[");      

    // expected particle markers (GeV/c^2)
    const Double_t mRho   = 0.775;  
    const Double_t mKstar = 0.892;  
    const Double_t mJpsi  = 3.0969; 
    const Double_t mPhi = 1.0195;
    const Double_t mPsi2S = 3.6861;

    TLine* lRho   = new TLine(mRho,   0, mRho,   1);
    TLine* lKstar = new TLine(mKstar, 0, mKstar, 1);
    TLine* lJpsi  = new TLine(mJpsi,  0, mJpsi,  1);
    TLine* lPhi   = new TLine(mPhi,   0, mPhi,   1);
    TLine* lPsi2S = new TLine(mPsi2S, 0, mPsi2S, 1);

    auto setupLine = [](TLine* l, Color_t col) {
        l->SetLineColor(col);
        l->SetLineStyle(7);
        l->SetLineWidth(2);
    };

    setupLine(lRho,   kBlack);
    setupLine(lKstar, kGreen+2);
    setupLine(lJpsi,  kRed);
    setupLine(lPhi,   kOrange+7);
    setupLine(lPsi2S, kViolet);
    
    TLegend* leg_lep = new TLegend(0.82, 0.79, 0.90, 0.89);
    leg_lep->AddEntry(lJpsi,  "J/#psi",    "l");
    leg_lep->AddEntry(lPsi2S, "#psi(2S)",  "l");
    leg_lep->SetBorderSize(0);
    leg_lep->SetFillStyle(0);
    leg_lep->SetTextSize(0.03);

    TLegend* leg_rho = new TLegend(0.80, 0.79, 0.88, 0.89);
    leg_rho->AddEntry(lRho, "#rho^{0}(770)", "l");
    leg_rho->SetBorderSize(0);
    leg_rho->SetFillStyle(0);
    leg_rho->SetTextSize(0.03);

    TLegend* leg_phi = new TLegend(0.80, 0.79, 0.88, 0.89);
    leg_phi->AddEntry(lPhi, "#varphi(1020)", "l");
    leg_phi->SetBorderSize(0);
    leg_phi->SetFillStyle(0);
    leg_phi->SetTextSize(0.03);

    TLegend* leg_kstar = new TLegend(0.80, 0.79, 0.88, 0.89);
    leg_kstar->AddEntry(lKstar, "K*(892)", "l");
    leg_kstar->SetFillStyle(0);
    leg_kstar->SetBorderSize(0);
    leg_kstar->SetTextSize(0.03);
    

    for (Int_t i = 0; i < 6; ++i) {
        c->Clear();
        c->SetLogy();

        Double_t y2 = hM[i]->GetMaximum()*1.65;
        lRho  ->SetY2(y2);
        lKstar->SetY2(y2);
        lJpsi ->SetY2(y2);
        lPhi  ->SetY2(y2);
        lPsi2S->SetY2(y2);

        hM[i]->Draw("E1");

        if (i == 0 || i == 1) {           // e+e- or mu+mu-
            lJpsi ->Draw("SAME");
            lPsi2S->Draw("SAME");
            leg_lep->Draw("SAME");
        } else if (i == 2) {              // pi+pi-
            lRho->Draw("SAME");
            leg_rho->Draw("SAME");
        } else if (i == 3) {              // K+K-
            lPhi->Draw("SAME");
            leg_phi->Draw("SAME");
        } else if (i == 4) {              // p+p- : no clear VM marker in this range
            // nothing extra
        } else if (i == 5) {              // K pi : K*
            lKstar->Draw("SAME");
            leg_kstar->Draw("SAME");
        }

        c->Print("InvariantMass.pdf");
    }

    c->Print("InvariantMass.pdf]");    
 
}


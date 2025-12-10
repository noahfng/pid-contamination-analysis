#include "TROOT.h"       
#include "TStyle.h"     
#include "TChain.h"
#include "TCanvas.h"
#include "TH2.h"          
#include "TMath.h"       
#include "TString.h"     
#include "TLorentzVector.h" 

#include <AddTrees.h>   // project-specific: file discovery/chain fill
#include <helpers.h>    // project-specific: masses, charges, colors, Bethe–Bloch

void Invariant_Mass_vs_pt() {
    auto help = new helper();
    gROOT->SetBatch(kTRUE);
    gStyle->SetPalette(kRainBow);

    // basic config
    const Double_t nEntriesLimit = 1e10;  
    const Bool_t applyTPCnSigmaFilter = true; // TPC PID gate
    const Float_t nSigmaTPC = 3.0;
    const Bool_t applyTOFEventfilter = true; // require valid TOF momentum for both tracks
    const Bool_t applyTOFnSigmaFilter = true; // TOF PID gate
    const Float_t nSigmaTOF = 3.0;

    // if both true, the 2D plot will show plotSysPt 
    const Bool_t plotSysPt = true; // pair pT: |pT(0)+pT(1)| → 1 entry/event
    const Bool_t plotTrackPt = false; // track pT: |pT(t)| for t=0,1 → 2 entries/event

    const Int_t nMassBins = 1000; const Float_t massMax = 10.0; // GeV/c^2
    const Int_t   nPtBins = 1000; const Float_t ptMax = 10.0; // GeV/c

    const Int_t nParts = helper::nParts;
    const Int_t NtrkMax = help->NtrkMax;

    TString yTitle = plotSysPt ? "p_{T,sys}" : "p_{T,track}";
    TString outName = plotSysPt ? "M_vs_PT_sys.pdf" : "M_vs_PT_track.pdf";

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
    Int_t colors[6] = { kBlue, kAzure+1, kCyan+1, kGreen+2, kMagenta+2, kOrange+1 };
    
    Long64_t nEntries = TMath::Min(chain.GetEntries(), static_cast<Long64_t>(nEntriesLimit));

    // output histograms
    TH2D* h2Mpt[6];
    for (Int_t i = 0; i < 6; ++i) {
        h2Mpt[i] = new TH2D(
            Form("M vs p_{T}%s", names[i]),
            Form("Mass vs %s %s; M (GeV/#it{c}^{2}); %s (GeV/#it{c})", yTitle.Data(), names[i], yTitle.Data()),
            nMassBins, 0.0, massMax,
            nPtBins,   0.0, ptMax
        );
        h2Mpt[i]->SetLineColor(colors[i]);
    }

    for (Long64_t ev = 0; ev < nEntries; ++ev) {
        chain.GetEntry(ev);

        // require valid TOF momentum for both tracks
        if (applyTOFEventfilter && (tofExpMom[0] < 0 || tofExpMom[1] < 0)) continue;

        for (Int_t j = 0; j < 5; ++j) {
            // PID gates per track
            if (applyTPCnSigmaFilter && (TMath::Abs(tpcNS[j][0]) > nSigmaTPC || TMath::Abs(tpcNS[j][1]) > nSigmaTPC)) continue;
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
            if (M2 <= 0) continue;

            Double_t mass = TMath::Sqrt(M2);
            if (plotSysPt) {
                Double_t ptSys = TMath::Hypot(pxsum, pysum);
                h2Mpt[j]->Fill(mass, ptSys);
            } else {
                for (Int_t t = 0; t < 2; ++t) {
                    Double_t ptTrk = TMath::Hypot(px[t], py[t]);
                    h2Mpt[j]->Fill(mass, ptTrk);
                }
            }
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

        // choose masses (MeV) -> energies (GeV)
        Int_t i1 = piK ? 2 : 3; // 2: π, 3: K
        Int_t i2 = piK ? 3 : 2;

        TLorentzVector l1, l2;
        Double_t p1  = std::hypot(px[0], py[0], pz[0]);
        Double_t p2  = std::hypot(px[1], py[1], pz[1]);
        Double_t e1  = TMath::Sqrt(p1*p1 + help->pMasses[i1]*help->pMasses[i1]/1e6);
        Double_t e2  = TMath::Sqrt(p2*p2 + help->pMasses[i2]*help->pMasses[i2]/1e6);
        l1.SetPxPyPzE(px[0], py[0], pz[0], e1);
        l2.SetPxPyPzE(px[1], py[1], pz[1], e2);
        Double_t ivm = (l1 + l2).M();
        if (ivm <= 0) continue;

        // fill histogram
        if (plotSysPt) {
            Double_t ptSys = TMath::Hypot(px[0]+px[1], py[0]+py[1]);
            h2Mpt[5]->Fill(ivm, ptSys);
        } 
        else {
            for (Int_t t = 0; t < 2; ++t) {
                Double_t ptTrk = TMath::Hypot(px[t], py[t]);
                h2Mpt[5]->Fill(ivm, ptTrk);
            }
        }
    }
    // Plotting 
    TCanvas* c = new TCanvas("c", "");
    c->Print(outName + "["); 
    for (Int_t i = 0; i < 6; ++i) {
        c->Clear();
        c->SetLogz();
        h2Mpt[i]->Draw("COLZ");
        c->Print(outName);
    }
    c->Print(outName + "]"); 
}
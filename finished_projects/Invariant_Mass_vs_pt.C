#include <TROOT.h>       
#include <TStyle.h>     
#include <TChain.h>
#include <TCanvas.h>
#include <TH2.h>          
#include <TMath.h>       
#include <TString.h>     
#include <TLorentzVector.h> 

#include "AddTrees.h"

void Invariant_Mass_vs_pt() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetPalette(kRainBow);
    const Char_t* baseDir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";

    const Bool_t    applyTPCnSigmaFilter = true;
    const Float_t nSigmaTPC            = 3.0;
    const Bool_t    applyTOFEventfilter  = false;
    const Bool_t    applyTOFnSigmaFilter = false;
    const Float_t nSigmaTOF            = 3.0;
    const Bool_t    plotSysPt            = true;
    const Bool_t    plotTrackPt          = false; // if both true, the 2D plot will show plotSysPt 
    const Double_t nEntriesLimit = 1e6;

    TString yTitle = plotSysPt
                     ? "p_{T,sys}"
                     : "p_{T,track}";
    TString outName = plotSysPt
                      ? "M_vs_PT_sys.pdf"
                      : "M_vs_PT_track.pdf";

    TChain chain("twotauchain");
    AddTrees(chain, baseDir);
    Long64_t nEntries = TMath::Min(chain.GetEntries(), static_cast<Long64_t>(nEntriesLimit));
    
    chain.SetBranchStatus("*",0);
    chain.SetBranchStatus("fTrkPx",1);
    chain.SetBranchStatus("fTrkPy",1);
    chain.SetBranchStatus("fTrkPz",1);
    chain.SetBranchStatus("fTrkTOFexpMom",1);
    const Char_t* subs[5]={"El","Mu","Pi","Ka","Pr"};
    for(Int_t i=0;i<5;++i){
        chain.SetBranchStatus(Form("fTrkTPCnSigma%s",subs[i]),1);
        chain.SetBranchStatus(Form("fTrkTOFnSigma%s",subs[i]),1);
    }

    Float_t px[2],py[2],pz[2],tofExp[2];
    Float_t tpcNS[5][2],tofNS[5][2];
    chain.SetBranchAddress("fTrkPx",px);
    chain.SetBranchAddress("fTrkPy",py);
    chain.SetBranchAddress("fTrkPz",pz);
    chain.SetBranchAddress("fTrkTOFexpMom",tofExp);
    for(Int_t i=0;i<5;++i){
        chain.SetBranchAddress(Form("fTrkTPCnSigma%s",subs[i]),tpcNS[i]);
        chain.SetBranchAddress(Form("fTrkTOFnSigma%s",subs[i]),tofNS[i]);
    }

    Double_t masses[5] = {
        0.00051099895,   // e
        0.1056583755,    // μ
        0.13957039,      // π
        0.493677,        // K
        0.93827208816    // p
    };
    const Char_t* names[6] = {
    "e^{+}e^{-}",
    "#mu^{+}#mu^{-}",
    "#pi^{+}#pi^{-}",
    "K^{+}K^{-}",
    "p^{+}p^{-}",
    "K#pi"};
    Int_t colors[6] = { kBlue, kAzure+1, kCyan+1, kGreen+2, kMagenta+2, kOrange+1 };
    const Int_t   nMassBins = 100;
    const Float_t massMax   =   5.0;
    const Int_t   nPtBins   = 100;
    const Float_t ptMax     =   5.0;
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
        if (applyTOFEventfilter && (tofExp[0] < 0 || tofExp[1] < 0))
            continue;

        for (Int_t j = 0; j < 5; ++j) {
            if (applyTPCnSigmaFilter &&
               (TMath::Abs(tpcNS[j][0]) > nSigmaTPC ||
                TMath::Abs(tpcNS[j][1]) > nSigmaTPC)) continue;
            if (applyTOFnSigmaFilter &&
               (TMath::Abs(tofNS[j][0]) > nSigmaTOF ||
                TMath::Abs(tofNS[j][1]) > nSigmaTOF)) continue;

            Float_t p1 = TMath::Sqrt(px[0]*px[0] + py[0]*py[0] + pz[0]*pz[0]);
            Float_t p2 = TMath::Sqrt(px[1]*px[1] + py[1]*py[1] + pz[1]*pz[1]);
            Float_t e1 = TMath::Sqrt(p1*p1 + masses[j]*masses[j]);
            Float_t e2 = TMath::Sqrt(p2*p2 + masses[j]*masses[j]);

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

        Bool_t piK = true, Kpi = true;
        if (applyTPCnSigmaFilter) {
            piK &= (TMath::Abs(tpcNS[2][0]) < nSigmaTPC &&
                    TMath::Abs(tpcNS[3][1]) < nSigmaTPC);
            Kpi &= (TMath::Abs(tpcNS[3][0]) < nSigmaTPC &&
                    TMath::Abs(tpcNS[2][1]) < nSigmaTPC);
        }
        if (applyTOFnSigmaFilter) {
            piK &= (TMath::Abs(tofNS[2][0]) < nSigmaTOF &&
                    TMath::Abs(tofNS[3][1]) < nSigmaTOF);
            Kpi &= (TMath::Abs(tofNS[3][0]) < nSigmaTOF &&
                    TMath::Abs(tofNS[2][1]) < nSigmaTOF);
        }
        if ((applyTPCnSigmaFilter || applyTOFnSigmaFilter) && (piK == Kpi))
            continue;

        Int_t i1 = piK ? 2 : 3;
        Int_t i2 = piK ? 3 : 2;

        TLorentzVector l1, l2;
        Double_t p1  = std::hypot(px[0], py[0], pz[0]);
        Double_t p2  = std::hypot(px[1], py[1], pz[1]);
        Double_t e1  = TMath::Sqrt(p1*p1 + masses[i1]*masses[i1]);
        Double_t e2  = TMath::Sqrt(p2*p2 + masses[i2]*masses[i2]);
        l1.SetPxPyPzE(px[0], py[0], pz[0], e1);
        l2.SetPxPyPzE(px[1], py[1], pz[1], e2);
        Double_t ivm = (l1 + l2).M();
        if (ivm <= 0) continue;

        if (plotSysPt) {
            Double_t ptSys = TMath::Hypot(px[0]+px[1], py[0]+py[1]);
            h2Mpt[5]->Fill(ivm, ptSys);
        } else {
            for (Int_t t = 0; t < 2; ++t) {
                Double_t ptTrk = TMath::Hypot(px[t], py[t]);
                h2Mpt[5]->Fill(ivm, ptTrk);
            }
        }
    }

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
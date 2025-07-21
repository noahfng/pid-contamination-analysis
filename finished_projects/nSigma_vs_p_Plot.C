#include <algorithm>
#include <vector>

#include "TROOT.h"       
#include "TStyle.h"       
#include "TMath.h"        
#include "TChain.h"      
#include "TCanvas.h"     
#include "TH2F.h"         
#include "TLegend.h"      
#include "TGraph.h"       
#include "TString.h" 
#include "TFile.h"    

#include <AddTrees.h>
#include <get_expected_signal.h>
#include <getReso.h>

void nSigma_vs_p_Plot() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);

    const Char_t* baseDir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";

    TChain chain("twotauchain");
    AddTrees(chain, baseDir);
    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);
    chain.SetBranchStatus("fTrkTOFexpMom", 1);
    const Char_t* subs[5] = {"El","Mu","Pi","Ka","Pr"};
    for (Int_t i = 0; i < 5; ++i) {
        chain.SetBranchStatus(Form("fTrkTPCnSigma%s", subs[i]), 1);
        chain.SetBranchStatus(Form("fTrkTOFnSigma%s", subs[i]), 1);
    }

    const Int_t NtrkMax = 2;
    Float_t inner[NtrkMax] = {0};
    Float_t tpcNS[5][NtrkMax] = {{0}};
    Float_t tofNS[5][NtrkMax] = {{0}};
    Float_t tpcSignal[NtrkMax] = {0};
    Float_t tofExpMom[NtrkMax] = {0};
    Bool_t plotTPC = true;
    Bool_t plotTOF = false;
    const Double_t nEntriesLimit = 1e6;
    
    chain.SetBranchAddress("fTrkTPCsignal", tpcSignal);
    chain.SetBranchAddress("fTrkTPCinnerParam", inner);
    chain.SetBranchAddress("fTrkTOFexpMom", tofExpMom);
    for (Int_t i = 0; i < 5; ++i) {
        chain.SetBranchAddress(Form("fTrkTPCnSigma%s", subs[i]), tpcNS[i]);
        chain.SetBranchAddress(Form("fTrkTOFnSigma%s", subs[i]), tofNS[i]);
    }
    const Char_t* names[5]   = {"e", "#mu", "#pi", "K", "p"};
    const Double_t masses[5] = {0.00051099895, 0.1056583755,  0.13957039, 0.493677, 0.93827208816};
    const Double_t resoTOF[5]   = {0.013, 0.013, 0.013, 0.019, 0.020};
    const Double_t resoTPC[5]   = {0.085, 0.072, 0.074, 0.09, 0.08}; 
    const Int_t colors[5]  = {kBlue, kGreen+2, kOrange+7, kMagenta+2, kCyan+1};
    const Int_t npoints = 200;
    const Double_t pMin = 0.01, pMax = 5.0;
    Double_t pgrid[npoints];
    for (Int_t i = 0; i < npoints; ++i) {
        pgrid[i] = pMin + i*(pMax - pMin)/(npoints - 1);
    }

    Long64_t nEntries = std::min(chain.GetEntries(), static_cast<Long64_t>(nEntriesLimit));
    TH2F* histTPC[5];
    TH2F* histTOF[5];
    for (Int_t i = 0; i < 5; ++i) {
        histTPC[i] = new TH2F(
          Form("tpc_%s", subs[i]),
          Form("n#sigma_{%s} vs p (TPC);p [GeV/c];n#sigma_{%s}", names[i], names[i]),
          100, pMin, pMax, 100, -50, 50
        );
        histTOF[i] = new TH2F(
          Form("tof_%s", subs[i]),
          Form("n#sigma_{%s} vs p (TOF);p [GeV/c];n#sigma_{%s}", names[i], names[i]),
          100, pMin, pMax, 100, -10, 10
        );
    }

    for (Long64_t ev = 0; ev < nTotal; ++ev) {
        chain.GetEntry(ev);
        for (Int_t tr = 0; tr < NtrkMax; ++tr) {
            Float_t p = inner[tr];
            if (p <= 0) continue;
            for (Int_t sp = 0; sp < 5; ++sp) {
                if(plotTPC) histTPC[sp]->Fill(p, tpcNS[sp][tr]);
                if (plotTOF && tofExpMom[tr] > 0)
                    histTOF[sp]->Fill(p, tofNS[sp][tr]);
            }
        }
    }

    TGraph* tpcCurves[5][5];
    TGraph* tofCurves[5][5];
    for (Int_t ref = 0; ref < 5; ++ref) {
        Double_t mRef = masses[ref];

        for (Int_t hyp = 0; hyp < 5; ++hyp) {
            Double_t mHyp = masses[hyp];

            std::vector<Double_t> xv_tpc, yv_tpc;
            std::vector<Double_t> xv_tof, yv_tof;
            xv_tpc.reserve(npoints);  yv_tpc.reserve(npoints);
            xv_tof.reserve(npoints);  yv_tof.reserve(npoints);

            for (Int_t ip = 0; ip < npoints; ++ip) {
                Double_t pg = pgrid[ip];
                if (pg < 0.1) continue;
                if (plotTPC){
                    Double_t dRef = get_expected_signal(pg*1000, mRef*1000, 1.0);
                    Double_t dHyp = get_expected_signal(pg*1000, mHyp*1000, 1.0);
                    Double_t resoRefAbs = getReso(kTPC, (Char_t*)subs[ref], pg);
                    Double_t fracRef = resoRefAbs / dRef;

                    xv_tpc.push_back(pg);
                    yv_tpc.push_back((dHyp/dRef - 1.0) / resoTPC[ref]);
                }
                if (plotTOF){
                    Double_t bRef = pg / TMath::Sqrt(mRef*mRef + pg*pg);
                    Double_t bHyp = pg / TMath::Sqrt(mHyp*mHyp + pg*pg);
                    Double_t resoHyp = getReso(kTOF, (Char_t*)subs[hyp], pg);
                    xv_tof.push_back(pg);
                    yv_tof.push_back((bRef - bHyp) / (bHyp*bHyp*resoHyp));
                }
            }
            
            if (plotTPC) {
                tpcCurves[ref][hyp] = new TGraph(xv_tpc.size(), xv_tpc.data(), yv_tpc.data());
                tpcCurves[ref][hyp]->SetLineColor(colors[hyp]);
                tpcCurves[ref][hyp]->SetLineWidth(2);
            }
            if (plotTOF) {
                tofCurves[ref][hyp] = new TGraph(xv_tof.size(), xv_tof.data(), yv_tof.data());
                tofCurves[ref][hyp]->SetLineColor(colors[hyp]);
                tofCurves[ref][hyp]->SetLineWidth(2);
            }
        }
    }
    
    TCanvas* c = new TCanvas("c", "", 800, 600);
    TLegend* leg = new TLegend(0, 0.10, 0.15, 0.30);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    const Char_t* pdfTPC = "nsigma_vs_p_tpc.pdf";
    const Char_t* pdfTOF = "nsigma_vs_p_tof.pdf";
    if (plotTPC) c->Print(Form("%s[", pdfTPC));
    if (plotTOF) c->Print(Form("%s[", pdfTOF));
    c->SetLogz();

    for (Int_t i = 0; i < 5; ++i) {
        if (plotTPC) {
            leg->Clear();
            for (Int_t hyp = 0; hyp < 5; ++hyp) leg->AddEntry(tpcCurves[i][hyp], subs[hyp], "l");
            histTPC[i]->Draw("COLZ");
            for (Int_t hyp = 0; hyp < 5; ++hyp) tpcCurves[i][hyp]->Draw("L SAME");
            leg->Draw();
            c->Print(pdfTPC);
        }
        if (plotTOF) {
            leg->Clear();
            for (Int_t hyp = 0; hyp < 5; ++hyp) leg->AddEntry(tofCurves[i][hyp], subs[hyp], "l");
            histTOF[i]->Draw("COLZ");
            for (Int_t hyp = 0; hyp < 5; ++hyp) tofCurves[i][hyp]->Draw("L SAME");
            leg->Draw();
            c->Print(pdfTOF);
        }
    }

    if (plotTPC) c->Print(Form("%s]", pdfTPC));
    if (plotTOF) c->Print(Form("%s]", pdfTOF));
}


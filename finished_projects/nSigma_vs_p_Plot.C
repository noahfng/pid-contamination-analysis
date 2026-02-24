#include <algorithm>
#include <vector>

#include "TROOT.h"       
#include "TStyle.h"       
#include "TMath.h"        
#include "TChain.h"      
#include "TCanvas.h"     
#include "TH2D.h"         
#include "TLegend.h"      
#include "TGraph.h"       
#include "TString.h" 
#include "TFile.h"    

#include <AddTrees.h>   // project-specific: file discovery/chain fill
#include <helpers.h>    // project-specific: masses, charges, colors, Bethe–Bloch

void nSigma_vs_p_Plot() {
    auto help = new helper();
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);
    gStyle->SetNumberContours(256);

    // basic config
    const Bool_t plotTPC = true; // draw TPC nσ vs p
    const Bool_t plotTOF = false; // draw TOF nσ vs p, for cleaner plots set pMin to 0.3 GeV/c
    const Bool_t KaExclusion = false; // TOF-based Kaon veto
    const Bool_t PrExclusion = false; // TOF-based Proton veto
    const Bool_t requireTOF = (KaExclusion || PrExclusion); // automatically require TOF info if any TOF-based veto is used
    const Double_t nEntriesLimit = 1e10;
    const Int_t npoints = 500;
    const Double_t yMin   = -20.0, yMax = 30.0;
    const Double_t pMin = 0.1, pMax = 10.0;

    const Int_t nParts = helper::nParts;
    auto NtrkMax = help->NtrkMax;
    // choose which PID reference species to plot data for (e, μ, π, K, p)
    const std::array<Bool_t, nParts> doPid = {{true, true, true, true, true}}; 

    // input data chain
    TChain chain("twotauchain");
    AddTrees(chain, help->base_dir);

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);
    chain.SetBranchStatus("fTrkTOFexpMom", 1);
    for (Int_t i = 0; i < nParts; ++i) {
        chain.SetBranchStatus(Form("fTrkTPCnSigma%s", help->pNames[i]), 1);
        chain.SetBranchStatus(Form("fTrkTOFnSigma%s", help->pNames[i]), 1);
    }

    std::vector<Float_t> inner(NtrkMax);
    std::vector<Float_t> tofExpMom(NtrkMax);
    std::vector<std::vector<Float_t>> tpcNS(nParts, std::vector<Float_t>(NtrkMax));
    std::vector<std::vector<Float_t>> tofNS(nParts, std::vector<Float_t>(NtrkMax));

    chain.SetBranchAddress("fTrkTPCinnerParam", inner.data());
    chain.SetBranchAddress("fTrkTOFexpMom", tofExpMom.data());

    for (Int_t i = 0; i < nParts; ++i) {
        chain.SetBranchAddress(Form("fTrkTPCnSigma%s", help->pNames[i]), tpcNS[i].data());
        chain.SetBranchAddress(Form("fTrkTOFnSigma%s", help->pNames[i]), tofNS[i].data());
    }
    
    // momentum grid for model curves (GeV/c)
    Double_t pgrid[npoints];
    for (Int_t i = 0; i < npoints; ++i) {
        pgrid[i] = pMin + i*(pMax - pMin)/(npoints - 1);
    }

    // histogram setup
    Long64_t nEntries = std::min(chain.GetEntries(), static_cast<Long64_t>(nEntriesLimit));
    TH2D* histTPC[nParts];
    TH2D* histTOF[nParts];
    for (Int_t i = 0; i < nParts; ++i) {
        if(!doPid[i]) continue;
        histTPC[i] = new TH2D(Form("tpc_%s", help->pNames[i]), Form("n#sigma_{%s} vs #it{p} (TPC);#it{p} (GeV/#it{c});n#sigma_{%s}", help->pCodes[i], help->pCodes[i]), 1000, pMin, pMax, 1000, yMin, yMax);
        histTPC[i]->GetXaxis()->SetTitleOffset(1.3);
        histTOF[i] = new TH2D(Form("tof_%s", help->pNames[i]), Form("n#sigma_{%s} vs #it{p} (TOF);#it{p} (GeV/#it{c});n#sigma_{%s}", help->pCodes[i], help->pCodes[i]), 1000, pMin, pMax, 1000, yMin, yMax);
        histTOF[i]->GetXaxis()->SetTitleOffset(1.3);
    }

    // fill histograms
    for (Long64_t ev = 0; ev < nEntries; ++ev) {
        chain.GetEntry(ev);
        for (Int_t tr = 0; tr < NtrkMax; ++tr) {
            Float_t p = inner[tr]; // GeV/c
            if (p <= 0) continue;
            if (KaExclusion && !TMath::IsNaN(tofNS[3][tr]) && TMath::Abs(tofNS[3][tr]) < 3.0) 
                continue;
            if (PrExclusion && !TMath::IsNaN(tofNS[4][tr]) && TMath::Abs(tofNS[4][tr]) < 3.0) 
                continue;
            if (requireTOF && tofExpMom[tr]< 0)
                continue;
            for (Int_t sp = 0; sp < nParts; ++sp) if(doPid[sp]) {
                if(plotTPC) histTPC[sp]->Fill(p, tpcNS[sp][tr]);
                if (plotTOF)
                    histTOF[sp]->Fill(p, tofNS[sp][tr]);
            }
        }
    }
    
    // Build model nσ curves (per ref species vs all hypotheses
    TGraph* tpcCurves[nParts][nParts];
    TGraph* tofCurves[nParts][nParts];
    for (Int_t ref = 0; ref < nParts; ++ref) {
        if (!doPid[ref]) continue;
        for (Int_t hyp = 0; hyp < nParts; ++hyp) {
            std::vector<Double_t> xv_tpc, yv_tpc;
            std::vector<Double_t> xv_tof, yv_tof;

            for (Int_t ip = 0; ip < npoints; ++ip) {
                Double_t pg = pgrid[ip];
                if (pg < 0.1) continue;

                // helper expects momentum in MeV/c → multiply by 1000
                if (plotTPC) {
                    Double_t ns = help->getnSigma(pg * 1000., help->dNames[0], ref, hyp);
                    xv_tpc.push_back(pg);
                    yv_tpc.push_back(ns);
                }
                if (plotTOF) {
                    Double_t ns = help->getnSigma(pg * 1000., help->dNames[1], ref, hyp);
                    xv_tof.push_back(pg);
                    yv_tof.push_back(ns);
                }
            }

            if (plotTPC) {
                tpcCurves[ref][hyp] = new TGraph(xv_tpc.size(), xv_tpc.data(), yv_tpc.data());
                tpcCurves[ref][hyp]->SetLineColor(help->colors[hyp]);
                tpcCurves[ref][hyp]->SetLineWidth(2);
            }
            if (plotTOF) {
                tofCurves[ref][hyp] = new TGraph(xv_tof.size(), xv_tof.data(), yv_tof.data());
                tofCurves[ref][hyp]->SetLineColor(help->colors[hyp]);
                tofCurves[ref][hyp]->SetLineWidth(2);
            }
        }
    }

    // Draw (one PDF per detector)
    TCanvas* c = new TCanvas("c", "", 800, 600);
    TLegend *leg = new TLegend(0.84, 0.70, 0.88, 0.90);
    leg->SetBorderSize(1); 
    leg->SetMargin(0.45);
    leg->SetFillColorAlpha(kWhite, 0.8);
    const Char_t* pdfTPC = "nsigma_vs_p_tpc.pdf";
    const Char_t* pdfTOF = "nsigma_vs_p_tof.pdf";
    if (plotTPC) c->Print(Form("%s[", pdfTPC));
    if (plotTOF) c->Print(Form("%s[", pdfTOF));
    c->SetRightMargin(0.12);
    c->SetLogz();
    c->SetLogx();

    for (Int_t i = 0; i < 5; ++i) {
        if (!doPid[i]) continue;
          
        if (plotTPC) {
            leg->Clear();
            for (Int_t hyp = 0; hyp < 5; ++hyp) leg->AddEntry(tpcCurves[i][hyp], help->pCodes[hyp], "l");
            histTPC[i]->Draw("COLZ");
            for (Int_t hyp = 0; hyp < 5; ++hyp) tpcCurves[i][hyp]->Draw("L SAME");
            leg->Draw();
            c->Print(pdfTPC);
        }
        if (plotTOF) {
            leg->Clear();
            for (Int_t hyp = 0; hyp < 5; ++hyp) leg->AddEntry(tofCurves[i][hyp], help->pCodes[hyp], "l");
            histTOF[i]->Draw("COLZ");
            for (Int_t hyp = 0; hyp < 5; ++hyp) tofCurves[i][hyp]->Draw("L SAME");
            leg->Draw();
            c->Print(pdfTOF);
        }
    }

    if (plotTPC) c->Print(Form("%s]", pdfTPC));
    if (plotTOF) c->Print(Form("%s]", pdfTOF));
}


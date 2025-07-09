#if defined(__CLING__)
  #pragma cling load("libHist")
  #pragma cling load("libCore")
  #pragma cling load("libRIO")
  #pragma cling load("libTree")
  #pragma cling load("libTreePlayer")
  #pragma cling load("libPhysics")
  #pragma link C++ class std::vector<float>+;
#endif

#include <TChain.h>
#include <TFile.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TIterator.h>
#include <TString.h>
#include <TMath.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <map>
#include "AddTrees.h"

Double_t bethe_bloch_aleph(Double_t bg, Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5) {
    Double_t beta = bg / TMath::Sqrt(1.0 + bg*bg);
    Double_t aa   = TMath::Power(beta, p4);
    Double_t bb   = TMath::Log(p3 + TMath::Power(1.0/bg, p5));
    return (p2 - aa - bb) * p1 / aa;
}

Double_t get_expected_signal(Double_t p, Double_t mass, Double_t charge) {
    const Double_t mMIP = 50.0;
    const Double_t params[5] = {0.19310481, 4.26696118, 0.00522579, 2.38124907, 0.98055396};
    const Double_t chFact = 2.3;

    Double_t bg = p / mass;
    if (bg < 0.05) return -999.;
    Double_t bethe = mMIP
                   * bethe_bloch_aleph(bg,
                                       params[0], params[1], params[2],
                                       params[3], params[4])
                   * TMath::Power(charge, chFact);
    return bethe >= 0 ? bethe : -999.;
}

void Sigma_plot() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);

    const char* baseDir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";

    TChain chain("twotauchain");
    AddTrees(chain, baseDir);
    Long64_t nTotal = chain.GetEntries();
    nTotal = std::min(nTotal, static_cast<Long64_t>(1e6));

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);
    const char* subs[5] = {"El","Mu","Pi","Ka","Pr"};
    for (int i = 0; i < 5; ++i) {
        chain.SetBranchStatus(Form("fTrkTPCnSigma%s", subs[i]), 1);
        chain.SetBranchStatus(Form("fTrkTOFnSigma%s", subs[i]), 1);
    }

    const int NtrkMax = 2;
    Float_t inner[NtrkMax] = {0};
    Float_t tpcNS[5][NtrkMax] = {{0}};
    Float_t tofNS[5][NtrkMax] = {{0}};
    Float_t tpcSignal[NtrkMax] = {0};
    chain.SetBranchAddress("fTrkTPCsignal", tpcSignal);
    chain.SetBranchAddress("fTrkTPCinnerParam", inner);
    for (int i = 0; i < 5; ++i) {
        chain.SetBranchAddress(Form("fTrkTPCnSigma%s", subs[i]), tpcNS[i]);
        chain.SetBranchAddress(Form("fTrkTOFnSigma%s", subs[i]), tofNS[i]);
    }

    Double_t masses[5] = {0.00051099895, 0.1056583755,  0.13957039, 0.493677, 0.93827208816};
    Double_t resoTOF[5]   = {0.013, 0.013, 0.013, 0.019, 0.020};
    Double_t resoTPC[5]   = {0.5, 0.075, 0.078, 0.09, 0.095}; 
    Int_t   colors[5]  = {kRed,kBlue,kMagenta,kOrange+1,kGreen+2};
    const int npoints = 200;
    const Double_t pMin = 0.01, pMax = 5.0;
    Double_t pgrid[npoints];
    for (int i = 0; i < npoints; ++i) {
        pgrid[i] = pMin + i*(pMax - pMin)/(npoints - 1);
    }

    std::map<TString, TH2F*> histTPC, histTOF;
    for (int i = 0; i < 5; ++i) {
        TString s = subs[i];
        histTPC[s] = new TH2F(Form("tpc_%s", s.Data()),
                              Form("TPC #sigma vs p (%s);p [GeV/c];n#sigma_{TPC}",s.Data()),
                              100, pMin, pMax, 100, -10, 10);
        histTOF[s] = new TH2F(Form("tof_%s", s.Data()),
                              Form("TOF #sigma vs p (%s);p [GeV/c];n#sigma_{TOF}",s.Data()),
                              100, pMin, pMax, 100, -10, 10);
    }

    for (Long64_t i = 0; i < nTotal; ++i) {
        chain.GetEntry(i);
        for (int j = 0; j < NtrkMax; ++j) {
            Float_t p = inner[j]; if (p<=0) continue;
            for (int pi = 0; pi < 5; ++pi) {
                histTPC[subs[pi]]->Fill(p, tpcNS[pi][j]);
                histTOF[subs[pi]]->Fill(p, tofNS[pi][j]);
            }
        }
    }

    TCanvas c("c","",800,600);
    const char* pdf = "nsigma_vs_p_tracks_with_curves.pdf";
    c.Print(Form("%s[",pdf));
    c.SetLogz();

    // TPC
    for (int i = 0; i < 5; ++i) {
        TString species = subs[i];
        TH2F* hTPC = histTPC[species];
        hTPC->SetTitle(Form("TPC - %s", species.Data()));
        hTPC->Draw("COLZ");
        int refIdx = -1;
        for (int i = 0; i < 5; ++i) {
            if (species == subs[i]) {
                refIdx = i; 
                break;
            }
        }
        if (refIdx < 0) continue;
        TLegend leg(0, 0.10, 0.15, 0.30, nullptr, "NDC");
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);

        for(int ih = 0; ih < 5; ++ih) {
            std::vector<double> xv, yv;
            xv.reserve(npoints);
            yv.reserve(npoints);

            for(int j = 0; j < npoints; ++j) {
                double pG = pgrid[j];
                double pMeV = pG * 1000.0;
                double dRef = get_expected_signal(pMeV, masses[refIdx]*1000.0, 1.0);
                double dHyp = get_expected_signal(pMeV, masses[ih]*1000.0, 1.0);
                if (dRef < 0 || dHyp < 0) continue;

                double nSigma = (dHyp/dRef - 1.0) / resoTPC[ih];

                xv.push_back(pG);
                yv.push_back(nSigma);

            }
            TGraph* g = new TGraph(xv.size(), xv.data(), yv.data());
            g->SetLineColor(colors[ih]);
            g->SetLineWidth(2);
            g->Draw("L SAME");
            leg.AddEntry(g, subs[ih], "l");
        }
        leg.Draw();
        c.Print(pdf);
    }

    // TOF
    for (int i = 0; i < 5; ++i) {
        TString species = subs[i];
        TH2F* hTOF = histTOF[species];
        hTOF->SetTitle(Form("TOF - %s", species.Data()));
        hTOF->Draw("COLZ");
        
        int idxRef = -1;
        for (int i = 0; i < 5; ++i) if (species.EqualTo(subs[i])) idxRef = i;
        if (idxRef<0) continue;


        TGraph* graphs[5];
        for (int ih = 0; ih < 5; ++ih) {
            Double_t y[npoints];
            for (int j=0; j<npoints; ++j) {
                Double_t p = pgrid[j];
                Double_t bref = p/TMath::Sqrt(masses[idxRef]*masses[idxRef]+p*p);
                Double_t bhyp = p/TMath::Sqrt(masses[ih]*masses[ih]+p*p);
                y[j] = (bref-bhyp)/(bhyp*bhyp*resoTOF[ih]);
            }
            graphs[ih] = new TGraph(npoints,pgrid,y);
            graphs[ih]->SetLineColor(colors[ih]);
            graphs[ih]->SetLineWidth(2);
            graphs[ih]->Draw("L SAME");
        }
        c.Update();
        TLegend* leg = new TLegend(0, 0.10, 0.15, 0.30, nullptr, "NDC");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);
        for (int ih = 0; ih < 5; ++ih) {
            leg->AddEntry(graphs[ih], subs[ih], "l");
        }
        leg->Draw();

        c.Print(pdf);
        for (int ih=0; ih<5; ++ih) {
            delete graphs[ih];}
    }

    c.Print(Form("%s]",pdf));
}


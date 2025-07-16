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

#include <AddTrees.h>

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

void nSigma_vs_P_Plot() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);

    const char* baseDir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";

    TChain chain("twotauchain");
    AddTrees(chain, baseDir);
    Long64_t nTotal = chain.GetEntries();
    nTotal = TMath::Min(nTotal, static_cast<Long64_t>(1e6));

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);
    chain.SetBranchStatus("fTrkTOFexpMom", 1);
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
    Float_t tofExpMom[NtrkMax] = {0};
    chain.SetBranchAddress("fTrkTPCsignal", tpcSignal);
    chain.SetBranchAddress("fTrkTPCinnerParam", inner);
    chain.SetBranchAddress("fTrkTOFexpMom", tofExpMom);
    for (int i = 0; i < 5; ++i) {
        chain.SetBranchAddress(Form("fTrkTPCnSigma%s", subs[i]), tpcNS[i]);
        chain.SetBranchAddress(Form("fTrkTOFnSigma%s", subs[i]), tofNS[i]);
    }

    Double_t masses[5] = {0.00051099895, 0.1056583755,  0.13957039, 0.493677, 0.93827208816};
    Double_t resoTOF[5]   = {0.013, 0.013, 0.013, 0.019, 0.020};
    Double_t resoTPC[5]   = {0.085, 0.072, 0.074, 0.09, 0.08}; 
    Int_t   colors[5]  = {kBlue, kGreen+2, kOrange+7, kMagenta+2, kCyan+1};
    const int npoints = 200;
    const Double_t pMin = 0.01, pMax = 5.0;
    Double_t pgrid[npoints];
    for (int i = 0; i < npoints; ++i) {
        pgrid[i] = pMin + i*(pMax - pMin)/(npoints - 1);
    }

    TH2F* histTPC[5];
    TH2F* histTOF[5];
    for (int i = 0; i < 5; ++i) {
        histTPC[i] = new TH2F(
          Form("tpc_%s", subs[i]),
          Form("TPC #sigma vs p (%s);p [GeV/c];n#sigma_{TPC}", subs[i]),
          100, pMin, pMax, 100, -10, 10
        );
        histTOF[i] = new TH2F(
          Form("tof_%s", subs[i]),
          Form("TOF #sigma vs p (%s);p [GeV/c];n#sigma_{TOF}", subs[i]),
          100, pMin, pMax, 100, -10, 10
        );
    }

    for (Long64_t ev = 0; ev < nTotal; ++ev) {
        chain.GetEntry(ev);
        for (int tr = 0; tr < NtrkMax; ++tr) {
            Float_t p = inner[tr];
            if (p <= 0) continue;
            for (int sp = 0; sp < 5; ++sp) {
                histTPC[sp]->Fill(p, tpcNS[sp][tr]);
                if (tofExpMom[tr] > 0)
                    histTOF[sp]->Fill(p, tofNS[sp][tr]);
            }
        }
    }

    TGraph* tpcCurves[5][5];
    TGraph* tofCurves[5][5];
    for (int ref = 0; ref < 5; ++ref) {
        double mRef = masses[ref];

        for (int hyp = 0; hyp < 5; ++hyp) {
            double mHyp = masses[hyp];

            std::vector<double> xv_tpc, yv_tpc;
            std::vector<double> xv_tof, yv_tof;
            xv_tpc.reserve(npoints);  yv_tpc.reserve(npoints);
            xv_tof.reserve(npoints);  yv_tof.reserve(npoints);

            for (int ip = 0; ip < npoints; ++ip) {
            double pg = pgrid[ip];
            if (pg < 0.1) continue;

            double dRef = get_expected_signal(pg*1000, mRef*1000, 1.0);
            double dHyp = get_expected_signal(pg*1000, mHyp*1000, 1.0);
            xv_tpc.push_back(pg);
            yv_tpc.push_back((dHyp/dRef - 1.0) / resoTPC[hyp]);

            double bRef = pg / TMath::Sqrt(mRef*mRef + pg*pg);
            double bHyp = pg / TMath::Sqrt(mHyp*mHyp + pg*pg);
            xv_tof.push_back(pg);
            yv_tof.push_back((bRef - bHyp) / (bHyp*bHyp * resoTOF[hyp]));
            }

            tpcCurves[ref][hyp] = new TGraph(xv_tpc.size(), xv_tpc.data(), yv_tpc.data());
            tpcCurves[ref][hyp]->SetLineColor(colors[hyp]);
            tpcCurves[ref][hyp]->SetLineWidth(2);

            tofCurves[ref][hyp] = new TGraph(xv_tof.size(), xv_tof.data(), yv_tof.data());
            tofCurves[ref][hyp]->SetLineColor(colors[hyp]);
            tofCurves[ref][hyp]->SetLineWidth(2);
        }
    }
    

    TCanvas* c = new TCanvas("c", "", 800, 600);
    TLegend* leg = new TLegend(0, 0.10, 0.15, 0.30);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    const char* pdf = "nsigma_vs_p_tracks_with_curves.pdf";
    c->Print(Form("%s[", pdf));
    c->SetLogz();

    for (int i = 0; i < 5; ++i) {
        leg->Clear();
        for (int hyp = 0; hyp < 5; ++hyp) {
            leg->AddEntry(tpcCurves[i][hyp], subs[hyp], "l");
        }
        histTPC[i]->SetTitle(Form("TPC #sigma vs p (%s)", subs[i]));
        histTPC[i]->Draw("COLZ");
        for (int hyp = 0; hyp < 5; ++hyp) {
            tpcCurves[i][hyp]->Draw("L SAME");
        }
        leg->Draw();
        c->Print(pdf);

        leg->Clear();
        for (int hyp = 0; hyp < 5; ++hyp) {
            leg->AddEntry(tofCurves[i][hyp], subs[hyp], "l");
        }
        histTOF[i]->SetTitle(Form("TOF #sigma vs p (%s)", subs[i]));
        histTOF[i]->Draw("COLZ");
        for (int hyp = 0; hyp < 5; ++hyp) {
            tofCurves[i][hyp]->Draw("L SAME");
        }
        leg->Draw();
        c->Print(pdf);
    }
    c->Print(Form("%s]", pdf));
}


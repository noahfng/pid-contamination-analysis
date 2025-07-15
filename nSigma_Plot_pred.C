#include <algorithm>
#include <vector>
#include <limits>
#include <cmath>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TF1.h"  
#include "TSpectrum.h"
#include "TSystem.h"
#include "RooFit.h"

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

std::vector<double> topBinCenters(TH1 *h, int nWanted)
{
   std::vector<std::pair<double,int>> bins;        
   for (int b = 1; b <= h->GetNbinsX(); ++b)
       bins.emplace_back(h->GetBinContent(b), b);

   std::partial_sort(bins.begin(), bins.begin()+nWanted, bins.end(),
                     std::greater<>());              

   std::vector<double> xc;
   for (int i = 0; i < nWanted; ++i)
       xc.push_back(h->GetXaxis()->GetBinCenter(bins[i].second));

   std::sort(xc.begin(), xc.end());                 
   return xc;
}

void nSigma_Plot_pred(){
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1);
    const char* baseDir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";
    TChain chain("twotauchain");
    AddTrees(chain, baseDir);

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);
    chain.SetBranchStatus("fTrkTPCsignal", 1);
    const char* subs[5] = {"El","Mu","Pi","Ka","Pr"};
    for (int i = 0; i < 5; ++i) {
        chain.SetBranchStatus(Form("fTrkTPCnSigma%s", subs[i]), 1);}
    Float_t inner[2], signal[2];
    Float_t tpcNS[5][2];
    chain.SetBranchAddress("fTrkTPCinnerParam", inner);
    chain.SetBranchAddress("fTrkTPCsignal", signal);
    for (int i = 0; i < 5; ++i) {
        chain.SetBranchAddress(Form("fTrkTPCnSigma%s", subs[i]), tpcNS[i]);}

    const Int_t nBins = 200;
    const Double_t pMin = 0.2, pMax = 0.3; 
    const Double_t xMin = -15.0, xMax = 15.0;
    const Int_t nParts = 5;
    const TString names[nParts] = {"e","#mu","#pi","K","p"};
    const Int_t colors[nParts] = {kBlue, kBlueYellow, kGreen+2, kOrange+7, kViolet};
    const Double_t resoTPC[5]   = {0.5, 0.075, 0.078, 0.09, 0.095}; 
    const Double_t masses[nParts] = {0.00051099895, 0.1056583755,  0.13957039, 0.493677, 0.93827208816};

    TH1F *hRes[nParts];
    for (int h = 0; h < nParts; ++h) {
        hRes[h] = new TH1F(Form("n#sigma_{%s}",names[h].Data()), 
        Form("n#sigma_{%s} for %.2f < p < %.2f GeV/c; n#sigma_{%s}; Counts", names[h].Data(),pMin, pMax, names[h].Data()),nBins, xMin, xMax);
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
    c->Print("nSigmaTPC_plot.pdf[");
    for (int h = 0; h < nParts; ++h) {
        c->Clear();
        hRes[h]->SetMarkerStyle(kFullCircle);
        hRes[h]->SetMarkerColor(kBlack);
        hRes[h]->SetLineColor(kRed);
        hRes[h]->SetMarkerSize(0.75); 
        hRes[h]->Draw("E1");

        const double pMid = 0.5*(pMin + pMax);            
        double mRef = masses[h]; 
        for (int hyp = 0; hyp < 5; ++hyp) {
            double mHyp = masses[hyp];
            double dRef = get_expected_signal(pMid*1000, mRef*1000, 1.0);
            double dHyp = get_expected_signal(pMid*1000, mHyp*1000, 1.0);                 

            double nSigmaPeak = (dHyp/dRef - 1.0)/ (resoTPC[hyp]);
            std::cout << "nSigmaPeak for " << subs[hyp] << " = " << nSigmaPeak << std::endl;
            double yMax = 1.05 * hRes[h]->GetMaximum();
            TLine *l = new TLine(nSigmaPeak, 0.0, nSigmaPeak, yMax);
            l->SetLineColor(colors[hyp]);
            l->SetLineStyle(1);
            l->Draw("same");

            TLatex txt;
            txt.SetTextSize(0.035);
            txt.SetTextAlign(22);                          
            txt.SetTextColor(colors[hyp]);
            txt.DrawLatex(nSigmaPeak, 0.9*yMax, subs[hyp]); 
        }

        const int nTryMax = 5;                         
        TF1  *bestFit = nullptr;
        double bestBIC = std::numeric_limits<double>::infinity();
        int    bestN   = 0;

        for (int nG = 1; nG <= nTryMax; ++nG) {

            std::ostringstream form;
            for (int i = 0; i < nG; ++i) {
                if (i) form << "+";
                form << "gaus(" << 3*i << ")";
            }
            TF1* fTmp = new TF1("fTmp", form.str().c_str(), xMin, xMax);

            auto mu = topBinCenters(hRes[h], nG);
            for (int i = 0; i < nG; ++i) {
                int    bin = hRes[h]->FindBin(mu[i]);
                double amp = hRes[h]->GetBinContent(bin);

                fTmp->SetParameter(3*i+0, amp);          // amplitude  (p = 3i)
                fTmp->SetParameter(3*i+1, mu[i]);        // mean       (p = 3i+1)
                fTmp->SetParameter(3*i+2, 0.30);         // sigma      (p = 3i+2)

                fTmp->SetParLimits(3*i+2, 0.05, 5.0);    
            }

            hRes[h]->Fit(fTmp, "QR0");                         
            double chi2 = fTmp->GetChisquare();
            int    k    = 3*nG;                                 
            double nObs = hRes[h]->GetEntries();
            double bic = chi2 +k * TMath::Log(nObs);
            const double deltaBICmin = 2;

            if (bestBIC - bic > deltaBICmin) {
                bestBIC = bic;
                bestN   = nG;
                delete bestFit;
                bestFit = (TF1*)fTmp->Clone(Form("best_%d",h));
            }
            else break;
        }

        bestFit->SetLineColor(kRed);
        bestFit->SetNpx(500);
        bestFit->SetLineWidth(2);
        bestFit->Draw("same");

        //for (int i = 0; i < bestN; ++i) {
        //    TF1* g = new TF1(Form("g_%d_%d",h,i),"gaus",xMin,xMax);
        //    for (int p = 0; p < 3; ++p)
        //        g->SetParameter(p, bestFit->GetParameter(3*i+p));
        //    g->SetLineStyle(2);
        //    g->SetLineColor(colors[i % nParts]);
        //    g->Draw("same");
        //}

        double chi2 = bestFit->GetChisquare();
        int    ndf  = bestFit->GetNDF();
        TPaveText* pt = new TPaveText(0.0,0.9,0.20,1.0,"NDC");
        pt->AddText(Form("N_{Gauss} = %d, #chi^{2}/NDF = %.2f, BIC = %.1f", bestN, chi2/ndf, bestBIC));
        pt->SetFillColorAlpha(0,0);
        pt->Draw("same");
        

        c->Print("nSigmaTPC_plot.pdf");} 
    
    c->Print("nSigmaTPC_plot.pdf]");
}
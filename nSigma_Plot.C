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
#include "TLine.h"
#include "TLatex.h"
#include "TPaveText.h"

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

void nSigma_Plot(double pStart = 0.1, double pEnd = 2.0, double step = 0.1, double muWindow = 1.0)   
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1);

    const char *baseDir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";
    TChain chain("twotauchain");
    AddTrees(chain, baseDir);

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);

    const char *subs[5] = {"El", "Mu", "Pi", "Ka", "Pr"};
    for (int i = 0; i < 5; ++i)
        chain.SetBranchStatus(Form("fTrkTPCnSigma%s", subs[i]), 1);

    Float_t inner[2];
    Float_t tpcNS[5][2];
    chain.SetBranchAddress("fTrkTPCinnerParam", inner);
    for (int i = 0; i < 5; ++i)
        chain.SetBranchAddress(Form("fTrkTPCnSigma%s", subs[i]), tpcNS[i]);

    const Int_t   nBins   = 200;
    const Double_t xMin   = -25.0, xMax = 150.0;
    const Int_t   nParts  = 5;
    const TString names[nParts]   = {"e", "#mu", "#pi", "K", "p"};
    const Int_t   colors[nParts] = {kBlue, kGreen+2, kOrange+7, kMagenta+2, kCyan+1};
    const Double_t resoTPC[nParts] = {0.085, 0.072, 0.074, 0.09, 0.08}; 
    const Double_t masses[nParts]  = {0.00051099895, 0.1056583755, 0.13957039, 0.493677, 0.93827208816};

    Long64_t nEntries = std::min(chain.GetEntries(), static_cast<Long64_t>(1e6));

    for (int ref = 0; ref < nParts; ++ref) {

        TString pdfName = Form("nSigmaTPC_%s.pdf", names[ref].Data());
        TCanvas *c = new TCanvas("c","n#sigma("+names[ref]+")",950,700);
        c->SetLeftMargin(0.15); 
        c->SetGrid(); 
        c->SetLogy();
        c->Print(pdfName+"[");

        for (double pMin = pStart; pMin < pEnd; pMin += step) {
            double pMax = pMin + step;
            TH1F *h = new TH1F(Form("h_%s_%g",names[ref].Data(),pMin),
                Form("n#sigma_{%s}, %.1f < p < %.1f GeV/c; n#sigma; Counts",
                     names[ref].Data(), pMin, pMax),
                nBins,xMin,xMax);
            h->SetMarkerStyle(kFullCircle); 
            h->SetMarkerSize(0.75);
            h->SetMarkerColor(kBlack); 
            h->SetLineColor(kBlack);
            TLegend* leg = new TLegend(0, 0.10, 0.15, 0.30);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);

            for(Long64_t ev=0; ev<nEntries; ++ev){
                chain.GetEntry(ev);
                for(int t=0;t<2;++t){
                    float pG=inner[t];
                    if(pG<pMin||pG>=pMax) continue;
                    if(!TMath::IsNaN(tpcNS[ref][t])) h->Fill(tpcNS[ref][t]);
                }
            }

            double pMid = 0.5*(pMin+pMax);
            double dRef = get_expected_signal(pMid*1000, masses[ref]*1000, 1.0);
            struct Peak {
                double A, mu, sigma; 
                int id; 
                std::vector<int> merged_ids;
            };
            std::vector<Peak> seeds;
            double yMax = 1.05*h->GetMaximum();

            for(int hyp=0; hyp<nParts; ++hyp){
                double dHyp = get_expected_signal(pMid*1000, masses[hyp]*1000, 1.0);
                double sigma0 = (resoTPC[hyp]/resoTPC[ref]) * (dHyp/dRef);
                double mu = (dHyp/dRef -1.0)/resoTPC[hyp];
                int    bin = h->FindBin(mu);        
                double amp = h->GetBinContent(bin);
                sigma0 = std::clamp(sigma0, 0.5, 15.0); 
                if(dRef<0||dHyp<0) continue;
                if(mu<xMin||mu>xMax) continue;

                seeds.push_back({amp, mu, sigma0, hyp, {hyp}});
            }
            std::vector<Peak> merged;
            std::sort(seeds.begin(), seeds.end(),
                    [](auto &a, auto &b){return a.mu<b.mu;});

            for (size_t i=0; i<seeds.size(); ++i) {
                if (!merged.empty() &&
                    std::abs(seeds[i].mu - merged.back().mu) <
                    1.0*std::max(seeds[i].sigma, merged.back().sigma))          
                {
                    auto &p = merged.back();
                    double I1 = p.A * p.sigma * sqrt(2*M_PI);
                    double I2 = seeds[i].A * seeds[i].sigma * sqrt(2*M_PI);
                    double It = I1 + I2;

                    double mu_eff = (I1*p.mu + I2*seeds[i].mu) / It;
                    double var_eff = (I1*(p.sigma*p.sigma + (p.mu-mu_eff)*(p.mu-mu_eff)) +
                                    I2*(seeds[i].sigma*seeds[i].sigma +
                                        (seeds[i].mu-mu_eff)*(seeds[i].mu-mu_eff))) / It;
                    double sig_eff = std::sqrt(var_eff);
                    double A_eff   = It / (std::sqrt(2*M_PI)*sig_eff);

                    std::vector<int> ids = p.merged_ids;
                    ids.insert(ids.end(), seeds[i].merged_ids.begin(), seeds[i].merged_ids.end());
                    p = {A_eff, mu_eff, sig_eff, -1, ids}; 
                } else {
                    merged.push_back(seeds[i]);
                }
            }
            size_t nG = merged.size();
            if (nG == 0) { c->Print(pdfName); delete h; continue; }

            std::ostringstream form;
            for (size_t i = 0; i < nG; ++i) {
                if (i) form << "+";
                form << "gaus(" << 3*i << ")";
            }
            TF1 *sum=new TF1("sum",form.str().c_str(),xMin,xMax);
            for (size_t i = 0; i < nG; ++i) {
                const auto &p = merged[i];

                sum->SetParLimits(3*i, 0.0, std::max(h->GetMaximum()*1.2, p.A*1.05));
                sum->SetParameter (3*i, p.A);

                double dMu = std::max(muWindow, 0.10*std::abs(p.mu));
                sum->SetParLimits(3*i+1, p.mu - dMu, p.mu + dMu);
                sum->SetParameter (3*i+1, p.mu);

                sum->SetParLimits(3*i+2, 0.5*p.sigma, 3.0*p.sigma);
                sum->SetParameter (3*i+2, p.sigma);
            }
            h->Fit(sum,"QR0");
            c->Clear(); 
            h->Draw("E1");
            sum->SetLineColor(kRed); 
            sum->SetLineWidth(2); 
            sum->SetNpx(500); 
            sum->Draw("same");

            for (const auto &pk : merged) {
                double mu = pk.mu;
                TLine *l = new TLine(mu, 0, mu, yMax);
                int col  = (pk.id >= 0) ? colors[pk.id] : kGray+2;
                l->SetLineColor(col);
                l->Draw();

                
                TString label;
                if (pk.merged_ids.size()==1)
                    label = names[pk.merged_ids[0]];
                else {
                    for (size_t j=0;j<pk.merged_ids.size();++j){
                        if (j) label += " + ";
                        label += names[pk.merged_ids[j]];
                    }
                }
            }

            for (size_t i = 0; i < nG; ++i) {
                if (sum->GetParameter(3*i) <= 0) continue;

                TF1 *g = new TF1(Form("g_%d_%zu", ref, i), "gaus", xMin, xMax);
                g->SetParameters(&sum->GetParameters()[3*i]);

                int col = (merged[i].id >= 0) ? colors[merged[i].id] : kGray+2;
                g->SetLineColor(col);
                g->SetLineStyle(2);
                g->Draw("same");

                TString label;
                if (merged[i].merged_ids.size() == 1) {
                    label = names[merged[i].merged_ids[0]];
                } else {
                    for (size_t j = 0; j < merged[i].merged_ids.size(); ++j) {
                        if (j > 0) label += " + ";
                        label += names[merged[i].merged_ids[j]];
                    }
                }
                leg->AddEntry(g, label, "l");
            }

            TPaveText *pt=new TPaveText(0.02,0.90,0.25,0.99,"NDC");
            pt->AddText(Form("#chi^{2}/NDF = %.2f", sum->GetChisquare()/sum->GetNDF()));
            pt->AddText(Form("N_{Gauss} = %zu", nG)); 
            pt->SetFillColorAlpha(0,0); pt->Draw("same");
            leg->Draw();
            c->Print(pdfName);
            delete sum; 
            delete h;
            delete leg;
        }
        c->Print(pdfName+"]"); 
        delete c;
    }
}


#include <algorithm>
#include <vector>
#include <cmath>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TF1.h"  
#include "TLine.h"
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

void nSigma_Plot(){
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1);

    const char *baseDir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";
    TChain chain("twotauchain");
    AddTrees(chain, baseDir);

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);
    chain.SetBranchStatus("fTrkTOFexpMom", 1);

    const char *subs[5] = {"El", "Mu", "Pi", "Ka", "Pr"};
    for (int i = 0; i < 5; ++i){
        chain.SetBranchStatus(Form("fTrkTPCnSigma%s", subs[i]), 1);
        chain.SetBranchStatus(Form("fTrkTOFnSigma%s", subs[i]), 1);
    }
    Float_t inner[2], tofExpMom[2];
    Float_t tpcNS[5][2], tofNS[5][2];
    chain.SetBranchAddress("fTrkTPCinnerParam", inner);
    chain.SetBranchAddress("fTrkTOFexpMom", tofExpMom);
    for (int i = 0; i < 5; ++i){
        chain.SetBranchAddress(Form("fTrkTPCnSigma%s", subs[i]), tpcNS[i]);
        chain.SetBranchAddress(Form("fTrkTOFnSigma%s", subs[i]), tofNS[i]);
    }
    const Int_t   nBins   = 500;
    const Double_t xMin   = -25.0, xMax = 50.0;
    const Int_t countThreshold = 0;
    const Double_t pStart = 0.2, pEnd = 2.0, step = 0.1;
    const Double_t muWindow = 1.0;
    const Double_t mergeDistanceFactor = 1.0;
    const bool useTOF = true;
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
        //c->SetGrid(); 
        c->SetLogy();
        c->Print(pdfName+"[");

        for (double pMin = pStart; pMin < pEnd; pMin += step) {
            double pMax = pMin + step;
            TH1F *h = new TH1F(Form("h_%s_%g",names[ref].Data(),pMin),
                Form("n#sigma_{%s}, %.1f < p < %.1f GeV/c; n#sigma; Counts",
                     names[ref].Data(), pMin, pMax),
                nBins,xMin,xMax);
            h->Sumw2(true);
            h->SetMarkerStyle(kFullCircle); 
            h->SetMarkerSize(0.75);
            h->SetMarkerColor(kBlack); 
            h->SetLineColor(kBlack);
            TLegend* leg = new TLegend(0, 0.10, 0.15, 0.30);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);

            for(Long64_t ev = 0; ev < nEntries; ++ev){
                chain.GetEntry(ev);
                for(int t = 0; t < 2; ++t){
                    if (useTOF && tofExpMom[t] < 0) continue;
                    float pG=inner[t];
                    if(pG<pMin||pG>=pMax) continue;
                    if(!TMath::IsNaN(tpcNS[ref][t])) h->Fill(tpcNS[ref][t]);
                }
            }
            for (int b = 1; b <= h->GetNbinsX(); ++b) {
                if (h->GetBinContent(b) < countThreshold){
                    h->SetBinContent(b, 0);
                }
            }

            double pMid = 0.5*(pMin+pMax);
            double dRef = get_expected_signal(pMid*1000, masses[ref]*1000, 1.0);
            struct Peak {
                double A, mu, sigma; 
                int id; 
                std::vector<int> merged_ids;
                bool alwaysSeparate = false;
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

                seeds.push_back({ amp, mu, sigma0, hyp, {hyp}, false });
            }

            Peak muPiCombined{0,0,0,-1,{},false};
            bool sawMu = false, sawPi = false;
            for (auto it = seeds.begin(); it != seeds.end(); ) {
                if (it->id == 1 || it->id == 2) {
                    if (!sawMu && it->id == 1) {
                        muPiCombined = *it;
                        sawMu = true;
                    }
                    else {
                        double I_old = muPiCombined.A * muPiCombined.sigma * sqrt(2*M_PI);
                        double I_new = it->A * it->sigma * sqrt(2*M_PI);
                        double I_tot = I_old + I_new;
                        double mu_eff = (I_old*muPiCombined.mu + I_new*it->mu) / I_tot;
                        double var_eff = (
                            I_old*(muPiCombined.sigma*muPiCombined.sigma 
                                + (muPiCombined.mu-mu_eff)*(muPiCombined.mu-mu_eff))
                        + I_new*(it->sigma*it->sigma 
                                + (it->mu-mu_eff)*(it->mu-mu_eff))
                        ) / I_tot;
                        muPiCombined.A     = I_tot / (muPiCombined.sigma * sqrt(2*M_PI));
                        muPiCombined.mu    = mu_eff;
                        muPiCombined.sigma = sqrt(var_eff);
                    }
                    muPiCombined.merged_ids.push_back(it->id);
                    sawPi = sawPi || (it->id == 2);
                    it = seeds.erase(it);
                }
                else {
                    ++it;
                }
            }
            if (sawMu || sawPi) {
                muPiCombined.id = 1;
                muPiCombined.alwaysSeparate = std::find(muPiCombined.merged_ids.begin(), muPiCombined.merged_ids.end(), ref) != muPiCombined.merged_ids.end();
                seeds.push_back(muPiCombined);
            }

            for (auto &s : seeds) {
                if (std::find(s.merged_ids.begin(), s.merged_ids.end(), ref) != s.merged_ids.end()) {
                    s.alwaysSeparate = true;
                }
            }

            std::sort(seeds.begin(), seeds.end(),
                    [](auto &a, auto &b){ return a.mu < b.mu; });

            std::vector<Peak> merged;
            for (auto &s : seeds) {
                if (!merged.empty()
                    && !s.alwaysSeparate
                    && !merged.back().alwaysSeparate
                    && std::abs(s.mu - merged.back().mu)
                        < mergeDistanceFactor * std::max(s.sigma, merged.back().sigma))
                {
                    Peak &p = merged.back();
                    double I1 = p.A * p.sigma * sqrt(2*M_PI);
                    double I2 = s.A * s.sigma * sqrt(2*M_PI);
                    double It = I1 + I2;
                    double mu_eff = (I1*p.mu + I2*s.mu) / It;
                    double var_eff = (I1*(p.sigma*p.sigma + (p.mu-mu_eff)*(p.mu-mu_eff))
                                    + I2*(s.sigma*s.sigma + (s.mu-mu_eff)*(s.mu-mu_eff)))
                                    / It;
                    p.sigma = sqrt(var_eff);
                    p.mu    = mu_eff;
                    p.A     = It / (sqrt(2*M_PI) * p.sigma);
                    p.merged_ids.insert(p.merged_ids.end(),
                                        s.merged_ids.begin(),
                                        s.merged_ids.end());
                }
                else {
                    merged.push_back(s);
                }
            }
            size_t nG = merged.size();
            if (nG == 0) { c->Print(pdfName); delete h; continue; }

            std::ostringstream form;
            for (size_t i = 0; i < nG; ++i) {
                if (i) form << "+";
                form << "gaus(" << 3 * i << ")";
            }
            size_t idxLine = 3 * nG;
            form << "+pol0(" << idxLine << ")";
            TF1 *sum=new TF1("sum",form.str().c_str(),xMin,xMax);
            for (size_t i = 0; i < nG; ++i) {
                const auto &p = merged[i];
                
                sum->SetParLimits(idxLine, 0, 10.0);
                sum->SetParameter(idxLine, 1.0);

                sum->SetParLimits(3 * i, 0.0, std::max(h->GetMaximum() * 1.2, p.A * 1.05));
                sum->SetParameter (3 * i, p.A);

                double dMu = std::max(muWindow, 0.10 * std::abs(p.mu));
                sum->SetParLimits(3 * i + 1, p.mu - dMu, p.mu + dMu);
                sum->SetParameter (3 * i + 1, p.mu);

                sum->SetParLimits(3 * i + 2, 0.5 * p.sigma, 3.0 * p.sigma);
                sum->SetParameter (3 * i + 2, p.sigma);

            }
            h->Fit(sum, "RLQ0S");
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
                
                auto &ids = merged[i].merged_ids;
                std::sort(ids.begin(), ids.end());
                ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

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
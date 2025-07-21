#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>

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
#include "TFile.h"
#include "TGraph.h"
#include "TSystem.h"


#include <AddTrees.h>
#include <get_expected_signal.h>
#include <getReso.h>

std::vector<Double_t> topBinCenters(TH1 *h, Int_t nWanted)
{
   std::vector<std::pair<Double_t, Int_t>> bins;        
   for (Int_t b = 1; b <= h->GetNbinsX(); ++b)
       bins.emplace_back(h->GetBinContent(b), b);

   std::partial_sort(bins.begin(), bins.begin()+nWanted, bins.end(),
                     std::greater<>());              

   std::vector<Double_t> xc;
   for (Int_t i = 0; i < nWanted; ++i)
       xc.push_back(h->GetXaxis()->GetBinCenter(bins[i].second));

   std::sort(xc.begin(), xc.end());                 
   return xc;
}

void nSigma_Plot(){

    const Int_t   nBins   = 500;
    const Double_t xMin   = -30.0, xMax = 70.0;
    const Double_t pStart = 0.6, pEnd = 0.8, step = 0.1;
    const Double_t muWindow = 1.0;
    const Double_t mergeDistanceFactor = 1.0;
    const Double_t nEntriesLimit = 1e6;
    const Bool_t TOFfilter = false;
    const Bool_t plotTPC = true;
    const Bool_t plotTOF = false;
    const Bool_t PeakZoom = true;
    const Bool_t manualPredictPeaks = false;

    if (manualPredictPeaks) gROOT->SetBatch(kFALSE);
    else gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1);

    const Char_t *baseDir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";
    TChain chain("twotauchain");
    AddTrees(chain, baseDir);

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);
    chain.SetBranchStatus("fTrkTOFexpMom", 1);

    const Char_t *subs[5] = {"El", "Mu", "Pi", "Ka", "Pr"};
    for (Int_t i = 0; i < 5; ++i){
        chain.SetBranchStatus(Form("fTrkTPCnSigma%s", subs[i]), 1);
        chain.SetBranchStatus(Form("fTrkTOFnSigma%s", subs[i]), 1);
    }
    Float_t inner[2], tofExpMom[2];
    Float_t tpcNS[5][2], tofNS[5][2];
    chain.SetBranchAddress("fTrkTPCinnerParam", inner);
    chain.SetBranchAddress("fTrkTOFexpMom", tofExpMom);
    for (Int_t i = 0; i < 5; ++i){
        chain.SetBranchAddress(Form("fTrkTPCnSigma%s", subs[i]), tpcNS[i]);
        chain.SetBranchAddress(Form("fTrkTOFnSigma%s", subs[i]), tofNS[i]);
    }
    
    const Int_t   nParts  = 5;
    const TString names[nParts]   = {"e", "#mu", "#pi", "K", "p"};
    const Int_t   colors[nParts] = {kBlue, kGreen+2, kOrange+7, kMagenta+2, kCyan+1};
    const Double_t resoTPC[nParts] = {0.085, 0.072, 0.074, 0.09, 0.08}; 
    const Double_t resoTOF[nParts]   = {0.013, 0.013, 0.013, 0.019, 0.020};
    const Double_t masses[nParts]  = {0.00051099895, 0.1056583755, 0.13957039, 0.493677, 0.93827208816};

    Long64_t nEntries = std::min(chain.GetEntries(), static_cast<Long64_t>(nEntriesLimit));

    auto drawNSigma = [&](Bool_t isTPCmode) {
        TString suffix = isTPCmode ? "TPC" : "TOF";
        Double_t pLoopStart = isTPCmode ? pStart : std::max(pStart, 0.4);
        for (Int_t ref = 0; ref < nParts; ++ref) {
            TString pdfName = Form("nSigma%s_%s.pdf",suffix.Data(), names[ref].Data());
            TCanvas *c = new TCanvas("c","n#sigma("+names[ref]+")", 950, 700);
            c->SetLeftMargin(0.15); 
            //c->SetGrid(); 
            c->SetLogy();
            c->Print(pdfName+"[");
            const Int_t nSteps = static_cast<Int_t>(std::floor((pEnd - pLoopStart) / step + 0.5));
            for (Int_t i = 0; i < nSteps; ++i) {
                Double_t pMin = pLoopStart + i * step;
                Double_t pMax = pMin + step;
                TH1F *h = new TH1F(Form("n#sigma_%s_%g < p < %g (%s)",names[ref].Data(),pMin, pMax, suffix.Data()),
                    Form("n#sigma_{%s} %g < p < %g GeV/c (%s); n#sigma_{%s}; Counts",
                        names[ref].Data(), pMin, pMax, suffix.Data(), names[ref].Data()),
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
                    for(Int_t t = 0; t < 2; ++t){
                        if (tofExpMom[t] < 0 && (!isTPCmode || TOFfilter)) continue;
                        Float_t pG=inner[t];
                        if(pG < pMin || pG >= pMax) continue;
                        Float_t val = isTPCmode ? tpcNS[ref][t] : tofNS[ref][t];
                        if (!TMath::IsNaN(val)) h->Fill(val);
                    }
                }

                struct Peak {Double_t A, mu, sigma; Int_t id; std::vector<Int_t> merged_ids;Bool_t alwaysSeparate = false;};
                std::vector<Peak> seeds;
                Double_t yMax = 1.05 * h->GetMaximum();
                Double_t pMid = 0.5 * (pMin+pMax);

                Double_t refMass = masses[ref];
                Double_t dRef = get_expected_signal(pMid * 1000, masses[ref] * 1000, 1.0);
                Double_t bRef = pMid / TMath::Sqrt(pMid * pMid + refMass * refMass);

                for(Int_t hyp=0; hyp < nParts; ++hyp){
                    Double_t hypMass = masses[hyp];
                    Double_t dHyp = get_expected_signal(pMid * 1000, hypMass * 1000, 1.0);
                    Double_t bHyp = pMid / TMath::Sqrt(pMid * pMid + hypMass * hypMass);
                    if(dRef < 0 || dHyp < 0) continue;

                    Double_t sigma0, mu;
                    if (isTPCmode) {
                        Double_t resoHypAbs = getReso(kTPC, (Char_t*)subs[hyp], pMid); 
                        Double_t resoRefAbs = getReso(kTPC, (Char_t*)subs[ref],  pMid); 

                        Double_t fracHyp = resoHypAbs / dHyp;
                        Double_t fracRef = resoRefAbs / dRef;

                        sigma0 = (fracHyp / fracRef) * (dHyp / dRef);
                        mu     = (dHyp/dRef - 1.0) / fracRef;
                    } else {
                        Double_t resoHypAbs = getReso(kTOF, (Char_t*)subs[hyp], pMid);
                        Double_t resoRefAbs = getReso(kTOF, (Char_t*)subs[ref], pMid);
                        Double_t fracHyp = resoHypAbs / dHyp;
                        Double_t fracRef = resoRefAbs / dRef;
                        sigma0 = (fracHyp / fracRef) * (1.0 / (bHyp * bHyp));
                        mu     = (dHyp/dRef - 1.0) / fracRef;
                    }
                    sigma0 = std::clamp(sigma0, 0.5, 15.0);
                    if (mu < xMin || mu > xMax) continue;

                    Int_t    bin = h->FindBin(mu);        
                    Double_t amp = h->GetBinContent(bin);

                    seeds.push_back({amp, mu, sigma0, hyp, {hyp}, false});
                }

                Peak muPiCombined{0,0,0,-1,{},false};
                Bool_t sawMu = false, sawPi = false;
                for (auto it = seeds.begin(); it != seeds.end();) {
                    if (it->id == 1 || it->id == 2) {
                        if (!sawMu && it->id == 1) {
                            muPiCombined = *it;
                            sawMu = true;
                        }
                        else {
                            Double_t I_old = muPiCombined.A * muPiCombined.sigma * TMath::Sqrt(2*TMath::Pi());
                            Double_t I_new = it->A * it->sigma * TMath::Sqrt(2*TMath::Pi());
                            Double_t I_tot = I_old + I_new;
                            Double_t mu_eff = (I_old*muPiCombined.mu + I_new*it->mu) / I_tot;
                            Double_t var_eff = (
                                I_old*(muPiCombined.sigma*muPiCombined.sigma 
                                    + (muPiCombined.mu-mu_eff)*(muPiCombined.mu-mu_eff))
                            + I_new*(it->sigma*it->sigma 
                                    + (it->mu-mu_eff)*(it->mu-mu_eff))
                            ) / I_tot;
                            muPiCombined.A     = I_tot / (muPiCombined.sigma * TMath::Sqrt(2*TMath::Pi()));
                            muPiCombined.mu    = mu_eff;
                            muPiCombined.sigma = TMath::Sqrt(var_eff);
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

                std::sort(seeds.begin(), seeds.end(), [](auto &a, auto &b){ return a.mu < b.mu; });
                std::vector<Peak> merged;
                for (auto &s : seeds) {
                    if (!merged.empty() && !s.alwaysSeparate && !merged.back().alwaysSeparate && std::abs(s.mu - merged.back().mu) < mergeDistanceFactor * std::max(s.sigma, merged.back().sigma))
                    {
                        Peak &p = merged.back();
                        Double_t I1 = p.A * p.sigma * TMath::Sqrt(2*TMath::Pi());
                        Double_t I2 = s.A * s.sigma * TMath::Sqrt(2*TMath::Pi());
                        Double_t It = I1 + I2;
                        Double_t mu_eff = (I1*p.mu + I2*s.mu) / It;
                        Double_t var_eff = (I1*(p.sigma*p.sigma + (p.mu-mu_eff)*(p.mu-mu_eff))
                                        + I2*(s.sigma*s.sigma + (s.mu-mu_eff)*(s.mu-mu_eff)))
                                        / It;
                        p.sigma = TMath::Sqrt(var_eff);
                        p.mu    = mu_eff;
                        p.A     = It / (TMath::Sqrt(2*TMath::Pi()) * p.sigma);
                        p.merged_ids.insert(p.merged_ids.end(),s.merged_ids.begin(), s.merged_ids.end());
                    } 
                    else {
                        merged.push_back(s);
                    }
                }
                
                if (PeakZoom) {
                    const Double_t Nsigma = 8.0;
                    Double_t x_low  = xMin;
                    Double_t x_high = xMax;

                    if (!merged.empty()) {
                    x_low  = merged[0].mu - Nsigma * merged[0].sigma;
                    x_high = merged[0].mu + Nsigma * merged[0].sigma;
                    for (size_t i = 1; i < merged.size(); ++i) {
                        Double_t this_lo = merged[i].mu - Nsigma * merged[i].sigma;
                        Double_t this_hi = merged[i].mu + Nsigma * merged[i].sigma;
                        x_low  = std::min(x_low,  this_lo);
                        x_high = std::max(x_high, this_hi);
                    }
                    x_low  = std::max(x_low,  xMin);
                    x_high = std::min(x_high, xMax);
                    }
                    h->GetXaxis()->SetRangeUser(x_low, x_high);
                }

                Int_t manualNGauss = 0;
                std::vector<double> manualMeans;
                if (manualPredictPeaks) {
                    c->cd();
                    c->Clear();
                    h->Draw("E1");
                    c->Modified();
                    c->Update();
                    c->RaiseWindow();
                    gSystem->ProcessEvents();
                    std::cout << "How many Gaussians? ";
                    std::cin >> manualNGauss;
                    if (manualNGauss > 0) {
                        manualMeans.resize(manualNGauss);
                        std::cout << "Enter " << manualNGauss << " mean positions space-separated: ";
                        for (int i = 0; i < manualNGauss; ++i) {
                            std::cin >> manualMeans[i];
                        }
                    } else {
                        std::cerr << "Invalid number of Gaussians; falling back to automatic mode." << std::endl;
                        manualNGauss = 0;
                    }
                }

                size_t nG;
                std::vector<double> means;
                if (manualPredictPeaks && manualNGauss > 0) {
                    nG = manualNGauss;
                    means = manualMeans;
                } else {
                    nG = merged.size();
                    means.reserve(nG);
                    for (auto &pk : merged) means.push_back(pk.mu);
                }
                if (nG < 1) { c->Print(pdfName); delete h; continue; }
                
                std::ostringstream form;
                for (size_t i = 0; i < nG; ++i) {
                    if (i) form << "+"; 
                    form << "gaus(" << 3 * i << ")";
                }
                size_t idxLine = 3 * nG;
                form << "+pol0(" << idxLine << ")";
                TF1 *sum = new TF1("sum", form.str().c_str(), xMin, xMax);
                sum->SetParLimits(idxLine, 0, 10.0);
                sum->SetParameter(idxLine, 1.0);
                for (size_t i = 0; i < nG; ++i) {
                    if (manualPredictPeaks && manualNGauss > 0) {
                        double mu0 = means[i];

                        sum->SetParLimits(3*i, 0.0, h->GetMaximum()*2);
                        sum->SetParameter(3*i, h->GetMaximum()/nG);

                        sum->SetParLimits(3*i+1, mu0 - muWindow, mu0 + muWindow);
                        sum->SetParameter(3*i+1, mu0);

                        sum->SetParLimits(3*i+2, 0.5, 10.0);
                        sum->SetParameter(3*i+2, 2.0);
                    }
                    else {
                        const auto &p = merged[i];

                        sum->SetParLimits(3*i, 0.0, std::max(h->GetMaximum()*1.2, p.A*1.05));
                        sum->SetParameter(3*i, p.A);

                        double dMu = std::max(muWindow, 0.1*std::abs(p.mu));
                        sum->SetParLimits(3*i+1, p.mu - dMu, p.mu + dMu);
                        sum->SetParameter(3*i+1, p.mu);

                        sum->SetParLimits(3*i+2, 0.5*p.sigma, 3.0*p.sigma);
                        sum->SetParameter(3*i+2, p.sigma);
                    }
                }
                
                h->Fit(sum, "RLQ0S");
                c->Clear();
                h->Draw("E1");
                sum->SetLineColor(kRed); 
                sum->SetLineWidth(2); 
                sum->SetNpx(500); 
                sum->Draw("same");
                if (!manualPredictPeaks) {
                    for (const auto &pk : merged) {
                        Double_t mu = pk.mu;
                        TLine *l = new TLine(mu, 0, mu, yMax);
                        Int_t col  = (pk.id >= 0) ? colors[pk.id] : kGray+2;
                        l->SetLineColor(col);
                        l->Draw();
                    }
                }
                else {
                    for (Int_t i = 0; i < nG; ++i) {
                        Double_t mu = means[i];
                        TLine *l = new TLine(mu, 0, mu, yMax);
                        l->SetLineColor(colors[i % nParts]);
                        l->Draw();
                    }
                }

                for (Int_t i = 0; i < nG; ++i) {
                    if (sum->GetParameter(3*i) <= 0) continue;
                    TF1 *g = new TF1(Form("g_%d_%d", ref, i), "gaus", xMin, xMax);
                    g->SetParameters(&sum->GetParameters()[3*i]);
                    if (!manualPredictPeaks) {
                        Int_t col = (merged[i].id >= 0) ? colors[merged[i].id] : kGray+2;
                        g->SetLineColor(col);
                    } else {
                        g->SetLineColor(colors[i % nParts]);   
                    }
                    g->SetLineStyle(2);
                    g->SetNpx(500);
                    g->Draw("same");

                    TString label;
                    if (!manualPredictPeaks) {
                        auto &ids = merged[i].merged_ids;
                        std::sort(ids.begin(), ids.end());
                        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
                        if (ids.size()==1) {
                            label = names[ids[0]];
                        } else {
                            for (size_t j=0; j<ids.size(); ++j) {
                                if (j) label += " + ";
                                label += names[ids[j]];
                            }
                        }
                    } else {
                        label = Form("Peak %d", i);
                    }
                    leg->AddEntry(g, label, "l");
                }
 
                TPaveText *pt=new TPaveText(0.02,0.90,0.25,0.99,"NDC");
                pt->AddText(Form("#chi^{2}/NDF = %.2f", sum->GetChisquare()/sum->GetNDF()));
                pt->SetFillColorAlpha(0,0); 
                pt->Draw("same");
                leg->Draw();
                c->Print(pdfName);
                delete sum; 
                delete h;
                delete leg;
            }
            c->Print(pdfName+"]"); 
        }
    };
    if (plotTPC)  drawNSigma(true);   
    if (plotTOF)  drawNSigma(false);  
}
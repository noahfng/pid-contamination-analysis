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
#include "TSystem.h"

#include <AddTrees.h>
#include <get_expected_signal.h>
#include <getReso.h>


void nSigmaDebug(){
    const Int_t   nBins   = 500;
    const Double_t xMin   = -30.0, xMax = 70.0;
    const Double_t pStart = 0.35, pEnd = 0.45, step = 0.1;
    const Double_t muWindow = 2.0;
    const Double_t mergeDistanceFactor = 1.0;
    const Double_t nEntriesLimit = 1e7;
    const Bool_t TOFfilter = false;
    const Bool_t plotTPC = true;
    const Bool_t plotTOF = false;

    gROOT->SetBatch(kFALSE);
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
        const Double_t pMin   = isTPCmode ? pStart : std::max(pStart, 0.4);
        const Int_t    nSteps = Int_t(std::floor((pEnd - pMin) / step + 0.5));
        std::vector<Double_t> pEdges(nSteps+1);
        for (int i = 0; i <= nSteps; ++i) {
            pEdges[i] = pMin + i * step;
        }
        TString suffix = isTPCmode ? "TPC" : "TOF";
        std::vector<std::vector<TH1F*>> hists(nParts, std::vector<TH1F*>(nSteps,nullptr));
        for (int pid = 0; pid < nParts; ++pid) {
            for (int i = 0; i < nSteps; ++i) {
                TString name1 = Form("n#sigma_{%s} %g < p < %g GeV/c (%s)", 
                        names[pid].Data(), pEdges[i], pEdges[i+1], suffix.Data());
                TString name2 = Form("n#sigma_{%s} %g < p < %g GeV/c (%s); n#sigma_{%s}; Counts",
                        names[pid].Data(), pEdges[i], pEdges[i+1], suffix.Data(), names[pid].Data());
                hists[pid][i] = new TH1F(name1, name2, nBins, xMin, xMax);
                hists[pid][i]->Sumw2(true);
                hists[pid][i]->SetMarkerStyle(kFullCircle);
                hists[pid][i]->SetMarkerSize(0.75);
                hists[pid][i]->SetMarkerColor(kBlack);
                hists[pid][i]->SetLineColor(kBlack);
            }
        }
        for (Long64_t ev = 0; ev < nEntries; ++ev) {
            chain.GetEntry(ev);
            for (int t = 0; t < 2; ++t) {
                if (tofExpMom[t] < 0 && (!isTPCmode || TOFfilter)) 
                    continue;
                Double_t pG = inner[t];
                int bin = std::lower_bound(pEdges.begin(), pEdges.end(), pG)
                        - pEdges.begin() - 1;
                if (bin < 0 || bin >= nSteps) 
                    continue;
                for (int pid = 0; pid < nParts; ++pid) {
                    Float_t val = isTPCmode ? tpcNS[pid][t] : tofNS[pid][t];
                    if (!TMath::IsNaN(val))
                        hists[pid][bin]->Fill(val);
                }
            }
        }        
        for (int ref = 0; ref < nParts; ++ref) {
            TString pdfName = Form("nSigma%s_%s.pdf", suffix.Data(), names[ref].Data());
            TCanvas* c = new TCanvas("c","", 950, 700);
            c->SetLeftMargin(0.15);
            c->SetLogy();
            c->Print(pdfName + "[");
            for (int i = 0; i < nSteps; ++i) {
                TH1F* h = hists[ref][i];

                struct Peak {Double_t A, mu, sigma; Int_t id; std::vector<Int_t> merged_ids;Bool_t alwaysSeparate = false;};
                std::vector<Peak> seeds;
                Double_t sliceMin = pEdges[i];
                Double_t sliceMax = pEdges[i+1];
                Double_t pMid     = 0.5 * (sliceMin + sliceMax);

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
                        sigma0 = (resoHypAbs / resoRefAbs) * (1.0 / (bHyp * bHyp));
                        mu     = (bRef - bHyp) / (bHyp * bHyp * resoHypAbs);
                    }
                    //sigma0 = std::clamp(sigma0, 0.1, 5.0);
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
                std::cout << ">> Automatic peaks for ref=" << names[ref] << "  p in [" 
                << sliceMin << "," << sliceMax << "] GeV/c:\n";
                for (size_t j = 0; j < merged.size(); ++j) {
                    const auto &pk = merged[j];
                    std::cout << Form("   peak %zu (ids:", j);
                    for (auto id : pk.merged_ids) std::cout << " " << names[id];
                    std::cout << Form(" ) -> μ = %.3f, σ = %.3f\n", pk.mu, pk.sigma);
                }

                Double_t x_low  = xMin;
                Double_t x_high = xMax;

                Int_t manualNGauss = 0;
                std::vector<double> manualMeans, manualSigmas;

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
                    manualSigmas.resize(manualNGauss);
                    std::cout << "Enter " << manualNGauss << " sigma (width) guesses for each peak (space.separated): ";
                    for (int i = 0; i < manualNGauss; ++i)
                        std::cin >> manualSigmas[i];
                } 
                

                size_t nG;
                std::vector<double> means;
                if (manualNGauss > 0) {
                    nG = manualNGauss;
                    means = manualMeans;
                } 
                if (nG < 1) { c->Print(pdfName); delete h; continue; }
                
                Double_t fit_lo = x_low;
                Double_t fit_hi = x_high;

                std::ostringstream form;
                for (size_t i = 0; i < nG; ++i) {
                    if (i) form << "+"; 
                    form << "gaus(" << 3 * i << ")";
                }
                size_t idxLine = 3 * nG;
                form << "+pol0(" << idxLine << ")";
                TF1 *sum = new TF1("sum", form.str().c_str(), fit_lo, fit_hi);
                sum->SetParLimits(idxLine, 0, 10.0);
                sum->SetParameter(idxLine, 1.0);
                for (size_t i = 0; i < nG; ++i) {
                    if (manualNGauss > 0) {
                        Double_t mu0 = means[i];
                        Double_t sig0 = manualSigmas[i]; 
                        Int_t    bin = h->FindBin(mu0);        
                        Double_t amp = h->GetBinContent(bin);

                        sum->SetParLimits(3*i+0, 0.0, std::max(h->GetMaximum()*1.2, 1.05*amp));
                        sum->SetParameter(3*i+0, amp);

                        sum->SetParLimits(3*i+1, mu0 - muWindow, mu0 + muWindow);
                        sum->SetParameter(3*i+1, mu0);

                        sum->SetParLimits(3*i+2, sig0*0.5, sig0*2.0);
                        sum->SetParameter(3*i+2, sig0);

                        Double_t lo, hi;
                        sum->GetParLimits(3*i+2, lo, hi);
                        std::cout << Form("  Peak %zu σ‑limits=[%.3f,%.3f], seed=%.3f\n", i, lo, hi, sig0);
                    }
                
                }
                h->Fit(sum, "RQ0S", "", fit_lo, fit_hi);
                h->Draw("E1");
                Double_t yMax = 1.25 * h->GetMaximum();
                TLegend* leg = new TLegend(0, 0.10, 0.15, 0.30);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                sum->SetLineColor(kRed); 
                sum->SetLineWidth(2); 
                sum->SetNpx(500); 
                sum->Draw("same");
                for (Int_t i = 0; i < nG; ++i) {
                    Double_t mu = means[i];
                    TLine *l = new TLine(mu, 0, mu, yMax);
                    l->SetLineColor(colors[i % nParts]);
                    l->Draw();
                }
                

                for (Int_t i = 0; i < nG; ++i) {
                    if (sum->GetParameter(3*i) <= 0) continue;
                    TF1 *g = new TF1(Form("g_%d_%d", ref, i), "gaus", x_low, x_high);
                    g->SetParameters(&sum->GetParameters()[3*i]);
                    Double_t A   = sum->GetParameter(3*i + 0);
                    Double_t mu  = sum->GetParameter(3*i + 1);
                    Double_t sig = sum->GetParameter(3*i + 2);
                    std::cout << Form("A = %.3f, μ = %.3f, σ = %.3f\n", A, mu, sig);
                    g->SetLineColor(colors[i % nParts]);   
                    g->SetLineStyle(2);
                    g->SetNpx(500);
                    g->Draw("same");

                    TString label;
                    label = Form("Peak %d", i);
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

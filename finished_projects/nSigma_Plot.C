#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>
#include <sstream>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TF1.h"  
#include "TLine.h"
#include "TPaveText.h"
#include "TSystem.h"

#include <AddTrees.h>   // project-specific: file discovery/chain fill
#include <helpers.h>    // project-specific: masses, charges, colors, Bethe–Bloch

void nSigma_Plot(){
    auto help = new helper();
    const Int_t nParts = helper::nParts;
    const Int_t NtrkMax = help->NtrkMax;

    // histogram & slicing config
    const Int_t nBins = 500;
    const Double_t xMin = -10.0, xMax = 10.0;
    const Double_t pStart = 0.45, pEnd = 0.55; // momentum slicing (GeV/c)
    const Double_t step = 0.1; // slice width (GeV/c)

    // peak finding & fitting config
    const Double_t muWindow = 0.5; // ±window for μ during fit init
    const Double_t mergeDistanceFactor = 1.0; // merge seeds if |Δμ| < f·max(σ)
    const Double_t nEntriesLimit = 1e7;

    // event/track filters
    const Bool_t TOFfilter = false; // require TOF info
    const Bool_t KaExclusion = false; // veto |nσ_K^TOF|<3
    const Bool_t PrExclusion = false; // veto |nσ_p^TOF|<3

    // what to draw
    const Bool_t plotTPC = true;
    const Bool_t plotTOF = false;
    const Bool_t PeakZoom = false; // autoset x-range around predicted peaks

    // seed strategy
    const Bool_t manualPredictPeaks = true; // interactive seeding (means/sigmas/amps); else auto from predicted model
    const std::array<Bool_t, nParts> doPid = {{true, false, false, false, false}}; // which ref species to plot
    using PeakPars = std::array<Double_t,4>;

    // batch only in auto mode; interactive needs windows
    gROOT->SetBatch(!manualPredictPeaks);
    gStyle->SetOptStat(1);

    // input data chain
    TChain chain("twotauchain");
    AddTrees(chain, help->base_dir);

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);
    chain.SetBranchStatus("fTrkTOFexpMom", 1);
    for (Int_t i = 0; i < nParts; ++i){
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

    Long64_t nEntries = std::min(chain.GetEntries(), static_cast<Long64_t>(nEntriesLimit));

    // lambda to draw nσ histograms (TPC or TOF)
    auto drawNSigma = [&](Bool_t isTPCmode) {
        // momentum bins for this detector mode
        const Double_t pMin   = isTPCmode ? pStart : std::max(pStart, 0.4);
        const Int_t    nSteps = Int_t(std::floor((pEnd - pMin) / step + 0.5));
        std::vector<Double_t> pEdges(nSteps+1);
        for (Int_t i = 0; i <= nSteps; ++i) pEdges[i] = pMin + i * step;
        
        TString suffix = isTPCmode ? "TPC" : "TOF";

        // book per-slice histograms of nσ (one set per reference species)
        std::vector<std::vector<TH1D*>> hists(nParts, std::vector<TH1D*>(nSteps,nullptr));
        for (Int_t pid = 0; pid < nParts; ++pid) {
            if (!doPid[pid]) continue;
            for (Int_t i = 0; i < nSteps; ++i) {
                TString name1 = Form("n#sigma_{%s} %g < p < %g GeV/c (%s)", 
                        help->pCodes[pid], pEdges[i], pEdges[i+1], suffix.Data());
                TString name2 = Form("n#sigma_{%s} %g < p < %g GeV/c (%s); n#sigma_{%s}; Counts",
                        help->pCodes[pid], pEdges[i], pEdges[i+1], suffix.Data(), help->pCodes[pid]);
                hists[pid][i] = new TH1D(name1, name2, nBins, xMin, xMax);
                hists[pid][i]->Sumw2(true);
                hists[pid][i]->SetMarkerStyle(kFullCircle);
                hists[pid][i]->SetMarkerSize(0.75);
                hists[pid][i]->SetMarkerColor(kBlack);
                hists[pid][i]->SetLineColor(kBlack);
            }
        }

        // fill histograms
        for (Long64_t ev = 0; ev < nEntries; ++ev) {
            chain.GetEntry(ev);
            for (Int_t t = 0; t < NtrkMax; ++t) {
                // TOF validity: in TOF mode always require; in TPC mode only if TOFfilter=true
                if (tofExpMom[t] < 0 && (!isTPCmode || TOFfilter)) continue;
                if (KaExclusion && !TMath::IsNaN(tofNS[3][t]) && TMath::Abs(tofNS[3][t]) < 3.0) continue;
                if (PrExclusion && !TMath::IsNaN(tofNS[4][t]) && TMath::Abs(tofNS[4][t]) < 3.0) continue;
                
                Double_t pG = inner[t]; // assumed GeV/c
                Int_t bin = std::lower_bound(pEdges.begin(), pEdges.end(), pG) - pEdges.begin() - 1;
                if (bin < 0 || bin >= nSteps) continue;

                for (Int_t pid = 0; pid < nParts; ++pid) {
                    if (!doPid[pid]) continue;
                    Float_t val = isTPCmode ? tpcNS[pid][t] : tofNS[pid][t];
                    if (!TMath::IsNaN(val)) hists[pid][bin]->Fill(val);
                }
            }
        }  

        // for each reference species, loop momentum slices and fit multi-Gaussian + constant background
        for (Int_t ref = 0; ref < nParts; ++ref) {
            if (!doPid[ref]) continue;
            TString pdfName = Form("nSigma%s_%s.pdf", suffix.Data(), help->pCodes[ref]);

            TCanvas* c = new TCanvas("c","", 950, 700);
            c->SetLeftMargin(0.15);
            c->SetLogy();
            c->Print(pdfName + "[");
            
            for (Int_t i = 0; i < nSteps; ++i) {
                TH1D* h1 = hists[ref][i];

                // seed peaks from physics model (μ,σ) and histogram amplitude at μ
                struct Peak {Double_t A, mu, sigma; Int_t id; std::vector<Int_t> merged_ids;Bool_t alwaysSeparate = false;};
                std::vector<Peak> seeds;
                Double_t sliceMin = pEdges[i];
                Double_t sliceMax = pEdges[i+1];
                Double_t pMid     = 0.5 * (sliceMin + sliceMax);

                for (Int_t hyp = 0; hyp < nParts; ++hyp) {
                    Double_t hypMass = help->pMasses[hyp];
                    Double_t sigma0, mu;
                    if (isTPCmode) {
                        // TPC: μ = (dHyp/dRef -1)/resoRef ; σ ∝ (resoHyp/resoRef)*(dHyp/dRef)
                        Double_t dRef = help->getTPCSignal(pMid * 1000, help->pMasses[ref], 1.0);
                        Double_t dHyp = help->getTPCSignal(pMid * 1000, hypMass, 1.0);
                        if (dRef < 0 || dHyp < 0) continue;
                        Double_t resoHyp = help->resoTPC[hyp][ref];
                        Double_t resoRef = help->resoTPC[ref][hyp];
                        sigma0 = (resoHyp / resoRef) * (dHyp / dRef);
                        mu     = (dHyp/dRef - 1.0) / resoRef;
                    }
                    else {
                        // TOF: μ = (β_ref - β_hyp)/(β_hyp^2 * resoRef) ; σ ∝ (resoHyp/resoRef)*(1/β_hyp^2)
                        Double_t bRef = help->beta(help->pMasses[ref], pMid*1000);
                        Double_t bHyp = help->beta(help->pMasses[hyp], pMid*1000);
                        Double_t resoHyp = help->resoTOF[hyp][ref];
                        Double_t resoRef = help->resoTOF[ref][hyp];
                        sigma0 = (resoHyp / resoRef) * (1.0 / (bHyp * bHyp));
                        mu     = (bRef - bHyp) / (bHyp * bHyp * resoRef);
                    }

                    if (mu < xMin || mu > xMax) continue;
                    Int_t    bin = h1->FindBin(mu);        
                    Double_t amp = h1->GetBinContent(bin);
                    seeds.push_back({amp, mu, sigma0, hyp, {hyp}, false});
                }

                // optionally merge μ and π seeds into a single seed (keeps IDs in merged_ids)
                Peak muPiCombined{0,0,0,-1,{},false};
                Bool_t sawMu = false, sawPi = false;
                for (auto it = seeds.begin(); it != seeds.end();) {
                    if (it->id == 1 || it->id == 2) {
                        if (!sawMu && it->id == 1) {
                            muPiCombined = *it;
                            sawMu = true;
                        }
                        else {
                            // merge by conserving integral and mean/variance
                            Double_t I_old = muPiCombined.A * muPiCombined.sigma * TMath::Sqrt(2*TMath::Pi());
                            Double_t I_new = it->A * it->sigma * TMath::Sqrt(2*TMath::Pi());
                            Double_t I_tot = I_old + I_new;
                            Double_t mu_eff = (I_old*muPiCombined.mu + I_new*it->mu) / I_tot;
                            Double_t var_eff = (I_old*(muPiCombined.sigma*muPiCombined.sigma + (muPiCombined.mu-mu_eff)*(muPiCombined.mu-mu_eff)) + I_new*(it->sigma*it->sigma + (it->mu-mu_eff)*(it->mu-mu_eff))) / I_tot;
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
                    muPiCombined.id = 1; // // tag as "μ/π"
                    // never merge away the reference species peak
                    muPiCombined.alwaysSeparate = std::find(muPiCombined.merged_ids.begin(), muPiCombined.merged_ids.end(), ref) != muPiCombined.merged_ids.end();
                    seeds.push_back(muPiCombined);
                }

                for (auto &s : seeds) {
                    if (std::find(s.merged_ids.begin(), s.merged_ids.end(), ref) != s.merged_ids.end()) {
                        s.alwaysSeparate = true;
                    }
                }

                // merge neighboring seeds if close in μ (unless flagged alwaysSeparate)
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

                 // optional x-range zoom around predicted peaks
                Double_t x_low  = xMin, x_high = xMax;
                if (PeakZoom) {
                    const Double_t Nsigma = 8.0;
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
                    h1->GetXaxis()->SetRangeUser(x_low, x_high);
                }

                // seeding mode: manual (interactive) or automatic (from merged)
                Int_t manualNGauss = 0;
                std::vector<Double_t> manualMeans, manualSigmas, manualAmps;
                if (manualPredictPeaks) {
                    c->cd();
                    c->Clear();
                    h1->Draw("E1");
                    c->Modified();
                    c->Update();
                    c->RaiseWindow();
                    gSystem->ProcessEvents();
                    std::cout << "How many Gaussians? ";
                    std::cin >> manualNGauss;
                    if (manualNGauss > 0) {
                        manualMeans.resize(manualNGauss);
                        std::cout << "Enter " << manualNGauss << " mean positions (space-separated): ";
                        for (Int_t i = 0; i < manualNGauss; ++i) {
                            std::cin >> manualMeans[i];
                        }
                        manualSigmas.resize(manualNGauss);
                        std::cout << "Enter " << manualNGauss << " sigma (width) guesses for each peak (space-separated): ";
                        for (Int_t i = 0; i < manualNGauss; ++i){
                            std::cin >> manualSigmas[i];
                        }
                        manualAmps.resize(manualNGauss);
                        std::cout << "Enter " << manualNGauss << " amplitude guesses for each peak (space-separated): ";
                        for (Int_t i = 0; i < manualNGauss; ++i){
                            std::cin >> manualAmps[i];
                        }
                    } else {
                        std::cerr << "Invalid number of Gaussians; falling back to automatic mode." << std::endl;
                        manualNGauss = 0;
                    }
                }

                size_t nG;
                std::vector<Double_t> means;
                if (manualPredictPeaks && manualNGauss > 0) {
                    nG = manualNGauss;
                    means = manualMeans;
                } else {
                    nG = merged.size();
                    means.reserve(nG);
                    for (auto &pk : merged) means.push_back(pk.mu);
                }
                if (nG < 1) { c->Print(pdfName); delete h1; continue; }
                
                // build sum-of-Gaussians + constant background and set parameter limits
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
                    if (manualPredictPeaks && manualNGauss > 0) {
                        Double_t mu0 = means[i];
                        Double_t sig0 = manualSigmas[i]; 
                        Double_t amp = manualAmps[i];

                        sum->SetParLimits(3*i+0, 0.0, 1e10);
                        sum->SetParameter(3*i+0, amp);

                        sum->SetParLimits(3*i+1, mu0 - muWindow, mu0 + muWindow);
                        sum->SetParameter(3*i+1, mu0);

                        sum->SetParLimits(3*i+2, sig0*0.5, sig0*2.0);
                        sum->SetParameter(3*i+2, sig0);

                        Double_t lo, hi;
                        sum->GetParLimits(3*i+2, lo, hi);
                    }
                    else {
                        const auto &p = merged[i];

                        sum->SetParLimits(3*i, 0.0, 1e10);
                        sum->SetParameter(3*i, p.A);

                        Double_t dMu = std::max(muWindow, 0.1*std::abs(p.mu));
                        sum->SetParLimits(3*i+1, p.mu - dMu, p.mu + dMu);
                        sum->SetParameter(3*i+1, p.mu);

                        sum->SetParLimits(3*i+2, 0.5*p.sigma, 2.0*p.sigma);
                        sum->SetParameter(3*i+2, p.sigma);
                    }
                }

                // χ² fit (custom helper uses bin-by-bin errors, respects TF1 limits)
                help->FitHistogramByChi2(h1, sum, x_low, x_high);

                // draw: data, total fit, components, labels
                h1->Draw("E1");
                Double_t yMax1 = 1.25 * h1->GetMaximum();

                TLegend* leg = new TLegend(0, 0.10, 0.15, 0.30);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                sum->SetLineColor(kRed); 
                sum->SetLineWidth(2); 
                sum->SetNpx(500); 
                sum->Draw("same");

                // vertical lines at predicted/entered μ
                if (!manualPredictPeaks) {
                    for (const auto &pk : merged) {
                        Double_t mu = pk.mu;
                        TLine *l = new TLine(mu, 0, mu, yMax1);
                        Int_t col  = (pk.id >= 0) ? help->colors[pk.id] : kGray+2;
                        l->SetLineColor(col);
                        l->Draw();
                    }
                }
                else {
                    for (Int_t i = 0; i < nG; ++i) {
                        Double_t mu = means[i];
                        TLine *l = new TLine(mu, 0, mu, yMax1);
                        l->SetLineColor(help->colors[i % nParts]);
                        l->Draw();
                    }
                }
                
                // individual Gaussian components + legend labels
                for (Int_t i = 0; i < nG; ++i) {
                    if (sum->GetParameter(3*i) <= 0) continue;
                    TF1 *g = new TF1(Form("g_%d_%d", ref, i), "gaus", x_low, x_high);
                    g->SetParameters(
                        sum->GetParameter(3*i),
                        sum->GetParameter(3*i+1),
                        sum->GetParameter(3*i+2)
                    );
                    if (!manualPredictPeaks) {
                        Int_t col = (merged[i].id >= 0) ? help->colors[merged[i].id] : kGray+2;
                        g->SetLineColor(col);
                    } else {
                        g->SetLineColor(help->colors[i % nParts]);   
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
                            label = help->pCodes[ids[0]];
                        } else {
                            for (size_t j=0; j<ids.size(); ++j) {
                                if (j) label += " + ";
                                label += help->pCodes[ids[j]];
                            }
                        }
                    } else {
                    label = Form("Peak %d", i+1);
                    }
                    leg->AddEntry(g, label, "l");
                }

                // χ²/ndf box
                TPaveText *pt=new TPaveText(0.02,0.90,0.25,0.99,"NDC");
                pt->AddText(Form("#chi^{2}/NDF = %.2f", sum->GetChisquare()/sum->GetNDF()));
                pt->SetFillColorAlpha(0,0); 
                pt->Draw("same");
                leg->Draw();
                c->Print(pdfName);
                delete sum;
                delete leg;
            }
            c->Print(pdfName+"]"); 
        }
    };

    // run for selected detectors
    if (plotTPC) drawNSigma(true);     
    if (plotTOF) drawNSigma(false);  
}

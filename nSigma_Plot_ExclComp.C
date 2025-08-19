#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>
#include <utility>
#include <fstream>

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
#include "TH1D.h"

#include <AddTrees.h>
#include <helpers.h>
#include <covarianceMatrix.h>

struct ndJsonLogger {
    std::ofstream f;

    void open(const std::string& path) {
        f.open(path, std::ios::out | std::ios::app);
        f.imbue(std::locale::classic());
        f << std::fixed << std::setprecision(4);
    }

    void write_config(Int_t nBins, Double_t xMin, Double_t xMax,
                      Double_t pStart, Double_t pEnd, Double_t step,
                      Double_t muWindow, Double_t mergeDistanceFactor,
                      Double_t nEntriesLimit, Bool_t FitKaonExclComp,
                      Bool_t FitProtonExclComp, Bool_t plotTPC, Bool_t plotTOF, Bool_t PeakZoom, Bool_t manualPredictPeaks)
    {
        f << "{\n"
          << "\"type\":\"config\",\n"
          << "\"nBins\":" << nBins << ",\n"
          << "\"xMin\":" << xMin << ",\n"
          << "\"xMax\":" << xMax << ",\n"
          << "\"pStart\":" << pStart << ",\n"
          << "\"pEnd\":" << pEnd << ",\n"
          << "\"step\":" << step << ",\n"
          << "\"muWindow\":" << muWindow << ",\n"
          << "\"mergeDistanceFactor\":" << mergeDistanceFactor << ",\n"
          << "\"nEntriesLimit\":" << nEntriesLimit << ",\n"
          << "\"FitKaonExclComp\":" << (FitKaonExclComp?"true":"false") << ",\n"
          << "\"FitProtonExclComp\":" << (FitProtonExclComp?"true":"false") << ",\n"
          << "\"plotTPC\":" << (plotTPC?"true":"false") << ",\n"
          << "\"plotTOF\":" << (plotTOF?"true":"false") << ",\n"
          << "\"PeakZoom\":" << (PeakZoom?"true":"false") << ",\n"
          << "\"manualPredictPeaks\":" << (manualPredictPeaks?"true":"false")
          << "\n}\n";
    }

    static void dump_vec(std::ostream& os, const std::vector<Double_t>& v) {
        os << "[";
        for (size_t i=0;i<v.size();++i){ if(i) os << ","; os << v[i]; }
        os << "]";
    }

    static inline void indent(std::ostream& os, int n) {
        for (int i = 0; i < n; ++i) os.put(' ');
    }

    static void dump_row_inline(std::ostream& os, const std::vector<Double_t>& v) {
        os << "[";
        for (size_t j = 0; j < v.size(); ++j) {
            if (j) os << ", ";
            if (std::isfinite(v[j])) os << v[j]; else os << "null";
        }
        os << "]";
    }

    static void dump_mat(std::ostream& os,
                                const std::vector<std::vector<Double_t>>& m,
                                int indent_level = 0,
                                int indent_width = 2) {
        os << "[\n";
        for (size_t i = 0; i < m.size(); ++i) {
            indent(os, indent_level + indent_width);
            dump_row_inline(os, m[i]);
            if (i + 1 < m.size()) os << ",";
            os << "\n";
        }
        indent(os, indent_level);
        os << "]";
    }

    void write_slice(const Char_t* mode,        
                     const Char_t* tag,         
                     Double_t sigmaExcl,
                     const Char_t* refCode,
                     Int_t slice,
                     Double_t pmin, Double_t pmax,
                     Double_t kLeft, Double_t kRight,
                     const std::vector<std::vector<Double_t>>& frac,
                     const std::vector<Double_t>& totCont,
                     const std::vector<Double_t>& area_noExcl,
                     const std::vector<Double_t>& err_noExcl,
                     const std::vector<Double_t>& area_wExcl,
                     const std::vector<Double_t>& err_wExcl,
                     Double_t D_over_N, Double_t chi2_over_ndf,
                     Bool_t manualPredictPeaks,
                     const std::vector<Double_t>& manualMeans,
                     const std::vector<Double_t>& manualSigmas,
                     const std::vector<Double_t>& manualAmps)
    {
        f << "{\n"
          << "\"tag\":\"" << tag << "\","
          << "\"sigmaExcl\":" << sigmaExcl << ",\n"
          << "\"ref\":\"" << refCode << "\",\n"
          << "\"slice\":" << slice << ",\n"
          << "\"pmin\":" << pmin << ",\n"
          << "\"pmax\":" << pmax << ",\n"
          << "\"band_window\":{\n"
              << "\"kLeft\":" << kLeft << ",\n"
              << "\"kRight\":" << kRight
          << "\n},\n";

        f << "\"frac\":";     dump_mat(f, frac);    f << ",\n";
        f << "\"totCont\":";  dump_vec(f, totCont); f << ",\n";

        f << "\"peak_areas\":{\n";
        f << "\"area_noExcl\":"; dump_vec(f, area_noExcl); f << ",\n";
        f << "\"err_noExcl\":";  dump_vec(f, err_noExcl);  f << ",\n";
        f << "\"area_wExcl\":";  dump_vec(f, area_wExcl);  f << ",\n";
        f << "\"err_wExcl\":";   dump_vec(f, err_wExcl);
        f << "\n},\n";

        f << "\"fit_quality\":{\n"
          << "\"D_over_N\":" << D_over_N << ",\n"
          << "\"chi2_over_ndf\":" << chi2_over_ndf
          << "\n},\n";

        f << "\"peak_seeding\":{\n";
        f << "\"mode\":\"" << (manualPredictPeaks ? "manual" : "auto") << ",\"";
        if (manualPredictPeaks) {
            f << "\"means\":";  dump_vec(f, manualMeans);  f << ",\n"
              << "\"sigmas\":"; dump_vec(f, manualSigmas); f << ",\n"
              << "\"amps\":";   dump_vec(f, manualAmps);
        }
        f << "\n}\n";

        f << "}\n";
    }
};

void nSigma_Plot_ExclComp(){
    auto help = new helper();
    const Int_t nParts = helper::nParts;
    const Int_t NtrkMax = help->NtrkMax;
    const Int_t   nBins   = 500;
    const Double_t xMin   = -10.0, xMax = 10.0;
    const Double_t pStart = 0.45, pEnd = 0.55, step = 0.1;
    const Double_t muWindow = 0.5;
    const Double_t mergeDistanceFactor = 1.0;
    const Double_t nEntriesLimit = 1e7;
    const Double_t eigenThr = 1e-4; 
    const Bool_t FitKaonExclComp = true;
    const Bool_t FitProtonExclComp = false;
    const Bool_t plotTPC = true;
    const Bool_t plotTOF = false;
    const Bool_t PeakZoom = false;
    const Bool_t manualPredictPeaks = true;
    const std::array<Bool_t, nParts> doPid = {{true, false, false, false, false}};
    using PeakPars = std::array<Double_t,4>;
    const std::vector<Double_t> sigmaExclList = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    ndJsonLogger ndlog;
    ndlog.open(Form("nSigmaEl-Plot-Data-%.2f<p%.2f.ndjson", pStart, pEnd));
    ndlog.write_config(nBins, xMin, xMax, pStart, pEnd, step, muWindow, mergeDistanceFactor, nEntriesLimit, FitKaonExclComp, FitProtonExclComp, plotTPC, plotTOF, PeakZoom, manualPredictPeaks);

    gROOT->SetBatch(!manualPredictPeaks);
    gStyle->SetOptStat(1);

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

    auto drawNSigma = [&](Bool_t isTPCmode, Double_t sigmaExcl) {
        const Double_t pMin   = isTPCmode ? pStart : std::max(pStart, 0.4);
        const Int_t    nSteps = Int_t(std::floor((pEnd - pMin) / step + 0.5));
        std::vector<Double_t> pEdges(nSteps+1);
        for (Int_t i = 0; i <= nSteps; ++i) pEdges[i] = pMin + i * step;
        TString suffix = isTPCmode ? "TPC" : "TOF";    
        enum ExclType{KaonExcl = 3, ProtExcl = 4};
        std::vector<std::pair<ExclType, Bool_t>> modes = {
        {KaonExcl,  FitKaonExclComp},
        {ProtExcl,  FitProtonExclComp}
        };
        for (auto [excl, doThis]: modes) {
            if(!doThis) continue;
            const Char_t* name = (excl==KaonExcl ? "KaonExcl" : "ProtonExcl");
            std::vector<std::vector<TH1F*>> h_no(nParts, std::vector<TH1F*>(nSteps,nullptr));
            std::vector<std::vector<TH1F*>> h_w(nParts, std::vector<TH1F*>(nSteps,nullptr));

            std::vector<std::vector<TH1D*>> hcm_no(nParts, std::vector<TH1D*>(nSteps,nullptr));
            std::vector<std::vector<TH1D*>> hcm_w (nParts, std::vector<TH1D*>(nSteps,nullptr));
            std::vector<std::vector<covarianceMatrix*>> cmObj(nParts, std::vector<covarianceMatrix*>(nSteps,nullptr));

        
            for (Int_t pid=0; pid<nParts; ++pid) if (doPid[pid]) {
                for (Int_t i=0; i<nSteps; ++i) {
                    auto name_no = Form("n#sigma_%s %g < p < %g GeV/c (%s-no%s)", help->pCodes[pid], pEdges[i], pEdges[i+1], suffix.Data(), name);
                    auto name_w = Form("n#sigma_%s %g < p < %g GeV/c (%s-%.2f-%s)", help->pCodes[pid], pEdges[i], pEdges[i+1], suffix.Data(), sigmaExcl, name);
                    auto title_no = Form("n#sigma_{%s} %g < p < %g GeV/c (%s-no%s); n#sigma_{%s}; Counts", help->pCodes[pid], pEdges[i], pEdges[i+1], suffix.Data(), name, help->pCodes[pid]);
                    auto title_w = Form("n#sigma_{%s} %g < p < %g GeV/c (%s-%.2f-%s); n#sigma_{%s}; Counts", help->pCodes[pid], pEdges[i], pEdges[i+1], suffix.Data(), sigmaExcl, name, help->pCodes[pid]);
                    h_no[pid][i] = new TH1F(name_no, title_no, nBins, xMin, xMax);
                    h_w [pid][i] = new TH1F(name_w, title_w, nBins, xMin, xMax);
                    h_no[pid][i]->SetMarkerStyle(kFullCircle);
                    h_no[pid][i]->SetMarkerSize(0.75);
                    h_no[pid][i]->SetMarkerColor(kBlack);
                    h_no[pid][i]->SetLineColor(kBlack);
                    h_w[pid][i]->SetMarkerStyle(kFullCircle);
                    h_w[pid][i]->SetMarkerSize(0.75);
                    h_w[pid][i]->SetMarkerColor(kBlack);
                    h_w[pid][i]->SetLineColor(kBlack);

                    
                    auto name_cm_no = Form("cm_nSTPC_%s_%g_%g_no%s", help->pCodes[pid], pEdges[i], pEdges[i+1], name);
                    auto name_cm_w  = Form("cm_nSTPC_%s_%g_%g_%s%.2f", help->pCodes[pid], pEdges[i], pEdges[i+1], name, sigmaExcl);

                    hcm_no[pid][i] = new TH1D(name_cm_no, h_no[pid][i]->GetTitle(), nBins, xMin, xMax);
                    hcm_w [pid][i] = new TH1D(name_cm_w , h_w [pid][i]->GetTitle(), nBins, xMin, xMax);

                    std::vector<TH1D*> cmHists{ hcm_no[pid][i], hcm_w[pid][i] };
                    std::vector<std::tuple<Double_t, Double_t>*> fitlims;
                    fitlims.push_back(new std::tuple<Double_t,Double_t>(xMin, xMax));   
                    fitlims.push_back(new std::tuple<Double_t,Double_t>(xMin, xMax)); 
                    cmObj[pid][i] = new covarianceMatrix(cmHists, fitlims);
                }
            }       

            for (Long64_t ev = 0; ev < nEntries; ++ev) {
                chain.GetEntry(ev);
                for (Int_t t = 0; t < NtrkMax; ++t) {
                    Double_t pG = inner[t];
                    if (tofExpMom[t] < 0) continue;
                    Int_t bin = std::lower_bound(pEdges.begin(), pEdges.end(), pG)
                            - pEdges.begin() - 1;
                    if (bin<0 || bin>=nSteps) continue;

                    for (Int_t pid=0; pid<nParts; ++pid) if (doPid[pid]) {
                        Float_t val = isTPCmode ? tpcNS[pid][t] : tofNS[pid][t];
                        if (TMath::IsNaN(val)) continue;

                        h_no[pid][bin]->Fill(val);

                        if (TMath::IsNaN(tofNS[excl][t]) || TMath::Abs(tofNS[excl][t])>=sigmaExcl) h_w[pid][bin]->Fill(val);

                        auto nsVals = new std::vector<Double_t>(); 
                        nsVals->reserve(2);
                        auto toadd  = new std::vector<Bool_t>();   
                        toadd->reserve(2);

                        nsVals->clear(); 
                        toadd->clear();
                        nsVals->push_back(val);  
                        toadd->push_back(true);
                        const Bool_t passExcl = (TMath::IsNaN(tofNS[excl][t]) || TMath::Abs(tofNS[excl][t]) >= sigmaExcl);
                        if (passExcl) {
                            nsVals->push_back(val); 
                            toadd->push_back(true);
                        }
                        else{
                            nsVals->push_back(-999.0); 
                            toadd->push_back(false);
                        }

                        cmObj[pid][bin]->addEvent(nsVals, toadd);
                    }
                }
            }
            
        
            for (Int_t ref = 0; ref < nParts; ++ref) {
                if (!doPid[ref]) continue;
                TString pdfName = Form("nSigma%s_%s_%sComp-%.2f.pdf", suffix.Data(), help->pCodes[ref], name, sigmaExcl);
                TCanvas* c = new TCanvas("c","", 950, 700);
                c->SetLeftMargin(0.15);
                c->SetLogy();
                c->Print(pdfName + "[");

                TString cmPdf = Form("CM_%s_%s_%s-%.2f.pdf",suffix.Data(), help->pCodes[ref], name, sigmaExcl);
                TCanvas* cCM = new TCanvas(Form("cCM_%s_%s", suffix.Data(), help->pCodes[ref]), "", 800, 700);
                cCM->SetLeftMargin(0.15);
                cCM->SetRightMargin(0.15);
                cCM->Print(cmPdf + "[");
                for (Int_t i = 0; i < nSteps; ++i) {
                    auto cm = cmObj[ref][i];
                    const Bool_t useOffDiag = kTRUE;
                    cm->make(useOffDiag, eigenThr);
                    cCM->cd();
                    cCM->Clear();
                    cm->plot(cCM, cmPdf);
                    c->cd();

                    TH1F* h1 = nullptr;
                    TH1F* h2 = nullptr;
                    h1 = h_no[ref][i];
                    h2 = h_w[ref][i];
                    struct Peak {Double_t A, mu, sigma; Int_t id; std::vector<Int_t> merged_ids;Bool_t alwaysSeparate = false;};
                    std::vector<Peak> seeds;
                    Double_t sliceMin = pEdges[i];
                    Double_t sliceMax = pEdges[i+1];
                    Double_t pMid     = 0.5 * (sliceMin + sliceMax);

                    for (Int_t hyp = 0; hyp < nParts; ++hyp) {
                        Double_t hypMass = help->pMasses[hyp];
                        Double_t sigma0, mu;
                        if (isTPCmode) {
                            Double_t dRef = help->getTPCSignal(pMid * 1000, help->pMasses[ref], 1.0);
                            Double_t dHyp = help->getTPCSignal(pMid * 1000, hypMass, 1.0);
                            if (dRef < 0 || dHyp < 0) 
                                continue;

                            Double_t resoHyp = help->resoTPC[hyp][ref];
                            Double_t resoRef = help->resoTPC[ref][hyp];
                            sigma0 = (resoHyp / resoRef) * (dHyp / dRef);
                            mu     = (dHyp/dRef - 1.0) / resoRef;
                        }
                        else {
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

                    Double_t x_low  = xMin;
                    Double_t x_high = xMax;
                    
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
                    
                    Double_t fit_lo = x_low;
                    Double_t fit_hi = x_high;

                    std::ostringstream form;
                    for (size_t i = 0; i < nG; ++i) {
                        if (i) form << "+";
                        form << "gaus(" << 3*i << ")+"
                        << "gaus(" << 3*(nG+i) << ")";
                    }
                    form << "+pol0(" << 4*nG << ")"  // background h1
                        << "+pol0(" << 4*nG+1 << ")"; // background h2

                    TF1* sum = new TF1("sum", form.str().c_str(), fit_lo, fit_hi);
                    sum->SetNpx(500); 
                    const Int_t offA1 = 0;
                    const Int_t offA2 = nG;
                    const Int_t offM  = 2*nG;
                    const Int_t offS  = 3*nG;
                    const Int_t offP1 = 4*nG;
                    const Int_t offP2 = 4*nG + 1;
                    sum->SetParLimits(offP1, 0, 10.0);
                    sum->SetParameter(offP1, 1.0);
                    sum->SetParLimits(offP2, 0, 10.0);
                    sum->SetParameter(offP2, 1.0);
                    for (size_t i = 0; i < nG; ++i) {
                        if (manualPredictPeaks && manualNGauss > 0) {
                            Double_t mu0 = means[i];
                            Double_t sig0 = manualSigmas[i]; 
                            Double_t amp = manualAmps[i]; 

                            sum->SetParLimits(offA1 + i, 0.0, std::max(h1->GetMaximum()*1.2, 1.05*amp));
                            sum->SetParameter(offA1 + i, amp);

                            sum->SetParLimits(offA2 + i, 0, std::max(h1->GetMaximum()*1.2, 1.05*amp));
                            sum->SetParameter (offA2 + i, amp);

                            sum->SetParLimits(offM + i, mu0 - muWindow, mu0 + muWindow);
                            sum->SetParameter(offM + i, mu0);

                            sum->SetParLimits(offS + i, sig0*0.5, sig0*2.0);
                            sum->SetParameter(offS + i, sig0);

                            Double_t lo, hi;
                            sum->GetParLimits(offS, lo, hi);
                        }
                        else {
                            const auto &p = merged[i];

                            sum->SetParLimits(offA1 + i, 0.0, std::max(h1->GetMaximum()*1.2, 1.05*p.A));
                            sum->SetParameter(offA1 + i, p.A);

                            sum->SetParLimits(offA2 + i, 0, std::max(h1->GetMaximum()*1.2, 1.05*p.A));
                            sum->SetParameter (offA2 + i, p.A);

                            sum->SetParLimits(offM + i, p.mu - muWindow, p.mu + muWindow);
                            sum->SetParameter(offM + i, p.mu);

                            sum->SetParLimits(offS + i, p.sigma*0.5, p.sigma*2.0);
                            sum->SetParameter(offS + i, p.sigma);
                        }
                    }
                    help->FitHistogramsByChi2(h1, h2, sum, nG, fit_lo, fit_hi, cmObj[ref][i], std::vector<TH1D*>{hcm_no[ref][i], hcm_w[ref][i]}, eigenThr, nBins);
                    c->Clear();
                    h1->Draw("E1");
                    Double_t yMax1 = 1.25 * h1->GetMaximum();
                    TLegend* leg1 = new TLegend(0, 0.10, 0.15, 0.30);
                    leg1->SetBorderSize(0);
                    leg1->SetFillStyle(0);
                    std::vector<Double_t> par(sum->GetNpar());
                    std::vector<Double_t> err(sum->GetNpar());
                    for(Int_t i=0; i<sum->GetNpar(); ++i)  {
                        par[i] = sum->GetParameter(i);
                        err[i] = sum->GetParError(i);
                    }                    

                    struct BandStats {
                        Double_t kLeft, kRight;
                        std::vector<std::vector<Double_t>> frac;
                        std::vector<Double_t> totCont;
                    };

                    auto computeBandFractions = [&](Int_t whichHist) -> BandStats {
                        const Int_t offA = (whichHist==0 ? offA1 : offA2);
                        const Double_t kLeft = 3.0;
                        const Double_t kRight = 3.0;

                        std::vector<std::vector<Double_t>> frac(nG, std::vector<Double_t>(nG, std::numeric_limits<Double_t>::quiet_NaN()));
                        std::vector<Double_t> totCont(nG, std::numeric_limits<Double_t>::quiet_NaN());

                        auto integrate_sum_others = [&](Int_t ii, Double_t a, Double_t b) -> Double_t {
                            const int N = 20480;
                            const Double_t dx = (b - a) / (N - 1);
                            Double_t acc = 0.0;
                            for (int t = 0; t < N; ++t) {
                                const Double_t x = a + t*dx;
                                const Double_t w = (t==0 || t==N-1) ? 0.5 : 1.0;
                                Double_t y = 0.0;
                                for (Int_t j = 0; j < nG; ++j) if (j != ii) {
                                    const Double_t Aj   = par[offA + j];
                                    const Double_t muj  = par[offM + j];
                                    const Double_t sigj = par[offS + j];
                                    y += Aj * TMath::Gaus(x, muj, sigj, kFALSE);
                                }
                                acc += w * y;
                            }
                            return acc * dx;
                        };

                        for (Int_t ii = 0; ii < nG; ++ii) {
                            const Double_t Ai   = par[offA + ii];
                            const Double_t mui  = par[offM + ii];
                            const Double_t sigi = par[offS + ii];
                            const Double_t a = mui - kLeft;
                            const Double_t b = mui + kRight;

                            const Double_t are_i = help->GaussIntegral(Ai, mui, sigi, a, b);

                            for (Int_t j = 0; j < nG; ++j) {
                                const Double_t Aj   = par[offA + j];
                                const Double_t muj  = par[offM + j];
                                const Double_t sigj = par[offS + j];
                                const Double_t num_j = help->GaussIntegral(Aj, muj, sigj, a, b);
                                frac[ii][j] = (are_i > 0.0) ? (num_j / are_i) : std::numeric_limits<Double_t>::quiet_NaN();
                            }

                            if (are_i > 0.0) {
                                const Double_t others_sum_int = integrate_sum_others(ii, a, b);
                                totCont[ii] = others_sum_int / are_i;
                            }
                        }

                        return BandStats{kLeft, kRight, std::move(frac), std::move(totCont)};
                    };

                    auto [D1,N1] = help->PoissonDeviance(h1, 0, par, nG, offA1, offA2, offM, offS, offP1, offP2);
                    auto [D2,N2] = help->PoissonDeviance(h2, 1, par, nG, offA1, offA2, offM, offS, offP1, offP2);
                    const Double_t D  = D1 + D2;
                    const Double_t N  = N1 + N2;
                    const Int_t    k  = 4*nG + 2; 

                    const Double_t chi2 = sum->GetChisquare();
                    const Int_t    ndf  = sum->GetNDF();
                    const Double_t pval = TMath::Prob(chi2, ndf);

                    std::vector<Double_t> areas_noK(nG), areas_wK(nG), err_areas_noK(nG), err_areas_wK(nG);
                    const Double_t factor = TMath::Sqrt(2*TMath::Pi());
                    for (Int_t ig = 0; ig < nG; ++ig) {
                        Double_t A1   = par[offA1 + ig];
                        Double_t A2   = par[offA2 + ig];
                        Double_t sig  = par[offS  + ig];
                        Double_t eA1 = err[offA1 + ig];
                        Double_t eA2 = err[offA2 + ig];
                        Double_t eS  = err[offS  + ig];

                        areas_noK[ig] = A1  * sig * factor;
                        areas_wK [ig] = A2  * sig * factor;

                        Double_t dFdA = sig * factor;
                        Double_t dFdS = A1 * factor;        
                        Double_t var_noK = dFdA*dFdA * eA1*eA1+ dFdS*dFdS * eS*eS; 
                        err_areas_noK[ig] = TMath::Sqrt(var_noK);

                        dFdS = A2 * factor;
                        Double_t var_wK = dFdA*dFdA * eA2*eA2 + dFdS*dFdS * eS*eS;
                        err_areas_wK[ig] = TMath::Sqrt(var_wK);
                    }
                    const Double_t D_over_N      = (N   > 0 ? D   / N   : std::numeric_limits<Double_t>::quiet_NaN());
                    const Double_t chi2_over_ndf = (ndf > 0 ? chi2/ ndf : std::numeric_limits<Double_t>::quiet_NaN());


                    auto band_no = computeBandFractions(0); 
                    auto band_w  = computeBandFractions(1); 
                    ndlog.write_slice(
                        suffix.Data(), "noExcl", sigmaExcl, help->pCodes[ref],
                        i, pEdges[i], pEdges[i+1],
                        band_no.kLeft, band_no.kRight, band_no.frac, band_no.totCont,
                        areas_noK, err_areas_noK,  
                        areas_wK,  err_areas_wK,  
                        D_over_N, chi2_over_ndf,
                        manualPredictPeaks, manualMeans, manualSigmas, manualAmps
                    );

                    ndlog.write_slice(
                        suffix.Data(), "wExcl", sigmaExcl, help->pCodes[ref],
                        i, pEdges[i], pEdges[i+1],
                        band_w.kLeft, band_w.kRight, band_w.frac, band_w.totCont,
                        areas_noK, err_areas_noK,
                        areas_wK,  err_areas_wK,
                        D_over_N, chi2_over_ndf,
                        manualPredictPeaks, manualMeans, manualSigmas, manualAmps
                    );
                    const Int_t nPoints = 500; 
                    Double_t xlo = fit_lo, xhi = fit_hi;
                    Double_t dx  = (xhi - xlo)/(nPoints-1);
                    std::vector<Double_t> xv(nPoints), y1(nPoints), y2(nPoints);

                    for(Int_t ip=0; ip<nPoints; ++ip) {
                        Double_t x = xlo + ip*dx;
                        xv[ip] = x;
                        Double_t yy1 = par[offP1];
                        Double_t yy2 = par[offP2];
                        for(Int_t ig=0; ig<nG; ++ig) {
                            Double_t A1  = par[offA1 + ig];
                            Double_t A2  = par[offA2 + ig];
                            Double_t mu  = par[offM  + ig];
                            Double_t sig = par[offS  + ig];
                            yy1 += A1 * TMath::Gaus(x, mu, sig, kFALSE);
                            yy2 += A2 * TMath::Gaus(x, mu, sig, kFALSE);
                        }
                        y1[ip] = yy1;
                        y2[ip] = yy2;
                    }

                    auto gSum1 = new TGraph(nPoints, xv.data(), y1.data());
                    gSum1->SetLineColor(kRed);
                    gSum1->SetLineWidth(2);
                    gSum1->Draw("L SAME");

                    auto gSum2 = new TGraph(nPoints, xv.data(), y2.data());
                    gSum2->SetLineColor(kRed);
                    gSum2->SetLineWidth(2);

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
                        
                    for (Int_t i = 0; i < nG; ++i) {
                        if (sum->GetParameter(offA1 + i) < 0) continue;
                        TF1 *g1 = new TF1(Form("g_%d_%d", ref, i), "gaus", x_low, x_high);
                        g1->SetParameters(
                            sum->GetParameter(offA1 + i),
                            sum->GetParameter(offM + i),
                            sum->GetParameter(offS + i)
                        );
                        if (!manualPredictPeaks) {
                            Int_t col = (merged[i].id >= 0) ? help->colors[merged[i].id] : kGray+2;
                            g1->SetLineColor(col);
                        } else {
                            g1->SetLineColor(help->colors[i % nParts]);   
                        }
                        g1->SetLineStyle(2);
                        g1->SetNpx(500);
                        g1->Draw("same");

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
                        leg1->AddEntry(g1, label, "l");
                    }
                    TPaveText *pt1=new TPaveText(0.02,0.90,0.15,0.99,"NDC");
                    pt1->AddText(Form("#chi^{2}/NDF = %.2f", sum->GetChisquare()/sum->GetNDF()));
                    pt1->SetFillColorAlpha(0,0); 
                    pt1->Draw("same");
                    leg1->Draw();
                    c->Print(pdfName);
                    delete leg1;
                    delete pt1;

                    c->Clear();
                    h2->Draw("E1");
                    gSum2->Draw("L SAME");
                    TLegend* leg2 = new TLegend(0, 0.10, 0.15, 0.30);
                    leg2->SetBorderSize(0);
                    leg2->SetFillStyle(0);
                    Double_t yMax2 = 1.25 * h2->GetMaximum();
                    if (!manualPredictPeaks) {
                        for (const auto &pk : merged) {
                            Double_t mu = pk.mu;
                            TLine *l = new TLine(mu, 0, mu, yMax2);
                            Int_t col  = (pk.id >= 0) ? help->colors[pk.id] : kGray+2;
                            l->SetLineColor(col);
                            l->Draw();
                        }
                    }
                    else {
                        for (Int_t i = 0; i < nG; ++i) {
                            Double_t mu = means[i];
                            TLine *l = new TLine(mu, 0, mu, yMax2);
                            l->SetLineColor(help->colors[i % nParts]);
                            l->Draw();
                        }
                    }

                    for (Int_t i = 0; i < nG; ++i) {
                        if (sum->GetParameter(offA2 + i) < 0) continue;
                        TF1 *g2 = new TF1(Form("g_%d_%d", ref, i), "gaus", x_low, x_high);
                        g2->SetParameters(
                            sum->GetParameter(offA2 + i),
                            sum->GetParameter(offM + i),
                            sum->GetParameter(offS + i)
                        );
                        if (!manualPredictPeaks) {
                            Int_t col = (merged[i].id >= 0) ? help->colors[merged[i].id] : kGray+2;
                            g2->SetLineColor(col);
                        } else {
                            g2->SetLineColor(help->colors[i % nParts]);   
                        }
                        g2->SetLineStyle(2);
                        g2->SetNpx(500);
                        g2->Draw("same");

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
                        leg2->AddEntry(g2, label, "l");
                    } 
                    TPaveText *pt2=new TPaveText(0.02,0.90,0.15,0.99,"NDC");
                    pt2->AddText(Form("#chi^{2}/NDF = %.2f", sum->GetChisquare()/sum->GetNDF()));
                    pt2->SetFillColorAlpha(0,0); 
                    pt2->Draw("same");
                    leg2->Draw();
                    c->Print(pdfName);
                    delete pt2;
                    delete leg2;
                    delete sum;
                }    
            c->Print(pdfName+"]"); 
            cCM->Print(cmPdf + "]");
            }
        }
    };
    for (Double_t s : sigmaExclList) {
        if (plotTPC) drawNSigma(true,  s);
        if (plotTOF) drawNSigma(false, s);
    } 
}

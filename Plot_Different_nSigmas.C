#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>
#include <utility>
#include <fstream>
#include <string>
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
#include "TH1D.h"

#include <AddTrees.h>
#include <helpers.h>
#include <covarianceMatrix.h>

enum PID {kEl=0, kMu=1, kPi=2, kKa=3, kPr=4};

struct ndJsonLogger {
    std::ofstream f;

    void open(const std::string& path) {
        f.open(path, std::ios::out | std::ios::app);
        f.imbue(std::locale::classic());
        f << std::fixed << std::setprecision(4);
    }

    void write_config(Int_t nBins, Double_t xMin1, Double_t xMax1, Double_t xMin2, Double_t xMax2,
                      Double_t pStart, Double_t pEnd, Double_t step,
                      Double_t muWindow, Double_t mergeDistanceFactor,
                      Double_t nEntries, Bool_t FitKaonExclComp,
                      Bool_t FitProtonExclComp, Bool_t plotTPC, Bool_t plotTOF, Bool_t PeakZoom, Bool_t manualPredictPeaks, Double_t eigenThr, Bool_t useOffDiag)
    {
        f << "{\n"
          << "\"type\":\"config\",\n"
          << "\"nBins\":" << nBins << ",\n"
          << "\"xMin1\":" << xMin1 << ",\n"
          << "\"xMax1\":" << xMax1 << ",\n"
          << "\"xMin2\":" << xMin2 << ",\n"
          << "\"xMax2\":" << xMax2 << ",\n"
          << "\"pStart\":" << pStart << ",\n"
          << "\"pEnd\":" << pEnd << ",\n"
          << "\"step\":" << step << ",\n"
          << "\"muWindow\":" << muWindow << ",\n"
          << "\"mergeDistanceFactor\":" << mergeDistanceFactor << ",\n"
          << "\"nEntries\":" << nEntries << ",\n"
          << "\"FitKaonExclComp\":" << (FitKaonExclComp?"true":"false") << ",\n"
          << "\"FitProtonExclComp\":" << (FitProtonExclComp?"true":"false") << ",\n"
          << "\"plotTPC\":" << (plotTPC?"true":"false") << ",\n"
          << "\"plotTOF\":" << (plotTOF?"true":"false") << ",\n"
          << "\"PeakZoom\":" << (PeakZoom?"true":"false") << ",\n"
          << "\"manualPredictPeaks\":" << (manualPredictPeaks?"true":"false")<< ",\n" 
          << "\"eigenThr\":" << eigenThr << ",\n"                
          << "\"useOffDiag\":" << (useOffDiag?"true":"false")
          << "\n}\n\n";
    }

    static void dump_vec(std::ostream& os, const std::vector<Double_t>& v) {
        os << "[";
        for (size_t i=0;i<v.size();++i){ if(i) os << ","; os << v[i]; }
        os << "]";
    }

    static inline void indent(std::ostream& os, Int_t n) {
        for (Int_t i = 0; i < n; ++i) os.put(' ');
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
                                Int_t indent_level = 0,
                                Int_t indent_width = 2) {
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
                     const std::vector<Double_t>& totContErr,
                     const std::vector<Double_t>& area,
                     const std::vector<Double_t>& area_err,
                     Double_t D_over_N, Double_t chi2_over_ndf,
                     Bool_t manualPredictPeaks,
                     const std::vector<Double_t>& seedMeans,
                     const std::vector<Double_t>& seedSigmas,
                     const std::vector<Double_t>& seedAmps,
                     const std::vector<Double_t>& fitMeans,
                     const std::vector<Double_t>& fitSigmas,
                     const std::vector<Double_t>& fitAmps,
                     const std::vector<Double_t>& err_fitMeans,
                     const std::vector<Double_t>& err_fitSigmas,
                     const std::vector<Double_t>& err_fitAmps,
                     Double_t bg1, Double_t err_bg1)
    {
        f << "{\n"
          << "\"tag\":\"" << tag << "\",\n"
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
        f << "\"totContErr\":"; dump_vec(f, totContErr); f << ",\n";

        f << "\"peak_areas\":{\n";
        f << "\"area\":";  dump_vec(f, area);  f << ",\n";
        f << "\"area_err\":";   dump_vec(f, area_err);
        f << "\n},\n";

        f << "\"fit_quality\":{"
          << "\"D_over_N\":" << D_over_N << ",\n"
          << "\"chi2_over_ndf\":" << chi2_over_ndf
          << "},\n";

        f << "\"peak_seeding\":{\n";
        f << "\"mode\":\"" << (manualPredictPeaks ? "manual" : "auto") << "\",\n";
        f << "\"seedMeans\":";  dump_vec(f, seedMeans);  f << ",\n";
        f << "\"seedSigmas\":"; dump_vec(f, seedSigmas); f << ",\n";
        f << "\"seedAmps\":";   dump_vec(f, seedAmps); f << "\n";
        f << "},\n";

        f << "\"fit_params\":{\n";
        f << "\"means\":"; dump_vec(f, fitMeans); f << ",\n";
        f << "\"means_err\":";  dump_vec(f, err_fitMeans); f << ",\n";
        f << "\"sigmas\":"; dump_vec(f, fitSigmas); f << ",\n";
        f << "\"sigmas_err\":"; dump_vec(f, err_fitSigmas); f << ",\n";
        f << "\"amps\":"; dump_vec(f, fitAmps); f << ",\n";
        f << "\"amps_err\":"; dump_vec(f, err_fitAmps); f << ",\n";
        f << "\"bg\":" << bg1 << ",\n";
        f << "\"bg_err\":" << err_bg1 << ",\n";
        f << "}\n";
        f << "}\n\n";
    }
};

static inline Double_t getenv_double(const Char_t* k, Double_t fallback){
  const Char_t* v = gSystem->Getenv(k);
  if(!v || !*v) return fallback;
  Char_t* end=nullptr;
  Double_t x = strtod(v, &end);
  return (end && *end==0) ? x : fallback;
}

static inline Int_t getenv_int(const Char_t* k, Int_t fallback){
  const Char_t* v = gSystem->Getenv(k);
  if(!v || !*v) return fallback;
  Char_t* end=nullptr;
  long x = strtol(v, &end, 10);
  return (end && *end==0) ? (Int_t)x : fallback;
}

static inline Bool_t getenv_bool(const Char_t* k, Bool_t fallback) {
  const Char_t* v = gSystem->Getenv(k);
  if (!v) return fallback;

  TString s(v);
  s = s.Strip(TString::kBoth); 
  s.ToLower();

  if (s=="1" || s=="true"  || s=="yes" || s=="on")  return kTRUE;
  if (s=="0" || s=="false" || s=="no"  || s=="off") return kFALSE;
  return fallback; 
}

static inline std::vector<Double_t> getenv_vecd(const Char_t* k, std::vector<Double_t> fallback){
  const Char_t* v = gSystem->Getenv(k);
  if(!v || !*v) return fallback;
  std::vector<Double_t> out;
  std::stringstream ss(v);
  std::string tok;
  while(std::getline(ss, tok, ',')){
    if(tok.empty()) continue;
    Char_t* end=nullptr;
    Double_t x = strtod(tok.c_str(), &end);
    if(!(end && *end==0)) { out.clear(); return fallback; }
    out.push_back(x);
  }
  return out.empty()? fallback : out;
}

static inline Int_t getenv_pid(const Char_t* key, Int_t fallback){
  const Char_t* v = gSystem->Getenv(key);
  if(!v || !*v) return fallback;

  Char_t* end = nullptr;
  long x = strtol(v, &end, 10);
  if (end && *end == 0 && x >= 0 && x <= 4) return (Int_t)x;

  TString s(v);
  s = s.Strip(TString::kBoth);
  s.ToLower();

  if (s=="el" || s=="e"  || s=="electron") return kEl;
  if (s=="mu" || s=="muon")                 return kMu;
  if (s=="pi" || s=="pion")                 return kPi;
  if (s=="ka" || s=="k"   || s=="kaon")     return kKa;
  if (s=="pr" || s=="p"   || s=="proton")   return kPr;

  return fallback;
}

void Plot_Different_nSigmas(){
    auto help = new helper();
    const Int_t nParts = helper::nParts;
    const Int_t NtrkMax = help->NtrkMax;
    const Int_t   nBins   = 500;
    const Double_t xMin1   = getenv_double("XMIN1", -12.0), xMax1 = getenv_double("XMAX1", 10.0);
    const Double_t xMin2   = getenv_double("XMIN2", -5), xMax2 = getenv_double("XMAX2", 15.0);
    const Double_t pStart = getenv_double("PSTART", 0.45), pEnd = getenv_double("PEND", 0.55), step = 0.1;
    const Double_t muWindow = 2.0;
    const Double_t mergeDistanceFactor = 1.0;
    const Double_t nEntriesLimit = 1e7; 
    const Bool_t FitKaonExclComp = getenv_bool("FITKAONEXCLCOMP", true);
    const Bool_t FitProtonExclComp = getenv_bool("FITPROTONEXCLCOMP", false);
    const Bool_t plotTPC = true;
    const Bool_t plotTOF = false;
    const Bool_t PeakZoom = false;
    const Bool_t manualPredictPeaks = true;
    const Double_t eigenThr = 0.0;
    const Bool_t useOffDiag = false;
    const Bool_t plotCM = false;
    const std::vector<Double_t> sigmaExclList = {0.0};
    const Bool_t KaonExcl = false;
    const Bool_t ProtonExcl = false;
    const Int_t primaryRef   = getenv_pid("PRIMARY_REF",   kEl);
    const Int_t secondaryRef = getenv_pid("SECONDARY_REF", kPi); 
    using PeakPars = std::array<Double_t,4>;

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

    ndJsonLogger ndlog;
    ndlog.open(Form("nSigmaEl-Plot-Data-%.2f<p%.2f.ndjson", pStart, pEnd));
    ndlog.write_config(nBins, xMin1, xMax1, xMin2, xMax2, pStart, pEnd, step, muWindow, mergeDistanceFactor, nEntries, FitKaonExclComp, FitProtonExclComp, plotTPC, plotTOF, PeakZoom, manualPredictPeaks, eigenThr, useOffDiag);

    auto drawNSigma = [&](Bool_t isTPCmode, Double_t sigmaExcl) {
        const Double_t pMin   = isTPCmode ? pStart : std::max(pStart, 0.4);
        const Int_t    nSteps = Int_t(std::floor((pEnd - pMin) / step + 0.5));
        std::vector<Double_t> pEdges(nSteps+1);
        for (Int_t i = 0; i <= nSteps; ++i) pEdges[i] = pMin + i * step;
        TString suffix = isTPCmode ? "TPC" : "TOF";    

        std::vector<TH1D*> h_pri(nSteps, nullptr);
        std::vector<TH1D*> h_sec(nSteps, nullptr);
        std::vector<TH1D*> h1cm(nSteps, nullptr);
        std::vector<TH1D*> h2cm(nSteps, nullptr);
        std::vector<covarianceMatrix*> cmObj(nSteps, nullptr);
    
        for (Int_t i=0; i<nSteps; ++i) {
            auto name1 = Form("n#sigma_%s %g < p < %g GeV/c (%s)", help->pCodes[primaryRef], pEdges[i], pEdges[i+1], suffix.Data());
            auto name2 = Form("n#sigma_%s %g < p < %g GeV/c (%s)", help->pCodes[secondaryRef], pEdges[i], pEdges[i+1], suffix.Data());
            auto title1 = Form("n#sigma_{%s} %g < p < %g GeV/c (%s); n#sigma_{%s}; Counts", help->pCodes[primaryRef], pEdges[i], pEdges[i+1], suffix.Data(), help->pCodes[primaryRef]);
            auto title2 = Form("n#sigma_{%s} %g < p < %g GeV/c (%s); n#sigma_{%s}; Counts", help->pCodes[secondaryRef], pEdges[i], pEdges[i+1], suffix.Data(), help->pCodes[secondaryRef]);
            h_pri[i] = new TH1D(name1, title1, nBins, xMin1, xMax1);
            h_sec[i] = new TH1D(name2, title2, nBins, xMin2, xMax2);
            h_pri[i]->SetMarkerStyle(kFullCircle);
            h_pri[i]->SetMarkerSize(0.75);
            h_pri[i]->SetMarkerColor(kBlack);
            h_pri[i]->SetLineColor(kBlack);
            h_sec[i]->SetMarkerStyle(kFullCircle);
            h_sec[i]->SetMarkerSize(0.75);
            h_sec[i]->SetMarkerColor(kBlack);
            h_sec[i]->SetLineColor(kBlack);

                    
            auto name_cm1 = Form("cm_nSTPC_%s_%g_%g", help->pCodes[primaryRef], pEdges[i], pEdges[i+1]);
            auto name_cm2 = Form("cm_nSTPC_%s_%g_%g", help->pCodes[secondaryRef], pEdges[i], pEdges[i+1]);

            h1cm[i] = new TH1D(name_cm1, h_pri[i]->GetTitle(), nBins, xMin1, xMax1);
            h2cm[i] = new TH1D(name_cm2 , h_sec[i]->GetTitle(), nBins, xMin2, xMax2);

            std::vector<TH1D*> cmHists{h1cm[i], h2cm[i]};
            std::vector<std::tuple<Double_t, Double_t>*> fitlims;
            fitlims.push_back(new std::tuple<Double_t,Double_t>(xMin1, xMax1));   
            fitlims.push_back(new std::tuple<Double_t,Double_t>(xMin2, xMax2)); 
            cmObj[i] = new covarianceMatrix(cmHists, fitlims);
        }
                   
        for (Long64_t ev = 0; ev < nEntries; ++ev) {
            chain.GetEntry(ev);
            for (Int_t t = 0; t < NtrkMax; ++t) {
                Double_t pG = inner[t];
                if (tofExpMom[t] < 0) continue;
                Int_t bin = std::lower_bound(pEdges.begin(), pEdges.end(), pG)
                        - pEdges.begin() - 1;
                if (bin<0 || bin>=nSteps) continue;

                Float_t val1 = isTPCmode ? tpcNS[primaryRef][t] : tofNS[primaryRef][t];
                Float_t val2 = isTPCmode ? tpcNS[secondaryRef][t] : tofNS[secondaryRef][t];
                if (TMath::IsNaN(val1)) continue;
                if (TMath::IsNaN(val2)) continue;
                
                const bool passK = TMath::IsNaN(tofNS[kKa][t]) || TMath::Abs(tofNS[kKa][t]) >= sigmaExcl;
                const bool passP = TMath::IsNaN(tofNS[kPr][t]) || TMath::Abs(tofNS[kPr][t]) >= sigmaExcl;

                if ((!KaonExcl && !ProtonExcl) || (KaonExcl && passK) || (ProtonExcl && passP)) {
                    h_pri[bin]->Fill(val1);
                    h_sec[bin]->Fill(val2);
                }
                
                auto nsVals = new std::vector<Double_t>(); 
                nsVals->reserve(2);
                nsVals->clear(); 
                auto toadd  = new std::vector<Bool_t>();   
                toadd->reserve(2);
                toadd->clear();

                const bool keep = (!KaonExcl && !ProtonExcl) || (KaonExcl && passK) || (ProtonExcl && passP);
                if (keep) {
                nsVals->push_back(val1); toadd->push_back(true);
                nsVals->push_back(val2); toadd->push_back(true);
                } else {
                nsVals->push_back(-999.0); toadd->push_back(false);
                nsVals->push_back(-999.0); toadd->push_back(false);
                }
                                
                cmObj[bin]->addEvent(nsVals, toadd);
            }
        }
        
        TString pdfName = Form("nSigma%s_%s-%s.pdf", suffix.Data(), help->pNames[primaryRef], help->pNames[secondaryRef]);
        TCanvas* c = new TCanvas("c","", 950, 700);
        c->SetLeftMargin(0.15);
        c->SetLogy();
        c->Print(pdfName + "[");

        TString  cmPdf;
        TCanvas* cCM = nullptr;
        if (plotCM) {
            cmPdf = Form("CM_%s_%s-%s.pdf",suffix.Data(), help->pNames[primaryRef], help->pNames[secondaryRef]);
            cCM = new TCanvas(Form("cCM_%s_%s-%s", suffix.Data(), help->pCodes[primaryRef], help->pCodes[secondaryRef]), "", 800, 700);
            cCM->SetLeftMargin(0.15);
            cCM->SetRightMargin(0.15);
            cCM->Print(cmPdf + "[");
        }

        for (Int_t i = 0; i < nSteps; ++i) {
            auto cm = cmObj[i];
            cm->make(useOffDiag, eigenThr);
            if (plotCM) {
                cCM->cd();
                cCM->Clear();
                cm->plot(cCM, cmPdf);
                c->cd();
            }
            TH1D* h1 = h_pri[i];
            TH1D* h2 = h_sec[i];

            struct Peak {Double_t A1, mu1, sigma1; Double_t A2, mu2, sigma2; Int_t id; std::vector<Int_t> merged_ids;Bool_t alwaysSeparate = false;};
            std::vector<Peak> seeds;
            Double_t sliceMin = pEdges[i];
            Double_t sliceMax = pEdges[i+1];
            Double_t pMid     = 0.5 * (sliceMin + sliceMax);

            auto predictMuSig = [&](Int_t ref, Int_t hyp)->std::pair<Double_t,Double_t>{
                Double_t mu=std::numeric_limits<Double_t>::quiet_NaN(), sig=std::numeric_limits<Double_t>::quiet_NaN();
                if (isTPCmode) {
                    Double_t dRef = help->getTPCSignal(pMid*1000, help->pMasses[ref], 1.0);
                    Double_t dHyp = help->getTPCSignal(pMid*1000, help->pMasses[hyp], 1.0);
                    if (dRef<0 || dHyp<0) return {mu,sig};
                    Double_t resoHyp = help->resoTPC[hyp][ref];
                    Double_t resoRef = help->resoTPC[ref][hyp];
                    sig = (resoHyp/resoRef) * (dHyp/dRef);
                    mu  = (dHyp/dRef - 1.0)/resoRef;
                } else {
                    Double_t bRef = help->beta(help->pMasses[ref], pMid*1000);
                    Double_t bHyp = help->beta(help->pMasses[hyp], pMid*1000);
                    Double_t resoHyp = help->resoTOF[hyp][ref];
                    Double_t resoRef = help->resoTOF[ref][hyp];
                    sig = (resoHyp/resoRef) * (1.0/(bHyp*bHyp));
                    mu  = (bRef - bHyp)/(bHyp*bHyp*resoRef);
                }
                return {mu,sig};
            };

            for (Int_t hyp=0; hyp<nParts; ++hyp) {
                auto [mu1,sig1] = predictMuSig(primaryRef,  hyp);
                auto [mu2,sig2] = predictMuSig(secondaryRef,hyp);
                if (!std::isfinite(mu1) || mu1<xMin1 || mu1>xMax1) continue; 

                Int_t bin1 = h1->FindBin(mu1);
                Int_t bin2 = h2->FindBin(mu2);
                Double_t A1 = h1->GetBinContent(bin1);
                Double_t A2 = h2->GetBinContent(bin2);

                seeds.push_back({A1,mu1,sig1, A2,mu2,sig2, hyp, {hyp}, false});
            }

            auto combine_one = [&](Double_t A_old,Double_t mu_old,Double_t sig_old, Double_t A_new,Double_t mu_new,Double_t sig_new) -> std::tuple<Double_t,Double_t,Double_t> {
                const Double_t I_old = A_old*sig_old*TMath::Sqrt(2*TMath::Pi());
                const Double_t I_new = A_new*sig_new*TMath::Sqrt(2*TMath::Pi());
                const Double_t It    = I_old + I_new;
                const Double_t mu_eff  = (It>0) ? (I_old*mu_old + I_new*mu_new)/It : mu_old;
                const Double_t var_eff = (It>0) ? (I_old*(sig_old*sig_old + (mu_old-mu_eff)*(mu_old-mu_eff)) + I_new*(sig_new*sig_new + (mu_new-mu_eff)*(mu_new-mu_eff)))/It : sig_old*sig_old;
                const Double_t sig_eff = std::sqrt(std::max(0.0, var_eff));
                const Double_t A_eff   = (sig_eff>0) ? (It/(TMath::Sqrt(2*TMath::Pi())*sig_eff)) : 0.0;
                return {A_eff, mu_eff, sig_eff};
            };

            Peak muPiCombined{};
            Bool_t sawMu=false, sawPi=false;
            for (auto it = seeds.begin(); it != seeds.end(); ) {
                if (it->id == kMu || it->id == kPi) {
                    if (!sawMu && !sawPi) {
                    muPiCombined = *it;
                    } else {
                    {
                        auto [Aeff,mueff,sigeff] =
                        combine_one(muPiCombined.A1, muPiCombined.mu1, muPiCombined.sigma1, it->A1, it->mu1, it->sigma1);
                        muPiCombined.A1=Aeff; muPiCombined.mu1=mueff; muPiCombined.sigma1=sigeff;
                    }
                    {
                        auto [Aeff,mueff,sigeff] =
                        combine_one(muPiCombined.A2, muPiCombined.mu2, muPiCombined.sigma2, it->A2, it->mu2, it->sigma2);
                        muPiCombined.A2=Aeff; muPiCombined.mu2=mueff; muPiCombined.sigma2=sigeff;
                    }
                    }
                    muPiCombined.merged_ids.push_back(it->id);
                    sawMu |= (it->id == kMu);
                    sawPi |= (it->id == kPi);
                    it = seeds.erase(it);
                } else {
                    ++it;
                }
            }

            if (sawMu || sawPi) {
                muPiCombined.id = kMu;
                muPiCombined.alwaysSeparate = std::find(muPiCombined.merged_ids.begin(), muPiCombined.merged_ids.end(), primaryRef)  != muPiCombined.merged_ids.end() || std::find(muPiCombined.merged_ids.begin(), muPiCombined.merged_ids.end(), secondaryRef) != muPiCombined.merged_ids.end();
                seeds.push_back(muPiCombined);
            }

            for (auto &s : seeds) {
                if (std::find(s.merged_ids.begin(), s.merged_ids.end(), primaryRef)  != s.merged_ids.end() ||
                    std::find(s.merged_ids.begin(), s.merged_ids.end(), secondaryRef) != s.merged_ids.end()) {
                    s.alwaysSeparate = true;
                }
}

            std::vector<Peak> merged;
            auto merge_into = [&](Peak &p, const Peak &s){
                {
                    auto [Aeff,mueff,sigeff] = combine_one(p.A1,p.mu1,p.sigma1, s.A1,s.mu1,s.sigma1);
                    p.A1=Aeff; p.mu1=mueff; p.sigma1=sigeff;
                }
                {
                    auto [Aeff,mueff,sigeff] = combine_one(p.A2,p.mu2,p.sigma2, s.A2,s.mu2,s.sigma2);
                    p.A2=Aeff; p.mu2=mueff; p.sigma2=sigeff;
                }
                p.merged_ids.insert(p.merged_ids.end(), s.merged_ids.begin(), s.merged_ids.end());
            };
            
            auto one_pass = [&](auto get_mu, auto get_sig, std::vector<Peak>& v){
                std::sort(v.begin(), v.end(),
                            [&](const Peak& a, const Peak& b){ return get_mu(a) < get_mu(b); });

                std::vector<Peak> out;
                out.reserve(v.size());

                for (const auto& s : v) {
                    if (!out.empty()) {
                    Peak& p = out.back();
                        const Bool_t sep = p.alwaysSeparate || s.alwaysSeparate;
                        const Bool_t fin = std::isfinite(get_mu(p)) && std::isfinite(get_mu(s)) && std::isfinite(get_sig(p)) && std::isfinite(get_sig(s));
                        const Bool_t close = fin && std::abs(get_mu(p) - get_mu(s)) < mergeDistanceFactor ;
                        if (!sep && close) {merge_into(p, s); continue;}
                    }
                    out.push_back(s);
                }
                v.swap(out);
            };
            merged = seeds;
            one_pass([](const Peak& p){return p.mu1;}, [](const Peak& p){return p.sigma1;}, merged);
            one_pass([](const Peak& p){return p.mu2;}, [](const Peak& p){return p.sigma2;}, merged);

            Double_t x_low1  = xMin1;
            Double_t x_high1 = xMax1;
            Double_t x_low2  = xMin2;
            Double_t x_high2 = xMax2;
            
            const Double_t Nsigma = 8.0;

            auto compute_zoom = [&](auto get_mu, auto get_sig, Double_t xMin, Double_t xMax) -> std::pair<Double_t,Double_t> {
                Double_t xl = xMin, xh = xMax;
                if (PeakZoom && !merged.empty()) {
                    Double_t lo = get_mu(merged[0]) - Nsigma * get_sig(merged[0]);
                    Double_t hi = get_mu(merged[0]) + Nsigma * get_sig(merged[0]);
                    for (size_t k = 1; k < merged.size(); ++k) {
                    Double_t this_lo = get_mu(merged[k]) - Nsigma * get_sig(merged[k]);
                    Double_t this_hi = get_mu(merged[k]) + Nsigma * get_sig(merged[k]);
                    lo = std::min(lo, this_lo);
                    hi = std::max(hi, this_hi);
                    }
                    xl = std::max(lo, xMin);
                    xh = std::min(hi, xMax);
                }
                return {xl, xh};
            };

            if (PeakZoom) {
                auto z1 = compute_zoom([](const Peak& p){return p.mu1;},[](const Peak& p){return p.sigma1;}, xMin1, xMax1);
                x_low1  = z1.first;  x_high1 = z1.second;
                auto z2 = compute_zoom([](const Peak& p){return p.mu2;},[](const Peak& p){return p.sigma2;}, xMin2, xMax2);
                x_low2  = z2.first;  x_high2 = z2.second;
            }

            Int_t manualNGauss = 0;
            std::vector<Double_t> manualMeans1, manualSigmas1, manualAmps1;
            std::vector<Double_t> manualMeans2, manualSigmas2, manualAmps2;
            if (manualPredictPeaks) {
                manualNGauss = getenv_int("MANUAL_NGAUSS", 4);
                auto defMeans1  = std::vector<Double_t>{-6.0, -5.0, 0.0, 1.5};
                auto defSigmas1 = std::vector<Double_t>{1.0, 1.0, 1.0, 1.0};
                auto defAmps1   = std::vector<Double_t>{1e5, 1e5, 1e3, 1e3};

                auto defMeans2  = std::vector<Double_t>{-1.0, 0.0, 7.0, 9.0};
                auto defSigmas2 = std::vector<Double_t>{1.0, 1.0, 1.0, 1.0};
                auto defAmps2   = std::vector<Double_t>{1e5, 1e5, 1e3, 1e3};

                auto mMeans1  = getenv_vecd("MANUAL_MEANS1",  defMeans1);
                auto mSigmas1 = getenv_vecd("MANUAL_SIGMAS1", defSigmas1);
                auto mAmps1   = getenv_vecd("MANUAL_AMPS1",   defAmps1);

                auto mMeans2  = getenv_vecd("MANUAL_MEANS2",  defMeans2);
                auto mSigmas2 = getenv_vecd("MANUAL_SIGMAS2", defSigmas2);
                auto mAmps2   = getenv_vecd("MANUAL_AMPS2",   defAmps2);

                Int_t L1 = std::min({manualNGauss, (Int_t)mMeans1.size(), (Int_t)mSigmas1.size(), (Int_t)mAmps1.size()});
                Int_t L2 = std::min({manualNGauss, (Int_t)mMeans2.size(), (Int_t) mSigmas2.size(), (Int_t) mAmps2.size()});
                Int_t L  = std::min(L1, L2);

                if (L > 0) {
                    manualMeans1.assign (mMeans1.begin(),  mMeans1.begin()+L);
                    manualSigmas1.assign(mSigmas1.begin(), mSigmas1.begin()+L);
                    manualAmps1.assign  (mAmps1.begin(),   mAmps1.begin()+L);

                    manualMeans2.assign (mMeans2.begin(),  mMeans2.begin()+L);
                    manualSigmas2.assign(mSigmas2.begin(), mSigmas2.begin()+L);
                    manualAmps2.assign  (mAmps2.begin(),   mAmps2.begin()+L);

                    manualNGauss = L;
                } else {
                    manualNGauss = 0;
                }
            }

            size_t nG = 0;
            std::vector<Double_t> seedMeans1, seedSigmas1, seedAmps1;
            std::vector<Double_t> seedMeans2, seedSigmas2, seedAmps2;
            if (manualPredictPeaks && manualNGauss > 0) {
                nG = manualNGauss;
                seedMeans1  = manualMeans1;  seedSigmas1 = manualSigmas1;  seedAmps1 = manualAmps1;
                seedMeans2  = manualMeans2;  seedSigmas2 = manualSigmas2;  seedAmps2 = manualAmps2;
            } else {
                nG = merged.size();
                seedMeans1.resize(nG); seedSigmas1.resize(nG); seedAmps1.resize(nG);
                seedMeans2.resize(nG); seedSigmas2.resize(nG); seedAmps2.resize(nG);
                for (size_t i = 0; i < nG; ++i) {
                    seedMeans1[i]  = merged[i].mu1;   seedSigmas1[i] = merged[i].sigma1;  seedAmps1[i] = merged[i].A1;
                    seedMeans2[i]  = merged[i].mu2;   seedSigmas2[i] = merged[i].sigma2;  seedAmps2[i] = merged[i].A2;
                }
            }
            if (nG < 1) { c->Print(pdfName); delete h1; continue; }

            h1->GetXaxis()->SetRangeUser(x_low1, x_high1);
            h2->GetXaxis()->SetRangeUser(x_low2, x_high2);

            const Int_t offA1 = 0; 
            const Int_t offA2 = nG; 
            const Int_t offM1 = 2*nG; 
            const Int_t offM2 = 3*nG; 
            const Int_t offS1 = 4*nG; 
            const Int_t offS2 = 5*nG; 
            const Int_t offP1 = 6*nG; 
            const Int_t offP2 = 6*nG + 1;

            std::ostringstream form;
            for (Int_t i = 0; i < nG; ++i) {
                if (i) form << "+";
                form << Form("[%d]*exp(-0.5*pow((x-[%d])/( [%d] ),2))",
                            offA1+i, offM1+i, offS1+i);
                form << "+";
                form << Form("[%d]*exp(-0.5*pow((x-[%d])/( [%d] ),2))",
                            offA2+i, offM2+i, offS2+i);
            }
            form << "+pol0(" << 6*nG << ")"  // background h1
                << "+pol0(" << 6*nG+1 << ")"; // background h2

            TF1* sum = new TF1("sum", form.str().c_str());
            sum->SetNpx(500); 
            sum->SetParLimits(offP1, 0, 10.0);
            sum->SetParameter(offP1, 1.0);
            sum->SetParLimits(offP2, 0, 10.0);
            sum->SetParameter(offP2, 1.0);
            for (size_t i = 0; i < nG; ++i) {
                if (manualPredictPeaks && manualNGauss > 0) {
                    Double_t mu1 = seedMeans1[i];
                    Double_t sig1 = seedSigmas1[i]; 
                    Double_t amp1 = seedAmps1[i]; 

                    Double_t mu2 = seedMeans2[i];
                    Double_t sig2 = seedSigmas2[i]; 
                    Double_t amp2 = seedAmps2[i];

                    sum->SetParLimits(offA1 + i, 0.0, 1e10);
                    sum->SetParameter(offA1 + i, amp1);

                    sum->SetParLimits(offA2 + i, 0, 1e10);
                    sum->SetParameter (offA2 + i, amp2);

                    sum->SetParLimits(offM1 + i, mu1 - muWindow, mu1 + muWindow);
                    sum->SetParameter(offM1 + i, mu1);

                    sum->SetParLimits(offM2 + i, mu2 - muWindow, mu2 + muWindow);
                    sum->SetParameter(offM2 + i, mu2);

                    sum->SetParLimits(offS1 + i, sig1*0.5, sig1*2.0);
                    sum->SetParameter(offS1 + i, sig1);

                    sum->SetParLimits(offS2 + i, sig2*0.5, sig2*2.0);
                    sum->SetParameter(offS2 + i, sig2);
                }
                else {
                    const auto &p = merged[i];

                    sum->SetParLimits(offA1 + i, 0.0, 1e10);
                    sum->SetParameter(offA1 + i, p.A1);

                    sum->SetParLimits(offA2 + i, 0.0, 1e10);
                    sum->SetParameter (offA2 + i, p.A2);

                    sum->SetParLimits(offM1 + i, p.mu1 - muWindow, p.mu1 + muWindow);
                    sum->SetParameter(offM1 + i, p.mu1);
                    
                    sum->SetParLimits(offM2 + i, p.mu2 - muWindow, p.mu2 + muWindow);
                    sum->SetParameter(offM2 + i, p.mu2);

                    sum->SetParLimits(offS1 + i, p.sigma1*0.5, p.sigma1*2.0);
                    sum->SetParameter(offS1 + i, p.sigma1);

                    sum->SetParLimits(offS2 + i, p.sigma2*0.5, p.sigma2*2.0);
                    sum->SetParameter(offS2 + i, p.sigma2);
                }
            }
            help->FitHistogramsByChi2(h1, h2, sum, nG, cmObj[i], std::vector<TH1D*>{h1cm[i], h2cm[i]}, eigenThr, nBins);
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
                std::vector<Double_t> totContErr;
            };

            struct IntDerivs { Double_t dA, dMu, dSig; };
            auto dGaussInt = [&](Double_t A, Double_t mu, Double_t sig, Double_t a, Double_t b)->IntDerivs {
                const Double_t ea = std::exp(-0.5 * std::pow((a-mu)/sig, 2.0));
                const Double_t eb = std::exp(-0.5 * std::pow((b-mu)/sig, 2.0));
                const Double_t I1 = help->GaussIntegral(1.0, mu, sig, a, b); 
                const Double_t dI_dMu  = A * (ea - eb);                      
                const Double_t I       = help->GaussIntegral(A,   mu, sig, a, b);
                const Double_t dI_dSig = (I/sig) + (A/sig)*((a-mu)*ea - (b-mu)*eb); 
                return { I1, dI_dMu, dI_dSig };
            };

            auto idxA = [&](Int_t which, Int_t i){ return (which==0) ? offA1 + i : offA2 + i; };
            auto idxM = [&](Int_t which, Int_t i){ return (which==0) ? offM1 + i : offM2 + i; };
            auto idxS = [&](Int_t which, Int_t i){ return (which==0) ? offS1 + i : offS2 + i; };
            auto idxBkg = [&](Int_t which){ return (which==0) ? offP1 : offP2; };

            auto computeBandFractions = [&](Int_t whichHist) -> BandStats {
                const Double_t kLeft = 3.0, kRight = 3.0;

                std::vector<std::vector<Double_t>> frac(nG, std::vector<Double_t>(nG, std::numeric_limits<Double_t>::quiet_NaN()));
                std::vector<Double_t> totCont(nG, std::numeric_limits<Double_t>::quiet_NaN());
                std::vector<Double_t> totContErr(nG, std::numeric_limits<Double_t>::quiet_NaN());

                auto integrate_sum_others = [&](Int_t ii, Double_t a, Double_t b) -> Double_t {
                    const Int_t N = 20480;
                    const Double_t dx = (b - a) / (N - 1);
                    Double_t acc = 0.0;
                    for (Int_t t = 0; t < N; ++t) {
                        const Double_t x = a + t*dx;
                        const Double_t w = (t==0 || t==N-1) ? 0.5 : 1.0;
                        Double_t y = 0.0;
                        for (Int_t j = 0; j < (Int_t)nG; ++j) if (j != ii) {
                            const Double_t Aj   = par[idxA(whichHist,j)];
                            const Double_t muj  = par[idxM(whichHist,j)];
                            const Double_t sigj = par[idxS(whichHist,j)];
                            y += Aj * TMath::Gaus(x, muj, sigj, kFALSE);
                        }
                        acc += w * y;
                    }
                    return acc * dx;
                };

                for (Int_t ii = 0; ii < (Int_t)nG; ++ii) {
                    const Double_t Ai   = par[idxA(whichHist,ii)];
                    const Double_t mui  = par[idxM(whichHist,ii)];
                    const Double_t sigi = par[idxS(whichHist,ii)];
                    if (!(std::isfinite(mui) && std::isfinite(sigi) && sigi>0)) continue;

                    const Double_t a = mui - kLeft  * sigi;
                    const Double_t b = mui + kRight * sigi;

                    const Double_t are_i = help->GaussIntegral(Ai, mui, sigi, a, b);

                    for (Int_t j = 0; j < (Int_t)nG; ++j) {
                        const Double_t Aj   = par[idxA(whichHist,j)];
                        const Double_t muj  = par[idxM(whichHist,j)];
                        const Double_t sigj = par[idxS(whichHist,j)];
                        const Double_t num_j = help->GaussIntegral(Aj, muj, sigj, a, b);
                        frac[ii][j] = (are_i > 0.0) ? (num_j / are_i) : std::numeric_limits<Double_t>::quiet_NaN();
                    }

                    if (are_i > 0.0) {
                        const Double_t S = integrate_sum_others(ii, a, b);
                        const Double_t R = S / are_i;
                        totCont[ii] = R;

                        Double_t varR = 0.0;
                        for (Int_t k = 0; k < (Int_t)nG; ++k) if (k != ii) {
                            const Double_t Ak   = par[idxA(whichHist,k)];
                            const Double_t muk  = par[idxM(whichHist,k)];
                            const Double_t sigk = par[idxS(whichHist,k)];
                            auto dK = dGaussInt(Ak, muk, sigk, a, b);

                            const Double_t gA  = dK.dA   / are_i;
                            const Double_t gMu = dK.dMu  / are_i;
                            const Double_t gS  = dK.dSig / are_i;

                            varR += gA  * gA  * err[idxA(whichHist,k)] * err[idxA(whichHist,k)];
                            varR += gMu * gMu * err[idxM(whichHist,k)] * err[idxM(whichHist,k)];
                            varR += gS  * gS  * err[idxS(whichHist,k)] * err[idxS(whichHist,k)];
                        }
                        {
                            auto dI = dGaussInt(Ai, mui, sigi, a, b);
                            const Double_t c = -(R / are_i);
                            const Double_t gA  = c * dI.dA;
                            const Double_t gMu = c * dI.dMu;
                            const Double_t gS  = c * dI.dSig;

                            varR += gA  * gA  * err[idxA(whichHist,ii)] * err[idxA(whichHist,ii)];
                            varR += gMu * gMu * err[idxM(whichHist,ii)] * err[idxM(whichHist,ii)];
                            varR += gS  * gS  * err[idxS(whichHist,ii)] * err[idxS(whichHist,ii)];
                        }
                        totContErr[ii] = (varR >= 0.0) ? std::sqrt(varR) : std::numeric_limits<Double_t>::quiet_NaN();
                    }
                }
                return BandStats{kLeft, kRight, std::move(frac), std::move(totCont), std::move(totContErr)};
            };

            auto [D1,N1] = help->PoissonDeviance(h1, 0, par, nG, offA1, offA2, offM1, offM2, offS1, offS2, offP1, offP2);
            auto [D2,N2] = help->PoissonDeviance(h2, 1, par, nG, offA1, offA2, offM1, offM2, offS1, offS2, offP1, offP2);
            const Double_t D  = D1 + D2;
            const Double_t N  = N1 + N2;
            const Int_t    k  = 4*nG + 2; 

            const Double_t chi2 = sum->GetChisquare();
            const Int_t    ndf  = sum->GetNDF();
            const Double_t pval = TMath::Prob(chi2, ndf);

            std::vector<Double_t> areas1(nG), areas2(nG), err_areas1(nG), err_areas2(nG);
            const Double_t factor = TMath::Sqrt(2*TMath::Pi());
            for (Int_t ig = 0; ig < nG; ++ig) {
                Double_t A1   = par[offA1 + ig];
                Double_t A2   = par[offA2 + ig];
                Double_t sig1  = par[offS1 + ig];
                Double_t sig2  = par[offS2 + ig];
                Double_t eA1 = err[offA1 + ig];
                Double_t eA2 = err[offA2 + ig];
                Double_t eS1  = err[offS1 + ig];
                Double_t eS2  = err[offS2 + ig];

                areas1[ig] = A1  * sig1 * factor;
                areas2[ig] = A2  * sig2 * factor;

                Double_t dFdA1 = sig1 * factor;
                Double_t dFdS1 = A1 * factor;        
                Double_t var1 = dFdA1*dFdA1 * eA1*eA1+ dFdS1*dFdS1 * eS1*eS1; 
                err_areas1[ig] = TMath::Sqrt(var1);

                Double_t dFdA2 = sig2 * factor;
                Double_t dFdS2 = A2 * factor;        
                Double_t var2 = dFdA2*dFdA2 * eA2*eA2+ dFdS2*dFdS2 * eS2*eS2; 
                err_areas2[ig] = TMath::Sqrt(var2);
            }
            const Double_t D_over_N      = (N   > 0 ? D   / N   : std::numeric_limits<Double_t>::quiet_NaN());
            const Double_t chi2_over_ndf = (ndf > 0 ? chi2/ ndf : std::numeric_limits<Double_t>::quiet_NaN());

            std::vector<Double_t> fitMeans1(nG), fitSigmas1(nG), fitAmps1(nG);
            std::vector<Double_t> fitMeans2(nG), fitSigmas2(nG), fitAmps2(nG);
            std::vector<Double_t> err_fitMeans1(nG), err_fitSigmas1(nG), err_fitAmps1(nG);
            std::vector<Double_t> err_fitMeans2(nG), err_fitSigmas2(nG), err_fitAmps2(nG);
            for (Int_t ig=0; ig<nG; ++ig) {
                fitAmps1[ig] = par[offA1 + ig];;
                fitMeans1[ig] = par[offM1  + ig];
                fitSigmas1[ig] = par[offS1  + ig];
                fitAmps2[ig] = par[offA2 + ig];;
                fitMeans2[ig] = par[offM2  + ig];
                fitSigmas2[ig] = par[offS2  + ig];

                err_fitAmps1[ig] = err[offA1 + ig];
                err_fitMeans1[ig] = err[offM1  + ig];
                err_fitSigmas1[ig] = err[offS1  + ig];
                err_fitAmps2[ig] = err[offA2 + ig];
                err_fitMeans2[ig] = err[offM2  + ig];
                err_fitSigmas2[ig] = err[offS2  + ig];
            }
            const Double_t bkg1   = par[offP1];   
            const Double_t bkg2   = par[offP2];
            const Double_t ebkg1  = err[offP1];
            const Double_t ebkg2   = err[offP2];

            auto band_no = computeBandFractions(0); 
            auto band_w  = computeBandFractions(1); 

            ndlog.write_slice(
                suffix.Data(), "h1", sigmaExcl, help->pNames[primaryRef],
                i, pEdges[i], pEdges[i+1],
                band_no.kLeft, band_no.kRight, band_no.frac, band_no.totCont, band_no.totContErr,
                areas1, err_areas1,  
                D_over_N, chi2_over_ndf,
                manualPredictPeaks, seedMeans1, seedSigmas1, seedAmps1,
                fitMeans1, fitSigmas1, fitAmps1, err_fitMeans1, err_fitSigmas1, err_fitAmps1,
                bkg1, ebkg1);

            ndlog.write_slice(
                suffix.Data(), "h2", sigmaExcl, help->pNames[secondaryRef],
                i, pEdges[i], pEdges[i+1],
                band_w.kLeft, band_w.kRight, band_w.frac, band_w.totCont, band_w.totContErr,
                areas2, err_areas2,  
                D_over_N, chi2_over_ndf,
                manualPredictPeaks, seedMeans2, seedSigmas2, seedAmps2,
                fitMeans2, fitSigmas2, fitAmps2, err_fitMeans2, err_fitSigmas2, err_fitAmps2,
                bkg2, ebkg2);

            const Int_t nPoints = 500; 
            Double_t dx1  = (x_high1 - x_low1)/(nPoints-1);
            Double_t dx2  = (x_high2 - x_low2)/(nPoints-1);
            std::vector<Double_t> xv1(nPoints), xv2(nPoints), y1(nPoints), y2(nPoints);
                
            const double bw1 = h1->GetXaxis()->GetBinWidth(1);
            const double bw2 = h2->GetXaxis()->GetBinWidth(1);
            for(Int_t ip=0; ip<nPoints; ++ip) {
                Double_t x1 = x_low1 + ip*dx1;
                Double_t x2 = x_low2 + ip*dx2;
                xv1[ip] = x1;
                xv2[ip] = x2;
                Double_t yy1 = par[offP1] * bw1;
                Double_t yy2 = par[offP2] * bw2;
                for(Int_t ig=0; ig<nG; ++ig) {
                    Double_t A1  = par[offA1 + ig];
                    Double_t A2  = par[offA2 + ig];
                    Double_t mu1  = par[offM1  + ig];
                    Double_t sig1 = par[offS1  + ig];
                    Double_t mu2  = par[offM2  + ig];
                    Double_t sig2 = par[offS2  + ig];
                    yy1 += A1 * TMath::Gaus(x1, mu1, sig1, kFALSE) * bw1;
                    yy2 += A2 * TMath::Gaus(x2, mu2, sig2, kFALSE) * bw2;
                }
                y1[ip] = yy1;
                y2[ip] = yy2;
            }

            auto gSum1 = new TGraph(nPoints, xv1.data(), y1.data());
            gSum1->SetLineColor(kRed);
            gSum1->SetLineWidth(2);
            gSum1->Draw("L SAME");

            auto gSum2 = new TGraph(nPoints, xv2.data(), y2.data());
            gSum2->SetLineColor(kRed);
            gSum2->SetLineWidth(2);

            if (!manualPredictPeaks) {
                for (const auto &pk : merged) {
                    Double_t mu = pk.mu1;
                    TLine *l = new TLine(mu, 0, mu, yMax1);
                    Int_t col  = (pk.id >= 0) ? help->colors[pk.id] : kGray+2;
                    l->SetLineColor(col);
                    l->Draw();
                }
            }
            else {
                for (Int_t i = 0; i < nG; ++i) {
                    Double_t mu = seedMeans1[i];
                    TLine *l = new TLine(mu, 0, mu, yMax1);
                    l->SetLineColor(help->colors[i % nParts]);
                    l->Draw();
                }
            }
                
            for (Int_t i = 0; i < nG; ++i) {
                if (sum->GetParameter(offA1 + i) < 0) continue;
                TF1 *g1 = new TF1(Form("g_%d", i), "gaus", x_low1, x_high1);
                g1->SetParameters(
                    sum->GetParameter(offA1 + i) * bw1,
                    sum->GetParameter(offM1 + i),
                    sum->GetParameter(offS1 + i)
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
                    Double_t mu = pk.mu2;
                    TLine *l = new TLine(mu, 0, mu, yMax2);
                    Int_t col  = (pk.id >= 0) ? help->colors[pk.id] : kGray+2;
                    l->SetLineColor(col);
                    l->Draw();
                }
            }
            else {
                for (Int_t i = 0; i < nG; ++i) {
                    Double_t mu = seedMeans2[i];
                    TLine *l = new TLine(mu, 0, mu, yMax2);
                    l->SetLineColor(help->colors[i % nParts]);
                    l->Draw();
                }
            }

            for (Int_t i = 0; i < nG; ++i) {
                if (sum->GetParameter(offA2 + i) < 0) continue;
                TF1 *g2 = new TF1(Form("g_%d", i), "gaus", x_low2, x_high2);
                g2->SetParameters(
                    sum->GetParameter(offA2 + i) * bw2,
                    sum->GetParameter(offM2 + i),
                    sum->GetParameter(offS2 + i)
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
        if (plotCM) {
            cCM->Print(cmPdf + "]"); 
        }
    
    };
    for (Double_t s : sigmaExclList) {
        if (plotTPC) drawNSigma(true,  s);
        if (plotTOF) drawNSigma(false, s);
    }  
}

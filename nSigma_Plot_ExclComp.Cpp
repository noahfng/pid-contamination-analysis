#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>
#include <utility>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

#include "TGraph.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"  
#include "TLine.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TH1D.h"

#include <AddTrees.h>   // project-specific: file discovery/chain fill
#include <helpers.h>    // project-specific: masses, charges, colors, Bethe–Bloch
#include <covarianceMatrix.h>

// small text-formatting utilities for logging and table generation
static std::string vec_to_txt(const std::vector<Double_t>& v)
{
    std::ostringstream os;
    os << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) os << ", ";
        if (std::isfinite(v[i])) os << v[i];
        else os << "null";
    }
    os << "]";
    return os.str();
}

static std::string mat_to_txt(const std::vector<std::vector<Double_t>>& m)
{
    std::ostringstream os;
    os << "[";
    for (size_t i = 0; i < m.size(); ++i) {
        if (i) os << ", ";
        os << vec_to_txt(m[i]);
    }
    os << "]";
    return os.str();
}

static std::string labels_to_txt(const std::vector<TString>& labels)
{
    std::ostringstream os;
    os << "[";
    for (size_t i = 0; i < labels.size(); ++i) {
        if (i) os << ", ";
        os << labels[i].Data();
    }
    os << "]";
    return os.str();
}

static std::string labeled_vec_to_txt(const std::vector<TString>& labels, const std::vector<Double_t>& vals, const std::vector<Double_t>* errs = nullptr)
{
    std::ostringstream os;
    for (size_t i = 0; i < vals.size(); ++i) {
        TString lab = (i < labels.size()) ? labels[i] : Form("peak_%zu", i);
        os << "  " << lab.Data() << " : " << vals[i];
        if (errs && i < errs->size()) os << " ± " << (*errs)[i];
        os << "\n";
    }
    return os.str();
}

static std::string labeled_mat_to_txt(const std::vector<TString>& labels, const std::vector<std::vector<Double_t>>& m)
{
    std::ostringstream os;
    for (size_t i = 0; i < m.size(); ++i) {
        TString rowlab = (i < labels.size()) ? labels[i] : Form("peak_%zu", i);
        os << "  from " << rowlab.Data() << " : [";
        for (size_t j = 0; j < m[i].size(); ++j) {
            if (j) os << ", ";
            TString collab = (j < labels.size()) ? labels[j] : Form("peak_%zu", j);
            os << m[i][j];
        }
        os << "]\n";
    }
    return os.str();
}

static std::string fmt_pm_aligned(Double_t v, Double_t e, Int_t prec_v = 6, Int_t prec_e = 6, Bool_t sci = kFALSE)
{
    std::ostringstream os;
    if (sci) os << std::scientific;
    else     os << std::fixed;

    os << std::setprecision(prec_v) << v
       << " ± "
       << std::setprecision(prec_e) << e;
    return os.str();
}

static std::string species_table_to_txt(const std::vector<TString>& labels, const std::vector<Double_t>& means, const std::vector<Double_t>& mean_errs, const std::vector<Double_t>& sigmas, const std::vector<Double_t>& sigma_errs, const std::vector<Double_t>& amps, const std::vector<Double_t>& amp_errs, const std::vector<Double_t>& areas, const std::vector<Double_t>& area_errs)
{
    std::ostringstream os;

    const int w_species = 12;
    const int w_mean    = 30;
    const int w_sigma   = 30;
    const int w_amp     = 34;
    const int w_area    = 34;
    const std::string sep = "  ";

    os << std::left
       << std::setw(w_species) << "species" << sep
       << std::setw(w_mean)    << "mean"    << sep
       << std::setw(w_sigma)   << "sigma"   << sep
       << std::setw(w_amp)     << "amp"     << sep
       << std::setw(w_area)    << "area"
       << "\n";

    os << std::string(
            w_species + w_mean + w_sigma + w_amp + w_area + 4 * sep.size(),
            '-'
         ) << "\n";

    for (size_t i = 0; i < labels.size(); ++i) {
        const std::string mean_s  = fmt_pm_aligned(means[i],  mean_errs[i],  6, 6, kFALSE);
        const std::string sigma_s = fmt_pm_aligned(sigmas[i], sigma_errs[i], 6, 6, kFALSE);
        const std::string amp_s   = fmt_pm_aligned(amps[i],   amp_errs[i],   6, 6, kTRUE);
        const std::string area_s  = fmt_pm_aligned(areas[i],  area_errs[i],  6, 6, kTRUE);

        os << std::left
           << std::setw(w_species) << labels[i].Data() << sep
           << std::setw(w_mean)    << mean_s           << sep
           << std::setw(w_sigma)   << sigma_s          << sep
           << std::setw(w_amp)     << amp_s            << sep
           << std::setw(w_area)    << area_s
           << "\n";
    }

    return os.str();
}

// Per-window output for contamination analysis around each fitted peak
struct WindowScanEntry {
    Double_t nsigma;
    std::vector<Double_t> target_yield;
    std::vector<Double_t> target_yield_err;

    std::vector<Double_t> other_yield;
    std::vector<Double_t> other_yield_err;

    std::vector<Double_t> bg_yield;
    std::vector<Double_t> bg_yield_err;

    std::vector<Double_t> total_yield;
    std::vector<Double_t> total_yield_err;

    std::vector<Double_t> contamination;
    std::vector<Double_t> contamination_err;
};

using WindowScan = std::vector<WindowScanEntry>;

static std::string scan_table_to_txt(const std::vector<TString>& labels, const WindowScan& scan)
{
    std::ostringstream os;
    os << std::fixed << std::setprecision(6);

    const int w_species = 12;
    const int w_yield   = 34;
    const int w_cont    = 26;
    const std::string sep = "  ";

    for (const auto& entry : scan) {
        os << "nSigmaWindow = " << entry.nsigma << "\n";

        os << std::left
           << std::setw(w_species) << "species"       << sep
           << std::setw(w_yield)   << "target_yield"  << sep
           << std::setw(w_yield)   << "other_yield"   << sep
           << std::setw(w_yield)   << "bg_yield"      << sep
           << std::setw(w_yield)   << "total_yield"   << sep
           << std::setw(w_cont)    << "contamination"
           << "\n";

        os << std::string(
                w_species + 4*w_yield + w_cont + 5*sep.size(),
                '-'
             ) << "\n";

        for (size_t i = 0; i < labels.size(); ++i) {
            TString lab = (i < labels.size()) ? labels[i] : Form("peak_%zu", i);

            os << std::left
               << std::setw(w_species) << lab.Data() << sep
               << std::setw(w_yield)   << fmt_pm_aligned(entry.target_yield[i], entry.target_yield_err[i], 6, 6, kFALSE) << sep
               << std::setw(w_yield)   << fmt_pm_aligned(entry.other_yield[i],  entry.other_yield_err[i],  6, 6, kFALSE) << sep
               << std::setw(w_yield)   << fmt_pm_aligned(entry.bg_yield[i],     entry.bg_yield_err[i],     6, 6, kFALSE) << sep
               << std::setw(w_yield)   << fmt_pm_aligned(entry.total_yield[i],  entry.total_yield_err[i],  6, 6, kFALSE) << sep
               << std::setw(w_cont)    << fmt_pm_aligned(entry.contamination[i], entry.contamination_err[i], 6, 6, kFALSE)
               << "\n";
        }

        os << "\n";
    }

    return os.str();
}

// Environment variable readers so the macor can be steered without editing code (will return fallback value if variable is not set or invalid)
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

static inline std::string trim_copy(std::string s){
  auto isws = [](unsigned char c){ return std::isspace(c); };
  while(!s.empty() && isws(s.front())) s.erase(s.begin());
  while(!s.empty() && isws(s.back()))  s.pop_back();
  return s;
}

static inline std::vector<TString> getenv_labels(const Char_t* k, const std::vector<TString>& fallback)
{
  const Char_t* v = gSystem->Getenv(k);
  if(!v || !*v) return fallback;

  std::string raw(v);
  // allow either ',' or '|'
  for (auto &ch : raw) if (ch=='|') ch=',';

  std::vector<TString> out;
  std::stringstream ss(raw);
  std::string tok;
  while(std::getline(ss, tok, ',')){
    tok = trim_copy(tok);
    if(tok.empty()) continue;
    out.push_back(TString(tok));
  }
  return out.empty() ? fallback : out;
}

void nSigma_Plot_ExclComp(){
    auto help = new helper();
    const Int_t nParts = helper::nParts;
    const Int_t NtrkMax = help->NtrkMax;

    // histogram & slicing 
    const Int_t nBins = 500;
    const Double_t xMin = getenv_double("XMIN", -10.0), xMax = getenv_double("XMAX", 10.0);
    const Double_t pStart = getenv_double("PSTART", 0.35), pEnd = getenv_double("PEND", 0.45), step = 0.1; // momentum slicing (GeV/c)

    // peak finding & fitting config
    const Double_t muWindow = 2.0; // ±window for μ during fit init
    const Double_t mergeDistanceFactor = 1.0; // merge seeds if |Δμ| < f·max(σ)
    const Double_t nEntriesLimit = 1e10;

    // Analysis modes
    const Bool_t FitKaonExclComp = getenv_bool("FITKAONEXCLCOMP", true);
    const Bool_t FitProtonExclComp = getenv_bool("FITPROTONEXCLCOMP", false);
    const Bool_t plotTPC = true;
    const Bool_t plotTOF = false;

    // drawing options
    const Bool_t drawPredictedLines = false; // vertical lines at predicted peak positions
    const Bool_t drawSum = true; // overlay total fit
    const Bool_t drawComponents = true; // overlay individual Gaussian components

    // manual peak control
    const Bool_t manualPredictPeaks = true;
    const Bool_t interactiveManualPeaks = getenv_bool("INTERACTIVE_MANUAL_PEAKS", manualPredictPeaks);
    
    // logging
    const Bool_t writeSummaryLog = getenv_bool("WRITE_SUMMARY_LOG", true);

    // covariance matrix config
    const Double_t eigenThr = 0.0;
    const Bool_t useOffDiag = false;
    const Bool_t plotCM = false;

    // contanimation analysis
    std::vector<double> sigmaScan = {0.5, 1.0, 1.5, 2.0, 3.0}; // width of the area around the peak to analyze for contamination (in units of fitted σ)

    // PID selection: which particle species to run the analysis for (in order: e, π, K, p, d)
    const std::array<Bool_t, nParts> doPid = {{true, false, false, false, false}};
    using PeakPars = std::array<Double_t,4>;

    // exclusion scan config
    const std::vector<Double_t> sigmaExclList = {3.0};

    gROOT->SetBatch(!(manualPredictPeaks && interactiveManualPeaks)); // if manual peak prediction is on, default to interactive mode (unless explicitly disabled)
    gStyle->SetOptStat(0);

    // build chain and set branch addresses
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

    std::ofstream summaryLog;
    if (writeSummaryLog) {
        TString logName = Form("nSigmaSummary_%.2f_%.2f.txt", pStart, pEnd);
        summaryLog.open(logName.Data(), std::ios::out);
        summaryLog.imbue(std::locale::classic());
        summaryLog << std::fixed << std::setprecision(6);
    }
    // Main worker lambda:
    //   - chooses TPC or TOF nSigma as the plotted observable
    //   - fills spectra with and without a TOF-based exclusion
    //   - fits the pair of histograms simultaneously
    //   - writes plots + optional log output
    auto drawNSigma = [&](Bool_t isTPCmode, Double_t sigmaExcl) {
        const Double_t pMin = isTPCmode ? pStart : std::max(pStart, 0.4);
        const Int_t nSteps = Int_t(std::floor((pEnd - pMin) / step + 0.5));
        
        std::vector<Double_t> pEdges(nSteps+1);
        for (Int_t i = 0; i <= nSteps; ++i) pEdges[i] = pMin + i * step;

        TString suffix = isTPCmode ? "TPC" : "TOF";  

        // species used only for the exclusion cut
        enum ExclType{KaonExcl = 3, ProtExcl = 4};
        std::vector<std::pair<ExclType, Bool_t>> modes = {{KaonExcl, FitKaonExclComp}, {ProtExcl, FitProtonExclComp}};
        for (auto [excl, doThis]: modes) {
            if(!doThis) continue;

            const Char_t* name = (excl==KaonExcl ? "KaonExcl" : "ProtonExcl");

            // h_no = inclusive histogram
            // h_w = same histogram after excluding tracks compatible with 'excl' in TOF
            std::vector<std::vector<TH1D*>> h_no(nParts, std::vector<TH1D*>(nSteps,nullptr));
            std::vector<std::vector<TH1D*>> h_w(nParts, std::vector<TH1D*>(nSteps,nullptr));

            std::vector<std::vector<TH1D*>> hcm_no(nParts, std::vector<TH1D*>(nSteps,nullptr));
            std::vector<std::vector<TH1D*>> hcm_w(nParts, std::vector<TH1D*>(nSteps,nullptr));
            std::vector<std::vector<covarianceMatrix*>> cmObj(nParts, std::vector<covarianceMatrix*>(nSteps,nullptr));
        
            for (Int_t pid=0; pid<nParts; ++pid) if (doPid[pid]) {
                for (Int_t i=0; i<nSteps; ++i) {
                    auto name_no = Form("n#sigma_%s %g < p < %g GeV/c (%s-no%s)", help->pCodes[pid], pEdges[i], pEdges[i+1], suffix.Data(), name);
                    auto name_w = Form("n#sigma_%s %g < p < %g GeV/c (%s-%.2f-%s)", help->pCodes[pid], pEdges[i], pEdges[i+1], suffix.Data(), sigmaExcl, name);
                    auto title_no = Form("n#sigma_{%s} %g < p < %g GeV/c (%s-no%s); n#sigma_{%s}; Counts", help->pCodes[pid], pEdges[i], pEdges[i+1], suffix.Data(), name, help->pCodes[pid]);
                    auto title_w = Form("n#sigma_{%s} %g < p < %g GeV/c (%s-%.2f-%s); n#sigma_{%s}; Counts", help->pCodes[pid], pEdges[i], pEdges[i+1], suffix.Data(), sigmaExcl, name, help->pCodes[pid]);
                    h_no[pid][i] = new TH1D(name_no, title_no, nBins, xMin, xMax);
                    h_w [pid][i] = new TH1D(name_w, title_w, nBins, xMin, xMax);
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

                    std::vector<TH1D*> cmHists{hcm_no[pid][i], hcm_w[pid][i] };
                    std::vector<std::tuple<Double_t, Double_t>*> fitlims;
                    fitlims.push_back(new std::tuple<Double_t,Double_t>(xMin, xMax));   
                    fitlims.push_back(new std::tuple<Double_t,Double_t>(xMin, xMax)); 
                    cmObj[pid][i] = new covarianceMatrix(cmHists, fitlims);
                }
            }       

            // Fill the inclusive / excluded histograms and covariance matrix objects
            for (Long64_t ev = 0; ev < nEntries; ++ev) {
                chain.GetEntry(ev);
                for (Int_t t = 0; t < NtrkMax; ++t) {
                    Double_t pG = inner[t];

                    // Require a valid TOF momentum estimate because the exclusion variable lives in TOF nSigma 
                    if (tofExpMom[t] < 0) continue;

                    Int_t bin = std::lower_bound(pEdges.begin(), pEdges.end(), pG) - pEdges.begin() - 1;
                    if (bin<0 || bin>=nSteps) continue;

                    for (Int_t pid=0; pid<nParts; ++pid) if (doPid[pid]) {
                        Float_t val = isTPCmode ? tpcNS[pid][t] : tofNS[pid][t];
                        if (TMath::IsNaN(val)) continue;

                        h_no[pid][bin]->Fill(val);

                        // Exclusion: keep track only if it is not within sigmaExcl of the excluded species in TOF.
                        if (TMath::IsNaN(tofNS[excl][t]) || TMath::Abs(tofNS[excl][t])>=sigmaExcl) h_w[pid][bin]->Fill(val);

                        auto nsVals = new std::vector<Double_t>(); 
                        nsVals->reserve(2);
                        auto toadd  = new std::vector<Bool_t>();   
                        toadd->reserve(2);

                        nsVals->clear(); 
                        toadd->clear();

                        // Entry 0 = inclusive histogram, always filled.
                        nsVals->push_back(val);  
                        toadd->push_back(true);

                        // Entry 1 = excluded histogram, filled only if the track passes the TOF nSigma exclusion cut.
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
                c->SetLogy();
                c->Print(pdfName + "[");

                TString  cmPdf;
                TCanvas* cCM = nullptr;
                if (plotCM) {
                    cmPdf = Form("CM_%s_%s_%s-%.2f.pdf",suffix.Data(), help->pCodes[ref], name, sigmaExcl);
                    cCM = new TCanvas(Form("cCM_%s_%s", suffix.Data(), help->pCodes[ref]), "", 800, 700);
                    cCM->SetLeftMargin(0.15);
                    cCM->SetRightMargin(0.15);
                    cCM->Print(cmPdf + "[");
                }

                for (Int_t i = 0; i < nSteps; ++i) {
                    auto cm = cmObj[ref][i];
                    cm->make(useOffDiag, eigenThr);
                    if (plotCM) {
                        cCM->cd();
                        cCM->Clear();
                        cm->plot(cCM, cmPdf);
                        c->cd();
                    }

                    TH1D* h1 = nullptr;
                    TH1D* h2 = nullptr;
                    h1 = h_no[ref][i];
                    h2 = h_w[ref][i];

                    // If no fit overlays are requested, just draw the raw data and continue to the next iteration
                    if (!drawSum && !drawComponents && !drawPredictedLines) {
                        c->cd();
                        c->Clear();
                        h1->Draw("E1");
                        c->Print(pdfName);

                        c->Clear();
                        h2->Draw("E1");
                        c->Print(pdfName);
                        continue;
                    }
                    
                    // Peak seed used to intialize the simulaneous fit
                    struct Peak {Double_t A, mu, sigma; Int_t id; std::vector<Int_t> merged_ids;Bool_t alwaysSeparate = false;};
                    std::vector<Peak> seeds;
                    Double_t sliceMin = pEdges[i];
                    Double_t sliceMax = pEdges[i+1];
                    Double_t pMid = 0.5 * (sliceMin + sliceMax);

                    // Buil expected peak positions from detector response
                    for (Int_t hyp = 0; hyp < nParts; ++hyp) {
                        Double_t hypMass = help->pMasses[hyp];
                        Double_t sigma0, mu;
                        if (isTPCmode) {
                            // For TPC, convert expected dE/dx separation into an nSigma position / width. 
                            Double_t dRef = help->getTPCSignal(pMid * 1000, help->pMasses[ref], 1.0);
                            Double_t dHyp = help->getTPCSignal(pMid * 1000, hypMass, 1.0);
                            
                            if (dRef < 0 || dHyp < 0) continue;

                            Double_t resoHyp = help->resoTPC[hyp][ref];
                            Double_t resoRef = help->resoTPC[ref][hyp];
                            sigma0 = (resoHyp / resoRef) * (dHyp / dRef);
                            mu     = (dHyp/dRef - 1.0) / resoRef;
                        }
                        else {
                            // For TOF, the expected nSigma separation comes from beta differences.
                            Double_t bRef = help->beta(help->pMasses[ref], pMid*1000);
                            Double_t bHyp = help->beta(help->pMasses[hyp], pMid*1000);

                            Double_t resoHyp = help->resoTOF[hyp][ref];
                            Double_t resoRef = help->resoTOF[ref][hyp];
                            sigma0 = (resoHyp / resoRef) * (1.0 / (bHyp * bHyp));
                            mu     = (bRef - bHyp) / (bHyp * bHyp * resoRef);
                        }

                        if (mu < xMin || mu > xMax) continue;

                        Int_t bin = h1->FindBin(mu);        
                        Double_t amp = h1->GetBinContent(bin);

                        seeds.push_back({amp, mu, sigma0, hyp, {hyp}, false});
                    }

                    // Explicitly merge mu and pi into one effective peak if both are present.
                    // This is separate from the generic proximity-based merging below.                    
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
                                Double_t mu_eff = (I_old * muPiCombined.mu + I_new * it->mu) / I_tot;
                                Double_t var_eff = (I_old * (muPiCombined.sigma*muPiCombined.sigma + (muPiCombined.mu-mu_eff)*(muPiCombined.mu-mu_eff)) + I_new *(it->sigma*it->sigma + (it->mu-mu_eff)*(it->mu-mu_eff))) / I_tot;
                                muPiCombined.A = I_tot / (muPiCombined.sigma * TMath::Sqrt(2*TMath::Pi()));
                                muPiCombined.mu = mu_eff;
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

                    // Prevent merging away the reference peak when relevant
                    for (auto &s : seeds) {
                        if (std::find(s.merged_ids.begin(), s.merged_ids.end(), ref) != s.merged_ids.end()) {
                            s.alwaysSeparate = true;
                        }
                    }

                    // Generic merging of nearby predicted peaks
                    std::sort(seeds.begin(), seeds.end(), [](auto &a, auto &b){ return a.mu < b.mu; });
                    std::vector<Peak> merged;
                    for (auto &s : seeds) {
                        if (!merged.empty() && !s.alwaysSeparate && !merged.back().alwaysSeparate && std::abs(s.mu - merged.back().mu) < mergeDistanceFactor * std::max(s.sigma, merged.back().sigma))
                        {
                            Peak &p = merged.back();
                            Double_t I1 = p.A * p.sigma * TMath::Sqrt(2*TMath::Pi());
                            Double_t I2 = s.A * s.sigma * TMath::Sqrt(2*TMath::Pi());
                            Double_t It = I1 + I2;
                            Double_t mu_eff = (I1 * p.mu + I2 * s.mu) / It;
                            Double_t var_eff = (I1 * (p.sigma * p.sigma + (p.mu-mu_eff) * (p.mu-mu_eff)) + I2 * (s.sigma * s.sigma + (s.mu-mu_eff) * (s.mu-mu_eff))) / It;
                            p.sigma = TMath::Sqrt(var_eff);
                            p.mu    = mu_eff;
                            p.A     = It / (TMath::Sqrt(2*TMath::Pi()) * p.sigma);
                            p.merged_ids.insert(p.merged_ids.end(),s.merged_ids.begin(), s.merged_ids.end());
                        } 
                        else {
                            merged.push_back(s);
                        }
                    }

                    // Optional manual override of the peak model
                    Int_t manualNGauss = 0;
                    std::vector<Double_t> manualMeans, manualSigmas, manualAmps;
                    std::vector<TString> manualLabels;
                    if (manualPredictPeaks) {
                        if (interactiveManualPeaks) {
                            c->cd();
                            c->Clear();
                            h1->Draw("E1");
                            c->Modified();
                            c->Update();
                            gSystem->ProcessEvents();

                            std::cout << "\nSlice " << i
                                    << " " << pEdges[i] << " < p < " << pEdges[i+1]
                                    << " ref=" << help->pCodes[ref] << "\n";

                            std::cout << "How many Gaussians? ";
                            std::cin >> manualNGauss;

                            if (manualNGauss > 0) {
                                manualMeans.resize(manualNGauss);
                                manualSigmas.resize(manualNGauss);
                                manualAmps.resize(manualNGauss);
                                manualLabels.resize(manualNGauss);

                                std::cout << "Enter " << manualNGauss << " mean positions: ";
                                for (Int_t ig = 0; ig < manualNGauss; ++ig) std::cin >> manualMeans[ig];

                                std::cout << "Enter " << manualNGauss << " sigma guesses: ";
                                for (Int_t ig = 0; ig < manualNGauss; ++ig) std::cin >> manualSigmas[ig];

                                std::cout << "Enter " << manualNGauss << " amplitude guesses: ";
                                for (Int_t ig = 0; ig < manualNGauss; ++ig) std::cin >> manualAmps[ig];

                                std::cout << "Enter " << manualNGauss << " labels (e, mu, pi, K, p, mu+pi, ...):\n";
                                for (Int_t ig = 0; ig < manualNGauss; ++ig) {
                                    std::cout << "label[" << ig << "] = ";
                                    std::cin >> manualLabels[ig];
                                }
                            } else {
                                manualNGauss = 0;
                            }
                        } else {
                            manualNGauss = getenv_int("MANUAL_NGAUSS", 4);
                            auto defMeans  = std::vector<Double_t>{-6.0, -5.0, 0.5, 5.0};
                            auto defSigmas = std::vector<Double_t>{1.0, 1.0, 1.0, 1.0};
                            auto defAmps   = std::vector<Double_t>{3000, 3000, 300, 10};
                            auto defLabels = std::vector<TString>{"#mu", "#pi", "e", "K"};

                            auto mMeans  = getenv_vecd("MANUAL_MEANS", defMeans);
                            auto mSigmas = getenv_vecd("MANUAL_SIGMAS", defSigmas);
                            auto mAmps   = getenv_vecd("MANUAL_AMPS", defAmps);
                            auto mLabels = getenv_labels("MANUAL_LABELS", defLabels);

                            Int_t L = std::min({manualNGauss, (Int_t)mMeans.size(), (Int_t)mSigmas.size(), (Int_t)mAmps.size()});
                            if (L <= 0) {
                                manualNGauss = 0;
                            } else {
                                manualNGauss = L;
                                manualMeans.assign(mMeans.begin(),  mMeans.begin()  + L);
                                manualSigmas.assign(mSigmas.begin(), mSigmas.begin() + L);
                                manualAmps.assign(mAmps.begin(),    mAmps.begin()    + L);

                                Int_t nLabels = std::min((Int_t)mLabels.size(), L);
                                manualLabels.assign(mLabels.begin(), mLabels.begin() + nLabels);
                            }
                        }
                    }

                    // Final seed list used for the fit: manual if requested and valid, otherwise the automatically merged peaks
                    size_t nG;
                    std::vector<Double_t> seedMeans, seedSigmas, seedAmps;
                    if (manualPredictPeaks && manualNGauss > 0) {
                        nG = manualNGauss;
                        seedMeans  = manualMeans;
                        seedSigmas = manualSigmas;
                        seedAmps   = manualAmps;
                    } else {
                        nG = merged.size();
                        seedMeans.resize(merged.size());
                        seedSigmas.resize(merged.size());
                        seedAmps.resize(merged.size());
                        for (size_t i = 0; i < merged.size(); ++i) {
                            seedMeans[i]  = merged[i].mu;
                            seedSigmas[i] = merged[i].sigma;
                            seedAmps[i]   = merged[i].A;   
                        }
                    }
                    if (nG < 1) { c->Print(pdfName); continue;}

                    // Human-readable lables for the fitted peaks
                    std::vector<TString> peakLabels(nG);
                    for (Int_t ig = 0; ig < (Int_t)nG; ++ig) {
                        if (!manualPredictPeaks) {
                            peakLabels[ig] = help->LegendFromPids(merged[ig].merged_ids, help);
                        } else {
                            if ((Int_t)manualLabels.size() > ig) {
                                auto pids = help->PidsFromLabel(manualLabels[ig]);
                                if (!pids.empty()) peakLabels[ig] = help->LegendFromPids(pids, help);
                                else               peakLabels[ig] = manualLabels[ig];
                            } else {
                                peakLabels[ig] = Form("Peak %d", ig+1);
                            }
                        }
                    }

                    // Optional vertical lines showin predicted/manuel seed positions
                    auto drawPredictedPeakLines = [&](Double_t yMax, TLegend* leg) {
                    if (!drawPredictedLines) return;
                    std::vector<TString> addedPredLabels;

                        if (!manualPredictPeaks) {
                            for (const auto &pk : merged) {
                                Double_t muPred = pk.mu;
                                if (muPred < xMin || muPred > xMax) continue;

                                TLine *l = new TLine(muPred, 0.0, muPred, yMax);
                                Int_t col = kGray+2;
                                TString label = "Pred.";

                                if (!pk.merged_ids.empty()) {
                                    std::vector<Int_t> ids = pk.merged_ids;
                                    std::sort(ids.begin(), ids.end());
                                    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

                                    const Bool_t hasMu = std::find(ids.begin(), ids.end(), 1) != ids.end();
                                    const Bool_t hasPi = std::find(ids.begin(), ids.end(), 2) != ids.end();

                                    if (hasMu && hasPi) col = help->colors[1];
                                    else col = help->colors[ids[0]];

                                    if ((int)ids.size() == 1) {
                                        label = help->pCodes[ids[0]];
                                    } else {
                                        label = "";
                                        for (size_t j = 0; j < ids.size(); ++j) {
                                            if (j) label += " + ";
                                            label += help->pCodes[ids[j]];
                                        }
                                    }
                                }

                                l->SetLineColor(col);
                                l->SetLineStyle(7);
                                l->SetLineWidth(2);
                                l->Draw("same");

                                if (leg && !drawComponents) {
                                    Bool_t alreadyAdded = kFALSE;
                                    for (const auto &lab : addedPredLabels) {
                                        if (lab == label) {
                                            alreadyAdded = kTRUE;
                                            break;
                                        }
                                    }
                                    if (!alreadyAdded) {
                                        leg->AddEntry(l, label, "l");
                                        addedPredLabels.push_back(label);
                                    }
                                }
                            }
                        } else {
                            for (Int_t ig = 0; ig < (Int_t)nG; ++ig) {
                                Double_t muPred = seedMeans[ig];
                                if (muPred < xMin || muPred > xMax) continue;

                                TLine *l = new TLine(muPred, 0.0, muPred, yMax);

                                Int_t col = kGray+2;
                                TString label = Form("Peak %d", ig+1);

                                if (ig < (Int_t)manualLabels.size() && manualLabels[ig].Length()) {
                                    auto pids = help->PidsFromLabel(manualLabels[ig]);

                                    if (!pids.empty()) {
                                        const Bool_t hasMu = std::find(pids.begin(), pids.end(), 1) != pids.end();
                                        const Bool_t hasPi = std::find(pids.begin(), pids.end(), 2) != pids.end();

                                        if (hasMu && hasPi) col = help->colors[1];
                                        else col = help->colors[pids[0]];

                                        label = help->LegendFromPids(pids, help);
                                    } else {
                                        label = manualLabels[ig];
                                    }
                                }

                                l->SetLineColor(col);
                                l->SetLineStyle(7);
                                l->SetLineWidth(2);
                                l->Draw("same");

                                if (leg && !drawComponents) {
                                    Bool_t alreadyAdded = kFALSE;
                                    for (const auto &lab : addedPredLabels) {
                                        if (lab == label) {
                                            alreadyAdded = kTRUE;
                                            break;
                                        }
                                    }
                                    if (!alreadyAdded) {
                                        leg->AddEntry(l, label, "l");
                                        addedPredLabels.push_back(label);
                                    }
                                }
                            }
                        }
                    };

                    // Build simultaneous model:
                    // h1 = sum_i G_i(A1_i, mu_i, sigma_i) + const_bg1
                    // h2 = sum_i G_i(A2_i, mu_i, sigma_i) + const_bg2
                    // Means and sigmas are shared; amplitudes/backgrounds differ between inclusive and excluded spectra.
                    std::ostringstream form;
                    for (size_t i = 0; i < nG; ++i) {
                        if (i) form << "+";
                        form << "gaus(" << 3*i << ")+"
                        << "gaus(" << 3*(nG+i) << ")";
                    }
                    form << "+pol0(" << 4*nG << ")"  // background h1
                        << "+pol0(" << 4*nG+1 << ")"; // background h2

                    TF1* sum = new TF1("sum", form.str().c_str());
                    sum->SetNpx(500); 
                    const Int_t offA1 = 0;
                    const Int_t offA2 = nG;
                    const Int_t offM = 2*nG;
                    const Int_t offS = 3*nG;
                    const Int_t offP1 = 4*nG;
                    const Int_t offP2 = 4*nG + 1;

                    sum->SetParLimits(offP1, 0, 10.0);
                    sum->SetParameter(offP1, 1.0);
                    sum->SetParLimits(offP2, 0, 10.0);
                    sum->SetParameter(offP2, 1.0);

                    // Initialize fit parameters with either the manual or automatically merged peak values
                    for (size_t i = 0; i < nG; ++i) {
                        if (manualPredictPeaks && manualNGauss > 0) {
                            Double_t mu0 = seedMeans[i];
                            Double_t sig0 = seedSigmas[i]; 
                            Double_t amp = seedAmps[i]; 

                            sum->SetParLimits(offA1 + i, 0.0, 1e10);
                            sum->SetParameter(offA1 + i, amp);

                            sum->SetParLimits(offA2 + i, 0.0, 1e10);
                            sum->SetParameter (offA2 + i, amp);

                            sum->SetParLimits(offM + i, mu0 - muWindow, mu0 + muWindow);
                            sum->SetParameter(offM + i, mu0);

                            sum->SetParLimits(offS + i, sig0 * 0.5, sig0 * 2.0);
                            sum->SetParameter(offS + i, sig0);

                            Double_t lo, hi;
                            sum->GetParLimits(offS, lo, hi);
                        }
                        else {
                            const auto &p = merged[i];

                            sum->SetParLimits(offA1 + i, 0.0, 1e10);
                            sum->SetParameter(offA1 + i, p.A);

                            sum->SetParLimits(offA2 + i, 0.0, 1e10);
                            sum->SetParameter (offA2 + i, p.A);

                            sum->SetParLimits(offM + i, p.mu - muWindow, p.mu + muWindow);
                            sum->SetParameter(offM + i, p.mu);

                            sum->SetParLimits(offS + i, p.sigma * 0.5, p.sigma * 2.0);
                            sum->SetParameter(offS + i, p.sigma);
                        }
                    }
                    help->FitHistogramsExclCompByChi2(h1, h2, sum, nG, cmObj[ref][i], std::vector<TH1D*>{hcm_no[ref][i], hcm_w[ref][i]}, eigenThr, nBins);
                    
                    c->Clear();
                    h1->Draw("E1");
                    Double_t yMax1 = 1.25 * h1->GetMaximum();
                    
                    TLegend *leg1 = nullptr;
                    if (drawSum || drawComponents || drawPredictedLines) {
                        leg1 = new TLegend(0.80, 0.70, 0.90, 0.90);
                        leg1->SetBorderSize(1);
                        leg1->SetMargin(0.45);
                        leg1->SetFillColorAlpha(kWhite, 0.8);
                        leg1->SetTextSize(0.030);
                        leg1->SetTextFont(42);
                    }
                    drawPredictedPeakLines(yMax1, leg1);

                    std::vector<Double_t> par(sum->GetNpar());
                    std::vector<Double_t> err(sum->GetNpar());
                    for(Int_t i = 0; i < sum->GetNpar(); ++i)  {
                        par[i] = sum->GetParameter(i);
                        err[i] = sum->GetParError(i);
                    }                    

                    // Compute yields and contamination in ±nSigma windows around each fitted peak
                    auto computeWindowScan = [&](Int_t whichHist) -> WindowScan {
                        const Int_t offA = (whichHist == 0 ? offA1 : offA2);
                        const Int_t offP = (whichHist == 0 ? offP1 : offP2);

                        WindowScan out;

                        struct WindowVals {
                            Double_t target = std::numeric_limits<Double_t>::quiet_NaN();
                            Double_t other = std::numeric_limits<Double_t>::quiet_NaN();
                            Double_t bg = std::numeric_limits<Double_t>::quiet_NaN();
                            Double_t total = std::numeric_limits<Double_t>::quiet_NaN();
                            Double_t cont = std::numeric_limits<Double_t>::quiet_NaN();
                        };

                        auto eval_window_vals = [&](const std::vector<Double_t>& pvar, Double_t nsigma, Int_t ii) -> WindowVals {
                            WindowVals outv;

                            const Double_t Ai = pvar[offA + ii];
                            const Double_t mui  = pvar[offM + ii];
                            const Double_t sigi = pvar[offS + ii];

                            if (!(Ai > 0.0) || !(sigi > 0.0)) return outv;

                            const Double_t a = mui - nsigma * sigi;
                            const Double_t b = mui + nsigma * sigi;

                            Double_t I_target = 0.0;
                            Double_t I_other = 0.0;

                            for (Int_t j = 0; j < nG; ++j) {
                                const Double_t Aj = pvar[offA + j];
                                const Double_t muj = pvar[offM + j];
                                const Double_t sigj = pvar[offS + j];

                                if (!(Aj > 0.0) || !(sigj > 0.0)) continue;

                                const Double_t Ij = help->GaussIntegral(Aj, muj, sigj, a, b);

                                if (j == ii) I_target += Ij;
                                else I_other += Ij;
                            }

                            const Double_t I_bg  = pvar[offP] * (b - a);
                            const Double_t I_tot = I_target + I_other + I_bg;

                            outv.target = I_target;
                            outv.other = I_other;
                            outv.bg = I_bg;
                            outv.total = I_tot;
                            outv.cont = (I_tot > 0.0) ? (I_other + I_bg) / I_tot : std::numeric_limits<Double_t>::quiet_NaN();
                            return outv;
                        };

                        // Numerical derivativ step size for error propagation 
                        auto step_for = [&](Int_t ipar, Double_t pval, Double_t perr) -> Double_t {
                            Double_t h = 0.0;
                            if (perr > 0.0 && std::isfinite(perr)) h = 0.25 * perr;
                            else h = std::max(1e-6, 1e-4 * std::abs(pval));

                            // avoid zero step
                            h = std::max(h, 1e-8);

                            // keep sigma variations positive
                            if (ipar >= offS && ipar < offS + nG) {
                                h = std::min(h, 0.45 * std::max(1e-8, pval));
                            }
                            return h;
                        };

                        for (Double_t nsigma : sigmaScan) {
                            std::vector<Double_t> target_yield(nG, std::numeric_limits<Double_t>::quiet_NaN());
                            std::vector<Double_t> target_yield_err(nG, std::numeric_limits<Double_t>::quiet_NaN());

                            std::vector<Double_t> other_yield(nG, std::numeric_limits<Double_t>::quiet_NaN());
                            std::vector<Double_t> other_yield_err(nG, std::numeric_limits<Double_t>::quiet_NaN());

                            std::vector<Double_t> bg_yield(nG, std::numeric_limits<Double_t>::quiet_NaN());
                            std::vector<Double_t> bg_yield_err(nG, std::numeric_limits<Double_t>::quiet_NaN());

                            std::vector<Double_t> total_yield(nG, std::numeric_limits<Double_t>::quiet_NaN());
                            std::vector<Double_t> total_yield_err(nG, std::numeric_limits<Double_t>::quiet_NaN());

                            std::vector<Double_t> contamination(nG, std::numeric_limits<Double_t>::quiet_NaN());
                            std::vector<Double_t> contamination_err(nG, std::numeric_limits<Double_t>::quiet_NaN());

                            for (Int_t ii = 0; ii < nG; ++ii) {
                                const WindowVals base = eval_window_vals(par, nsigma, ii);

                                target_yield[ii] = base.target;
                                other_yield[ii] = base.other;
                                bg_yield[ii] = base.bg;
                                total_yield[ii] = base.total;
                                contamination[ii] = base.cont;

                                if (!std::isfinite(base.total)) continue;

                                Double_t var_target = 0.0;
                                Double_t var_other = 0.0;
                                Double_t var_bg = 0.0;
                                Double_t var_total = 0.0;
                                Double_t var_cont = 0.0;

                                // Active parameters for one histogram:
                                // amplitudes of that histogram + shared means/sigmas + its own background
                                std::vector<Int_t> activePars;
                                for (Int_t j = 0; j < nG; ++j) activePars.push_back(offA + j);
                                for (Int_t j = 0; j < nG; ++j) activePars.push_back(offM + j);
                                for (Int_t j = 0; j < nG; ++j) activePars.push_back(offS + j);
                                activePars.push_back(offP);

                                for (Int_t ipar : activePars) {
                                    const Double_t pk = par[ipar];
                                    const Double_t ek = err[ipar];

                                    if (!(ek > 0.0) || !std::isfinite(ek)) continue;

                                    const Double_t h = step_for(ipar, pk, ek);
                                    if (!(h > 0.0)) continue;

                                    std::vector<Double_t> pplus  = par;
                                    std::vector<Double_t> pminus = par;

                                    pplus[ipar] = pk + h;
                                    pminus[ipar] = pk - h;

                                    // guard against nonphysical sigma / negative background in derivative eval
                                    if (ipar >= offS && ipar < offS + nG) {
                                        pplus[ipar] = std::max(pplus[ipar],  1e-8);
                                        pminus[ipar] = std::max(pminus[ipar], 1e-8);
                                    }

                                    const WindowVals vp = eval_window_vals(pplus,  nsigma, ii);
                                    const WindowVals vm = eval_window_vals(pminus, nsigma, ii);

                                    auto deriv = [&](Double_t fp, Double_t fm) -> Double_t {
                                        if (!std::isfinite(fp) || !std::isfinite(fm)) return 0.0;
                                        return (fp - fm) / (2.0 * h);
                                    };

                                    const Double_t d_target = deriv(vp.target, vm.target);
                                    const Double_t d_other = deriv(vp.other,  vm.other);
                                    const Double_t d_bg = deriv(vp.bg,     vm.bg);
                                    const Double_t d_total = deriv(vp.total,  vm.total);
                                    const Double_t d_cont = deriv(vp.cont,   vm.cont);

                                    var_target += d_target * d_target * ek * ek;
                                    var_other += d_other * d_other * ek * ek;
                                    var_bg += d_bg * d_bg * ek * ek;
                                    var_total += d_total * d_total * ek * ek;
                                    var_cont += d_cont * d_cont * ek * ek;
                                }

                                target_yield_err[ii] = std::sqrt(std::max(0.0, var_target));
                                other_yield_err[ii] = std::sqrt(std::max(0.0, var_other));
                                bg_yield_err[ii] = std::sqrt(std::max(0.0, var_bg));
                                total_yield_err[ii] = std::sqrt(std::max(0.0, var_total));
                                contamination_err[ii] = std::sqrt(std::max(0.0, var_cont));
                            }

                            out.push_back(WindowScanEntry{
                                nsigma,
                                std::move(target_yield),
                                std::move(target_yield_err),
                                std::move(other_yield),
                                std::move(other_yield_err),
                                std::move(bg_yield),
                                std::move(bg_yield_err),
                                std::move(total_yield),
                                std::move(total_yield_err),
                                std::move(contamination),
                                std::move(contamination_err)
                            });
                        }

                        return out;
                    };

                    // Diagnostic of the final fit
                    auto [D1,N1] = help->PoissonDeviance(h1, 0, par, nG, offA1, offA2, offM, offM, offS, offS, offP1, offP2);
                    auto [D2,N2] = help->PoissonDeviance(h2, 1, par, nG, offA1, offA2, offM, offM, offS, offS, offP1, offP2);
                    const Double_t D = D1 + D2;
                    const Double_t N = N1 + N2;
                    const Int_t k = 4 * nG + 2; 

                    const Double_t chi2 = sum->GetChisquare();
                    const Int_t ndf  = sum->GetNDF();
                    const Double_t pval = TMath::Prob(chi2, ndf);

                    // Total Gaussian areas in the two histograms
                    std::vector<Double_t> areas_noK(nG), areas_wK(nG), err_areas_noK(nG), err_areas_wK(nG);
                    const Double_t factor = TMath::Sqrt(2*TMath::Pi());
                    for (Int_t ig = 0; ig < nG; ++ig) {
                        Double_t A1 = par[offA1 + ig];
                        Double_t A2 = par[offA2 + ig];
                        Double_t sig = par[offS + ig];
                        Double_t eA1 = err[offA1 + ig];
                        Double_t eA2 = err[offA2 + ig];
                        Double_t eS = err[offS + ig];

                        areas_noK[ig] = A1 * sig * factor;
                        areas_wK [ig] = A2 * sig * factor;

                        Double_t dFdA = sig * factor;
                        Double_t dFdS = A1 * factor;        
                        Double_t var_noK = dFdA*dFdA * eA1*eA1 + dFdS*dFdS * eS*eS; 
                        err_areas_noK[ig] = TMath::Sqrt(var_noK);

                        dFdS = A2 * factor;
                        Double_t var_wK = dFdA*dFdA * eA2*eA2 + dFdS*dFdS * eS*eS;
                        err_areas_wK[ig] = TMath::Sqrt(var_wK);
                    }
                    const Double_t D_over_N = (N > 0 ? D / N : std::numeric_limits<Double_t>::quiet_NaN());
                    const Double_t chi2_over_ndf = (ndf > 0 ? chi2 / ndf : std::numeric_limits<Double_t>::quiet_NaN());

                    std::vector<Double_t> fitMeans(nG), fitSigmas(nG);
                    std::vector<Double_t> fitAmps_noExcl(nG), fitAmps_wExcl(nG);
                    std::vector<Double_t> e_fitMeans(nG), e_fitSigmas(nG);
                    std::vector<Double_t> e_fitAmps_noExcl(nG), e_fitAmps_wExcl(nG);

                    for (Int_t ig = 0; ig < nG; ++ig) {
                        fitAmps_noExcl[ig] = par[offA1 + ig];
                        fitAmps_wExcl[ig] = par[offA2 + ig];
                        fitMeans[ig] = par[offM + ig];
                        fitSigmas[ig] = par[offS + ig];

                        e_fitAmps_noExcl[ig] = err[offA1 + ig];
                        e_fitAmps_wExcl[ig] = err[offA2 + ig];
                        e_fitMeans[ig] = err[offM + ig];
                        e_fitSigmas[ig] = err[offS + ig];
                    }
                    const Double_t bkg_noExcl = par[offP1];   
                    const Double_t bkg_wExcl = par[offP2];
                    const Double_t ebkg_noExcl = err[offP1];
                    const Double_t ebkg_wExcl = err[offP2];

                    auto scan_no = computeWindowScan(0); 
                    auto scan_w = computeWindowScan(1); 

                    // optional logging of the fit results and diagnostics in txt format for later analysis
                    if (writeSummaryLog && summaryLog.is_open()) {
                        summaryLog
                            << "mode: " << suffix.Data() << "\n"
                            << "tag: noExcl\n"
                            << "ref: " << help->pNames[ref] << "\n"
                            << "slice: " << i << "\n"
                            << "pmin: " << pEdges[i] << "\n"
                            << "pmax: " << pEdges[i+1] << "\n"
                            << "chi2_ndf: " << chi2_over_ndf << "\n"
                            << "labels: " << labels_to_txt(peakLabels) << "\n\n"
                            << "fit_table:\n"
                            << species_table_to_txt(peakLabels,
                                                    fitMeans, e_fitMeans,
                                                    fitSigmas, e_fitSigmas,
                                                    fitAmps_noExcl, e_fitAmps_noExcl,
                                                    areas_noK, err_areas_noK)
                            << "\n"
                            << "background: " << bkg_noExcl << " ± " << ebkg_noExcl << "\n\n"
                            << "window_yield_scan:\n"
                            << scan_table_to_txt(peakLabels, scan_no)
                            << "----------------------------------------\n\n";

                        summaryLog
                            << "mode: " << suffix.Data() << "\n"
                            << "tag: wExcl\n"
                            << "KaonExcl: " << (FitKaonExclComp ? "yes" : "no") << "\n"
                            << "ProtonExcl: " << (FitProtonExclComp ? "yes" : "no") << "\n"
                            << "ExclSigma: " << sigmaExcl << "\n"
                            << "ref: " << help->pNames[ref] << "\n"
                            << "slice: " << i << "\n"
                            << "pmin: " << pEdges[i] << "\n"
                            << "pmax: " << pEdges[i+1] << "\n"
                            << "chi2_ndf: " << chi2_over_ndf << "\n"
                            << "labels: " << labels_to_txt(peakLabels) << "\n\n"
                            << "fit_table:\n"
                            << species_table_to_txt(peakLabels,
                                                    fitMeans, e_fitMeans,
                                                    fitSigmas, e_fitSigmas,
                                                    fitAmps_wExcl, e_fitAmps_wExcl,
                                                    areas_wK, err_areas_wK)
                            << "\n"
                            << "background: " << bkg_wExcl << " ± " << ebkg_wExcl << "\n\n"
                            << "window_yield_scan:\n"
                            << scan_table_to_txt(peakLabels, scan_w)
                            << "----------------------------------------\n\n";
                    }

                    // Build curves for the total fit and each individual component
                    const Int_t nPoints = 500; 
                    Double_t dx = (xMax - xMin)/(nPoints-1);
                    std::vector<Double_t> xv(nPoints), y1(nPoints), y2(nPoints);

                    const double bw1 = h1->GetXaxis()->GetBinWidth(1);
                    const double bw2 = h2->GetXaxis()->GetBinWidth(1);

                    for(Int_t ip = 0; ip < nPoints; ++ip) {
                        Double_t x = xMin + ip * dx;
                        xv[ip] = x;
                        Double_t yy1 = par[offP1] * bw1;
                        Double_t yy2 = par[offP2] * bw2;
                        for(Int_t ig=0; ig<nG; ++ig) {
                            Double_t A1 = par[offA1 + ig];
                            Double_t A2 = par[offA2 + ig];
                            Double_t mu = par[offM + ig];
                            Double_t sig = par[offS + ig];
                            yy1 += A1 * TMath::Gaus(x, mu, sig, kFALSE) * bw1;
                            yy2 += A2 * TMath::Gaus(x, mu, sig, kFALSE) * bw2;
                        }
                        y1[ip] = yy1;
                        y2[ip] = yy2;
                    }

                    TGraph* gSum1 = nullptr;
                    if (drawSum) {
                        gSum1 = new TGraph(nPoints, xv.data(), y1.data());
                        gSum1->SetLineColor(kRed);
                        gSum1->SetLineWidth(2);
                        gSum1->Draw("L SAME");
                        if (leg1) leg1->AddEntry(gSum1, "Sum", "l");
                    }

                    TGraph* gSum2 = nullptr;
                    if (drawSum) {
                        gSum2 = new TGraph(nPoints, xv.data(), y2.data());
                        gSum2->SetLineColor(kRed);
                        gSum2->SetLineWidth(2);
                    }

                        
                    for (Int_t i = 0; i < nG; ++i) {
                        if (sum->GetParameter(offA1 + i) < 0) continue;

                        TF1 *g1 = new TF1(Form("g1_%d_%d", ref, i), "gaus", xMin, xMax);
                        g1->SetParameters(
                            sum->GetParameter(offA1 + i) * bw1,
                            sum->GetParameter(offM + i),
                            sum->GetParameter(offS + i)
                        );
                        TString label;
                        Int_t col = kGray+2;
                        if (!manualPredictPeaks) {
                            // automatic mode: use merged IDs
                            col = (merged[i].id >= 0) ? help->colors[merged[i].id] : kGray+2;
                            g1->SetLineColor(col);

                            auto ids = merged[i].merged_ids;
                            std::sort(ids.begin(), ids.end());
                            ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

                            if ((int)ids.size() == 1) {
                                label = help->pCodes[ids[0]];
                            } else {
                                for (size_t j = 0; j < ids.size(); ++j) {
                                    if (j) label += " + ";
                                    label += help->pCodes[ids[j]];
                                }
                            }
                        } else {
                            TString raw = ((Int_t)manualLabels.size() > i) ? manualLabels[i] : "";

                            auto pids = help->PidsFromLabel(raw);

                            const Bool_t hasMu = (std::find(pids.begin(), pids.end(), 1) != pids.end());
                            const Bool_t hasPi = (std::find(pids.begin(), pids.end(), 2) != pids.end());

                            Int_t colorPid = -1;
                            if (hasMu && hasPi) colorPid = 1;
                            else if ((Int_t)pids.size() == 1) colorPid = pids[0];
                            else if (!pids.empty()) colorPid = pids[0];

                            col = (colorPid >= 0) ? help->colors[colorPid] : kGray+2;
                            label = !pids.empty() ? help->LegendFromPids(pids, help) : (raw.Length() ? raw : Form("Peak %d", i+1));
                            g1->SetLineColor(col);
                        }
                        g1->SetLineStyle(2);
                        g1->SetNpx(500);

                        if (drawComponents) {
                            g1->Draw("same");
                            if (leg1) leg1->AddEntry(g1, label, "l");
                        }

                    }
                    if(leg1) leg1->Draw();
                    c->Print(pdfName);
                    if(leg1) delete leg1;

                    c->Clear();
                    h2->Draw("E1");
                    if (drawSum && gSum2) gSum2->Draw("L SAME");
                    Double_t yMax2 = 1.25 * h2->GetMaximum();

                    TLegend *leg2 = nullptr;
                    if (drawSum || drawComponents || drawPredictedLines) {
                        leg2 = new TLegend(0.80, 0.70, 0.90, 0.90);
                        leg2->SetBorderSize(1); 
                        leg2->SetMargin(0.45);
                        leg2->SetFillColorAlpha(kWhite, 0.8);
                        leg2->SetTextSize(0.030);
                        leg2->SetTextFont(42);
                        if (drawSum && gSum2) leg2->AddEntry(gSum2, "Sum", "l");
                    }
                    drawPredictedPeakLines(yMax2, leg2);

                    for (Int_t i = 0; i < nG; ++i) {
                        if (sum->GetParameter(offA2 + i) < 0) continue;
                        TF1 *g2 = new TF1(Form("g2_%d_%d", ref, i), "gaus", xMin, xMax);
                        g2->SetParameters(
                            sum->GetParameter(offA2 + i) * bw2,
                            sum->GetParameter(offM + i),
                            sum->GetParameter(offS + i)
                        );

                        TString label;
                        Int_t col = kGray+2;
                        if (!manualPredictPeaks) {
                            // automatic mode: use merged IDs
                            col = (merged[i].id >= 0) ? help->colors[merged[i].id] : kGray+2;
                            g2->SetLineColor(col);

                            auto ids = merged[i].merged_ids;
                            std::sort(ids.begin(), ids.end());
                            ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

                            if ((int)ids.size() == 1) {
                                label = help->pCodes[ids[0]];
                            } else {
                                for (size_t j = 0; j < ids.size(); ++j) {
                                    if (j) label += " + ";
                                    label += help->pCodes[ids[j]];
                                }
                            }
                        } else {
                            TString raw = ((Int_t)manualLabels.size() > i) ? manualLabels[i] : "";

                            auto pids = help->PidsFromLabel(raw);

                            const Bool_t hasMu = (std::find(pids.begin(), pids.end(), 1) != pids.end());
                            const Bool_t hasPi = (std::find(pids.begin(), pids.end(), 2) != pids.end());

                            Int_t colorPid = -1;
                            if (hasMu && hasPi) colorPid = 1;
                            else if ((Int_t)pids.size() == 1) colorPid = pids[0];
                            else if (!pids.empty()) colorPid = pids[0];

                            col = (colorPid >= 0) ? help->colors[colorPid] : kGray+2;
                            label = !pids.empty() ? help->LegendFromPids(pids, help): (raw.Length() ? raw : Form("Peak %d", i+1));
                            g2->SetLineColor(col);
                        }
                        g2->SetLineStyle(2);
                        g2->SetNpx(500);

                        if (drawComponents) {
                            g2->Draw("same");
                            if (leg2) leg2->AddEntry(g2, label, "l");
                        }
                    }      
                    if (leg2) leg2->Draw();   
                    c->Print(pdfName);
                    if (leg2) delete leg2;
                    delete sum;
                }
                c->Print(pdfName+"]"); 
                if (plotCM) {
                    cCM->Print(cmPdf + "]"); 
                }
            }
        }
    };
    // Run over the requested sigma exclusion values
    for (Double_t s : sigmaExclList) {
        if (plotTPC) drawNSigma(true, s);
        if (plotTOF) drawNSigma(false, s);
    } 
}

#if defined(__CLING__)
  #pragma cling load("libHist")
  #pragma cling load("libCore")
  #pragma cling load("libRIO")
  #pragma cling load("libTree")
  #pragma cling load("libTreePlayer")
  #pragma cling load("libPhysics")
  #pragma link C++ class std::vector<float>+;
#endif

#include <TChain.h>
#include <TFile.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TIterator.h>
#include <TString.h>
#include <TMath.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <TSpectrum.h>  


void AddTrees(TChain &chain, const char* baseDir) {
    TSystemDirectory dir("base", baseDir);
    TList *subdirs = dir.GetListOfFiles();
    TSystemFile *sysfile;
    TIterator *itSub = subdirs->MakeIterator();
    while ((sysfile = (TSystemFile*)itSub->Next())) {
        TString dname = sysfile->GetName();
        if (!sysfile->IsDirectory() || !dname.BeginsWith("hy_")) continue;
        TString fullDir = TString(baseDir) + "/" + dname;

        TSystemDirectory subdir(dname, fullDir);
        TList *files = subdir.GetListOfFiles();
        TSystemFile *f2;
        TIterator *itFile = files->MakeIterator();
        while ((f2 = (TSystemFile*)itFile->Next())) {
            TString fname = f2->GetName();
            if (!fname.BeginsWith("RLAnalysisTree") || !fname.EndsWith(".root")) continue;
            TString path = fullDir + "/" + fname;

            TFile tf(path, "READ");
            if (tf.IsZombie()) { tf.Close(); continue; }
            TIterator *itKey = tf.GetListOfKeys()->MakeIterator();
            TKey *key;
            while ((key = (TKey*)itKey->Next())) {
                TString keyName = key->GetName();
                if (keyName.BeginsWith("DF_")) {
                    chain.Add(path + "/" + keyName + "/O2tautwotrack");
                    break;
                }
            }
            delete itKey;
            tf.Close();
        }
        delete itFile;
    }
    delete itSub;
}

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

void DrawResidual(){
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1);
    const char* baseDir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";
    TChain chain("twotauchain");
    AddTrees(chain, baseDir);

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);
    chain.SetBranchStatus("fTrkTPCsignal", 1);
    Float_t inner[2], signal[2];
    chain.SetBranchAddress("fTrkTPCinnerParam", inner);
    chain.SetBranchAddress("fTrkTPCsignal", signal);

    const Int_t nBins = 200;
    const Double_t pMin = 0.45, pMax = 0.55; 

    const Int_t nParts = 5;
    const TString names[nParts]   = {"e","mu","pi","K","p"};
    const Double_t masses[nParts] = {0.51099895, 105.6583755,  139.57039, 493.677, 938.27208816};
    const Double_t charges[nParts]= {1,1,1,1,1};
    const Int_t colors[nParts]    = {kBlue, kRed, kGreen+2, kOrange+7, kViolet};

    TH1F *hRes = new TH1F("hRes", Form("Residuals for %.2f < p < %.2f GeV/c;(dE/dx_{meas}-dE/dx_{calc})/dE/dx_{calc};Counts",pMin, pMax),nBins, -1.0, +1.0);
    Long64_t nEntries = std::min(chain.GetEntries(), (Long64_t)1e6);
    for(Long64_t i=0; i<nEntries; ++i){
      chain.GetEntry(i);
      for(int j=0; j<2; ++j){
        Float_t pG = inner[j]; 
        if(pG < pMin || pG > pMax)   continue;

        Float_t s  = signal[j]; 

        Double_t pMeV = pG * 1000.;

        for(int h=0; h < nParts; ++h){
          Double_t exp = get_expected_signal(pMeV, masses[h], charges[h]);
          if(exp < 0) continue;
          Double_t R = (s - exp) / exp;
          hRes->Fill(R);
        }
      }
    }

    TSpectrum spec(5);
    Int_t nFound = spec.Search(hRes, 1, "nodraw goff", 0.05);
    Double_t* xPeaks = spec.GetPositionX();
    std::vector<Double_t> peakPositions(xPeaks, xPeaks + nFound);
    std::sort(peakPositions.begin(), peakPositions.end());

    struct Window { const char* name; Double_t lo, hi; Int_t col; };
    std::vector<Window> fitWindows = {
        {"kaon",    -0.9, -0.6, kBlue   },
        {"pion",    -0.6, -0.2, kRed    },
        {"electron",-0.2,  0.4, kMagenta}
    };

    for (auto &w : fitWindows) {
        Int_t binLo = hRes->GetXaxis()->FindBin(w.lo);
        Int_t binHi = hRes->GetXaxis()->FindBin(w.hi);

        Double_t amp = 0;
        Int_t    iMax = binLo;
        for (Int_t b = binLo; b <= binHi; ++b) {
            Double_t c = hRes->GetBinContent(b);
            if (c > amp) { amp = c; iMax = b; }
        }
        Double_t mean  = hRes->GetBinCenter(iMax);
        Double_t sigma = (w.hi - w.lo) / 4.;  

        TF1 *f = new TF1(Form("f_%s", w.name),
                        "gaus", w.lo, w.hi);
        f->SetLineColor(w.col);
        f->SetLineWidth(2);
        f->SetParameters(amp, mean, sigma);

        hRes->Fit(f, "R+");  // R = restrict to [lo,hi], + = keep earlier fits
        Double_t name = f->GetParameter(0);  // update parameters after fit
    }

    TCanvas *c = new TCanvas("c", "Residuals", 800, 600);
    c->SetLeftMargin(0.15);
    c->SetGrid();
    hRes->SetMarkerStyle(kFullCircle);
    hRes->SetMarkerColor(kBlack);
    hRes->SetLineColor(kRed);
    hRes->SetMarkerSize(0.75);
    hRes->Draw("E1"); 
    
    c->Print("Residuals_of_Bethe_Bloch.pdf");
}
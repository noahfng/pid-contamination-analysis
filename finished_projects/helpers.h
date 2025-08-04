#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

class helper{
 public:

  // constants
  const Char_t* base_dir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";
  const Int_t NtrkMax = 2;

  enum { nParts = 5 };
  const Int_t colors[nParts]  = {kBlue, kGreen+2, kOrange+7, kMagenta+2, kCyan+1};
  const Char_t* dNames[2] = {"TPC", "TOF"};
  const Char_t* pNames[nParts]   = {"El", "Mu", "Pi", "Ka", "Pr"};
  const Char_t* pCodes[nParts]   = {"e", "#mu", "#pi", "K", "p"};
  const Double_t pMasses[nParts] = {0.51099895, 105.6583755,  139.57039, 493.677, 938.27208816};
  const Double_t pCharges[nParts] = {1, 1, 1, 1, 1};

  const Double_t resoTPC[nParts] = {0.085, 0.072, 0.074, 0.09, 0.08}; 
  const Double_t resoTOF[nParts]   = {0.013, 0.013, 0.013, 0.019, 0.020};

  // functions
  // relativ velocity beta
  Double_t beta(Float_t mass, Float_t mom)
  {
    return (mom / TMath::Sqrt(mass*mass + mom*mom));
  };

  // compute expected TPC signal at given particle momentum, mass and charge
  Double_t getTPCSignal(Double_t p, Double_t mass, Double_t charge) {
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
  };
  
  Double_t getnSigma(Double_t mom, TString det, Int_t ref, Int_t hyp)
  {
    Double_t val = 0.;
    if (det == dNames[0])
    {
      auto dRef = getTPCSignal(mom, pMasses[ref], pCharges[ref]);
      auto dHyp = getTPCSignal(mom, pMasses[hyp], pCharges[hyp]);
      auto rRef = resoTPC[hyp];
      val = (dHyp/dRef - 1.0) / rRef;
    } else {
      auto bRef = beta(pMasses[ref], mom);
      auto bHyp = beta(pMasses[hyp], mom);
      auto rRef = resoTOF[hyp];
      val = (bRef - bHyp) / (bHyp*bHyp * rRef);
    }
    return val;
  }

  enum ResoMode { kTPC, kTOF };

  Double_t getReso(ResoMode mode, const Char_t* hypo, Float_t mom) {
      // open the same file just once
      static TFile* file = TFile::Open(
          "UD_LHC23_pass4_SingleGap/nSigma_Resolution.root", "READ"
      );
      if (!file || file->IsZombie()) {
          std::cerr << "Error opening resolution file\n";
          return -1;
      }

      // build name: "nSigmaTPCresEl" or "nSigmaTOFresPr", etc.
      const Char_t* prefix = (mode == kTPC ? "nSigmaTPCres" : "nSigmaTOFres");
      TString gname = Form("%s%s", prefix, hypo);

      auto gr = dynamic_cast<TGraph*>(file->Get(gname));
      if (!gr) {
          std::cerr << "Graph " << gname.Data() << " not found\n";
          return -1;
      }

      return gr->Eval(mom);
  };

  inline void FitHistogramByChi2(TH1* hist, TF1* func, Double_t xlo, Double_t xhi) {
    func->SetRange(xlo, xhi);
    std::unique_ptr<ROOT::Math::Minimizer> minimizer(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", ""));
    const Int_t nPar = func->GetNpar();
    auto chi2_fcn = [&](const Double_t *par) {
        for (Int_t i = 0; i < nPar; ++i) {
            func->SetParameter(i, par[i]);
        }
        Double_t chi2 = 0;
        Int_t usedBins = 0;
        Int_t nbins = hist->GetNbinsX();
        for (Int_t ib = 1; ib <= nbins; ++ib) {
            Double_t x   = hist->GetBinCenter(ib);
            if (x < xlo || x > xhi) continue;
            Double_t y   = hist->GetBinContent(ib);
            Double_t err = hist->GetBinError(ib);
            if (err <= 0) continue;
            Double_t yfit = func->Eval(x);
            chi2 += (y - yfit)*(y - yfit)/(err*err);
            ++usedBins;
        }
        return chi2;
    };

    ROOT::Math::Functor fcn(chi2_fcn, nPar);
    minimizer->SetFunction(fcn);

    for (Int_t i = 0; i < nPar; ++i) {
        const Char_t* name = func->GetParName(i);
        Double_t      val  = func->GetParameter(i);
        Double_t      err  = func->GetParError(i);
        Double_t      step = (err > 0 ? err : (fabs(val)*0.1 + 1e-3));

        Double_t low = 0, up = 0;
        func->GetParLimits(i, low, up);

        if (up > low) {
            minimizer->SetLimitedVariable(i, name, val, step, low, up);
        } else {
            minimizer->SetVariable(i, name, val, step);
        }
    }

    minimizer->Minimize();

    const Double_t* xs = minimizer->X();
    for (Int_t i = 0; i < nPar; ++i) {
        func->SetParameter(i, xs[i]);
    }

    Double_t chi2 = 0;
    Int_t usedBins = 0;
    Int_t nbins = hist->GetNbinsX();
    for (Int_t ib = 1; ib <= nbins; ++ib) {
        Double_t x   = hist->GetBinCenter(ib);
        if (x < xlo || x > xhi) continue;
        Double_t y   = hist->GetBinContent(ib);
        Double_t err = hist->GetBinError(ib);
        if (err <= 0) continue;
        Double_t yfit = func->Eval(x);
        chi2 += (y - yfit)*(y - yfit)/(err*err);
        ++usedBins;
    }
    Int_t ndf = usedBins - nPar;
    func->SetChisquare(chi2);
    func->SetNDF(ndf);

}
 private:
  // bethe-bloch according to ALICE parametrisation 
  Double_t bethe_bloch_aleph(Double_t bg, Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5) {
    Double_t beta = bg / TMath::Sqrt(1.0 + bg*bg);
    Double_t aa   = TMath::Power(beta, p4);
    Double_t bb   = TMath::Log(p3 + TMath::Power(1.0/bg, p5));
    return (p2 - aa - bb) * p1 / aa;
  };

};

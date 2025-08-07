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
  const Char_t* base_dir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";
  const Int_t NtrkMax = 2;

  enum { nParts = 5 };
  const Int_t colors[nParts]  = {kBlue, kGreen+2, kOrange+7, kMagenta+2, kCyan+1};
  const Char_t* dNames[2] = {"TPC", "TOF"};
  const Char_t* pNames[nParts]   = {"El", "Mu", "Pi", "Ka", "Pr"};
  const Char_t* pCodes[nParts]   = {"e", "#mu", "#pi", "K", "p"};
  const Double_t pMasses[nParts] = {0.51099895, 105.6583755,  139.57039, 493.677, 938.27208816};
  const Double_t pCharges[nParts] = {1, 1, 1, 1, 1};

  const Double_t resoTPC[nParts][nParts] = {
    {1.000, 0.080, 0.080, 0.085, 0.075},
    {0.090, 1.000, 0.080, 0.090, 0.085},
    {0.090, 0.080, 1.000, 0.090, 0.085},
    {0.090, 0.080, 0.080, 1.000, 0.085},
    {0.090, 0.080, 0.080, 0.090, 1.000}};

  const Double_t resoTOF[nParts][nParts] = {
    {1.0000, 0.009, 0.009, 0.019, 0.021},
    {0.0075, 1.000, 0.009, 0.019, 0.021},
    {0.0075, 0.009, 1.000, 0.019, 0.021},
    {0.0073, 0.012, 0.012, 1.000, 0.021},
    {0.0073, 0.012, 0.012, 0.014, 1.000}};

  Double_t beta(Float_t mass, Float_t mom)
  {
    return (mom / TMath::Sqrt(mass*mass + mom*mom));
  };

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
      auto rRef = resoTPC[ref][hyp];
      val = (dHyp/dRef - 1.0) / rRef;
    } else {
      auto bRef = beta(pMasses[ref], mom);
      auto bHyp = beta(pMasses[hyp], mom);
      auto rRef = resoTOF[ref][hyp];
      val = (bRef - bHyp) / (bHyp*bHyp * rRef);
    }
    return val;
  }

  enum ResoMode { kTPC, kTOF };

  Double_t getReso(ResoMode mode, const Char_t* hypo, Float_t mom) {
      static TFile* file = TFile::Open(
          "UD_LHC23_pass4_SingleGap/nSigma_Resolution.root", "READ"
      );
      if (!file || file->IsZombie()) {
          std::cerr << "Error opening resolution file\n";
          return -1;
      }

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
            chi2 += (y - yfit) * (y - yfit) / (err * err);
            ++usedBins;
        }
        return chi2;
    };
    std::unique_ptr<ROOT::Math::Minimizer> minimizer(ROOT::Math::Factory::CreateMinimizer("Minuit2", ""));

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
    minimizer->SetMaxFunctionCalls(50000);   
    minimizer->SetMaxIterations(50000);
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
    for (Int_t i = 0; i < nPar; ++i) {
        func->SetParError(i, TMath::Sqrt(minimizer->CovMatrix(i, i)));
    }
    func->SetChisquare(chi2);
    func->SetNDF(ndf);
}

  inline void FitHistogramsByChi2(TH1* h1, TH1* h2, TF1* func, Int_t nG, Double_t xlo, Double_t xhi) {
    func->SetRange(xlo, xhi);
    const Int_t nPar = 4*nG +2;
    const Int_t offA1 = 0;         // Amplitudes for h1
    const Int_t offA2 = nG;        // Amplitudes for h2
    const Int_t offM  = 2*nG;      // means
    const Int_t offS  = 3*nG;      // sigmas
    const Int_t offP1 = 4*nG;      // constant background h1
    const Int_t offP2 = 4*nG + 1;  // constant background h2

    auto chi2_fcn = [&](const Double_t* par0) {
        std::vector<Double_t> par(par0, par0 + nPar);
        //par[offM + 2] = -0.272804;
        //par[offS + 2] = 1.08;
        //par[offA1 + 2] = 77.9813;
        //par[offA2 + 2] = par[offA1 + 2];
        //par[offA2 + 3] = 0;
        for (Int_t i = 0; i < nPar; ++i) func->SetParameter(i, par[i]);
        auto calc = [&](TH1* h, Int_t offA, Int_t offP) {
            Double_t chi2 = 0; 
            for (Int_t ib = 1; ib <= h->GetNbinsX(); ++ib) {
                Double_t x   = h->GetBinCenter(ib);
                Double_t y   = h->GetBinContent(ib);
                Double_t err = h->GetBinError(ib);
                if (err <= 0 || x < xlo || x > xhi) continue;
                Double_t yfit = par[offP];
                for (Int_t ig = 0; ig < nG; ++ig) {
                    Double_t A  = par[offA + ig];
                    Double_t mu = par[offM + ig];
                    Double_t s = par[offS + ig];
                    yfit += A * TMath::Gaus(x, mu, s, kFALSE);
                }

                chi2 += (y - yfit) * (y - yfit) / (err * err);
            }
            return chi2;
        };
        return calc(h1, offA1, offP1) + calc(h2, offA2, offP2);
    };

    auto minimizer = std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2",""));

    ROOT::Math::Functor fcn(chi2_fcn, nPar);
    minimizer->SetFunction(fcn);

    for (Int_t i = 0; i < nPar; ++i) {
        const char* name = func->GetParName(i);
        Double_t    val  = func->GetParameter(i);
        Double_t    err  = func->GetParError(i);
        Double_t    step = (err>0?err:(fabs(val)*0.1+1e-3));
        Double_t    lo, up;
        func->GetParLimits(i, lo, up);

        if (up > lo) {
            minimizer->SetLimitedVariable(i, name, val, step, lo, up);
        } else {
            minimizer->SetVariable(i, name, val, step);
        }
    }
    minimizer->SetMaxFunctionCalls(50000);   
    minimizer->SetMaxIterations(50000);
    minimizer->Minimize();

    const Double_t* xs = minimizer->X();
    for (Int_t i = 0; i < nPar; ++i) func->SetParameter(i, xs[i]);

    Double_t chi2 = chi2_fcn(xs);
    Int_t usedBins = 0;
    auto countBins = [&](TH1* h) {
        for (Int_t ib = 1; ib <= h->GetNbinsX(); ++ib) {
            Double_t x   = h->GetBinCenter(ib);
            Double_t err = h->GetBinError(ib);
            if (err > 0 && x >= xlo && x <= xhi) ++usedBins;
        }
    };
    countBins(h1); 
    countBins(h2);
    Double_t low = 0;
    Double_t high = 0;
    for (Int_t i = 0; i < nPar; ++i) {
        func->SetParError(i, TMath::Sqrt(minimizer->CovMatrix(i, i)));
    }
    Int_t ndf = usedBins - nPar;
    func->SetChisquare(chi2);
    func->SetNDF(ndf);
}
private:
  Double_t bethe_bloch_aleph(Double_t bg, Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5) {
    Double_t beta = bg / TMath::Sqrt(1.0 + bg*bg);
    Double_t aa   = TMath::Power(beta, p4);
    Double_t bb   = TMath::Log(p3 + TMath::Power(1.0/bg, p5));
    return (p2 - aa - bb) * p1 / aa;
  };
};
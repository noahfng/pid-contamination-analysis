#pragma once
#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "TH1.h"
#include "TF1.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "covarianceMatrix.h"

class helper{
public:
    // I/O
    const Char_t* base_dir = "/Users/noahfng/SMI/UD_LHC23_pass4_SingleGap/0106/B";

    // amount of tracks to consider per event
    const Int_t NtrkMax = 2;

    // Species bookkeeping
    enum {nParts = 5};
    const Int_t colors[nParts]  = {kOrange+7, kBlue, kViolet+1, kBlack, kGreen+2};
    const Char_t* dNames[2] = {"TPC", "TOF"};
    const Char_t* pNames[nParts] = {"El", "Mu", "Pi", "Ka", "Pr"};
    const Char_t* pCodes[nParts] = {"e", "#mu", "#pi", "K", "p"};

    // Masses in MeV/c^2; charges in |e|
    const Double_t pMasses[nParts] = {0.51099895, 105.6583755,  139.57039, 493.677, 938.27208816};
    const Double_t pCharges[nParts] = {1, 1, 1, 1, 1};

    // Resolution matrices r_{ref,hyp} used by getnSigma() (TPC/TOF)
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

    // beta = p / E  
    Double_t beta(Float_t mass, Float_t mom){
        return (mom / TMath::Sqrt(mass*mass + mom*mom));
    };

    // ALEPH-like Bethe–Bloch response. p in MeV/c.
    Double_t getTPCSignal(Double_t p, Double_t mass, Double_t charge) {
        const Double_t mMIP = 50.0; // overall scale
        const Double_t params[5] = {0.19310481, 4.26696118, 0.00522579, 2.38124907, 0.98055396};
        const Double_t chFact = 2.3; // dE/dx ~ |q|^{chFact}

        Double_t bg = p / mass;
        if (bg < 0.05) return -999.; // guard against unphysical region
        Double_t bethe = mMIP
                   * bethe_bloch_aleph(bg,
                                       params[0], params[1], params[2],
                                       params[3], params[4])
                   * TMath::Power(charge, chFact);
        return bethe >= 0 ? bethe : -999.;
    };

    // Compute nσ separation between ref and hyp for TPC or TOF at momentum mom (MeV/c).
    Double_t getnSigma(Double_t mom, TString det, Int_t ref, Int_t hyp){
        Double_t val = 0.;
        if (det == dNames[0]){ // TPC
            auto dRef = getTPCSignal(mom, pMasses[ref], pCharges[ref]);
            auto dHyp = getTPCSignal(mom, pMasses[hyp], pCharges[hyp]);
            auto rRef = resoTPC[ref][hyp];
            val = (dHyp/dRef - 1.0) / rRef;
        } 
        else { // TOF (beta-based)
            auto bRef = beta(pMasses[ref], mom);
            auto bHyp = beta(pMasses[hyp], mom);
            auto rRef = resoTOF[ref][hyp];
            val = (bRef - bHyp) / (bHyp*bHyp * rRef);
        }
        return val;
    };

    // Gaussian integral between a and b
    Double_t GaussIntegral(Double_t A, Double_t mu, Double_t sig, Double_t a, Double_t b) {
        const Double_t invsqrt2 = 1.0 / TMath::Sqrt(2.0);
        const Double_t t1 = (b - mu) * invsqrt2 / sig;
        const Double_t t0 = (a - mu) * invsqrt2 / sig;
        return A * sig * TMath::Sqrt(TMath::Pi()/2.0) * (TMath::Erf(t1) - TMath::Erf(t0));
    };

     // Piecewise-integrated model: nG Gaussians + constant per histogram 
    inline Double_t ModelIntegral(Int_t h, Double_t a, Double_t b, const std::vector<Double_t>& par, Int_t nG, Int_t offA1, Int_t offA2, Int_t offM1, Int_t offM2, Int_t offS1, Int_t offS2, Int_t offP1, Int_t offP2) {
        const Int_t offA = (h == 0 ? offA1 : offA2);
        const Int_t offM = (h == 0 ? offM1 : offM2);
        const Int_t offS = (h == 0 ? offS1 : offS2);
        const Int_t offP = (h == 0 ? offP1 : offP2);

        Double_t I = par[offP] * (b - a); // constant background
        for (Int_t j = 0; j < nG; ++j) {
            const Double_t A   = par[offA + j];
            const Double_t mu  = par[offM + j];
            const Double_t sig = par[offS + j];
            I += GaussIntegral(A, mu, sig, a, b);
        }
        return I;
    }

    // Poisson deviance (2× log-likelihood ratio vs saturated model) and total counts
    std::pair<Double_t,Double_t> PoissonDeviance(TH1* h, Int_t whichHist, const std::vector<Double_t>& par, Int_t nG, Int_t offA1, Int_t offA2, Int_t offM1, Int_t offM2, Int_t offS1, Int_t offS2, Int_t offP1, Int_t offP2){
        Double_t D = 0.0, N = 0.0;
        const Int_t nb = h->GetNbinsX();
        for (Int_t ib=1; ib<=nb; ++ib) {
            const Double_t n  = h->GetBinContent(ib);
            const Double_t a  = h->GetBinLowEdge(ib);
            const Double_t b  = h->GetBinLowEdge(ib+1);
            const Double_t mu = ModelIntegral(whichHist, a, b, par, nG, offA1, offA2, offM1, offM2, offS1, offS2, offP1, offP2);
            N += n;
            if (n > 0.0) D += 2.0 * (mu - n + n * std::log(n/mu));
            else         D += 2.0 * mu;
        }
        return {D, N};
    };

    enum ResoMode {kTPC, kTOF};
    // Resolution vs momentum from file (graph name: nSigma{TPC|TOF}res<hypo>)
    Double_t getReso(ResoMode mode, const Char_t* hypo, Float_t mom) {
        static TFile* file = TFile::Open("UD_LHC23_pass4_SingleGap/nSigma_Resolution.root", "READ");
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

    // Assemble expected vector for CM fit by integrating the model per CM bin
    TVectorD BuildExpectedVector(TH1* h1, TH1* h2, covarianceMatrix* cm, const std::vector<TH1D*>& histos, const std::vector<Double_t>& par, Int_t nG, Int_t offA1, Int_t offA2, Int_t offM1, Int_t offM2, Int_t offS1, Int_t offS2, Int_t offP1, Int_t offP2, Double_t nBins) {
        const auto bins2c4f = cm->bins2c4f(); 
        const Int_t nbins   = cm->observations().GetNrows();
        TVectorD expected(nbins);
        Int_t k = 0;
        for (Int_t h = 0; h < (Int_t)histos.size(); ++h) {
            auto* ax = histos[h]->GetXaxis();
            for (Int_t binAbs : bins2c4f[h]) {
            const Double_t a = ax->GetBinLowEdge(binAbs);
            const Double_t b = ax->GetBinUpEdge(binAbs);
            expected[k++] = ModelIntegral(h, a, b, par, nG, offA1, offA2, offM1, offM2, offS1, offS2, offP1, offP2);
            }
        }
        return expected;
    };

    // χ² = rᵀPr  with P = precomputed inverse covariance (from cm)
    Double_t Chi2_withCM(covarianceMatrix* cm, const TVectorD& expected) {
        const TVectorD y = cm->observations();
        TVectorD r = y; 
        r -= expected;
        TMatrixDSym* Pinv = cm->preMatrix();
        if (!Pinv) return std::numeric_limits<double>::infinity();
        TVectorD tmp = (*Pinv) * r;
        return r * tmp; 
    };

    // Standalone χ² minimization for a single histogram and TF1
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

        // Seed parameters (+ respect TF1 limits)
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

        // Export result + errors
        const Double_t* xs = minimizer->X();
        for (Int_t i = 0; i < nPar; ++i) {
            func->SetParameter(i, xs[i]);
        }
        
        // Store χ²/ndf on TF1
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
    };

    // Joint fit of two hists with common μ,σ; A2 enforced via A2 = A1*(σ1/σ2)
    inline void FitHistogramsExclCompByChi2(TH1* h1, TH1* h2, TF1* func, Int_t nG, covarianceMatrix* cm, const std::vector<TH1D*>& cmHists, const Double_t eigenThr, const Double_t nBins) {
        const Int_t nPar = 4*nG +2;
        const Int_t offA1 = 0;         // A for h1
        const Int_t offA2 = nG;        // A for h2
        const Int_t offM  = 2*nG;      // common μ
        const Int_t offS  = 3*nG;      // common σ
        const Int_t offP1 = 4*nG;      // const bg for h1
        const Int_t offP2 = 4*nG + 1;  // const bg for h2

        auto chi2_fcn = [&](const Double_t* par0) {
            std::vector<Double_t> par(par0, par0 + nPar);
            for (Int_t i = 0; i < nPar; ++i) func->SetParameter(i, par[i]);
            TVectorD expected = BuildExpectedVector(h1, h2, cm, cmHists, par, nG, offA1, offA2, offM, offM, offS, offS, offP1, offP2, nBins);
            const Double_t chi2_cov = Chi2_withCM(cm, expected);
            return chi2_cov;
        };

        auto minimizer = std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2",""));
        ROOT::Math::Functor fcn(chi2_fcn, nPar);
        minimizer->SetFunction(fcn);

        // Seed + limits
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

        for (Int_t i = 0; i < nPar; ++i) {
            func->SetParameter(i, minimizer->X()[i]);
            func->SetParError(i, TMath::Sqrt(minimizer->CovMatrix(i, i)));
        }
        const Double_t chi2 = minimizer->MinValue();
        const Int_t nb = cm->observations().GetNrows();
        const Int_t ndf = std::max(1, nb - nPar); 
        func->SetChisquare(chi2);
        func->SetNDF(ndf);
    }

    // Joint fit with per-hist μ,σ; link A2 to A1 via σ ratio and fix A2 during minimization
    inline void FitHistogramsByChi2(TH1* h1, TH1* h2, TF1* func, Int_t nG, covarianceMatrix* cm, const std::vector<TH1D*>& cmHists, const Double_t eigenThr, const Double_t nBins) {
        const Int_t nPar = 6*nG +2;
        const Int_t offA1 = 0;         // A for h1
        const Int_t offA2 = nG;        // A for h2
        const Int_t offM1 = 2*nG;      // μ for h1
        const Int_t offM2 = 3*nG;      // μ for h2
        const Int_t offS1 = 4*nG;      // σ for h1
        const Int_t offS2 = 5*nG;      // σ for h2
        const Int_t offP1 = 6*nG;      // const bg for h1
        const Int_t offP2 = 6*nG + 1;  // const bg for h2

        auto chi2_fcn = [&](const Double_t* par0) {
            std::vector<Double_t> par(par0, par0 + nPar);
            // enforce A2 from A1 and σ ratio during evaluation
            for (Int_t j = 0; j < nG; ++j) {
                const Double_t s1 = std::max(par[offS1 + j], 1e-12);
                const Double_t s2 = std::max(par[offS2 + j], 1e-12);
                par[offA2 + j] = par[offA1 + j] * (s1 / s2);
            }
            for (Int_t i = 0; i < nPar; ++i) func->SetParameter(i, par[i]);
            TVectorD expected = BuildExpectedVector(h1, h2, cm, cmHists, par, nG, offA1, offA2, offM1, offM2, offS1, offS2, offP1, offP2, nBins);
            const Double_t chi2_cov = Chi2_withCM(cm, expected);
            return chi2_cov;
        };

        auto minimizer = std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2",""));
        ROOT::Math::Functor fcn(chi2_fcn, nPar);
        minimizer->SetFunction(fcn);

         // Seed; fix A2 parameters (derived)
        for (Int_t i = 0; i < nPar; ++i) {
            const char* name = func->GetParName(i);
            Double_t    val  = func->GetParameter(i);
            Double_t    err  = func->GetParError(i);
            Double_t    step = (err>0?err:(fabs(val)*0.1+1e-3));
            Double_t    lo, up;
            func->GetParLimits(i, lo, up);
            if (i >= offA2 && i < offA2 + nG) {
                minimizer->SetFixedVariable(i, name, val);
                continue;
            }

            if (up > lo) {
                minimizer->SetLimitedVariable(i, name, val, step, lo, up);
            } else {
                minimizer->SetVariable(i, name, val, step);
            }
        }

        minimizer->SetMaxFunctionCalls(50000);
        minimizer->SetMaxIterations(50000);
        minimizer->Minimize();

        // Export best-fit + errors
        for (Int_t i = 0; i < nPar; ++i) {
            func->SetParameter(i, minimizer->X()[i]);
            func->SetParError(i, TMath::Sqrt(std::max(0.0, minimizer->CovMatrix(i, i))));
        }

        // Fill derived A2 and zero its error (since it’s constrained)
        for (Int_t j = 0; j < nG; ++j) {
            const double s1 = std::max(func->GetParameter(offS1 + j), 1e-12);
            const double s2 = std::max(func->GetParameter(offS2 + j), 1e-12);
            const double A2 = func->GetParameter(offA1 + j) * (s1 / s2);
            func->SetParameter(offA2 + j, A2);
            func->SetParError(offA2 + j, 0.0); 
        }

        const Double_t chi2 = minimizer->MinValue();
        const Int_t nb = cm->observations().GetNrows();
        const Int_t ndf = std::max(1, nb - nPar);
        func->SetChisquare(chi2);
        func->SetNDF(ndf);
    }
private:
    // ALEPH parameterization; input γβ, returns normalized dE/dx shape
    Double_t bethe_bloch_aleph(Double_t bg, Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5) {
        Double_t beta = bg / TMath::Sqrt(1.0 + bg*bg);
        Double_t aa   = TMath::Power(beta, p4);
        Double_t bb   = TMath::Log(p3 + TMath::Power(1.0/bg, p5));
        return (p2 - aa - bb) * p1 / aa;
    };
};
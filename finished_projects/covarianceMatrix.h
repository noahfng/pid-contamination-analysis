#include "TMath.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "TLatex.h"
#include "TDecompSVD.h"
#include "TMatrixD.h"
#include "TH2.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TMatrixDSymEigen.h"


class covarianceMatrix
{
  private:
    // options
    Bool_t moffDiagonal;
    Double_t meigenThreshold;

    // histogram parameters
    Int_t mnHistos;
    std::vector<TH1D*> mHistos;
    std::vector<std::tuple<Double_t, Double_t>*> mFitLims;
    std::vector<std::tuple<Int_t, Int_t>> mFitBinLims;
    std::vector<Int_t> mnBins;
    std::vector<std::vector<Int_t>> mBinsInLims; 
    std::vector<Int_t> mnBins2c4f;                    // 2c4f = to consider for fit
    std::vector<std::vector<Int_t>> mBinsInLims2c4f; 

    // structures which are filled in the event loop
    Long64_t mnEvents;
    Int_t mDim;
    std::vector<Int_t> mOffsets;
    std::vector<Double_t>* mCMVec;
    TMatrixDSym* mCMOuter;

    // covariance parameters
    Int_t mDim2c4f;
    TVectorD mObservations;
    TMatrixDSym* mCM;
    TMatrixDSym* mCMInverted;
    Bool_t mCMready;

    // for plotting
    TLatex* mtTxt;

    // this should be called when the filling of the histogram is finished
    // calculation of the bins with !0 entries and within the limits - 2c4f = to consider for fit
    void updateBins2c4f()
    {
      std::vector<Int_t> binsInLims2c4f;
      mBinsInLims2c4f.clear();
      mDim2c4f = 0;
      for (Int_t h=0; h<mnHistos; ++h)
      {
        binsInLims2c4f.clear();
        for (auto bin : mBinsInLims[h])
        {
          if (mHistos[h]->GetBinContent(bin) > 0.) binsInLims2c4f.push_back(bin);
        }
        mBinsInLims2c4f.push_back(binsInLims2c4f);
        mnBins2c4f.push_back(binsInLims2c4f.size());
        mDim2c4f += mnBins2c4f.back();
      }
    }

    // make a pseudo inversion
    Bool_t doCMInversion()
    {
      // do SVD decomposition
      TMatrixD mat(*mCM);
      TDecompSVD svd(mat);

      // get singular values
      TVectorD singVals = svd.GetSig();
    
      // truncate small singular values
      TMatrixD invS(mDim2c4f, mDim2c4f);
      for (int i = 0; i < mDim2c4f; i++) {
        if (singVals[i]/singVals[0] < meigenThreshold) {
          invS(i,i) = 0;
        } else {
          invS(i,i) = 1.0 / singVals[i];
        }
      }

      // Compute pseudo-inverse with truncated singular values
      // V * invS * U^T
      TMatrixD V = svd.GetV();
      TMatrixD U = svd.GetU();
      TMatrixD pseudoInv = V * invS * U.T();
      // make sure that it is symmetric
      TMatrixD pseudoInvSym = (pseudoInv + pseudoInv.T()) * 0.5;

      // create mCMInverted
      //delete mCMInverted;
      mCMInverted = new TMatrixDSym(mDim2c4f, pseudoInvSym.GetMatrixArray());

      return true;
    }
  
  public:
    // constructor
    //  hists: list of histograms of type TH1D*
    //  fitlims: a list of tuples with min/max nSigma values to consider (fit limits per histogram)
    covarianceMatrix(std::vector<TH1D*> hists, std::vector<std::tuple<Double_t, Double_t>*> fitlims)
    {
      // number of histograms
      mnHistos = hists.size();

      // update list of histograms and related items
      mDim = 0;
      for (Int_t h=0; h<mnHistos; ++h)
      {
        // update mHistos and mFitLims
        mHistos.push_back(hists[h]);
        mFitLims.push_back(fitlims[h]);
        
        // update mnBins, mFitBinLims, and mBinsInLims
        auto bmin = hists[h]->GetXaxis()->FindBin(std::get<0>(*(fitlims[h])));
        auto bmax = hists[h]->GetXaxis()->FindBin(std::get<1>(*(fitlims[h])));
        mnBins.push_back(bmax - bmin + 1);
        mFitBinLims.push_back(std::make_tuple(bmin, bmax));
        
        std::vector<Int_t> binsInLims;
        for (Int_t b=bmin; b<=bmax; ++b) binsInLims.push_back(b);
        mBinsInLims.push_back(binsInLims);

        // update mOffsetsn and mDim
        mOffsets.push_back(mDim);
        mDim += mnBins.back();
      }

      // create the CMVec and CMOuter
      mCMVec = new std::vector<Double_t>(mDim, 0.0);
      mCMOuter = new TMatrixDSym(mDim);

      // reset number of events and CM readiness
      mnEvents = 0;
      mCMready = false;
      
      // prepare the label
      mtTxt = new TLatex();
      mtTxt->SetTextFont(42);
      mtTxt->SetTextSize(0.02);
      mtTxt->SetTextColor(kBlack);
      mtTxt->SetNDC();
    };

    // add an event, update histograms and CM
    //  nsVals: list of nSigma values, one per histogram
    //  toadd: indicates if a nSigma value is added to a histogram
    Bool_t addEvent(std::vector<Double_t>* nsVals, std::vector<Bool_t>* toadd)
    {
      // update histogram and sparse list of activated global bins for this event
      std::vector<Int_t> activeBins;
      activeBins.reserve(mnHistos);
      for (Int_t h = 0; h < mnHistos; ++h)
      {
        if ((*toadd)[h])
        {
          mHistos[h]->Fill((*nsVals)[h], 1.);
        
          auto blims = mFitBinLims[h];
          auto b = mHistos[h]->GetXaxis()->FindBin((*nsVals)[h]) - std::get<0>(blims);
          if ( b >= 0 && b < mnBins[h] )
          {
            activeBins.push_back(mOffsets[h] + b);
          }
        }
      }
      Int_t nAct = activeBins.size();
      if (nAct != 0)
      {
        // Update sums
        for (Int_t idx : activeBins) {
          (*mCMVec)[idx] += 1.0;
          // Diagonal term
          (*mCMOuter)(idx, idx) += 1.0;
        }
        // Cross terms
        for (Int_t a = 0; a < nAct; ++a)
        {
          for (Int_t b = a + 1; b < nAct; ++b)
          {
            Int_t r = activeBins[a];
            Int_t c = activeBins[b];
            (*mCMOuter)(r, c) += 1.0;
            (*mCMOuter)(c, r) += 1.0; // symmetric
          }
        }
        
        // update the total number of events
        ++mnEvents;
      }
      
      return true;
    }

    // compute the covariance matrix and create a TVectorD with the histogram bin values to fit
    void make(Bool_t offDiagonal = true, Double_t eigenThreshold = 1.E-5)
    {
      // update the list of bins within limits and with entries
      updateBins2c4f();
      
      // prepare the TVector with the histogram bins to fit
      mObservations.ResizeTo(mDim2c4f);
      Int_t cnt = 0;
      for (Int_t h=0; h<mnHistos; ++h)
      {
        for (auto i : mBinsInLims2c4f[h])
        {
          mObservations[cnt] = mHistos[h]->GetBinContent(i);
          ++cnt;
        }
      }
      
      // prepare the final covariance matrix
      //delete mCM;
      mCM = new TMatrixDSym(mDim2c4f);
            
      // Compute covariance: C = (mCMOuter - (mCMVec*mCMVec^T)/mnEvents) * (mnEvents / (mnEvents-1))
      // skip empty elements
      moffDiagonal = offDiagonal;
      Int_t d1 = 0, d2;
      for (Int_t r = 0; r < mDim; ++r) {
        if ((*mCMVec)[r] == 0.) continue;
        d2 = -1;
        for (Int_t c = 0; c <= r; ++c) {
          if ((*mCMVec)[c] == 0.) continue;
          ++d2;
          if (!moffDiagonal && d1 != d2) continue;
          
          Double_t term = (*mCMOuter)(r, c) - (*mCMVec)[r] * (*mCMVec)[c] / Double_t(mnEvents);
          term *= Double_t(mnEvents) / Double_t(mnEvents - 1);
          (*mCM)(d1, d2) = term;
          (*mCM)(d2, d1) = term; // symmetric
        }
        ++d1;
      }
      
      // computer inverted CM
      meigenThreshold = eigenThreshold;
      if (doCMInversion()) mCMready = true;
    }
    
    // some getters
    std::vector<TH1D*> histograms() { return mHistos; }
    std::vector<std::vector<Int_t>> bins2c4f() { return mBinsInLims2c4f; }; 
    TVectorD observations() { return mObservations; }
    TMatrixDSym* covMatrix() { return mCM; }
    TMatrixDSym* preMatrix() { return mCMInverted; }

    void plot(TCanvas* c, TString pdfName)
    {
      // make must be run first
      if (!mCMready) make();
      
      // convert mCM and mCMInverted to TH2D
      auto hCM = new TH2D("CM", "Covariance matrix;;;", mDim2c4f, 0, mDim2c4f, mDim2c4f, 0, mDim2c4f);
      hCM->SetStats(0);
      auto hCMInverted = new TH2D("CMInverted", "Inverted covariance matrix;;;", mDim2c4f, 0, mDim2c4f, mDim2c4f, 0, mDim2c4f);
      hCMInverted->SetStats(0);
      
      // fill histograms with matrix values
      for (int i = 0; i < mDim2c4f; ++i) {
        for (int j = 0; j < mDim2c4f; ++j) {
          hCM->SetBinContent(i + 1, j + 1, (*mCM)(i, j));
          hCMInverted->SetBinContent(i + 1, j + 1, (*mCMInverted)(i, j));
        }
      }
      
      // draw the CM
      c->SetLogx(false);
      c->SetLogy(false);
      c->SetLogz(true);
      hCM->Draw("COLZ");
      c->Print(pdfName);
      
      // draw the CMInverted
      hCMInverted->Draw("COLZ");
      mtTxt->SetTextAlign(11);
      mtTxt->DrawLatex(0.16, 0.91, Form("eigen threshold: %.2e", meigenThreshold));
      c->Print(pdfName); 

      // add plot of eigen values 
      TMatrixDSymEigen eigen(*mCM);
      TVectorD eigenValues = eigen.GetEigenValues();
      
      TGraph* gr = new TGraph(mDim2c4f);
      for (Int_t i = 0; i < mDim2c4f; ++i) {
        gr->SetPoint(i, i, eigenValues[i]/eigenValues[0]);
      }
      gr->SetTitle(";Index;Relative eigen value");
      c->SetLogy(true);
      gr->Draw("ALP");
      if (meigenThreshold > 0.)
      {
        TLine* lethr = new TLine(0, meigenThreshold, mDim2c4f, meigenThreshold);
        lethr->SetLineStyle(3); lethr->SetLineWidth(1);
        lethr->Draw("SAME");
      }
      mtTxt->DrawLatex(0.16, 0.91, Form("eigen threshold: %.2e", meigenThreshold));
      c->Print(pdfName);
      
      // clean up
      delete hCM;
      delete hCMInverted;
      delete gr;
    }
};

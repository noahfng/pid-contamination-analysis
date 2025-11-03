#include <stdio.h>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iomanip>

#include "TChain.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"

#include <AddTrees.h>
#include <helpers.h>

void countTOF() {
  auto help = new helper();
  TChain chain("twotauchain");
  AddTrees(chain, help->base_dir);
  const Int_t NtrkMax = help->NtrkMax;

  Long64_t total = chain.GetEntries();
  Long64_t nEntries = std::min(total, static_cast<Long64_t>(1e7));
  printf("Chain has %lld entries (processing %lld)\n", total, nEntries);

  chain.SetBranchStatus("*", 0);
  chain.SetBranchStatus("fRunNumber", 1);
  chain.SetBranchStatus("fTrkTOFexpMom", 1);
  chain.SetBranchStatus("fTrkTPCinnerParam", 1);

  std::vector<Int_t> runNumber(NtrkMax);
  std::vector<Float_t> TOFexpMom(NtrkMax);
  std::vector<Float_t> innerParam(NtrkMax);

  chain.SetBranchAddress("fRunNumber",        runNumber.data());
  chain.SetBranchAddress("fTrkTOFexpMom",     TOFexpMom.data());
  chain.SetBranchAddress("fTrkTPCinnerParam", innerParam.data());

  std::map<Int_t, Long64_t> withTOF;
  std::map<Int_t, Long64_t> withoutTOF;
  std::map<Int_t, Long64_t> withP03;

  for (Long64_t i = 0; i < nEntries; ++i) {
    chain.GetEntry(i);

    Int_t rn = runNumber[0]; 
    Bool_t hasTOF = false;
    Bool_t hasP03 = false;

    
    for (Int_t j = 0; j < NtrkMax; ++j) {
      if (TOFexpMom[j] >= 0) {
        hasTOF = true;
        break;
      }
    }

    
    for (Int_t j = 0; j < NtrkMax; ++j) {
      if (innerParam[j] >= 0.3) { 
        hasP03 = true;
        break;
      }
    }

    if (hasTOF) ++withTOF[rn];
    else        ++withoutTOF[rn];

    if (hasP03) ++withP03[rn];
  }

  std::set<Int_t> allRuns;
  for (auto& p : withTOF)    allRuns.insert(p.first);
  for (auto& p : withoutTOF) allRuns.insert(p.first);
  for (auto& p : withP03)    allRuns.insert(p.first);


  std::ofstream jsonFile("run_counts.json");
  jsonFile << std::fixed << std::setprecision(2);
  jsonFile << "{\n";
  Bool_t first = true;
  for (auto rn : allRuns) {
    if (!first) jsonFile << ",\n";
    first = false;
    Long64_t totalRunEntries = withTOF[rn] + withoutTOF[rn];
    jsonFile << "  \"" << rn << "\": {\n"
            << "    \"withTOF\": " << (100.0 * withTOF[rn] / totalRunEntries) << ",\n"
            << "    \"withoutTOF\": " << (100.0 * withoutTOF[rn] / totalRunEntries) << ",\n"
            << "    \"withP03\": " << (100.0 * withP03[rn] / totalRunEntries) << "\n"
            << "  }";
  }
  jsonFile << "\n}\n";
  jsonFile.close();
  printf("\nSaved run statistics to run_counts.json\n");
}
#include <stdio.h>
#include <map>
#include <set>
#include <vector>
#include <algorithm>          

#include "TChain.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "AddTrees.h"

void countTOF() {
  TChain chain("twotauchain");
  const Char_t* base_dir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";
  AddTrees(chain, base_dir);

  Long64_t total = chain.GetEntries();
  Long64_t nEntries = std::min(total, static_cast<Long64_t>(1e4));
  printf("Chain has %lld entries (processing %lld)\n", total, nEntries);

  chain.SetBranchStatus("*",            0);
  chain.SetBranchStatus("fRunNumber",    1);
  chain.SetBranchStatus("fTrkTOFexpMom", 1);

  TTreeReader reader(&chain);
  TTreeReaderValue<Int_t>    run   (reader, "fRunNumber");
  TTreeReaderArray<Float_t>  tofArr(reader, "fTrkTOFexpMom");

  std::map<Int_t, Long64_t> withTOF, withoutTOF;

  Long64_t i = 0;
  while (i < nEntries && reader.Next()) {
    Int_t rn = *run;
    bool has = false;
    for (auto v : tofArr) {
      if (v >= 0) {
        has = true;
        break;
      }
    }
    if (has)       ++withTOF   [rn];
    else           ++withoutTOF[rn];
    ++i;
  }

  std::set<Int_t> allRuns;
  for (auto &p : withTOF)    allRuns.insert(p.first);
  for (auto &p : withoutTOF) allRuns.insert(p.first);

  printf("\n%10s | %12s | %14s\n",
         "RunNumber", "With TOF data", "Without TOF data");
  printf("---------------------------------------------\n");
  for (auto rn : allRuns) {
    printf("%10d | %12lld | %14lld\n",
           rn,
           withTOF   [rn],
           withoutTOF[rn]);
  }
}

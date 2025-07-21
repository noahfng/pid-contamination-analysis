enum ResoMode { kTPC, kTOF };

Float_t getReso(ResoMode mode, const Char_t* hypo, Float_t mom) {
    // open the same file just once
    static TFile* file = TFile::Open(
        "UD_LHC23_pass4_SingleGap/nSigma_Resolution.root", "READ"
    );
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening resolution file\n";
        return -1;
    }

    // build name: "nSigmaTPCresEl" or "nSigmaTOFresPr", etc.
    const char* prefix = (mode == kTPC ? "nSigmaTPCres" : "nSigmaTOFres");
    TString gname = Form("%s%s", prefix, hypo);

    auto gr = dynamic_cast<TGraph*>(file->Get(gname));
    if (!gr) {
        std::cerr << "Graph " << gname.Data() << " not found\n";
        return -1;
    }

    return gr->Eval(mom);
}
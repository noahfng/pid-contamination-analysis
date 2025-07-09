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

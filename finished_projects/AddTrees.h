void AddTrees(TChain &chain, const Char_t* baseDir) {
    // walk baseDir/* subdirs and add O2tautwotrack trees found under DF_* keys
    TSystemDirectory dir("base", baseDir);
    TList *subdirs = dir.GetListOfFiles();
    TSystemFile *sysfile;
    TIterator *itSub = subdirs->MakeIterator();

    while ((sysfile = (TSystemFile*)itSub->Next())) {
        TString dname = sysfile->GetName();
        if (!sysfile->IsDirectory() || !dname.BeginsWith("hy_")) continue;  // only hy_* dirs
        TString fullDir = TString(baseDir) + "/" + dname;

        TSystemDirectory subdir(dname, fullDir);
        TList *files = subdir.GetListOfFiles();
        TSystemFile *f2;
        TIterator *itFile = files->MakeIterator();

        while ((f2 = (TSystemFile*)itFile->Next())) {
            TString fname = f2->GetName();
            if (!fname.BeginsWith("RLAnalysisTree") || !fname.EndsWith(".root")) continue; // target files
            TString path = fullDir + "/" + fname;

            TFile tf(path, "READ");
            if (tf.IsZombie()) { tf.Close(); continue; } // skip broken files

            // find first top-level key like DF_* and add its O2tautwotrack subdir as a tree
            TIterator *itKey = tf.GetListOfKeys()->MakeIterator();
            TKey *key;
            while ((key = (TKey*)itKey->Next())) {
                TString keyName = key->GetName();
                if (keyName.BeginsWith("DF_")) {
                    chain.Add(path + "/" + keyName + "/O2tautwotrack");
                    break; // one DF_* per file is enough
                }
            }
            delete itKey;
            tf.Close();
        }
        delete itFile;
    }
    delete itSub;
}

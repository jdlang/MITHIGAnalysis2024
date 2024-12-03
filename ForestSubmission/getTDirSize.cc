#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TObject.h>
#include <iostream>

long long GetDirectorySize(TDirectory* dir, bool printSize) {
    long long size = 0;
    TIter next(dir->GetListOfKeys());
    TKey* key;

    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        long long dirSize = 0;
        if (obj) {
            TDirectory* subdir = dynamic_cast<TDirectory*>(obj);
            if (subdir) {
                // Recursive call for subdirectories
                size += GetDirectorySize(subdir, 0);
                dirSize += GetDirectorySize(subdir, 0);
            } else {
                // Get the size of the object
                size += key->GetNbytes();
                dirSize += key->GetNbytes();
            }
        }
        if (printSize) std::cout<<obj->GetName()<<" - "<<dirSize<<std::endl;
    }
    return size;
}

int getTDirSize() {
    // Open the ROOT file
    TFile* file = TFile::Open("HiForestMiniAOD_1.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file." << std::endl;
        return 1;
    }

    // Get the size of the top-level directory
    long long totalSize = GetDirectorySize(file, 1);
    std::cout << "Total size of TDirectories: " << totalSize << " bytes" << std::endl;

    // Clean up
    file->Close();
    delete file;

    return 0;
}

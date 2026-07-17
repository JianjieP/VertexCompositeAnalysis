#include <iostream>

#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

using namespace std;

void merge_chained_ParticleTrees() {
    TString outDir = "/eos/user/j/jianjie/fourpi/HIRun2025A_FourPi_0716_MINIAOD";
    TString mergedFile = outDir + "/chained_ParticleTree_merged.root";

    if (gSystem->mkdir(outDir, true) != 0 && gSystem->AccessPathName(outDir)) {
        cout << "Failed to create output directory: " << outDir << endl;
        return;
    }

    TString findCmd = Form("find %s -maxdepth 1 -type f -name 'chained_ParticleTree_*.root' | sort || true", outDir.Data());
    TString fileText = gSystem->GetFromPipe(findCmd);
    TObjArray* fileList = fileText.Tokenize("\n");

    Int_t totalFiles = fileList->GetEntries();
    cout << "Found " << totalFiles << " per-i chained ROOT files." << endl;

    if (totalFiles == 0) {
        cout << "No input files found for merging." << endl;
        return;
    }

    TChain chain("ParticleTree");
    Int_t addedFiles = 0;

    for (Int_t i = 0; i < totalFiles; ++i) {
        TString fileName = ((TObjString*)fileList->At(i))->GetString();
        if (fileName.IsNull()) {
            continue;
        }

        Long64_t result = chain.Add(fileName);
        if (result > 0) {
            ++addedFiles;
            cout << "Added " << fileName << endl;
        } else {
            cout << "Skipped " << fileName << endl;
        }
    }

    cout << "Merging " << addedFiles << " files." << endl;
    cout << "Merged entries available in chain: " << chain.GetEntries() << endl;

    if (chain.GetEntries() == 0) {
        cout << "No entries found to merge." << endl;
        return;
    }

    TFile* outFile = TFile::Open(mergedFile, "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        cout << "Failed to create " << mergedFile << endl;
        return;
    }

    outFile->cd();
    TTree* outTree = chain.CloneTree(0);
    if (!outTree) {
        cout << "CloneTree failed." << endl;
        outFile->Close();
        return;
    }

    TStopwatch writeWatch;
    writeWatch.Start();

    Long64_t totalEntries = chain.GetEntries();
    Long64_t writeStep = totalEntries / 100;
    if (writeStep < 1) {
        writeStep = 1;
    }

    bool writeFailed = false;
    for (Long64_t entry = 0; entry < totalEntries; ++entry) {
        if (chain.GetEntry(entry) <= 0) {
            continue;
        }
        Int_t filled = outTree->Fill();
        if (filled <= 0 || outFile->TestBit(TFile::kWriteError)) {
            cout << "Write failed around entry " << entry
                 << ". Check output filesystem/quota and free space." << endl;
            writeFailed = true;
            break;
        }

        if ((entry + 1) % writeStep == 0 || entry + 1 == totalEntries) {
            double elapsed = writeWatch.RealTime();
            writeWatch.Continue();
            double fraction = static_cast<double>(entry + 1) / static_cast<double>(totalEntries);
            double eta = (fraction > 0.0) ? elapsed * (1.0 - fraction) / fraction : 0.0;
            cout << "Merging entries: " << (entry + 1) << "/" << totalEntries
                 << " (" << TMath::Nint(fraction * 100.0) << "%)"
                 << ", elapsed " << Form("%.1f", elapsed) << " s"
                 << ", ETA " << Form("%.1f", eta) << " s" << endl;
        }
    }

    if (writeFailed) {
        outFile->Close();
        return;
    }

    outTree->SetAutoSave(0);
    Long64_t written = outTree->Write("ParticleTree", TObject::kOverwrite);
    outFile->Close();

    cout << "Write() returned " << written << endl;
    cout << "Saved merged file " << mergedFile << endl;
}
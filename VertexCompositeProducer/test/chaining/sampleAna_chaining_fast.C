#include <iostream>

#include "TChain.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

using namespace std;

void sampleAna_chaining_fast() {
    TChain chain("fourPiAna/ParticleTree");

    cout << "Searching for EOS files..." << endl;

    TString findCmd = "find /eos/user/j/jianjie/fourpi -type f "
                      "-path '*/HIForward*/HIForward*/HIForward*_HIRun2025A_FourPi_0410/*' "
                      "-name 'fourpi_ana_*.root' 2>/dev/null || true";
    TString fileText = gSystem->GetFromPipe(findCmd);
    TObjArray* fileList = fileText.Tokenize("\n");

    Int_t totalFiles = fileList->GetEntries();
    Int_t addedFiles = 0;
    Int_t skippedFiles = 0;
    Int_t progressStep = totalFiles / 100;
    if (progressStep < 1) {
        progressStep = 1;
    }

    TStopwatch checkWatch;
    checkWatch.Start();

    cout << "Found " << totalFiles << " candidate files. Checking keys and adding valid trees..." << endl;

    for (Int_t i = 0; i < totalFiles; ++i) {
        TString fileName = ((TObjString*)fileList->At(i))->GetString();
        if (fileName.IsNull()) {
            continue;
        }

        TFile* testFile = TFile::Open(fileName, "READ");
        if (!testFile || testFile->IsZombie()) {
            ++skippedFiles;
            if (testFile) {
                testFile->Close();
            }
            continue;
        }

        TDirectoryFile* dir = (TDirectoryFile*)testFile->Get("fourPiAna");
        if (dir && dir->GetListOfKeys() && dir->GetListOfKeys()->FindObject("ParticleTree")) {
            Long64_t result = chain.Add(fileName);
            if (result > 0) {
                ++addedFiles;
            }
        } else {
            ++skippedFiles;
        }

        testFile->Close();

        if ((i + 1) % progressStep == 0 || i + 1 == totalFiles) {
            double elapsed = checkWatch.RealTime();
            checkWatch.Continue();
            double fraction = static_cast<double>(i + 1) / static_cast<double>(totalFiles);
            double eta = (fraction > 0.0) ? elapsed * (1.0 - fraction) / fraction : 0.0;
            cout << "Checking files: " << (i + 1) << "/" << totalFiles
                 << " (" << TMath::Nint(fraction * 100.0) << "%)"
                 << ", accepted " << addedFiles
                 << ", skipped " << skippedFiles
                 << ", elapsed " << Form("%.1f", elapsed) << " s"
                 << ", ETA " << Form("%.1f", eta) << " s" << endl;
        }
    }

    cout << "Direct add complete." << endl;
    cout << "Accepted " << addedFiles << " files out of " << totalFiles << endl;
    cout << "Skipped " << skippedFiles << " files without fourPiAna/ParticleTree" << endl;
    cout << "Number of entries in chain: " << chain.GetEntries() << endl;

    if (chain.GetEntries() == 0) {
        cout << "No entries found in chain!" << endl;
        return;
    }

    cout << "Writing chained tree to chained_ParticleTree.root ..." << endl;
    TFile* outFile = TFile::Open("chained_ParticleTree.root", "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        cout << "Failed to create chained_ParticleTree.root" << endl;
        return;
    }

    outFile->cd();
    TStopwatch writeWatch;
    writeWatch.Start();

    TTree* outTree = chain.CloneTree(0);
    if (!outTree) {
        cout << "CloneTree failed." << endl;
        outFile->Close();
        return;
    }

    Long64_t totalEntries = chain.GetEntries();
    Long64_t writeStep = totalEntries / 50;
    if (writeStep < 1) {
        writeStep = 1;
    }

    for (Long64_t entry = 0; entry < totalEntries; ++entry) {
        if (chain.GetEntry(entry) <= 0) {
            continue;
        }
        outTree->Fill();

        if ((entry + 1) % writeStep == 0 || entry + 1 == totalEntries) {
            double elapsed = writeWatch.RealTime();
            writeWatch.Continue();
            double fraction = static_cast<double>(entry + 1) / static_cast<double>(totalEntries);
            double eta = (fraction > 0.0) ? elapsed * (1.0 - fraction) / fraction : 0.0;
            cout << "Writing output: " << (entry + 1) << "/" << totalEntries
                 << " (" << TMath::Nint(fraction * 100.0) << "%)"
                 << ", elapsed " << Form("%.1f", elapsed) << " s"
                 << ", ETA " << Form("%.1f", eta) << " s" << endl;
        }
    }

    Long64_t written = outTree->Write("ParticleTree");
    outFile->Close();

    cout << "Write() returned " << written << endl;
    cout << "Saved chained tree with " << chain.GetEntries() << " entries" << endl;
}
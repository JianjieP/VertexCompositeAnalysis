#include <iostream>

#include "TChain.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

using namespace std;

void sampleAna_chaining() {
    // Create a TChain and add your ROOT files
    // Note: Tree is inside fourPiAna directory, so we specify the full path
    TChain* chain = new TChain("fourPiAna/ParticleTree");
    
    cout << "Searching for files and checking for valid trees..." << endl;
    
    // Search from a stable root path; avoid shell-glob expansion failures in GetFromPipe.
    TString findCmd = "find /eos/user/j/jianjie/fourpi -type f "
                      "-path '*/HIForward*/HIForward*/HIForward*_HIRun2025A_FourPi_0410/*' "
                      "-name 'fourpi_ana_*.root' 2>/dev/null || true";
    TObjArray* fileList = gSystem->GetFromPipe(findCmd).Tokenize("\n");
    
    Int_t validFiles = 0;
    Int_t totalFiles = fileList->GetEntries();
    Int_t checkedFiles = 0;
    
    cout << "Found " << totalFiles << " candidate files in HIForward*_HIRun2025A_FourPi_0410 directories. Checking for valid trees..." << endl;
    
    for (Int_t i = 0; i < totalFiles; i++) {
        TString fileName = ((TObjString*)fileList->At(i))->GetString();
        if (fileName.IsNull()) continue;
        
        checkedFiles++;
        if (checkedFiles % 100 == 0) {  // Progress indicator
            cout << "Checked " << checkedFiles << "/" << totalFiles << " files, found " << validFiles << " valid files so far..." << endl;
        }
        
        // Check if this file has the tree we want
        TFile* testFile = TFile::Open(fileName, "READ");
        if (testFile && !testFile->IsZombie()) {
            TDirectoryFile* dir = (TDirectoryFile*)testFile->Get("fourPiAna");
            if (dir) {
                TTree* tree = (TTree*)dir->Get("ParticleTree");
                if (tree && tree->GetEntries() > 0) {
                    chain->Add(fileName);
                    validFiles++;
                    if (validFiles <= 3) {  // Show first few valid files
                        cout << "  Added: " << fileName << " (entries: " << tree->GetEntries() << ")" << endl;
                    }
                }
            }
            testFile->Close();
        }
    }
    
    cout << "Processing complete!" << endl;
    cout << "Added " << validFiles << " valid files out of " << totalFiles << " total files" << endl;
    
    cout << "Number of entries in chain: " << chain->GetEntries() << endl;
    
    // If no entries, show some diagnostic info
    if (chain->GetEntries() == 0) {
        cout << "No entries found in chain!" << endl;
        return;
    }
    
    // Save the chained tree to a new ROOT file for future analysis
    cout << "Saving chained tree to file..." << endl;
    cout << "Current working directory: " << gSystem->pwd() << endl;

    TString outputFileName = "chained_ParticleTree.root";
    TFile* chainOutFile = TFile::Open(outputFileName, "RECREATE");
    if (!chainOutFile || chainOutFile->IsZombie()) {
        cout << "Failed to create output file: " << outputFileName << endl;
        return;
    }

    chainOutFile->cd();
    TTree* outTree = chain->CloneTree(0); // Clone structure only; fill entries with progress below
    if (!outTree) {
        cout << "CloneTree failed; output file will not be written." << endl;
        chainOutFile->Close();
        return;
    }

    Long64_t nEntries = chain->GetEntries();
    TStopwatch stopwatch;
    stopwatch.Start();

    Long64_t progressStep = nEntries / 100;
    if (progressStep < 1) {
        progressStep = 1;
    }
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        if (chain->GetEntry(entry) <= 0) {
            continue;
        }

        outTree->Fill();

        if ((entry + 1) % progressStep == 0 || entry + 1 == nEntries) {
            double elapsed = stopwatch.RealTime();
            stopwatch.Continue();
            double doneFraction = static_cast<double>(entry + 1) / static_cast<double>(nEntries);
            double etaSeconds = (doneFraction > 0.0) ? elapsed * (1.0 - doneFraction) / doneFraction : 0.0;
            cout << "Writing output: " << (entry + 1) << "/" << nEntries
                 << " entries (" << TMath::Nint(doneFraction * 100.0) << "%)"
                 << ", elapsed " << Form("%.1f", elapsed) << " s"
                 << ", ETA " << Form("%.1f", etaSeconds) << " s" << endl;
        }
    }

    Long64_t nWritten = outTree->Write("ParticleTree");
    chainOutFile->Close();
    cout << "Write() returned " << nWritten << endl;
    cout << "Chained tree saved to " << outputFileName << " with " << nEntries << " entries" << endl;
}


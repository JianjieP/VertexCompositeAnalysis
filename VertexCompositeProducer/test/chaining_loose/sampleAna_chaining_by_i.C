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

void sampleAna_chaining_by_i(Int_t i) {
    TChain chain("fourPiAna/ParticleTree");

    TString outDir = "/eos/user/j/jianjie/fourpi/HIRun2023A-14Feb2025-v1_FourPi_0520_loose";
    TString baseDir = Form("/eos/user/j/jianjie/fourpi/HIForward%d/HIForward%d/HIForward%d_HIRun2023A-14Feb2025-v1_FourPi_0520-loose", i, i, i);
    TString outputFileName = Form("%s/chained_ParticleTree_%d.root", outDir.Data(), i);

    if (gSystem->mkdir(outDir, true) != 0 && gSystem->AccessPathName(outDir)) {
        cout << "[i=" << i << "] Failed to create output directory: " << outDir << endl;
        return;
    }

    cout << "[i=" << i << "] Searching in " << baseDir << endl;

    TString findCmd = Form("find %s -type f -name 'fourpi_ana_*.root' 2>/dev/null || true", baseDir.Data());
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

    cout << "[i=" << i << "] Found " << totalFiles << " candidate files." << endl;

    for (Int_t index = 0; index < totalFiles; ++index) {
        TString fileName = ((TObjString*)fileList->At(index))->GetString();
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

        if ((index + 1) % progressStep == 0 || index + 1 == totalFiles) {
            double elapsed = checkWatch.RealTime();
            checkWatch.Continue();
            double fraction = static_cast<double>(index + 1) / static_cast<double>(totalFiles);
            double eta = (fraction > 0.0) ? elapsed * (1.0 - fraction) / fraction : 0.0;
            cout << "[i=" << i << "] Checking files: " << (index + 1) << "/" << totalFiles
                 << " (" << TMath::Nint(fraction * 100.0) << "%)"
                 << ", accepted " << addedFiles
                 << ", skipped " << skippedFiles
                 << ", elapsed " << Form("%.1f", elapsed) << " s"
                 << ", ETA " << Form("%.1f", eta) << " s" << endl;
        }
    }

    cout << "[i=" << i << "] Accepted " << addedFiles << " files out of " << totalFiles << endl;
    cout << "[i=" << i << "] Skipped " << skippedFiles << " files without fourPiAna/ParticleTree" << endl;
    cout << "[i=" << i << "] Number of entries in chain: " << chain.GetEntries() << endl;

    if (chain.GetEntries() == 0) {
        cout << "[i=" << i << "] No entries found in chain!" << endl;
        return;
    }

    cout << "[i=" << i << "] Writing to " << outputFileName << endl;
    TFile* outFile = TFile::Open(outputFileName, "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        cout << "[i=" << i << "] Failed to create " << outputFileName << endl;
        return;
    }

    outFile->cd();
    TStopwatch writeWatch;
    writeWatch.Start();

    TTree* outTree = chain.CloneTree(0);
    if (!outTree) {
        cout << "[i=" << i << "] CloneTree failed." << endl;
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
              cout << "[i=" << i << "] Writing output: " << (entry + 1) << "/" << totalEntries
                  << " (" << TMath::Nint(fraction * 100.0) << "%)"
                 << ", elapsed " << Form("%.1f", elapsed) << " s"
                 << ", ETA " << Form("%.1f", eta) << " s" << endl;
        }
    }

    Long64_t written = outTree->Write("ParticleTree");
    outFile->Close();

    cout << "[i=" << i << "] Write() returned " << written << endl;
    cout << "[i=" << i << "] Saved " << outputFileName << " with " << chain.GetEntries() << " entries" << endl;
}
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"

namespace {

double computeQuantile(std::vector<double>& values, double p) {
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (p <= 0.0) {
        return *std::min_element(values.begin(), values.end());
    }
    if (p >= 1.0) {
        return *std::max_element(values.begin(), values.end());
    }

    std::sort(values.begin(), values.end());
    if (values.size() == 1) {
        return values[0];
    }

    const double index = (values.size() - 1) * p;
    const size_t lo = static_cast<size_t>(std::floor(index));
    const size_t hi = static_cast<size_t>(std::ceil(index));
    const double frac = index - static_cast<double>(lo);
    return values[lo] * (1.0 - frac) + values[hi] * frac;
}

TTree* findParticleTree(TFile* file) {
    if (!file) {
        return nullptr;
    }

    TTree* tree = dynamic_cast<TTree*>(file->Get("ParticleTree"));
    if (tree) {
        return tree;
    }

    TDirectory* dir = dynamic_cast<TDirectory*>(file->Get("fourPiAna"));
    if (!dir) {
        return nullptr;
    }

    return dynamic_cast<TTree*>(dir->Get("ParticleTree"));
}

} // namespace

void compute_zdc_empty_noise_quantile(
    const char* inputFile = "/eos/user/j/jianjie/fourpi/HIRun2025A_FourPi_0410/chained_ParticleTree_EmptyBX.root",
    double quantile = 0.9998
) {
    TFile* file = TFile::Open(inputFile, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: cannot open input file: " << inputFile << std::endl;
        return;
    }

    TTree* tree = findParticleTree(file);
    if (!tree) {
        std::cerr << "Error: cannot find ParticleTree (tried ParticleTree and fourPiAna/ParticleTree)." << std::endl;
        file->Close();
        return;
    }

    float zdcMinus = 0.0f;
    float zdcPlus = 0.0f;

    if (tree->SetBranchAddress("ZDCMinus", &zdcMinus) < 0) {
        std::cerr << "Error: branch ZDCMinus/F not found." << std::endl;
        file->Close();
        return;
    }
    if (tree->SetBranchAddress("ZDCPlus", &zdcPlus) < 0) {
        std::cerr << "Error: branch ZDCPlus/F not found." << std::endl;
        file->Close();
        return;
    }

    const Long64_t nEntries = tree->GetEntries();
    std::vector<double> minusVals;
    std::vector<double> plusVals;
    minusVals.reserve(nEntries);
    plusVals.reserve(nEntries);

    for (Long64_t i = 0; i < nEntries; ++i) {
        if (tree->GetEntry(i) <= 0) {
            continue;
        }
        if (std::isfinite(zdcMinus)) {
            minusVals.push_back(static_cast<double>(zdcMinus));
        }
        if (std::isfinite(zdcPlus)) {
            plusVals.push_back(static_cast<double>(zdcPlus));
        }
    }

    const double qMinus = computeQuantile(minusVals, quantile);
    const double qPlus = computeQuantile(plusVals, quantile);

    std::cout.setf(std::ios::fixed);
    std::cout.precision(12);

    std::cout << "Input file: " << inputFile << std::endl;
    std::cout << "Tree: " << tree->GetName() << std::endl;
    std::cout << "Entries total: " << nEntries
              << ", finite ZDCMinus: " << minusVals.size()
              << ", finite ZDCPlus: " << plusVals.size() << std::endl;
    std::cout << std::endl;

    std::cout << "Requirements for empty ZDC <" << (1.0 - quantile) * 100.0 << "% noise" << std::endl;
    std::cout << "0n+ => ADC<" << qPlus << std::endl;
    std::cout << "0n- => ADC<" << qMinus << std::endl;

    file->Close();
}

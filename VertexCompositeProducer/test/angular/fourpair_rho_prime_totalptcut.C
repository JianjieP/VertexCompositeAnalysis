// Standalone extraction of the rho' (4pi) and pi+ azimuthal-angle difference in the pair rest frame.
//
// What it does:
// - Builds the four possible pion pairings: (13,24), (14,23)
// - Chooses the pair whose invariant mass is closest to the rho mass
// - Boosts the selected pi+ and the 4pi system into the selected pair rest frame
// - Measures the azimuthal-angle difference Delta phi between the 4pi system
//   direction in the lab frame and the selected pi+ direction in the pair rest frame
// - Applies a strict total-pT cut: total pT < 0.1 GeV
// - Fits the Delta phi distribution in bins of the selected pair pT
// - Draws the four cosine-fit parameters A1..A4 versus pair pT

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMath.h"
#include "TLatex.h"
#include "TROOT.h"

namespace {

constexpr double kRhoMass = 0.770;
constexpr double kTotalPtCut = 0.1;
constexpr int kFitOrder = 4;
constexpr bool kOnlyEvenCosine = false;
constexpr int kPhiBins = 24;
constexpr int kFinePtBins = 20;
constexpr double kFinePtMin = 0.0;
constexpr double kFinePtMax = 1.0;
constexpr double kFinePtStep = (kFinePtMax - kFinePtMin) / kFinePtBins;

std::string buildCosineFormula(int fitOrder, bool onlyEvenCosine, std::vector<int> &terms) {
    terms.clear();
    terms.push_back(0);
    for (int n = 1; n <= fitOrder; ++n) {
        if (!onlyEvenCosine || (n % 2 == 0)) {
            terms.push_back(n);
        }
    }

    std::string formula = "[0]";
    for (std::size_t i = 1; i < terms.size(); ++i) {
        formula += "+";
        formula += Form("[0]*[%d]*cos(%d*x)", static_cast<int>(i), terms[i]);
    }
    return formula;
}

TLorentzVector makePion(const std::vector<float> &cand_p,
                        const std::vector<float> &cand_pT,
                        const std::vector<float> &cand_phi,
                        const std::vector<float> &cand_mass,
                        const std::vector<float> &cand_eta,
                        int index) {
    float pT = cand_pT[index];
    float phi = cand_phi[index];
    float px = pT * std::cos(phi);
    float py = pT * std::sin(phi);
    float p = cand_p[index];
    float pzSign = (cand_eta[index] >= 0.f) ? 1.f : -1.f;
    float pz = pzSign * std::sqrt(std::max(0.f, p * p - pT * pT));
    float mass = cand_mass[index];
    float energy = std::sqrt(mass * mass + p * p);
    return TLorentzVector(px, py, pz, energy);
}

TLorentzVector boostToRestFrame(const TLorentzVector &vec, const TLorentzVector &frame) {
    TLorentzVector boosted = vec;
    boosted.Boost(-frame.BoostVector());
    return boosted;
}

double wrapDeltaPhi(double deltaPhi) {
    while (deltaPhi > TMath::Pi()) deltaPhi -= 2.0 * TMath::Pi();
    while (deltaPhi < -TMath::Pi()) deltaPhi += 2.0 * TMath::Pi();
    return deltaPhi;
}

struct PhiFitResult {
    std::array<double, 5> params{};
    std::array<double, 5> errors{};
    double chi2 = std::numeric_limits<double>::quiet_NaN();
    int ndf = 0;
    long long entries = 0;
};

PhiFitResult fitPhiHistogram(TH1F *hist, int fitOrder, bool onlyEvenCosine) {
    PhiFitResult result;
    if (!hist) {
        return result;
    }

    result.entries = static_cast<long long>(hist->GetEntries());
    std::vector<int> terms;
    std::string formula = buildCosineFormula(fitOrder, onlyEvenCosine, terms);
    TF1 fit("fitPhiHistogram", formula.c_str(), -TMath::Pi(), TMath::Pi());
    fit.SetNpx(500);
    for (std::size_t i = 0; i < terms.size(); ++i) {
        fit.SetParName(static_cast<int>(i), Form("A%d", terms[i]));
        fit.SetParameter(static_cast<int>(i), (i == 0) ? std::max(1.0, static_cast<double>(hist->GetMaximum())) : 0.05);
    }
    if (hist->GetEntries() > 0) {
        hist->Fit(&fit, "Q0");
    }
    for (std::size_t i = 0; i < terms.size() && i < result.params.size(); ++i) {
        result.params[i] = fit.GetParameter(static_cast<int>(i));
        result.errors[i] = fit.GetParError(static_cast<int>(i));
    }
    result.chi2 = fit.GetChisquare();
    result.ndf = fit.GetNDF();
    return result;
}

void styleHist(TH1F *hist, int color) {
    if (!hist) return;
    hist->SetStats(0);
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(0.85);
}

struct SelectedPairing {
    int pairIdx = 0;
    int positivePionIdx = 3;
    const char *label = "13/24";
};

SelectedPairing chooseBestPairing(const std::array<double, 4> &masses) {
    int bestIdx = 0;
    double bestDelta = std::abs(masses[0] - kRhoMass);
    for (int idx = 1; idx < 4; ++idx) {
        double delta = std::abs(masses[idx] - kRhoMass);
        if (delta < bestDelta) {
            bestDelta = delta;
            bestIdx = idx;
        }
    }

    switch (bestIdx) {
        case 0: return SelectedPairing{0, 3, "13/24"};
        case 1: return SelectedPairing{1, 4, "13/24"};
        case 2: return SelectedPairing{2, 4, "14/23"};
        case 3: return SelectedPairing{3, 3, "14/23"};
        default: return SelectedPairing{0, 3, "13/24"};
    }
}

TLorentzVector getPairVector(int pairIdx,
                             const TLorentzVector &p1,
                             const TLorentzVector &p2,
                             const TLorentzVector &p3,
                             const TLorentzVector &p4) {
    switch (pairIdx) {
        case 0: return p1 + p3;
        case 1: return p2 + p4;
        case 2: return p1 + p4;
        case 3: return p2 + p3;
        default: return p1 + p3;
    }
}

void drawFitPanel(TPad *pad,
                  TH1F *hist,
                  const char *title,
                  int color,
                  int fitOrder,
                  bool onlyEvenCosine,
                  PhiFitResult &fitOut) {
    pad->cd();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.08);

    styleHist(hist, color);
    hist->SetTitle(title);
    hist->Draw("E1");

    std::vector<int> terms;
    std::string formula = buildCosineFormula(fitOrder, onlyEvenCosine, terms);
    fitOut = fitPhiHistogram(hist, fitOrder, onlyEvenCosine);
    TF1 fit(Form("f_panel_%s", hist->GetName()), formula.c_str(), -TMath::Pi(), TMath::Pi());
    fit.SetNpx(500);
    for (std::size_t i = 0; i < terms.size(); ++i) {
        fit.SetParName(static_cast<int>(i), Form("A%d", terms[i]));
        fit.SetParameter(static_cast<int>(i), fitOut.params[i]);
    }
    fit.SetLineColor(kRed + 1);
    fit.SetLineWidth(2);
    fit.Draw("SAME");

    hist->GetListOfFunctions()->Add((TF1 *)fit.Clone(Form("fit_%s_saved", hist->GetName())));

    auto legend = new TLegend(0.12, 0.53, 0.52, 0.88);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.028);
    legend->AddEntry(hist, title, "lep");
    legend->AddEntry(&fit, "Cosine fit", "l");
    for (std::size_t i = 1; i < terms.size(); ++i) {
        legend->AddEntry((TObject *)nullptr,
                         Form("A%d = %.4f #pm %.4f", terms[i], fitOut.params[i], fitOut.errors[i]),
                         "");
    }
    legend->AddEntry((TObject *)nullptr,
                     Form("#chi^{2}/NDF = %.1f/%d", fitOut.chi2, fitOut.ndf),
                     "");
    legend->Draw();

    pad->Modified();
    pad->Update();
}

void saveDeltaPhiPlot(TH1F *hist,
                      const std::string &pngName,
                      const std::string &pdfName,
                      bool generatePng,
                      bool generatePdf,
                      PhiFitResult &fitOut,
                      const char *title,
                      int color,
                      int fitOrder,
                      bool onlyEvenCosine) {
    if (!hist) return;
    TCanvas canvas("c_deltaPhi_rhoPrime", title, 900, 700);
    drawFitPanel(static_cast<TPad *>(&canvas), hist, title, color, fitOrder, onlyEvenCosine, fitOut);
    canvas.Modified();
    canvas.Update();
    if (generatePng) canvas.SaveAs(pngName.c_str());
    if (generatePdf) canvas.SaveAs(pdfName.c_str());
}

void saveA1234VsPtPlots(const std::vector<double> &ptCenters,
                        const std::vector<double> &ptErr,
                        const std::array<std::vector<double>, 5> &aVals,
                        const std::array<std::vector<double>, 5> &aErrs,
                        const std::string &pngDir,
                        bool generatePng,
                        bool generatePdf) {
    if (ptCenters.empty()) return;

    const std::array<int, 4> colors = {kBlue + 2, kMagenta + 2, kGreen + 2, kRed + 1};
    const std::array<const char *, 4> labels = {"A1", "A2", "A3", "A4"};

    for (int harmonic = 1; harmonic <= 4; ++harmonic) {
        TCanvas canvas(Form("c_a%d_vs_pairPt", harmonic), Form("A%d versus pair pT", harmonic), 900, 700);
        canvas.SetLeftMargin(0.12);
        canvas.SetRightMargin(0.05);
        canvas.SetBottomMargin(0.12);
        canvas.SetTopMargin(0.08);
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.12);
        gPad->SetTopMargin(0.08);

        double yMin = 1e9;
        double yMax = -1e9;
        for (std::size_t i = 0; i < ptCenters.size(); ++i) {
            yMin = std::min(yMin, aVals[harmonic][i] - aErrs[harmonic][i]);
            yMax = std::max(yMax, aVals[harmonic][i] + aErrs[harmonic][i]);
        }
        if (yMin > yMax) {
            yMin = -0.1;
            yMax = 0.1;
        }
        double margin = 0.2 * std::max(1e-6, yMax - yMin);

        TH1F frame(Form("h_frame_a%d", harmonic), Form(";pair p_{T} (GeV);A_{%d}", harmonic), 100, kFinePtMin, kFinePtMax);
        frame.SetStats(0);
        frame.SetMinimum(yMin - margin);
        frame.SetMaximum(yMax + margin);
        frame.Draw();

        TGraphErrors graph(static_cast<int>(ptCenters.size()), ptCenters.data(), aVals[harmonic].data(), ptErr.data(), aErrs[harmonic].data());
        graph.SetMarkerStyle(20 + harmonic);
        graph.SetMarkerColor(colors[harmonic - 1]);
        graph.SetLineColor(colors[harmonic - 1]);
        graph.SetLineWidth(2);
        graph.Draw("P SAME");

        auto legend = new TLegend(0.66, 0.78, 0.90, 0.90);
        legend->SetFillStyle(0);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.03);
        legend->AddEntry(&graph, labels[harmonic - 1], "lp");
        legend->Draw();

        canvas.Modified();
        canvas.Update();

        if (generatePng) canvas.SaveAs((pngDir + Form("/A%d_vs_pairPt.png", harmonic)).c_str());
        if (generatePdf) canvas.SaveAs((pngDir + Form("/A%d_vs_pairPt.pdf", harmonic)).c_str());
    }
}

} // namespace

void fourpair_rho_prime_totalptcut(const char *input_path = "/eos/user/j/jianjie/fourpi/HIRun2025A_FourPi_0410/filtered_chained_ParticleTree.root",
                                   const char *output_root = "fourpair_rho_prime_totalptcut.root",
                                   bool generate_png = true,
                                   bool generate_pdf = false) {
    gROOT->SetBatch(kTRUE);

    const std::string png_dir = "fourpair_rho_prime_totalptcut_PNG";
    std::filesystem::create_directories(png_dir);

    TFile *fin = TFile::Open(input_path);
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: cannot open input file: " << input_path << std::endl;
        return;
    }

    TTree *tin = static_cast<TTree *>(fin->Get("ParticleTree"));
    if (!tin) {
        std::cerr << "Error: TTree 'ParticleTree' not found in input file" << std::endl;
        return;
    }

    std::vector<float> *cand_p = nullptr;
    std::vector<float> *cand_pT = nullptr;
    std::vector<float> *cand_phi = nullptr;
    std::vector<float> *cand_mass = nullptr;
    std::vector<float> *cand_eta = nullptr;

    tin->SetBranchAddress("cand_p", &cand_p);
    tin->SetBranchAddress("cand_pT", &cand_pT);
    tin->SetBranchAddress("cand_phi", &cand_phi);
    tin->SetBranchAddress("cand_mass", &cand_mass);
    tin->SetBranchAddress("cand_eta", &cand_eta);

    TFile *fout = new TFile(output_root, "RECREATE");

    TH1F *h_deltaPhi_all = new TH1F("h_deltaPhi_rhoPrime_all", ";#Delta#phi(#rho^{\\prime}_{lab},#pi^{+}_{pairRF});Counts", kPhiBins, -TMath::Pi(), TMath::Pi());
    h_deltaPhi_all->SetStats(0);

    std::vector<TH1F *> h_deltaPhi_ptBins;
    h_deltaPhi_ptBins.reserve(kFinePtBins);
    for (int bin = 0; bin < kFinePtBins; ++bin) {
        double low = kFinePtMin + bin * kFinePtStep;
        double high = low + kFinePtStep;
        h_deltaPhi_ptBins.push_back(new TH1F(Form("h_deltaPhi_rhoPrime_pt_%02d", bin),
                     Form(";#Delta#phi(#rho^{\\prime}_{lab},#pi^{+}_{pairRF});Counts (%.2f<=p_{T}<%.2f)", low, high),
                                             kPhiBins, -TMath::Pi(), TMath::Pi()));
        h_deltaPhi_ptBins.back()->SetStats(0);
    }

    Long64_t nentries = tin->GetEntries();
    std::cout << "Processing " << nentries << " entries..." << std::endl;

    for (Long64_t i = 0; i < nentries; ++i) {
        tin->GetEntry(i);
        if (!cand_p || !cand_pT || !cand_phi || !cand_mass || !cand_eta ||
            cand_p->size() != 5 || cand_pT->size() != 5 || cand_phi->size() != 5 ||
            cand_mass->size() != 5 || cand_eta->size() != 5) {
            continue;
        }

        TLorentzVector p1 = makePion(*cand_p, *cand_pT, *cand_phi, *cand_mass, *cand_eta, 1);
        TLorentzVector p2 = makePion(*cand_p, *cand_pT, *cand_phi, *cand_mass, *cand_eta, 2);
        TLorentzVector p3 = makePion(*cand_p, *cand_pT, *cand_phi, *cand_mass, *cand_eta, 3);
        TLorentzVector p4 = makePion(*cand_p, *cand_pT, *cand_phi, *cand_mass, *cand_eta, 4);

        TLorentzVector total = p1 + p2 + p3 + p4;
        if (total.Pt() >= kTotalPtCut) {
            continue;
        }

        TLorentzVector pair13 = p1 + p3;
        TLorentzVector pair24 = p2 + p4;
        TLorentzVector pair14 = p1 + p4;
        TLorentzVector pair23 = p2 + p3;

        std::array<double, 4> masses = {pair13.M(), pair24.M(), pair14.M(), pair23.M()};
        std::array<double, 4> pairPts = {pair13.Pt(), pair24.Pt(), pair14.Pt(), pair23.Pt()};
        SelectedPairing selected = chooseBestPairing(masses);
        TLorentzVector selectedPair = getPairVector(selected.pairIdx, p1, p2, p3, p4);
        TLorentzVector selectedPiPlus = (selected.positivePionIdx == 3) ? p3 : p4;

        TLorentzVector rhoPrimeLab = total;
        TLorentzVector piPlusRF = boostToRestFrame(selectedPiPlus, selectedPair);
        double deltaPhi = wrapDeltaPhi(piPlusRF.Phi() - rhoPrimeLab.Phi());

        h_deltaPhi_all->Fill(deltaPhi);

        double bestPairPt = pairPts[selected.pairIdx];
        if (bestPairPt >= kFinePtMin && bestPairPt < kFinePtMax) {
            int fineBin = static_cast<int>((bestPairPt - kFinePtMin) / kFinePtStep);
            if (fineBin >= 0 && fineBin < kFinePtBins) {
                h_deltaPhi_ptBins[fineBin]->Fill(deltaPhi);
            }
        }
    }

    PhiFitResult fitAll = fitPhiHistogram(h_deltaPhi_all, kFitOrder, kOnlyEvenCosine);

    saveDeltaPhiPlot(h_deltaPhi_all,
                     png_dir + "/deltaPhi_rhoPrime_all.png",
                     png_dir + "/deltaPhi_rhoPrime_all.pdf",
                     generate_png,
                     generate_pdf,
                     fitAll,
                     "#Delta#phi(#rho^{\\prime}_{lab},#pi^{+}_{pairRF}) with total p_{T}<0.1",
                     kBlue + 2,
                     kFitOrder,
                     kOnlyEvenCosine);

    std::ofstream txt(png_dir + "/A1234_vs_pairPt_fit_results.txt");
    txt << "pt_low pt_high a1 a2 a3 a4 chi2 ndf entries\n";

    std::vector<double> ptCenters;
    std::vector<double> ptErr;
    std::array<std::vector<double>, 5> aVals;
    std::array<std::vector<double>, 5> aErrs;

    for (int bin = 0; bin < kFinePtBins; ++bin) {
        double low = kFinePtMin + bin * kFinePtStep;
        double high = low + kFinePtStep;
        PhiFitResult fit = fitPhiHistogram(h_deltaPhi_ptBins[bin], kFitOrder, kOnlyEvenCosine);
        txt << low << " " << high << " "
            << fit.params[1] << " " << fit.params[2] << " " << fit.params[3] << " " << fit.params[4] << " "
            << fit.chi2 << " " << fit.ndf << " " << fit.entries << "\n";

        if (fit.entries > 0) {
            ptCenters.push_back(0.5 * (low + high));
            ptErr.push_back(0.5 * (high - low));
            for (int harmonic = 1; harmonic <= 4; ++harmonic) {
                aVals[harmonic].push_back(fit.params[harmonic]);
                aErrs[harmonic].push_back(fit.errors[harmonic]);
            }
        }
    }

    saveA1234VsPtPlots(ptCenters,
                       ptErr,
                       aVals,
                       aErrs,
                       png_dir,
                       generate_png,
                       generate_pdf);

    fout->cd();
    h_deltaPhi_all->Write();
    for (auto *hist : h_deltaPhi_ptBins) {
        hist->Write();
    }
    fout->Write();
    fout->Close();
    fin->Close();

    std::cout << "Done. Output written to " << output_root << std::endl;
}

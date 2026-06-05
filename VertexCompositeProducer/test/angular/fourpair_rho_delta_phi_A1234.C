// Standalone rho-pair Delta phi analysis for four-pion events.
//
// What it does:
// - Builds the two possible rho pairings: (13,24) and (14,23)
// - Picks the pairing whose member pair mass is closest to m_rho
// - Measures the azimuthal angle of the two pi+ in the corresponding rho rest frames
// - Fills four histograms for the 2x2 combinations of total pT and pair pT cuts
// - Fits each histogram with a cosine series up to 4th harmonic and reports A1..A4
//
// Output:
// - PNG/PDF images in a new folder
// - A ROOT file with the histograms and fitted summaries

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
constexpr double kPairPtCut = 0.1;
constexpr int kFitOrder = 4;
constexpr bool kOnlyEvenCosine = false;

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

struct FitSummary {
    std::array<double, 5> params{};
    std::array<double, 5> errors{};
    double chi2 = std::numeric_limits<double>::quiet_NaN();
    int ndf = 0;
};

FitSummary fitHistogram(TH1F *hist, const char *fitTag, int fitOrder, bool onlyEvenCosine) {
    FitSummary summary;
    if (!hist) {
        return summary;
    }

    std::vector<int> terms;
    std::string formula = buildCosineFormula(fitOrder, onlyEvenCosine, terms);
    TF1 fit(Form("f_%s", fitTag), formula.c_str(), -TMath::Pi(), TMath::Pi());
    fit.SetNpx(500);
    for (std::size_t i = 0; i < terms.size(); ++i) {
        fit.SetParName(static_cast<int>(i), Form("A%d", terms[i]));
        fit.SetParameter(static_cast<int>(i), (i == 0) ? std::max(1.0, static_cast<double>(hist->GetMaximum())) : 0.05);
    }
    if (hist->GetEntries() > 0) {
        hist->Fit(&fit, "Q0");
    }

    for (std::size_t i = 0; i < terms.size() && i < summary.params.size(); ++i) {
        summary.params[i] = fit.GetParameter(static_cast<int>(i));
        summary.errors[i] = fit.GetParError(static_cast<int>(i));
    }
    summary.chi2 = fit.GetChisquare();
    summary.ndf = fit.GetNDF();
    return summary;
}

void stylePhiHist(TH1F *hist, int color) {
    if (!hist) {
        return;
    }
    hist->SetStats(0);
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(0.85);
}

TLorentzVector boostToPairRestFrame(const TLorentzVector &pion, const TLorentzVector &pair) {
    TLorentzVector boosted = pion;
    boosted.Boost(-pair.BoostVector());
    return boosted;
}

double wrapDeltaPhi(double deltaPhi) {
    while (deltaPhi > TMath::Pi()) deltaPhi -= 2.0 * TMath::Pi();
    while (deltaPhi < -TMath::Pi()) deltaPhi += 2.0 * TMath::Pi();
    return deltaPhi;
}

struct SelectedPairing {
    int pairA = 0;
    int pairB = 1;
    int positivePiA = 3;
    int positivePiB = 4;
    const char *label = "13/24";
};

SelectedPairing choosePairing(const std::array<double, 4> &masses) {
    const double score13_24 = std::min(std::abs(masses[0] - kRhoMass), std::abs(masses[1] - kRhoMass));
    const double score14_23 = std::min(std::abs(masses[2] - kRhoMass), std::abs(masses[3] - kRhoMass));

    if (score14_23 < score13_24) {
        return SelectedPairing{2, 3, 4, 3, "14/23"};
    }
    return SelectedPairing{0, 1, 3, 4, "13/24"};
}

void savePanel(TPad *pad,
               TH1F *hist,
               const char *drawTitle,
               int color,
               FitSummary &fitOut,
               const char *fitTag,
               int fitOrder,
               bool onlyEvenCosine) {
    pad->cd();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.08);

    stylePhiHist(hist, color);
    hist->SetTitle(drawTitle);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->Draw("E1");

    fitOut = fitHistogram(hist, fitTag, fitOrder, onlyEvenCosine);
    std::vector<int> terms;
    std::string formula = buildCosineFormula(fitOrder, onlyEvenCosine, terms);
    TF1 fit(Form("f_draw_%s", fitTag), formula.c_str(), -TMath::Pi(), TMath::Pi());
    fit.SetNpx(500);
    for (std::size_t i = 0; i < terms.size(); ++i) {
        fit.SetParName(static_cast<int>(i), Form("A%d", terms[i]));
        fit.SetParameter(static_cast<int>(i), fitOut.params[i]);
    }
    fit.SetLineColor(kRed + 1);
    fit.SetLineWidth(2);
    fit.Draw("SAME");
    hist->GetListOfFunctions()->Add((TF1 *)fit.Clone(Form("fit_%s_saved", fitTag)));

    auto legend = new TLegend(0.12, 0.53, 0.52, 0.88);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.028);
    legend->AddEntry(hist, drawTitle, "lep");
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

void saveSinglePlot(TH1F *hist,
                    const char *canvasName,
                    const char *canvasTitle,
                    const char *drawTitle,
                    int color,
                    const std::string &pngName,
                    const std::string &pdfName,
                    bool generatePng,
                    bool generatePdf,
                    FitSummary &fitOut,
                    const char *fitTag,
                    int fitOrder,
                    bool onlyEvenCosine) {
    if (!hist) {
        return;
    }

    TCanvas canvas(canvasName, canvasTitle, 900, 700);
    savePanel(static_cast<TPad *>(&canvas), hist, drawTitle, color, fitOut, fitTag, fitOrder, onlyEvenCosine);
    canvas.Modified();
    canvas.Update();
    if (generatePng) {
        canvas.SaveAs(pngName.c_str());
    }
    if (generatePdf) {
        canvas.SaveAs(pdfName.c_str());
    }
}

void saveComparisonPlot(const std::array<FitSummary, 4> &fits,
                        const std::string &pngName,
                        const std::string &pdfName,
                        bool generatePng,
                        bool generatePdf) {
    TCanvas canvas("c_deltaPhi_A1234_fit_comparison", "Comparison of fit results", 1000, 700);
    canvas.SetLeftMargin(0.12);
    canvas.SetRightMargin(0.05);
    canvas.SetBottomMargin(0.12);
    canvas.SetTopMargin(0.08);

    TH1F frame("h_deltaPhi_A1234_fit_frame", ";Harmonic;Fit value", 4, 0.5, 4.5);
    frame.SetStats(0);
    frame.GetXaxis()->SetBinLabel(1, "A1");
    frame.GetXaxis()->SetBinLabel(2, "A2");
    frame.GetXaxis()->SetBinLabel(3, "A3");
    frame.GetXaxis()->SetBinLabel(4, "A4");

    double yMin = 0.0;
    double yMax = 0.0;
    bool firstValue = true;
    for (const auto &fit : fits) {
        for (int n = 1; n <= 4; ++n) {
            double low = fit.params[n] - fit.errors[n];
            double high = fit.params[n] + fit.errors[n];
            if (firstValue) {
                yMin = low;
                yMax = high;
                firstValue = false;
            } else {
                yMin = std::min(yMin, low);
                yMax = std::max(yMax, high);
            }
        }
    }
    if (firstValue) {
        yMin = -0.1;
        yMax = 0.1;
    }
    const double margin = 0.2 * std::max(1e-6, yMax - yMin);
    frame.SetMinimum(yMin - margin);
    frame.SetMaximum(yMax + margin);
    frame.Draw();

    std::array<std::array<double, 4>, 4> xVals{};
    std::array<std::array<double, 4>, 4> yVals{};
    std::array<std::array<double, 4>, 4> xErrs{};
    std::array<std::array<double, 4>, 4> yErrs{};
    for (int cat = 0; cat < 4; ++cat) {
        for (int n = 0; n < 4; ++n) {
            xVals[cat][n] = static_cast<double>(n + 1);
            yVals[cat][n] = fits[cat].params[n + 1];
            xErrs[cat][n] = 0.0;
            yErrs[cat][n] = fits[cat].errors[n + 1];
        }
    }

    const std::array<int, 4> colors = {kBlack, kBlue + 2, kMagenta + 2, kRed + 1};
    const std::array<int, 4> markers = {20, 21, 22, 23};
    const std::array<const char *, 4> labels = {"all", "total pT < 0.1", "pair pT < 0.1", "both cuts"};

    auto legend = new TLegend(0.68, 0.62, 0.93, 0.88);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.03);

    for (int cat = 0; cat < 4; ++cat) {
        auto graph = new TGraphErrors(4, xVals[cat].data(), yVals[cat].data(), xErrs[cat].data(), yErrs[cat].data());
        graph->SetMarkerStyle(markers[cat]);
        graph->SetMarkerColor(colors[cat]);
        graph->SetLineColor(colors[cat]);
        graph->SetLineWidth(2);
        graph->Draw(cat == 0 ? "P SAME" : "P SAME");
        legend->AddEntry(graph, labels[cat], "lp");
    }

    legend->Draw();

    TLatex note;
    note.SetNDC();
    note.SetTextSize(0.03);
    note.DrawLatex(0.14, 0.86, "Errors are fit parameter uncertainties");

    if (generatePng) {
        canvas.SaveAs(pngName.c_str());
    }
    if (generatePdf) {
        canvas.SaveAs(pdfName.c_str());
    }
}

} // namespace

void fourpair_rho_delta_phi_A1234(const char *input_path = "/eos/user/j/jianjie/fourpi/HIRun2025A_FourPi_0410/filtered_chained_ParticleTree.root",
                                  const char *output_root = "fourpair_rho_delta_phi_A1234.root",
                                  bool generate_png = true,
                                  bool generate_pdf = false) {
    gROOT->SetBatch(kTRUE);

    const std::string outputDir = "fourpair_rho_delta_phi_A1234_PNG";
    std::filesystem::create_directories(outputDir);

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

    const int phiBins = 24;
    const double phiMin = -TMath::Pi();
    const double phiMax = TMath::Pi();

    TH1F *h_deltaPhi_all = new TH1F("h_deltaPhi_all", ";#Delta#phi;Counts", phiBins, phiMin, phiMax);
    TH1F *h_deltaPhi_totalPt = new TH1F("h_deltaPhi_totalPt", ";#Delta#phi;Counts", phiBins, phiMin, phiMax);
    TH1F *h_deltaPhi_pairPt = new TH1F("h_deltaPhi_pairPt", ";#Delta#phi;Counts", phiBins, phiMin, phiMax);
    TH1F *h_deltaPhi_both = new TH1F("h_deltaPhi_both", ";#Delta#phi;Counts", phiBins, phiMin, phiMax);

    h_deltaPhi_all->SetStats(0);
    h_deltaPhi_totalPt->SetStats(0);
    h_deltaPhi_pairPt->SetStats(0);
    h_deltaPhi_both->SetStats(0);

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

        TLorentzVector pair13 = p1 + p3;
        TLorentzVector pair24 = p2 + p4;
        TLorentzVector pair14 = p1 + p4;
        TLorentzVector pair23 = p2 + p3;
        TLorentzVector total = p1 + p2 + p3 + p4;

        std::array<double, 4> invMasses = {pair13.M(), pair24.M(), pair14.M(), pair23.M()};
        std::array<double, 4> pairPts = {pair13.Pt(), pair24.Pt(), pair14.Pt(), pair23.Pt()};

        SelectedPairing pairing = choosePairing(invMasses);
        TLorentzVector selectedPairA = (pairing.pairA == 0) ? pair13 : pair14;
        TLorentzVector selectedPairB = (pairing.pairB == 1) ? pair24 : pair23;
        TLorentzVector selectedPiA = (pairing.positivePiA == 3) ? p3 : p4;
        TLorentzVector selectedPiB = (pairing.positivePiB == 3) ? p3 : p4;

        TLorentzVector piA_rest = boostToPairRestFrame(selectedPiA, selectedPairA);
        TLorentzVector piB_rest = boostToPairRestFrame(selectedPiB, selectedPairB);
        double deltaPhi = wrapDeltaPhi(piA_rest.Phi() - piB_rest.Phi());

        const bool passTotalPt = total.Pt() < kTotalPtCut;
        const bool passPairPt = std::max(pairPts[pairing.pairA], pairPts[pairing.pairB]) < kPairPtCut;

        h_deltaPhi_all->Fill(deltaPhi);
        if (passTotalPt) {
            h_deltaPhi_totalPt->Fill(deltaPhi);
        }
        if (passPairPt) {
            h_deltaPhi_pairPt->Fill(deltaPhi);
        }
        if (passTotalPt && passPairPt) {
            h_deltaPhi_both->Fill(deltaPhi);
        }
    }

    FitSummary fitAll;
    FitSummary fitTotal;
    FitSummary fitPair;
    FitSummary fitBoth;

    saveSinglePlot(h_deltaPhi_all,
                   "c_deltaPhi_A1234_all",
                   "Delta phi in rho rest frame - all events",
                   "All events",
                   kBlack,
                   outputDir + "/delta_phi_A1234_all.png",
                   outputDir + "/delta_phi_A1234_all.pdf",
                   generate_png,
                   generate_pdf,
                   fitAll,
                   "all",
                   kFitOrder,
                   kOnlyEvenCosine);

    saveSinglePlot(h_deltaPhi_totalPt,
                   "c_deltaPhi_A1234_totalPt",
                   "Delta phi in rho rest frame - total pT cut",
                   Form("total p_{T} < %.1f GeV", kTotalPtCut),
                   kBlue + 2,
                   outputDir + "/delta_phi_A1234_totalPt.png",
                   outputDir + "/delta_phi_A1234_totalPt.pdf",
                   generate_png,
                   generate_pdf,
                   fitTotal,
                   "totalPt",
                   kFitOrder,
                   kOnlyEvenCosine);

    saveSinglePlot(h_deltaPhi_pairPt,
                   "c_deltaPhi_A1234_pairPt",
                   "Delta phi in rho rest frame - pair pT cut",
                   Form("pair p_{T} < %.1f GeV", kPairPtCut),
                   kMagenta + 2,
                   outputDir + "/delta_phi_A1234_pairPt.png",
                   outputDir + "/delta_phi_A1234_pairPt.pdf",
                   generate_png,
                   generate_pdf,
                   fitPair,
                   "pairPt",
                   kFitOrder,
                   kOnlyEvenCosine);

    saveSinglePlot(h_deltaPhi_both,
                   "c_deltaPhi_A1234_both",
                   "Delta phi in rho rest frame - both cuts",
                   Form("total p_{T} < %.1f GeV, pair p_{T} < %.1f GeV", kTotalPtCut, kPairPtCut),
                   kRed + 1,
                   outputDir + "/delta_phi_A1234_both.png",
                   outputDir + "/delta_phi_A1234_both.pdf",
                   generate_png,
                   generate_pdf,
                   fitBoth,
                   "bothCuts",
                   kFitOrder,
                   kOnlyEvenCosine);

    std::array<FitSummary, 4> fitSummaries = {fitAll, fitTotal, fitPair, fitBoth};
    saveComparisonPlot(fitSummaries,
                       outputDir + "/delta_phi_A1234_fit_comparison.png",
                       outputDir + "/delta_phi_A1234_fit_comparison.pdf",
                       generate_png,
                       generate_pdf);

    std::ofstream report(outputDir + "/A1234_fit_results.txt");
    report << "category entries a1 a2 a3 a4 chi2 ndf\n";
    report << "all " << h_deltaPhi_all->GetEntries() << " "
           << fitAll.params[1] << " " << fitAll.params[2] << " " << fitAll.params[3] << " " << fitAll.params[4] << " "
           << fitAll.chi2 << " " << fitAll.ndf << "\n";
    report << "totalPt " << h_deltaPhi_totalPt->GetEntries() << " "
           << fitTotal.params[1] << " " << fitTotal.params[2] << " " << fitTotal.params[3] << " " << fitTotal.params[4] << " "
           << fitTotal.chi2 << " " << fitTotal.ndf << "\n";
    report << "pairPt " << h_deltaPhi_pairPt->GetEntries() << " "
           << fitPair.params[1] << " " << fitPair.params[2] << " " << fitPair.params[3] << " " << fitPair.params[4] << " "
           << fitPair.chi2 << " " << fitPair.ndf << "\n";
    report << "bothCuts " << h_deltaPhi_both->GetEntries() << " "
           << fitBoth.params[1] << " " << fitBoth.params[2] << " " << fitBoth.params[3] << " " << fitBoth.params[4] << " "
           << fitBoth.chi2 << " " << fitBoth.ndf << "\n";

    fout->cd();
    h_deltaPhi_all->Write();
    h_deltaPhi_totalPt->Write();
    h_deltaPhi_pairPt->Write();
    h_deltaPhi_both->Write();
    fout->Write();
    fout->Close();
    fin->Close();

    std::cout << "Finished. Output written to " << outputDir << std::endl;
}

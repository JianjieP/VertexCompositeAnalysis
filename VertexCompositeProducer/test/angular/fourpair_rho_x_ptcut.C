// Standalone extraction of the rho-closest x-frame pion angle plots with pair-pT cuts.
// Input: ParticleTree with cand_p, cand_pT, cand_phi, cand_mass, cand_eta branches.
// Output: ROOT file with the histograms and optional PNG/PDF plots.

#include <array>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <limits>
#include <vector>

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMath.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TROOT.h"

namespace {

std::string buildCosineFormula(int fitOrder, bool onlyEvenCosine, std::vector<int> &terms);

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
    float pz_sign = (cand_eta[index] >= 0.f) ? 1.f : -1.f;
    float pz = pz_sign * std::sqrt(std::max(0.f, p * p - pT * pT));
    float mass = cand_mass[index];
    float energy = std::sqrt(mass * mass + p * p);
    return TLorentzVector(px, py, pz, energy);
}

enum NeutronState {
    STATE_0n0n = 0,
    STATE_0nXn = 1,
    STATE_XnXn = 2
};

constexpr double kZdcPlusThreshold = 21163.969742187866;
constexpr double kZdcMinusThreshold = 2006.606927734500;

NeutronState classifyNeutronState(float zdcPlus, float zdcMinus) {
    bool plusBreak = zdcPlus >= kZdcPlusThreshold;
    bool minusBreak = zdcMinus >= kZdcMinusThreshold;
    if (!plusBreak && !minusBreak) return STATE_0n0n;
    if (plusBreak && minusBreak) return STATE_XnXn;
    return STATE_0nXn;
}

const char* getNeutronStateLabel(NeutronState state) {
    switch (state) {
        case STATE_0n0n: return "0n0n";
        case STATE_0nXn: return "0nXn";
        case STATE_XnXn: return "XnXn";
        default: return "unknown";
    }
}

int getNeutronStateColor(NeutronState state, bool isBest) {
    if (isBest) {
        if (state == STATE_0n0n) return kBlue + 2;
        if (state == STATE_0nXn) return kGreen + 2;
        return kRed + 1;
    }
    if (state == STATE_0n0n) return kMagenta + 2;
    if (state == STATE_0nXn) return kCyan + 2;
    return kOrange + 7;
}

struct PhiFitResult {
    std::array<double, 5> params{};
    std::array<double, 5> paramErrs{};
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
        result.paramErrs[i] = fit.GetParError(static_cast<int>(i));
    }
    result.chi2 = fit.GetChisquare();
    result.ndf = fit.GetNDF();
    return result;
}

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

void stylePhiHist(TH1F *hist, int color) {
    hist->SetStats(0);
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(0.9);
}

void savePhiPlot(TH1F *hist,
                 const char *canvasName,
                 const char *canvasTitle,
                 const char *legendTitle,
                 const std::string &pngName,
                 const std::string &pdfName,
                 bool generatePng,
                 bool generatePdf,
                 int color,
                 int fitOrder,
                 bool onlyEvenCosine,
                 const char *noteText) {
    if (!hist) return;
    TCanvas canvas(canvasName, canvasTitle, 900, 700);
    stylePhiHist(hist, color);
    hist->Draw("E1");

    std::vector<int> terms;
    std::string formula = buildCosineFormula(fitOrder, onlyEvenCosine, terms);
    TF1 fit(Form("f_%s", canvasName), formula.c_str(), -TMath::Pi(), TMath::Pi());
    fit.SetNpx(500);
    for (std::size_t i = 0; i < terms.size(); ++i) {
        fit.SetParName(static_cast<int>(i), Form("A%d", terms[i]));
        fit.SetParameter(static_cast<int>(i), (i == 0) ? std::max(1.0, static_cast<double>(hist->GetMaximum())) : 0.05);
    }
    if (hist->GetEntries() > 0) {
        hist->Fit(&fit, "Q0");
    }
    fit.SetLineColor(kRed);
    fit.SetLineWidth(2);
    fit.Draw("SAME");

    auto legend = new TLegend(0.13, 0.58, 0.43, 0.88);
    legend->SetTextSize(0.03);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->AddEntry(hist, legendTitle, "p");
    legend->AddEntry(&fit, "Cosine fit", "l");
    for (std::size_t i = 0; i < terms.size(); ++i) {
        legend->AddEntry((TObject *)nullptr,
                         Form("A%d = %.3f #pm %.3f", terms[i], fit.GetParameter(static_cast<int>(i)), fit.GetParError(static_cast<int>(i))),
                         "");
    }
    legend->AddEntry((TObject *)nullptr,
                     Form("#chi^{2}/NDF = %.1f/%d", fit.GetChisquare(), fit.GetNDF()),
                     "");
    legend->Draw();

    if (noteText && noteText[0] != '\0') {
        TLatex note;
        note.SetNDC();
        note.SetTextSize(0.03);
        note.DrawLatex(0.60, 0.86, noteText);
    }

    if (generatePng) canvas.SaveAs(pngName.c_str());
    if (generatePdf) canvas.SaveAs(pdfName.c_str());
}

void savePhiPlotWithNeutronStates(std::array<TH1F*, 3> hists,
                                  const char *canvasName,
                                  const char *canvasTitle,
                                  const std::string &pngName,
                                  const std::string &pdfName,
                                  bool generatePng,
                                  bool generatePdf,
                                  bool isBest,
                                  int fitOrder,
                                  bool onlyEvenCosine,
                                  const char *noteText) {
    TCanvas canvas(canvasName, canvasTitle, 1000, 700);

    double maxVal = 0.0;
    for (int i = 0; i < 3; ++i) {
        if (hists[i]) maxVal = std::max(maxVal, static_cast<double>(hists[i]->GetMaximum()));
    }

    auto legend = new TLegend(0.70, 0.55, 0.95, 0.88);
    legend->SetTextSize(0.03);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);

    std::vector<int> terms;
    std::string formula = buildCosineFormula(fitOrder, onlyEvenCosine, terms);
    bool firstDrawn = false;

    for (int i = 0; i < 3; ++i) {
        if (!hists[i] || hists[i]->GetEntries() == 0) continue;

        NeutronState state = static_cast<NeutronState>(i);
        int color = getNeutronStateColor(state, isBest);
        stylePhiHist(hists[i], color);
        hists[i]->SetMarkerStyle(20 + i);
        hists[i]->SetLineWidth(2);

        if (!firstDrawn) {
            hists[i]->SetTitle(canvasTitle); // avoid using the first neutron state's title (e.g. 0n0n) as the plot title
            hists[i]->SetMinimum(0.0);
            hists[i]->SetMaximum(std::max(1.0, maxVal * 1.30));
            hists[i]->Draw("E1");
            firstDrawn = true;
        } else {
            hists[i]->Draw("E1 SAME");
        }

        TF1* fit = new TF1(Form("f_%d_%s", i, canvasName), formula.c_str(), -TMath::Pi(), TMath::Pi());
        fit->SetNpx(500);
        for (std::size_t j = 0; j < terms.size(); ++j) {
            fit->SetParameter(static_cast<int>(j), (j == 0) ? std::max(1.0, static_cast<double>(hists[i]->GetMaximum())) : 0.05);
        }
        hists[i]->Fit(fit, "Q0");
        fit->SetLineColor(color);
        fit->SetLineWidth(2);
        fit->Draw("SAME");

        legend->AddEntry(hists[i], Form("%s (entries=%lld)", getNeutronStateLabel(state), static_cast<long long>(hists[i]->GetEntries())), "p");
    }

    legend->Draw();

    if (noteText && noteText[0] != '\0') {
        TLatex note;
        note.SetNDC();
        note.SetTextSize(0.03);
        note.DrawLatex(0.15, 0.86, noteText);
    }

    if (generatePng) canvas.SaveAs(pngName.c_str());
    if (generatePdf) canvas.SaveAs(pdfName.c_str());
}

void saveAnVsPtPlot(const std::vector<double> &ptCenters,
                    const std::vector<double> &ptErr,
                    const std::vector<double> &anBest,
                    const std::vector<double> &anBestErr,
                    const std::vector<double> &anOther,
                    const std::vector<double> &anOtherErr,
                    int harmonic,
                    const std::string &pngName,
                    const std::string &pdfName,
                    bool generatePng,
                    bool generatePdf) {
    if (ptCenters.empty()) return;

    TCanvas canvas(Form("c_a%d_vs_pt", harmonic), Form("A_{%d} vs p_{T}", harmonic), 1000, 700);

    double yMin = 1e9;
    double yMax = -1e9;
    for (std::size_t i = 0; i < ptCenters.size(); ++i) {
        yMin = std::min(yMin, std::min(anBest[i] - anBestErr[i], anOther[i] - anOtherErr[i]));
        yMax = std::max(yMax, std::max(anBest[i] + anBestErr[i], anOther[i] + anOtherErr[i]));
    }
    if (yMin > yMax) { yMin = -0.1; yMax = 0.1; }
    double margin = 0.2 * std::max(1e-6, yMax - yMin);

    TH1F frame(Form("h_frame_a%d", harmonic), Form(";pair p_{T};A_{%d}", harmonic), 100, 0.0, 1.00);
    frame.SetMinimum(yMin - margin);
    frame.SetMaximum(yMax + margin);
    frame.SetStats(0);
    frame.Draw();

    TGraphErrors gBest(static_cast<int>(ptCenters.size()), ptCenters.data(), anBest.data(), ptErr.data(), anBestErr.data());
    gBest.SetMarkerStyle(20);
    gBest.SetMarkerColor(kBlue + 2);
    gBest.SetLineColor(kBlue + 2);
    gBest.SetLineWidth(2);
    gBest.Draw("P SAME");

    TGraphErrors gOther(static_cast<int>(ptCenters.size()), ptCenters.data(), anOther.data(), ptErr.data(), anOtherErr.data());
    gOther.SetMarkerStyle(21);
    gOther.SetMarkerColor(kMagenta + 2);
    gOther.SetLineColor(kMagenta + 2);
    gOther.SetLineWidth(2);
    gOther.Draw("P SAME");

    auto legend = new TLegend(0.65, 0.72, 0.92, 0.88);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.032);
    legend->SetEntrySeparation(0.005);
    legend->AddEntry(&gBest, Form("best (A_{%d})", harmonic), "lp");
    legend->AddEntry(&gOther, Form("other (A_{%d})", harmonic), "lp");
    legend->Draw();

    if (generatePng) canvas.SaveAs(pngName.c_str());
    if (generatePdf) canvas.SaveAs(pdfName.c_str());
}

void saveAnVsPtNeutronPlot(const std::vector<double> &ptCenters,
                           const std::vector<double> &ptErr,
                           const std::array<std::vector<double>, 3> &an,
                           const std::array<std::vector<double>, 3> &anErr,
                           int harmonic,
                           bool isBest,
                           const std::string &pngName,
                           const std::string &pdfName,
                           bool generatePng,
                           bool generatePdf) {
    if (ptCenters.empty()) return;

    double yMin = 1e9;
    double yMax = -1e9;
    bool hasPoint = false;
    for (int state = 0; state < 3; ++state) {
        for (std::size_t i = 0; i < an[state].size(); ++i) {
            yMin = std::min(yMin, an[state][i] - anErr[state][i]);
            yMax = std::max(yMax, an[state][i] + anErr[state][i]);
            hasPoint = true;
        }
    }
    if (!hasPoint) return;
    if (yMin > yMax) { yMin = -0.1; yMax = 0.1; }
    double margin = 0.2 * std::max(1e-6, yMax - yMin);

    const char *which = isBest ? "closest #rho" : "counterpart";
    const char *noteWhich = isBest ? "best" : "other";
    TCanvas canvas(Form("c_a%d_vs_pt_neutron_%s", harmonic, isBest ? "closest" : "counterpart"),
                   Form("A_{%d} vs p_{T}, %s, neutron selection", harmonic, which), 1000, 700);
    TH1F frame(Form("h_frame_a%d_neutron_%s", harmonic, isBest ? "closest" : "counterpart"),
               Form(";pair p_{T};A_{%d}", harmonic), 100, 0.0, 1.00);
    frame.SetMinimum(yMin - margin);
    frame.SetMaximum(yMax + margin);
    frame.SetStats(0);
    frame.Draw();

    auto legend = new TLegend(0.65, 0.68, 0.92, 0.88);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);

    std::vector<TGraphErrors*> graphs;
    for (int state = 0; state < 3; ++state) {
        if (an[state].empty()) continue;
        auto g = new TGraphErrors(static_cast<int>(ptCenters.size()), ptCenters.data(), an[state].data(), ptErr.data(), anErr[state].data());
        g->SetMarkerStyle(20 + state);
        g->SetMarkerColor(getNeutronStateColor(static_cast<NeutronState>(state), isBest));
        g->SetLineColor(getNeutronStateColor(static_cast<NeutronState>(state), isBest));
        g->SetLineWidth(2);
        g->Draw("P SAME");
        legend->AddEntry(g, Form("%s", getNeutronStateLabel(static_cast<NeutronState>(state))), "lp");
        graphs.push_back(g);
    }
    legend->Draw();

    TLatex note;
    note.SetNDC();
    note.SetTextSize(0.035);
    note.DrawLatex(0.15, 0.86, noteWhich);

    if (generatePng) canvas.SaveAs(pngName.c_str());
    if (generatePdf) canvas.SaveAs(pdfName.c_str());
}

} // namespace

void fourpair_rho_x_ptcut(const char *input_path = "/eos/user/j/jianjie/fourpi/HIRun2025A_FourPi_0410/filtered_chained_ParticleTree.root",
                          const char *output_root = "fourpair_rho_x_ptcut.root",
                          bool generate_png = true,
                          bool generate_pdf = false) {
    const bool export_pion_angles_rho_x_ptcut = false;
    const bool export_pion_angles_rho_x_ptcut_mid_high = false;
    const double rhoMass = 0.770;
    const int piPhiBins = 20;
    const double phiMin = -TMath::Pi();
    const double phiMax = TMath::Pi();
    const int fitOrder = 4;
    const bool onlyEvenCosine = false;
    constexpr int finePtBins = 20;
    constexpr double finePtMin = 0.0;
    constexpr double finePtMax = 1.00;
    constexpr double finePtStep = 0.05;
    const std::string png_dir = "fourpair_rho_x_ptcut_PNG";
    const std::string a1234_dir = png_dir + "/a1234_closest_counterpart";
    const std::string a1234_neutron_dir = png_dir + "/a1234_neutron_selection";

    std::filesystem::create_directories(png_dir);
    std::filesystem::create_directories(a1234_dir);
    std::filesystem::create_directories(a1234_neutron_dir);

    gROOT->SetBatch(kTRUE);

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
    float zdcPlus = 0.f;
    float zdcMinus = 0.f;

    tin->SetBranchAddress("cand_p", &cand_p);
    tin->SetBranchAddress("cand_pT", &cand_pT);
    tin->SetBranchAddress("cand_phi", &cand_phi);
    tin->SetBranchAddress("cand_mass", &cand_mass);
    tin->SetBranchAddress("cand_eta", &cand_eta);
    tin->SetBranchAddress("ZDCPlus", &zdcPlus);
    tin->SetBranchAddress("ZDCMinus", &zdcMinus);

    TFile *fout = new TFile(output_root, "RECREATE");

    TH1F *h_pi_phi_x_rho_best_pt = new TH1F("h_pi_phi_x_rho_best_pt", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_x_rho_best_pt_mid = new TH1F("h_pi_phi_x_rho_best_pt_mid", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_x_rho_best_pt_high = new TH1F("h_pi_phi_x_rho_best_pt_high", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_x_rho_other_pt = new TH1F("h_pi_phi_x_rho_other_pt", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_x_rho_other_pt_mid = new TH1F("h_pi_phi_x_rho_other_pt_mid", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_x_rho_other_pt_high = new TH1F("h_pi_phi_x_rho_other_pt_high", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    std::vector<TH1F *> h_pi_phi_x_rho_best_pt_fine;
    std::vector<TH1F *> h_pi_phi_x_rho_other_pt_fine;
    std::vector<std::array<TH1F *, 3>> h_pi_phi_x_rho_best_pt_fine_by_neutron;
    std::vector<std::array<TH1F *, 3>> h_pi_phi_x_rho_other_pt_fine_by_neutron;
    h_pi_phi_x_rho_best_pt_fine.reserve(finePtBins);
    h_pi_phi_x_rho_other_pt_fine.reserve(finePtBins);
    h_pi_phi_x_rho_best_pt_fine_by_neutron.reserve(finePtBins);
    h_pi_phi_x_rho_other_pt_fine_by_neutron.reserve(finePtBins);

    for (int bin = 0; bin < finePtBins; ++bin) {
        double low = finePtMin + bin * finePtStep;
        double high = low + finePtStep;
        h_pi_phi_x_rho_best_pt_fine.push_back(new TH1F(Form("h_pi_phi_x_rho_best_pt_fine_%02d", bin),
                                                       Form("#pi^{+} #Delta#varphi (closest #rho, %.2f<=p_{T}<%.2f);#Delta#varphi;Counts", low, high),
                                                       piPhiBins, phiMin, phiMax));
        h_pi_phi_x_rho_other_pt_fine.push_back(new TH1F(Form("h_pi_phi_x_rho_other_pt_fine_%02d", bin),
                                                        Form("#pi^{+} #Delta#varphi (counterpart, %.2f<=p_{T}<%.2f);#Delta#varphi;Counts", low, high),
                                                        piPhiBins, phiMin, phiMax));
        stylePhiHist(h_pi_phi_x_rho_best_pt_fine.back(), kBlue + 2);
        stylePhiHist(h_pi_phi_x_rho_other_pt_fine.back(), kMagenta + 2);

        std::array<TH1F *, 3> bestByNeutron{};
        std::array<TH1F *, 3> otherByNeutron{};
        for (int state = 0; state < 3; ++state) {
            bestByNeutron[state] = new TH1F(Form("h_pi_phi_x_rho_best_pt_fine_%02d_%s", bin, getNeutronStateLabel(static_cast<NeutronState>(state))),
                                            Form("#pi^{+} #Delta#varphi (closest #rho, %s, %.2f<=p_{T}<%.2f);#Delta#varphi;Counts",
                                                 getNeutronStateLabel(static_cast<NeutronState>(state)), low, high),
                                            piPhiBins, phiMin, phiMax);
            otherByNeutron[state] = new TH1F(Form("h_pi_phi_x_rho_other_pt_fine_%02d_%s", bin, getNeutronStateLabel(static_cast<NeutronState>(state))),
                                             Form("#pi^{+} #Delta#varphi (counterpart, %s, %.2f<=p_{T}<%.2f);#Delta#varphi;Counts",
                                                  getNeutronStateLabel(static_cast<NeutronState>(state)), low, high),
                                             piPhiBins, phiMin, phiMax);
            stylePhiHist(bestByNeutron[state], getNeutronStateColor(static_cast<NeutronState>(state), true));
            stylePhiHist(otherByNeutron[state], getNeutronStateColor(static_cast<NeutronState>(state), false));
        }
        h_pi_phi_x_rho_best_pt_fine_by_neutron.push_back(bestByNeutron);
        h_pi_phi_x_rho_other_pt_fine_by_neutron.push_back(otherByNeutron);
    }

    stylePhiHist(h_pi_phi_x_rho_best_pt, kBlue + 2);
    stylePhiHist(h_pi_phi_x_rho_best_pt_mid, kBlue + 2);
    stylePhiHist(h_pi_phi_x_rho_best_pt_high, kBlue + 2);
    stylePhiHist(h_pi_phi_x_rho_other_pt, kMagenta + 2);
    stylePhiHist(h_pi_phi_x_rho_other_pt_mid, kMagenta + 2);
    stylePhiHist(h_pi_phi_x_rho_other_pt_high, kMagenta + 2);

    Long64_t nentries = tin->GetEntries();
    std::cout << "Processing " << nentries << " entries..." << std::endl;

    for (Long64_t i = 0; i < nentries; ++i) {
        tin->GetEntry(i);
        if (!cand_p || !cand_pT || !cand_phi || !cand_mass || !cand_eta ||
            cand_p->size() != 5 || cand_pT->size() != 5 || cand_phi->size() != 5 || cand_mass->size() != 5 || cand_eta->size() != 5) {
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

        TLorentzVector p1x = p1, p2x = p2, p3x = p3, p4x = p4;
        double rotPhix1 = -pair13.Phi();
        double rotPhix2 = -pair24.Phi();
        double rotPhix3 = -pair14.Phi();
        double rotPhix4 = -pair23.Phi();
        p1x = p3x;
        p2x = p4x;
        p1x.RotateZ(rotPhix1);
        p2x.RotateZ(rotPhix2);
        p4x.RotateZ(rotPhix3);
        p3x.RotateZ(rotPhix4);
        pair13.RotateZ(rotPhix1);
        pair24.RotateZ(rotPhix2);
        pair14.RotateZ(rotPhix3);
        pair23.RotateZ(rotPhix4);

        TVector3 boostVec13x = -pair13.BoostVector();
        TVector3 boostVec24x = -pair24.BoostVector();
        TVector3 boostVec14x = -pair14.BoostVector();
        TVector3 boostVec23x = -pair23.BoostVector();
        p1x.Boost(boostVec13x);
        p2x.Boost(boostVec24x);
        p4x.Boost(boostVec14x);
        p3x.Boost(boostVec23x);

        std::array<double, 4> piPhiX = {p1x.Phi(), p2x.Phi(), p4x.Phi(), p3x.Phi()};

        int bestIdx = 0;
        double bestDelta = std::abs(invMasses[0] - rhoMass);
        for (int idx = 1; idx < 4; ++idx) {
            double delta = std::abs(invMasses[idx] - rhoMass);
            if (delta < bestDelta) {
                bestDelta = delta;
                bestIdx = idx;
            }
        }

        int otherPair = (bestIdx == 0) ? 1 : (bestIdx == 1) ? 0 : (bestIdx == 2) ? 3 : 2;
        double bestPairPt = pairPts[bestIdx];
        double otherPairPt = pairPts[otherPair];
        //bool passTotalPtCut = total.Pt() < 0.1;
        bool passPtCut = export_pion_angles_rho_x_ptcut  && bestPairPt >= 0.0 && bestPairPt < 0.1 && otherPairPt >= 0.0 && otherPairPt < 0.1;
        bool passPtCutMid = export_pion_angles_rho_x_ptcut && export_pion_angles_rho_x_ptcut_mid_high  && bestPairPt >= 0.1 && bestPairPt < 0.2 && otherPairPt >= 0.1 && otherPairPt < 0.2;
        bool passPtCutHigh = export_pion_angles_rho_x_ptcut && export_pion_angles_rho_x_ptcut_mid_high  && bestPairPt >= 0.2 && bestPairPt < 1.0 && otherPairPt >= 0.2 && otherPairPt < 1.0;

        if (passPtCut) {
            h_pi_phi_x_rho_best_pt->Fill(piPhiX[bestIdx]);
            h_pi_phi_x_rho_other_pt->Fill(piPhiX[otherPair]);
        }
        if (bestPairPt >= finePtMin && bestPairPt < finePtMax && otherPairPt >= finePtMin && otherPairPt < finePtMax) {
            int fineBin = static_cast<int>((bestPairPt - finePtMin) / finePtStep);
            if (fineBin >= 0 && fineBin < finePtBins) {
                h_pi_phi_x_rho_best_pt_fine[fineBin]->Fill(piPhiX[bestIdx]);
                h_pi_phi_x_rho_other_pt_fine[fineBin]->Fill(piPhiX[otherPair]);
                NeutronState state = classifyNeutronState(zdcPlus, zdcMinus);
                h_pi_phi_x_rho_best_pt_fine_by_neutron[fineBin][static_cast<int>(state)]->Fill(piPhiX[bestIdx]);
                h_pi_phi_x_rho_other_pt_fine_by_neutron[fineBin][static_cast<int>(state)]->Fill(piPhiX[otherPair]);
            }
        }
        if (passPtCutMid) {
            h_pi_phi_x_rho_best_pt_mid->Fill(piPhiX[bestIdx]);
            h_pi_phi_x_rho_other_pt_mid->Fill(piPhiX[otherPair]);
        }
        if (passPtCutHigh) {
            h_pi_phi_x_rho_best_pt_high->Fill(piPhiX[bestIdx]);
            h_pi_phi_x_rho_other_pt_high->Fill(piPhiX[otherPair]);
        }
    }

    if (export_pion_angles_rho_x_ptcut) {
        savePhiPlot(h_pi_phi_x_rho_best_pt,
                    "c_pi_phi_x_rho_best_pt",
                    "Pion #Delta#varphi (closest #rho, pair p_{T}<0.1)",
                    "#pi^{+} #Delta#varphi (closest #rho, pair p_{T}<0.1)",
                    png_dir + "/pi_phi_x_rho_best_pt.png",
                    "pi_phi_x_rho_best_pt.pdf",
                    generate_png,
                    generate_pdf,
                    kBlue + 2,
                    fitOrder,
                    onlyEvenCosine,
                    "pair p_{T}<0.1");

        savePhiPlot(h_pi_phi_x_rho_other_pt,
                    "c_pi_phi_x_rho_other_pt",
                    "Pion #Delta#varphi (counterpart, pair p_{T}<0.1)",
                    "#pi^{+} #Delta#varphi (counterpart, pair p_{T}<0.1)",
                    png_dir + "/pi_phi_x_rho_other_pt.png",
                    "pi_phi_x_rho_other_pt.pdf",
                    generate_png,
                    generate_pdf,
                    kMagenta + 2,
                    fitOrder,
                    onlyEvenCosine,
                    "pair p_{T}<0.1");
    }

    if (export_pion_angles_rho_x_ptcut_mid_high) {
        savePhiPlot(h_pi_phi_x_rho_best_pt_mid,
                    "c_pi_phi_x_rho_best_pt_mid",
                    "Pion #Delta#varphi (closest #rho, 0.1<pair p_{T}<0.2)",
                    "#pi^{+} #Delta#varphi (closest #rho, 0.1<pair p_{T}<0.2)",
                    png_dir + "/pi_phi_x_rho_best_pt_0p1_0p2.png",
                    "pi_phi_x_rho_best_pt_0p1_0p2.pdf",
                    generate_png,
                    generate_pdf,
                    kBlue + 2,
                    fitOrder,
                    onlyEvenCosine,
                    "0.1<pair p_{T}<0.2");

        savePhiPlot(h_pi_phi_x_rho_best_pt_high,
                    "c_pi_phi_x_rho_best_pt_high",
                    "Pion #Delta#varphi (closest #rho, pair p_{T}>0.2)",
                    "#pi^{+} #Delta#varphi (closest #rho, pair p_{T}>0.2)",
                    png_dir + "/pi_phi_x_rho_best_pt_gt_0p2.png",
                    "pi_phi_x_rho_best_pt_gt_0p2.pdf",
                    generate_png,
                    generate_pdf,
                    kBlue + 2,
                    fitOrder,
                    onlyEvenCosine,
                    "pair p_{T}>0.2");

        savePhiPlot(h_pi_phi_x_rho_other_pt_mid,
                    "c_pi_phi_x_rho_other_pt_mid",
                    "Pion #Delta#varphi (counterpart, 0.1<pair p_{T}<0.2)",
                    "#pi^{+} #Delta#varphi (counterpart, 0.1<pair p_{T}<0.2)",
                    png_dir + "/pi_phi_x_rho_other_pt_0p1_0p2.png",
                    "pi_phi_x_rho_other_pt_0p1_0p2.pdf",
                    generate_png,
                    generate_pdf,
                    kMagenta + 2,
                    fitOrder,
                    onlyEvenCosine,
                    "0.1<pair p_{T}<0.2");

        savePhiPlot(h_pi_phi_x_rho_other_pt_high,
                    "c_pi_phi_x_rho_other_pt_high",
                    "Pion #Delta#varphi (counterpart, pair p_{T}>0.2)",
                    "#pi^{+} #Delta#varphi (counterpart, pair p_{T}>0.2)",
                    png_dir + "/pi_phi_x_rho_other_pt_gt_0p2.png",
                    "pi_phi_x_rho_other_pt_gt_0p2.pdf",
                    generate_png,
                    generate_pdf,
                    kMagenta + 2,
                    fitOrder,
                    onlyEvenCosine,
                    "pair p_{T}>0.2");
    }

    std::ofstream bestTxt(png_dir + "/pi_phi_x_rho_best_pt_fits.txt");
    std::ofstream otherTxt(png_dir + "/pi_phi_x_rho_other_pt_fits.txt");
    bestTxt << "pt_low pt_high a1 a2 a3 a4 chi2 ndf entries\n";
    otherTxt << "pt_low pt_high a1 a2 a3 a4 chi2 ndf entries\n";

    std::vector<double> ptCenters;
    std::vector<double> ptErr;
    std::array<std::vector<double>, 5> anBest;
    std::array<std::vector<double>, 5> anBestErr;
    std::array<std::vector<double>, 5> anOther;
    std::array<std::vector<double>, 5> anOtherErr;
    std::array<std::array<std::vector<double>, 3>, 5> anBestNeutron;
    std::array<std::array<std::vector<double>, 3>, 5> anBestNeutronErr;
    std::array<std::array<std::vector<double>, 3>, 5> anOtherNeutron;
    std::array<std::array<std::vector<double>, 3>, 5> anOtherNeutronErr;

    for (int bin = 0; bin < finePtBins; ++bin) {
        double low = finePtMin + bin * finePtStep;
        double high = low + finePtStep;
        std::string binTag = Form("%02d", bin);
        std::string lowHighTag = Form("%.2f_%.2f", low, high);

        std::string bestCanvas = Form("c_pi_phi_x_rho_best_pt_fine_%s", binTag.c_str());
        std::string otherCanvas = Form("c_pi_phi_x_rho_other_pt_fine_%s", binTag.c_str());
        std::string bestTitle = Form("Pion #Delta#varphi (closest #rho, %.2f<=p_{T}<%.2f)", low, high);
        std::string otherTitle = Form("Pion #Delta#varphi (counterpart, %.2f<=p_{T}<%.2f)", low, high);
        std::string bestLegend = Form("#pi^{+} #Delta#varphi (closest #rho, %.2f<=p_{T}<%.2f)", low, high);
        std::string otherLegend = Form("#pi^{+} #Delta#varphi (counterpart, %.2f<=p_{T}<%.2f)", low, high);
        std::string note = Form("%.2f<=p_{T}<%.2f", low, high);

        savePhiPlot(h_pi_phi_x_rho_best_pt_fine[bin],
                    bestCanvas.c_str(),
                    bestTitle.c_str(),
                    bestLegend.c_str(),
                    png_dir + "/pi_phi_x_rho_best_pt_" + lowHighTag + ".png",
                    "pi_phi_x_rho_best_pt_" + lowHighTag + ".pdf",
                    generate_png,
                    generate_pdf,
                    kBlue + 2,
                    fitOrder,
                    onlyEvenCosine,
                    note.c_str());

        savePhiPlot(h_pi_phi_x_rho_other_pt_fine[bin],
                    otherCanvas.c_str(),
                    otherTitle.c_str(),
                    otherLegend.c_str(),
                    png_dir + "/pi_phi_x_rho_other_pt_" + lowHighTag + ".png",
                    "pi_phi_x_rho_other_pt_" + lowHighTag + ".pdf",
                    generate_png,
                    generate_pdf,
                    kMagenta + 2,
                    fitOrder,
                    onlyEvenCosine,
                    note.c_str());

        savePhiPlotWithNeutronStates(h_pi_phi_x_rho_best_pt_fine_by_neutron[bin],
                                     Form("c_pi_phi_x_rho_best_pt_neutron_%s", binTag.c_str()),
                                     Form("Pion #Delta#varphi (closest #rho, neutron split, %.2f<=p_{T}<%.2f)", low, high),
                                     png_dir + "/pi_phi_x_rho_best_pt_neutron_" + lowHighTag + ".png",
                                     "pi_phi_x_rho_best_pt_neutron_" + lowHighTag + ".pdf",
                                     generate_png,
                                     generate_pdf,
                                     true,
                                     fitOrder,
                                     onlyEvenCosine,
                                     note.c_str());

        savePhiPlotWithNeutronStates(h_pi_phi_x_rho_other_pt_fine_by_neutron[bin],
                                     Form("c_pi_phi_x_rho_other_pt_neutron_%s", binTag.c_str()),
                                     Form("Pion #Delta#varphi (counterpart, neutron split, %.2f<=p_{T}<%.2f)", low, high),
                                     png_dir + "/pi_phi_x_rho_other_pt_neutron_" + lowHighTag + ".png",
                                     "pi_phi_x_rho_other_pt_neutron_" + lowHighTag + ".pdf",
                                     generate_png,
                                     generate_pdf,
                                     false,
                                     fitOrder,
                                     onlyEvenCosine,
                                     note.c_str());

        PhiFitResult bestFit = fitPhiHistogram(h_pi_phi_x_rho_best_pt_fine[bin], fitOrder, onlyEvenCosine);
        PhiFitResult otherFit = fitPhiHistogram(h_pi_phi_x_rho_other_pt_fine[bin], fitOrder, onlyEvenCosine);
        bestTxt << low << " " << high << " "
                << bestFit.params[1] << " " << bestFit.params[2] << " " << bestFit.params[3] << " " << bestFit.params[4] << " "
                << bestFit.chi2 << " " << bestFit.ndf << " " << bestFit.entries << "\n";
        otherTxt << low << " " << high << " "
                 << otherFit.params[1] << " " << otherFit.params[2] << " " << otherFit.params[3] << " " << otherFit.params[4] << " "
                 << otherFit.chi2 << " " << otherFit.ndf << " " << otherFit.entries << "\n";

        if (bestFit.entries > 0 && otherFit.entries > 0) {
            ptCenters.push_back(0.5 * (low + high));
            ptErr.push_back(0.5 * (high - low));
            for (int harmonic = 1; harmonic <= 4; ++harmonic) {
                anBest[harmonic].push_back(bestFit.params[harmonic]);
                anBestErr[harmonic].push_back(bestFit.paramErrs[harmonic]);
                anOther[harmonic].push_back(otherFit.params[harmonic]);
                anOtherErr[harmonic].push_back(otherFit.paramErrs[harmonic]);
            }

            for (int state = 0; state < 3; ++state) {
                PhiFitResult bestNeutronFit = fitPhiHistogram(h_pi_phi_x_rho_best_pt_fine_by_neutron[bin][state], fitOrder, onlyEvenCosine);
                PhiFitResult otherNeutronFit = fitPhiHistogram(h_pi_phi_x_rho_other_pt_fine_by_neutron[bin][state], fitOrder, onlyEvenCosine);
                for (int harmonic = 1; harmonic <= 4; ++harmonic) {
                    anBestNeutron[harmonic][state].push_back(bestNeutronFit.params[harmonic]);
                    anBestNeutronErr[harmonic][state].push_back(bestNeutronFit.paramErrs[harmonic]);
                    anOtherNeutron[harmonic][state].push_back(otherNeutronFit.params[harmonic]);
                    anOtherNeutronErr[harmonic][state].push_back(otherNeutronFit.paramErrs[harmonic]);
                }
            }
        }
    }

    bestTxt.close();
    otherTxt.close();

    for (int harmonic = 1; harmonic <= 4; ++harmonic) {
        saveAnVsPtPlot(ptCenters,
                       ptErr,
                       anBest[harmonic],
                       anBestErr[harmonic],
                       anOther[harmonic],
                       anOtherErr[harmonic],
                       harmonic,
                       a1234_dir + Form("/a%d_vs_pt_closest_counterpart.png", harmonic),
                       Form("a%d_vs_pt_closest_counterpart.pdf", harmonic),
                       generate_png,
                       generate_pdf);

        saveAnVsPtNeutronPlot(ptCenters,
                              ptErr,
                              anBestNeutron[harmonic],
                              anBestNeutronErr[harmonic],
                              harmonic,
                              true,
                              a1234_neutron_dir + Form("/a%d_vs_pt_closest_neutron.png", harmonic),
                              Form("a%d_vs_pt_closest_neutron.pdf", harmonic),
                              generate_png,
                              generate_pdf);

        saveAnVsPtNeutronPlot(ptCenters,
                              ptErr,
                              anOtherNeutron[harmonic],
                              anOtherNeutronErr[harmonic],
                              harmonic,
                              false,
                              a1234_neutron_dir + Form("/a%d_vs_pt_counterpart_neutron.png", harmonic),
                              Form("a%d_vs_pt_counterpart_neutron.pdf", harmonic),
                              generate_png,
                              generate_pdf);
    }

    fout->cd();
    h_pi_phi_x_rho_best_pt->Write();
    h_pi_phi_x_rho_best_pt_mid->Write();
    h_pi_phi_x_rho_best_pt_high->Write();
    h_pi_phi_x_rho_other_pt->Write();
    h_pi_phi_x_rho_other_pt_mid->Write();
    h_pi_phi_x_rho_other_pt_high->Write();
    for (int bin = 0; bin < finePtBins; ++bin) {
        h_pi_phi_x_rho_best_pt_fine[bin]->Write();
        h_pi_phi_x_rho_other_pt_fine[bin]->Write();
        for (int state = 0; state < 3; ++state) {
            h_pi_phi_x_rho_best_pt_fine_by_neutron[bin][state]->Write();
            h_pi_phi_x_rho_other_pt_fine_by_neutron[bin][state]->Write();
        }
    }
    fout->Write();
    fout->Close();
    fin->Close();

    std::cout << "Done. Output written to " << output_root << std::endl;
}
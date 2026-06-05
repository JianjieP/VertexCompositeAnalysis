#include <iostream>
#include <filesystem>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TROOT.h"

namespace {

bool fillHistFromFile(const char* input_path,
                      const char* tree_name,
                      TH1F* h_zdc_plus,
                      TH1F* h_zdc_minus,
                      const char* tag) {
    TFile* fin = TFile::Open(input_path, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: cannot open input file (" << tag << "): " << input_path << std::endl;
        return false;
    }

    TTree* tree = static_cast<TTree*>(fin->Get(tree_name));
    if (!tree) {
        std::cerr << "Error: tree not found in " << tag << ": " << tree_name << std::endl;
        fin->Close();
        return false;
    }

    float zdcPlus = 0.f;
    float zdcMinus = 0.f;
    tree->SetBranchAddress("ZDCPlus", &zdcPlus);
    tree->SetBranchAddress("ZDCMinus", &zdcMinus);

    Long64_t nentries = tree->GetEntries();
    std::cout << "Processing " << tag << " with " << nentries << " entries..." << std::endl;
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        h_zdc_plus->Fill(zdcPlus);
        h_zdc_minus->Fill(zdcMinus);

        if ((i + 1) % 100000 == 0) {
            std::cout << "  [" << tag << "] processed " << (i + 1) << " entries..." << std::endl;
        }
    }

    fin->Close();
    return true;
}


} // namespace

void plot_zdc_distributions(
    const char* input_filtered = "/eos/user/j/jianjie/fourpi/HIRun2025A_FourPi_0410/filtered_chained_ParticleTree.root",
    const char* input_emptybx = "/eos/user/j/jianjie/fourpi/HIRun2025A_FourPi_0410/chained_ParticleTree_EmptyBX.root",
    const char* output_dir = "zdc_distributions",
    const char* tree_name = "ParticleTree") {
    
    std::filesystem::create_directories(output_dir);
    
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    
    constexpr double kZdcPlusThreshold = 21163.969742187866;
    constexpr double kZdcMinusThreshold = 2006.606927734500;

    // Histograms for filtered sample - wider range for full distribution
    TH1F* h_plus_filtered = new TH1F("h_plus_filtered", "ZDC Plus (filtered);ADC;Counts", 500, -100000.f, 1000000.f);
    TH1F* h_minus_filtered = new TH1F("h_minus_filtered", "ZDC Minus (filtered);ADC;Counts", 500, -100000.f, 800000.f);

    // Histograms for EmptyBX sample
    TH1F* h_plus_emptybx = new TH1F("h_plus_emptybx", "ZDC Plus (EmptyBX);ADC;Counts", 500, -100000.f, 1000000.f);
    TH1F* h_minus_emptybx = new TH1F("h_minus_emptybx", "ZDC Minus (EmptyBX);ADC;Counts", 500, -100000.f, 800000.f);

    // Style setup
    h_plus_filtered->SetLineColor(kBlue + 1);
    h_plus_filtered->SetLineWidth(2);
    h_minus_filtered->SetLineColor(kRed + 1);
    h_minus_filtered->SetLineWidth(2);

    h_plus_emptybx->SetLineColor(kAzure + 2);
    h_plus_emptybx->SetLineWidth(2);
    h_minus_emptybx->SetLineColor(kOrange + 7);
    h_minus_emptybx->SetLineWidth(2);

    // Fill histograms from files
    if (!fillHistFromFile(input_filtered, tree_name,
                          h_plus_filtered, h_minus_filtered,
                          "filtered")) {
        return;
    }
    if (!fillHistFromFile(input_emptybx, tree_name,
                          h_plus_emptybx, h_minus_emptybx,
                          "emptybx")) {
        return;
    }

    // ========== ZDC Plus ==========
    
    // Plot 1: Only filtered (Plus)
    {
        TCanvas c_plus_f("c_plus_f", "ZDCPlus filtered only", 900, 700);
        c_plus_f.SetLogy();
        h_plus_filtered->SetMinimum(0.5);
        h_plus_filtered->Draw("HIST");
        c_plus_f.SaveAs(Form("%s/zdc_plus_filtered_only.png", output_dir));
    }

    // Plot 2: Only EmptyBX (Plus)
    {
        TCanvas c_plus_e("c_plus_e", "ZDCPlus EmptyBX only", 900, 700);
        c_plus_e.SetLogy();
        h_plus_emptybx->SetMinimum(0.5);
        h_plus_emptybx->Draw("HIST");
        c_plus_e.SaveAs(Form("%s/zdc_plus_emptybx_only.png", output_dir));
    }

    // Plot 3: Filtered + EmptyBX overlay (Plus)
    {
        TCanvas c_plus_both("c_plus_both", "ZDCPlus filtered vs EmptyBX", 1000, 700);
        c_plus_both.SetLogy();
        h_plus_filtered->SetTitle("ZDCPlus: Filtered vs EmptyBX;ADC;Counts");
        double ymax_plus = h_plus_filtered->GetMaximum();
        if (h_plus_emptybx->GetMaximum() > ymax_plus) ymax_plus = h_plus_emptybx->GetMaximum();
        h_plus_filtered->SetMaximum(ymax_plus * 1.35);
        h_plus_filtered->SetMinimum(0.5);
        h_plus_filtered->Draw("HIST");
        h_plus_emptybx->Draw("HIST SAME");
        auto leg = new TLegend(0.65, 0.70, 0.92, 0.88);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.035);
        leg->AddEntry(h_plus_filtered, "Filtered", "l");
        leg->AddEntry(h_plus_emptybx, "EmptyBX", "l");
        leg->Draw();
        TLine thr_plus(kZdcPlusThreshold, 0.5, kZdcPlusThreshold, ymax_plus * 1.35);
        thr_plus.SetLineStyle(2);
        thr_plus.SetLineWidth(2);
        thr_plus.SetLineColor(kBlack);
        thr_plus.Draw("SAME");
        c_plus_both.SaveAs(Form("%s/zdc_plus_comparison.png", output_dir));
    }

    // ========== ZDC Minus ==========

    // Plot 4: Only filtered (Minus)
    {
        TCanvas c_minus_f("c_minus_f", "ZDCMinus filtered only", 900, 700);
        c_minus_f.SetLogy();
        h_minus_filtered->SetMinimum(0.5);
        h_minus_filtered->Draw("HIST");
        c_minus_f.SaveAs(Form("%s/zdc_minus_filtered_only.png", output_dir));
    }

    // Plot 5: Only EmptyBX (Minus)
    {
        TCanvas c_minus_e("c_minus_e", "ZDCMinus EmptyBX only", 900, 700);
        c_minus_e.SetLogy();
        h_minus_emptybx->SetMinimum(0.5);
        h_minus_emptybx->Draw("HIST");
        c_minus_e.SaveAs(Form("%s/zdc_minus_emptybx_only.png", output_dir));
    }

    // Plot 6: Filtered + EmptyBX overlay (Minus)
    {
        TCanvas c_minus_both("c_minus_both", "ZDCMinus filtered vs EmptyBX", 1000, 700);
        c_minus_both.SetLogy();
        h_minus_filtered->SetTitle("ZDCMinus: Filtered vs EmptyBX;ADC;Counts");
        double ymax_minus = h_minus_filtered->GetMaximum();
        if (h_minus_emptybx->GetMaximum() > ymax_minus) ymax_minus = h_minus_emptybx->GetMaximum();
        h_minus_filtered->SetMaximum(ymax_minus * 1.35);
        h_minus_filtered->SetMinimum(0.5);
        h_minus_filtered->Draw("HIST");
        h_minus_emptybx->Draw("HIST SAME");
        auto leg = new TLegend(0.65, 0.70, 0.92, 0.88);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.035);
        leg->AddEntry(h_minus_filtered, "Filtered", "l");
        leg->AddEntry(h_minus_emptybx, "EmptyBX", "l");
        leg->Draw();
        TLine thr_minus(kZdcMinusThreshold, 0.5, kZdcMinusThreshold, ymax_minus * 1.35);
        thr_minus.SetLineStyle(2);
        thr_minus.SetLineWidth(2);
        thr_minus.SetLineColor(kBlack);
        thr_minus.Draw("SAME");
        c_minus_both.SaveAs(Form("%s/zdc_minus_comparison.png", output_dir));
    }
    
    // ========== Plus vs Minus in filtered ==========
    {
        TCanvas c_filtered_both("c_filtered_both", "ZDCPlus vs ZDCMinus (filtered)", 1000, 700);
        c_filtered_both.SetLogy();
        h_plus_filtered->SetTitle("ZDC: Plus vs Minus (filtered);ADC;Counts");
        h_plus_filtered->SetMinimum(0.5);
        h_plus_filtered->Draw("HIST");
        h_minus_filtered->Draw("HIST SAME");
        auto leg = new TLegend(0.65, 0.70, 0.92, 0.88);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.035);
        leg->AddEntry(h_plus_filtered, "ZDC Plus", "l");
        leg->AddEntry(h_minus_filtered, "ZDC Minus", "l");
        leg->Draw();
        c_filtered_both.SaveAs(Form("%s/zdc_plus_vs_minus_filtered.png", output_dir));
    }
    
    // ========== Plus vs Minus in EmptyBX ==========
    {
        TCanvas c_emptybx_both("c_emptybx_both", "ZDCPlus vs ZDCMinus (EmptyBX)", 1000, 700);
        c_emptybx_both.SetLogy();
        h_plus_emptybx->SetTitle("ZDC: Plus vs Minus (EmptyBX);ADC;Counts");
        h_plus_emptybx->SetMinimum(0.5);
        h_plus_emptybx->Draw("HIST");
        h_minus_emptybx->Draw("HIST SAME");
        double ymax_emptybx = std::max(h_plus_emptybx->GetMaximum(), h_minus_emptybx->GetMaximum());
        TLine thr_plus(kZdcPlusThreshold, 0.5, kZdcPlusThreshold, ymax_emptybx * 1.2);
        thr_plus.SetLineStyle(2);
        thr_plus.SetLineWidth(2);
        thr_plus.SetLineColor(kBlack);
        thr_plus.Draw("SAME");
        TLine thr_minus(kZdcMinusThreshold, 0.5, kZdcMinusThreshold, ymax_emptybx * 1.2);
        thr_minus.SetLineStyle(2);
        thr_minus.SetLineWidth(2);
        thr_minus.SetLineColor(kGray + 2);
        thr_minus.Draw("SAME");
        auto leg = new TLegend(0.65, 0.70, 0.92, 0.88);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.035);
        leg->AddEntry(h_plus_emptybx, "ZDC Plus", "l");
        leg->AddEntry(h_minus_emptybx, "ZDC Minus", "l");
        leg->Draw();
        c_emptybx_both.SaveAs(Form("%s/zdc_plus_vs_minus_emptybx.png", output_dir));
    }
    
    std::cout << "ZDC plots saved to: " << output_dir << std::endl;
    std::cout << "  1. zdc_plus_filtered_only.png" << std::endl;
    std::cout << "  2. zdc_plus_emptybx_only.png" << std::endl;
    std::cout << "  3. zdc_plus_comparison.png" << std::endl;
    std::cout << "  4. zdc_minus_filtered_only.png" << std::endl;
    std::cout << "  5. zdc_minus_emptybx_only.png" << std::endl;
    std::cout << "  6. zdc_minus_comparison.png" << std::endl;
    std::cout << "  7. zdc_plus_vs_minus_filtered.png" << std::endl;
    std::cout << "  8. zdc_plus_vs_minus_emptybx.png" << std::endl;
}

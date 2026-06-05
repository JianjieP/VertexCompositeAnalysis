#include <TComplex.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"

const double kFitMin = 0.8;
const double kFitMax = 4.0;
const double mass_bins = 160;

// Fit initial parameters and limits (edit these values as needed)
struct ParLimit {
  double min;
  double max;
};

const double kF1Init[3] = {1e3, 1.385, 0.431};
const double kF2Init[8] = {1.0, 1.0, 1.385, 0.431, 0.7, 1.663, 0.357, 0.0};
const double kF3Init[12] = {
  1.0,
  1.0, 1.385, 0.431,
  0.7, 1.663, 0.357, -1.22,
  0.05, 2.1, 0.2, -0.5
};

const ParLimit kF2LimM01  = {1.35, 1.55};
const ParLimit kF2LimG1   = {0.00, 1.0};
const ParLimit kF2LimM02  = {1.6, 1.8};
const ParLimit kF2LimG2   = {0.00, 1.0};
const ParLimit kF2LimPhi  = {-TMath::Pi(), TMath::Pi()};

const ParLimit kF3LimM01  = {1.35, 1.55};
const ParLimit kF3LimG1   = {0.00, 1.0};
const ParLimit kF3LimM02  = {1.6, 1.8};
const ParLimit kF3LimG2   = {0.00, 1.0};
const ParLimit kF3LimPhi2 = {-TMath::Pi(), TMath::Pi()};
const ParLimit kF3LimM03  = {2.08, 2.22};
const ParLimit kF3LimG3   = {0.00, 1.0};
const ParLimit kF3LimPhi3 = {-TMath::Pi(), TMath::Pi()};

// Plot axis controls (edit these values as needed)
const double kFit2YMin = -8000.0;
const double kFit2YMax = 6000.0;
const double kFit3CompYMin = -8000.0;
const double kFit3CompYMax = 6000.0;

// 最简单的相对论 BW（常宽）： BW = 1 / (m0^2 - m^2 - i m0 Gamma)
static TComplex BW_rel(double m, double m0, double g0)
{
    const double re = m0*m0 - m*m;
    const double im = - m0*g0;
    return TComplex(1.0, 0.0) / TComplex(re, im);
}

double model_single(double *x, double *p)
{
    const double m = x[0];

// 参数：p0=scale, p1=m0, p2=Gamma0
    const double scale = p[0];
    const double m0    = p[1];
    const double g0    = p[2];

    TComplex amp = BW_rel(m, m0, g0);
    const double val = scale * amp.Rho2(); // |amp|^2

    return val;
}

double model_interf(double *x, double *p)
{
  const double m = x[0];

  // 参数：
  // p0=scale
  // p1=A, p2=m01, p3=G1
  // p4=B, p5=m02, p6=G2
  // p7=phi  (相位，单位：弧度)
  const double scale = p[0];

  const double A  = p[1], m01 = p[2], g1 = p[3];
  const double B  = p[4], m02 = p[5], g2 = p[6];
  const double phi = p[7];

  TComplex bw1 = BW_rel(m, m01, g1);
  TComplex bw2 = BW_rel(m, m02, g2);

  TComplex phase(TMath::Cos(phi), -TMath::Sin(phi)); // e^{-i phi}

  TComplex amp = A*bw1 + phase*(B*bw2);

  return scale * amp.Rho2(); // |amp|^2
}

double comp_bw1(double *x, double *p){
  double m=x[0];
  double scale=p[0], A=p[1], m01=p[2], g1=p[3];
  TComplex bw1 = BW_rel(m,m01,g1);
  return scale * (A*A) * bw1.Rho2();
}

double comp_bw2(double *x, double *p){
  double m=x[0];
  double scale=p[0], B=p[4], m02=p[5], g2=p[6];
  TComplex bw2 = BW_rel(m,m02,g2);
  return scale * (B*B) * bw2.Rho2();
}

double comp_int(double *x, double *p){
  double m=x[0];
  double scale=p[0];
  double A=p[1], m01=p[2], g1=p[3];
  double B=p[4], m02=p[5], g2=p[6];
  double phi=p[7];

  TComplex bw1 = BW_rel(m,m01,g1);
  TComplex bw2 = BW_rel(m,m02,g2);

  // 2AB Re[ bw1 * (e^{-i phi} bw2)^* ]
  TComplex phase(TMath::Cos(phi), -TMath::Sin(phi)); // e^{-i phi}
  TComplex term = bw1 * TComplex::Conjugate(phase*bw2);

  return scale * 2.0*A*B * term.Re();
}

double model_interf3(double *x, double *p)
{
  const double m = x[0];

  // 参数：
  // p0=scale
  // p1=A1, p2=m01, p3=G1
  // p4=A2, p5=m02, p6=G2, p7=phi2
  // p8=A3, p9=m03, p10=G3, p11=phi3
  const double scale = p[0];

  const double A1 = p[1], m01 = p[2], g1 = p[3];
  const double A2 = p[4], m02 = p[5], g2 = p[6], phi2 = p[7];
  const double A3 = p[8], m03 = p[9], g3 = p[10], phi3 = p[11];

  TComplex bw1 = BW_rel(m, m01, g1);
  TComplex bw2 = BW_rel(m, m02, g2);
  TComplex bw3 = BW_rel(m, m03, g3);

  TComplex phase2(TMath::Cos(phi2), -TMath::Sin(phi2));
  TComplex phase3(TMath::Cos(phi3), -TMath::Sin(phi3));

  TComplex amp = A1*bw1 + phase2*(A2*bw2) + phase3*(A3*bw3);
  return scale * amp.Rho2();
}

double comp3_bw1(double *x, double *p){
  const double m = x[0];
  const double scale = p[0], A1 = p[1], m01 = p[2], g1 = p[3];
  TComplex bw1 = BW_rel(m, m01, g1);
  return scale * (A1*A1) * bw1.Rho2();
}

double comp3_bw2(double *x, double *p){
  const double m = x[0];
  const double scale = p[0], A2 = p[4], m02 = p[5], g2 = p[6];
  TComplex bw2 = BW_rel(m, m02, g2);
  return scale * (A2*A2) * bw2.Rho2();
}

double comp3_bw3(double *x, double *p){
  const double m = x[0];
  const double scale = p[0], A3 = p[8], m03 = p[9], g3 = p[10];
  TComplex bw3 = BW_rel(m, m03, g3);
  return scale * (A3*A3) * bw3.Rho2();
}

double comp3_int12(double *x, double *p){
  const double m = x[0];
  const double scale = p[0];
  const double A1 = p[1], m01 = p[2], g1 = p[3];
  const double A2 = p[4], m02 = p[5], g2 = p[6], phi2 = p[7];

  TComplex bw1 = BW_rel(m, m01, g1);
  TComplex bw2 = BW_rel(m, m02, g2);
  TComplex phase2(TMath::Cos(phi2), -TMath::Sin(phi2));
  TComplex term = bw1 * TComplex::Conjugate(phase2*bw2);

  return scale * 2.0*A1*A2 * term.Re();
}

double comp3_int13(double *x, double *p){
  const double m = x[0];
  const double scale = p[0];
  const double A1 = p[1], m01 = p[2], g1 = p[3];
  const double A3 = p[8], m03 = p[9], g3 = p[10], phi3 = p[11];

  TComplex bw1 = BW_rel(m, m01, g1);
  TComplex bw3 = BW_rel(m, m03, g3);
  TComplex phase3(TMath::Cos(phi3), -TMath::Sin(phi3));
  TComplex term = bw1 * TComplex::Conjugate(phase3*bw3);

  return scale * 2.0*A1*A3 * term.Re();
}

double comp3_int23(double *x, double *p){
  const double m = x[0];
  const double scale = p[0];
  const double A2 = p[4], m02 = p[5], g2 = p[6], phi2 = p[7];
  const double A3 = p[8], m03 = p[9], g3 = p[10], phi3 = p[11];

  TComplex bw2 = BW_rel(m, m02, g2);
  TComplex bw3 = BW_rel(m, m03, g3);
  TComplex phase2(TMath::Cos(phi2), -TMath::Sin(phi2));
  TComplex phase3(TMath::Cos(phi3), -TMath::Sin(phi3));
  TComplex term = (phase2*bw2) * TComplex::Conjugate(phase3*bw3);

  return scale * 2.0*A2*A3 * term.Re();
}

  static void write_fit_params(std::ofstream &out, const char *title, TF1 *func)
  {
    out << "[" << title << "]\n";
    out << "chi2 = " << func->GetChisquare()
      << ", ndf = " << func->GetNDF()
      << ", chi2/ndf = "
      << (func->GetNDF() > 0 ? func->GetChisquare()/func->GetNDF() : 0.0)
      << "\n";

    const int npar = func->GetNpar();
    for (int i = 0; i < npar; ++i) {
      out << "p" << i << " (" << func->GetParName(i) << ")"
        << " = " << func->GetParameter(i)
        << " +- " << func->GetParError(i) << "\n";
    }
    out << "\n";
  }

void fit(){
  gStyle->SetOptStat(0); // Hide the default Entries/Mean/RMS stats box

  std::filesystem::path output_dir = "fit_output";
  try {
    std::filesystem::create_directories(output_dir);
  } catch (const std::exception& ex) {
    std::cerr << "Error: cannot create output directory '" << output_dir << "': " << ex.what() << std::endl;
    return;
  }

  auto out_path = [&](const char* filename) {
    return (output_dir / filename).string();
  };

  auto draw_chi2_ndf = [&](TF1* f) {
    if (!f) return;
    const double ndf = f->GetNDF();
    const double chi2 = f->GetChisquare();
    TLatex txt;
    txt.SetNDC();
    txt.SetTextSize(0.03);
    txt.SetTextAlign(31); // right-aligned
    if (ndf > 0.0) {
      txt.DrawLatex(0.88, 0.14, Form("#chi^{2}/NDF = %.2f/%d = %.3f", chi2, (int)ndf, chi2/ndf));
    } else {
      txt.DrawLatex(0.88, 0.14, Form("#chi^{2}/NDF = %.2f/0", chi2));
    }
  };

  const char* input_path = "/eos/user/j/jianjie/fourpi/HIRun2025A_FourPi_0410/filtered_chained_ParticleTree.root";
  TFile *file = TFile::Open(input_path);
  if (!file || file->IsZombie()) {
    std::cerr << "Error: cannot open input file: " << input_path << std::endl;
    return;
  }
    TTree *tree = (TTree*)file->Get("ParticleTree");
  if (!tree) {
    std::cerr << "Error: TTree 'ParticleTree' not found in input file" << std::endl;
    file->Close();
    return;
  }
    //std::vector<float> *cand_mass = nullptr;
    //tree->SetBranchAddress("cand_mass", &cand_mass);
    
    TCanvas *c1 = new TCanvas("c1", "Mass Fit", 800, 600);
    TH1D *h = new TH1D("h", "Invariant Mass", mass_bins, kFitMin, kFitMax);
    tree->Draw("cand_mass[0]>>h", "cand_pT[0]<0.15"); // 仅当 cand_pT[0]<0.15 时填充 cand_mass[0]
    h->Draw("E");
    c1->SaveAs(out_path("h.png").c_str());

    TF1 *f1 = new TF1("f1", model_single, kFitMin, kFitMax, 3);
    f1->SetParNames("scale","m0","Gamma0");
    f1->SetParameters(kF1Init[0], kF1Init[1], kF1Init[2]);

    h->Fit(f1, "R I S");  // R:用范围, I:bin积分(推荐), S:返回结果
    h->Draw("E");
    draw_chi2_ndf(f1);
    c1->SaveAs(out_path("h_fit1.png").c_str());
    
    TF1 *f2 = new TF1("f2", model_interf, kFitMin, kFitMax, 8);
    f2->SetParNames("scale","A","m01","G1","B","m02","G2","phi");
    f2->SetParameters(kF2Init[0], kF2Init[1], kF2Init[2], kF2Init[3],
              kF2Init[4], kF2Init[5], kF2Init[6], kF2Init[7]);

    //f2->SetParLimits(0, 1e-6, 1e6);
    //f2->SetParLimits(1, -100.0, 100.0);
    f2->SetParLimits(2, kF2LimM01.min, kF2LimM01.max);
    f2->SetParLimits(3, kF2LimG1.min, kF2LimG1.max);
    //f2->SetParLimits(4, -100.0, 100.0);
    f2->SetParLimits(5, kF2LimM02.min, kF2LimM02.max);
    f2->SetParLimits(6, kF2LimG2.min, kF2LimG2.max);
    f2->SetParLimits(7, kF2LimPhi.min, kF2LimPhi.max);
    h->Fit(f2, "R I S");

    TF1 *fbw1 = new TF1("fbw1", comp_bw1, kFitMin, kFitMax, 8);
    TF1 *fbw2 = new TF1("fbw2", comp_bw2, kFitMin, kFitMax, 8);
    TF1 *fint = new TF1("fint", comp_int, kFitMin, kFitMax, 8);

    for(int i=0;i<8;i++){
        fbw1->SetParameter(i, f2->GetParameter(i));
        fbw2->SetParameter(i, f2->GetParameter(i));
        fint->SetParameter(i, f2->GetParameter(i));
    }

    f2->SetLineColor(kBlack);
    f2->SetLineWidth(2);

    fbw1->SetLineColor(kRed + 1);
    fbw1->SetLineWidth(2);

    fbw2->SetLineColor(kBlue + 1);
    fbw2->SetLineWidth(2);

    fint->SetLineColor(kGreen + 2);
    fint->SetLineStyle(2);
    fint->SetLineWidth(2);

    // total 就是 f2
    h->SetMinimum(kFit2YMin);
    h->SetMaximum(kFit2YMax);
    h->GetXaxis()->SetRangeUser(kFitMin, kFitMax);
    h->Draw("E");
    f2->Draw("SAME");
    fbw1->Draw("SAME");
    fbw2->Draw("SAME");
    fint->Draw("SAME");

    TLegend *leg = new TLegend(0.60, 0.62, 0.88, 0.88);
    leg->AddEntry(h, "Data", "lep");
    leg->AddEntry(f2, "Total fit", "l");
    leg->AddEntry(fbw1, "BW1", "l");
    leg->AddEntry(fbw2, "BW2", "l");
    leg->AddEntry(fint, "Interference", "l");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();
    draw_chi2_ndf(f2);

    c1->SaveAs(out_path("h_fit2.png").c_str());

    TF1 *f3 = new TF1("f3", model_interf3, kFitMin, kFitMax, 12);
    const char *f3ParNames[12] = {
      "scale","A1","m01","G1","A2","m02","G2","phi2","A3","m03","G3","phi3"
    };
    for(int i=0;i<12;i++){
      f3->SetParName(i, f3ParNames[i]);
      f3->SetParameter(i, kF3Init[i]);
    }
    f3->SetParLimits(2, kF3LimM01.min, kF3LimM01.max);
    f3->SetParLimits(3, kF3LimG1.min, kF3LimG1.max);
    f3->SetParLimits(5, kF3LimM02.min, kF3LimM02.max);
    f3->SetParLimits(6, kF3LimG2.min, kF3LimG2.max);
    f3->SetParLimits(7, kF3LimPhi2.min, kF3LimPhi2.max);
    f3->SetParLimits(9, kF3LimM03.min, kF3LimM03.max);
    f3->SetParLimits(10, kF3LimG3.min, kF3LimG3.max);
    f3->SetParLimits(11, kF3LimPhi3.min, kF3LimPhi3.max);
    h->Fit(f3, "R I S");

    TCanvas *c2 = new TCanvas("c2", "Mass Fit 3BW", 800, 600);
    h->SetMaximum(305);
    h->Draw("E");
    f3->SetLineColor(kMagenta + 1);
    f3->SetLineWidth(2);
    f3->Draw("SAME");

    TLegend *leg3 = new TLegend(0.60, 0.75, 0.88, 0.88);
    leg3->AddEntry(h, "Data", "lep");
    leg3->AddEntry(f3, "3-BW coherent fit", "l");
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->Draw();
    draw_chi2_ndf(f3);
    c2->SaveAs(out_path("h_fit3.png").c_str());

    TF1 *f3bw1 = new TF1("f3bw1", comp3_bw1, kFitMin, kFitMax, 12);
    TF1 *f3bw2 = new TF1("f3bw2", comp3_bw2, kFitMin, kFitMax, 12);
    TF1 *f3bw3 = new TF1("f3bw3", comp3_bw3, kFitMin, kFitMax, 12);
    TF1 *f3int12 = new TF1("f3int12", comp3_int12, kFitMin, kFitMax, 12);
    TF1 *f3int13 = new TF1("f3int13", comp3_int13, kFitMin, kFitMax, 12);
    TF1 *f3int23 = new TF1("f3int23", comp3_int23, kFitMin, kFitMax, 12);

    for(int i=0;i<12;i++){
      f3bw1->SetParameter(i, f3->GetParameter(i));
      f3bw2->SetParameter(i, f3->GetParameter(i));
      f3bw3->SetParameter(i, f3->GetParameter(i));
      f3int12->SetParameter(i, f3->GetParameter(i));
      f3int13->SetParameter(i, f3->GetParameter(i));
      f3int23->SetParameter(i, f3->GetParameter(i));
    }

    f3->SetLineColor(kBlack);
    f3->SetLineWidth(2);
    f3bw1->SetLineColor(kRed + 1);
    f3bw1->SetLineWidth(2);
    f3bw2->SetLineColor(kBlue + 1);
    f3bw2->SetLineWidth(2);
    f3bw3->SetLineColor(kMagenta + 2);
    f3bw3->SetLineWidth(2);

    f3int12->SetLineColor(kGreen + 2);
    f3int12->SetLineStyle(2);
    f3int12->SetLineWidth(2);
    f3int13->SetLineColor(kOrange + 7);
    f3int13->SetLineStyle(2);
    f3int13->SetLineWidth(2);
    f3int23->SetLineColor(kCyan + 2);
    f3int23->SetLineStyle(2);
    f3int23->SetLineWidth(2);

    TCanvas *c3 = new TCanvas("c3", "Mass Fit 3BW Components", 900, 700);
    //c3->SetLogy();
    TH1D *hlog = (TH1D*)h->Clone("hlog");
    hlog->SetMinimum(kFit3CompYMin);
    hlog->SetMaximum(kFit3CompYMax);
    hlog->GetXaxis()->SetRangeUser(kFitMin, kFitMax);
    hlog->Draw("E");
    f3->Draw("SAME");
    f3bw1->Draw("SAME");
    f3bw2->Draw("SAME");
    f3bw3->Draw("SAME");
    f3int12->Draw("SAME");
    f3int13->Draw("SAME");
    f3int23->Draw("SAME");

    TLegend *legc = new TLegend(0.56, 0.52, 0.88, 0.88);
    legc->AddEntry(hlog, "Data", "lep");
    legc->AddEntry(f3, "Total 3-BW", "l");
    legc->AddEntry(f3bw1, "BW1", "l");
    legc->AddEntry(f3bw2, "BW2", "l");
    legc->AddEntry(f3bw3, "BW3", "l");
    legc->AddEntry(f3int12, "Interf 12", "l");
    legc->AddEntry(f3int13, "Interf 13", "l");
    legc->AddEntry(f3int23, "Interf 23", "l");
    legc->SetBorderSize(0);
    legc->SetFillStyle(0);
    legc->Draw();
    draw_chi2_ndf(f3);

    c3->SaveAs(out_path("h_fit3_components.png").c_str());

    std::ofstream fout(out_path("fit_parameters.txt"));
    fout << std::setprecision(10);
    fout << "Fit range: [" << kFitMin << ", " << kFitMax << "]\n\n";
    write_fit_params(fout, "Single BW (f1)", f1);
    write_fit_params(fout, "Interference total (f2)", f2);
    write_fit_params(fout, "Three-BW coherent total (f3)", f3);
    fout.close();

    
}
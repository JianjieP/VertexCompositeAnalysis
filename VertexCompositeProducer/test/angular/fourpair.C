// Four-Pion Fixed Pairing Mass Plot (13, 24, 14, 23)
// This macro computes invariant masses for the fixed pairs:
//  - 1&3, 2&4, 1&4, 2&3
// and produces a single overlay plot with all four mass distributions.
// phi is the angle between the directions of the mother in lab frame and one of the daughter in the pair rest frame.


// Required includes
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <string>
#include <functional>
#include <filesystem>
#include <fstream>
#include "TLorentzVector.h"
#include <cstdlib>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TLine.h"
#include "TLatex.h"
#include <stdlib.h>
#include <time.h>

void fourpair() {
    // Output generation flags
    bool generate_png = true;      // Save PNG overlays
        bool generate_pdf = false;     // Save PDF overlays
    bool export_mass = false;   // Control saving mass PNGs
    bool export_mass_2d = false; // Control saving pair-mass 2D colz plots
    bool export_pT = false;     // Control saving pT-related PNGs
    bool export_pT_sum = true;    // Control saving sum-pT plot
        bool export_sum_filtered = false;  // Control saving filtered-sum plot
        bool export_sum_minpair = false;   // Control saving min-pair plot
    bool export_pair_angles = false;        // Control saving pair theta/phi plots
    bool export_pion_angles = false;   // Control saving per-pion theta/phi plots
    bool export_theta_output = false;   // Control saving theta-related outputs
    bool export_phi_output = true;     // Control saving phi-related outputs
        bool export_pion_angles_win = false; // Control saving mass-windowed per-pion angle plots
        bool export_pion_angles_ab = false; // Control saving A/B-rule filtered per-pion angles(A1 is pion look at A2)
        bool export_pion_angles_min = false; // Control saving min-pair filtered per-pion angles(|pTA|+|pTB|)
        bool export_pion_angles_rho = false; // Control saving rho-closest pair pion angles (best/counterpart)
        bool export_pion_angles_rho_rest = false; // Control saving rho-closest remaining-two pion angles
        bool export_pion_angles_rho_x = false; // Control saving rho-closest pion angles in pair-rest x-frame x aligns with total pT
        bool export_pion_angles_rho_x_rest = false; // Control saving rho-closest remaining-two pion angles in pair-rest x-frame
        bool export_pion_angles_rho_x_ptcut = true; // Control saving x-frame rho pion angles with pair pT < 0.1
        bool export_pion_angles_rho_x_ptcut_mid_high = true; // Control saving mid/high pT x-frame rho theta plots
    bool export_pairmass_phi = false; // Control saving pair-mass-category pi-phi and pair-pT plots
    bool export_pairmass_phi_angles = false; // Control saving pair-mass-category pi-phi plots
    bool export_pairmass_phi_combined = false; // Save pi-phi graphs (A, C, B+D low, B+D mid)
    bool export_pairmass_phi_a_combined = false; // Save pi_a_phi graphs (A, C, B+D low, B+D mid) x aligns with pair pT
    bool export_pairmass_pair_pt = false; // Control saving pair-mass-category pair pT plots
    bool export_pair_phi_xphi_selected = false; // Pair phi/xphi with 0.6<M<0.9 && pT<0.1, prefer pairs 13/24

    bool export_pion_pion_angles = tr; // Control saving angles between pairs
    
    bool randomization = true;   // Enable randomization of pair assignments

    int pionPhiFitOrder = 4;      // Cosine order for pi_phi sum fit (e.g., 1->cos(phi), 2->cos(2phi))
    bool only_even_cosine = false; // Use only even cosine terms in pion phi fit
    std::cout << "Generating single overlay plot for pairs 13, 24, 14, 23" << std::endl;

    // Histogram binning (adjust here)
    const int massBins = 300;   // number of bins for invariant mass
    const double massMin = 0.;  // GeV
    const double massMax = 4.5;  // GeV
    const int pTBins = 200;     // number of bins for pair pT
    const double pTMin = 0.;    // GeV
    const double pTMax = 2.;    // GeV
    const int totalPtBins = 200;      // number of bins for 4pi total pT plot
    const double totalPtPlotMin = 0.; // GeV
    const double totalPtPlotMax = 1.; // GeV
    const double totalPtMax = 0.1; // GeV, event-level cut on total four-pion pT
    bool applyTotalPtCut = true;
    if (const char* envCut = std::getenv("FOURPAIR_TOTALPT_CUT")) {
        if (std::string(envCut) == "0" || std::string(envCut) == "false" || std::string(envCut) == "FALSE") {
            applyTotalPtCut = false;
        }
    }
    const int thetaBins = 30;   // number of bins for pair-theta histograms
    const double thetaMin = 0.0;
    const double thetaMax = TMath::Pi();
    const int phiBins = 22;     // number of bins for pair-phi histograms
    const double phiMin = -TMath::Pi();
    const double phiMax = TMath::Pi();
    // Separate binning for per-pion angle histograms
    const int piThetaBins = 36; // per-pion theta bins
    const int piPhiBins   = 20
    
    
    
    ; // per-pion phi bins
    // Mass window for filtered pT plot
    
    const double massWinMin = 0.6; // GeV
    const double massWinMax = 0.9; // GeV
    
    // Read filtered input from EOS
    const char* input_path = "/eos/user/j/jianjie/fourpi/HIRun2025A_FourPi_0410/filtered_chained_ParticleTree.root";
    TFile *fin = TFile::Open(input_path);
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: cannot open input file: " << input_path << std::endl;
        return;
    }
    TTree *tin = (TTree*)fin->Get("ParticleTree");
    if (!tin) {
        std::cerr << "Error: TTree 'ParticleTree' not found in filtered_chained_ParticleTree.root" << std::endl;
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

    TFile *fout = new TFile("new_tree.root", "RECREATE");
    TTree *tout = new TTree("NewTree", "Tree with invariant masses");
    // Optimize ROOT I/O performance
    tout->SetAutoFlush(10000);
    tout->SetAutoSave(100000);
    fout->SetCompressionLevel(1);

    std::vector<float> inv_masses(4, 0);
    tout->Branch("inv_masses", &inv_masses);
    std::vector<float> pair_pT(4, 0);
    tout->Branch("pair_pT", &pair_pT);
    std::vector<float> pair_pT_scalar_sum(4, 0);
    tout->Branch("pair_pT_scalar_sum", &pair_pT_scalar_sum);
    std::vector<float> pair_theta(4, 0);
    tout->Branch("pair_theta", &pair_theta);
    std::vector<float> pair_phi(4, 0);
    tout->Branch("pair_phi", &pair_phi);
    TH1F *h_total_pT = new TH1F("h_total_pT", "Total p_{T};p_{T} (GeV);Counts", totalPtBins, totalPtPlotMin, totalPtPlotMax);
    h_total_pT->SetStats(0);
    // Per-pion angle branches
    std::vector<float> pi_theta(4, 0);
    std::vector<float> pi_phi(4, 0);
    std::vector<float> pi_theta_x(4, 0);
    std::vector<float> pi_phi_x(4, 0);
    std::vector<float> pi_theta_a(4, 0);
    std::vector<float> pi_phi_a(4, 0);
    tout->Branch("pi_theta", &pi_theta);
    tout->Branch("pi_phi", &pi_phi);
    tout->Branch("pi_theta_x", &pi_theta_x);
    tout->Branch("pi_phi_x", &pi_phi_x);
    tout->Branch("pi_theta_a", &pi_theta_a);
    tout->Branch("pi_phi_a", &pi_phi_a);
    float total_pT = 0.f;
    tout->Branch("total_pT", &total_pT);
    float total_mass = 0.f;
    tout->Branch("total_mass", &total_mass);
    // (Optional) Add minimal four-pion info if needed later
    // float four_pi_mass = 0; float four_pi_pT = 0;

    auto get_pion = [&](int i) {
        float pT = (*cand_pT)[i];
        float phi = (*cand_phi)[i];
        float px = pT * cos(phi);
        float py = pT * sin(phi);
        float p = (*cand_p)[i];
        float pz_sign = ((*cand_eta)[i] >= 0) ? 1.0 : -1.0;
        float pz = pz_sign * sqrt(p*p - pT*pT);
        float mass = (*cand_mass)[i];
        float E = sqrt(mass*mass + p*p);
        TLorentzVector v(px, py, pz, E);
        return v;
    };

    Long64_t nentries = tin->GetEntries();
    std::cout << "Processing " << nentries << " entries..." << std::endl;
    
    // Progress reporting
    Long64_t progress_interval = nentries / 20;  // Report every 5%
    if (progress_interval == 0) progress_interval = 1000;
    
    for (Long64_t i=0; i<nentries; ++i) {
        if (i % progress_interval == 0) {
            std::cout << "Processed " << i << " / " << nentries << " entries (" 
                      << (100.0 * i / nentries) << "%)" << std::endl;
        }
        
        tin->GetEntry(i);

        // Protect against malformed events with insufficient candidates
        if (!cand_p || cand_p->size() != 5 || !cand_pT || cand_pT->size() != 5 ||
            !cand_phi || cand_phi->size() != 5 || !cand_mass || cand_mass->size() != 5 ||
            !cand_eta || cand_eta->size() != 5) {
            continue; // skip this entry to avoid out-of-range access
        }

        TLorentzVector p1 = get_pion(1);
        TLorentzVector p2 = get_pion(2);
        TLorentzVector p3 = get_pion(3);
        TLorentzVector p4 = get_pion(4);

        if(randomization) {
            // Seed random number generator
            static bool seeded = false;
            if (!seeded) {
                srand(static_cast<unsigned int>(time(0)));
                seeded = true;
            }
            int r1 = rand();
            int r2 = rand();
            if(r1 % 2 == 0) {
                std::swap(p1, p2);
            }
            if(r2 % 2 == 0) {
                std::swap(p3, p4);
            }
        }
        
        // Event-level total pT cut
        TLorentzVector total = p1 + p2 + p3 + p4;
        h_total_pT->Fill(total.Pt());
        if (applyTotalPtCut && total.Pt() >= totalPtMax) continue;

        // Pre-calculate all pair combinations once
        TLorentzVector pair13 = p1 + p3;
        TLorentzVector pair24 = p2 + p4;
        TLorentzVector pair14 = p1 + p4;
        TLorentzVector pair23 = p2 + p3;


        // Calculate masses once
        float m13 = pair13.M(); // 1&3
        float m24 = pair24.M(); // 2&4
        float m14 = pair14.M(); // 1&4
        float m23 = pair23.M(); // 2&3

        inv_masses[0] = m13;
        inv_masses[1] = m24;
        inv_masses[2] = m14;
        inv_masses[3] = m23;

        pair_pT[0] = pair13.Pt();
        pair_pT[1] = pair24.Pt();
        pair_pT[2] = pair14.Pt();
        pair_pT[3] = pair23.Pt();
        pair_pT_scalar_sum[0] = pair13.Pt();// pair pT for scalar sum
        pair_pT_scalar_sum[1] = pair24.Pt();
        pair_pT_scalar_sum[2] = pair14.Pt();
        pair_pT_scalar_sum[3] = pair23.Pt();

        // Rotate so x aligns with total pT, then boost pir pairs to 4pi rest frame to get angles, pir pix pia to pair restframe, 
        // pir x aligns with total pT, pix pia x aligns with pair pT
        TLorentzVector p1r = p1, p2r = p2, p3r = p3, p4r = p4;
        TLorentzVector p1x = p1, p2x = p2, p3x = p3, p4x = p4;
        TLorentzVector p1a = p1, p2a = p2, p3a = p3, p4a = p4;
        TLorentzVector total_rot = total;
        double rotPhi = -total.Phi();
        p1r.RotateZ(rotPhi); p2r.RotateZ(rotPhi); p3r.RotateZ(rotPhi); p4r.RotateZ(rotPhi);
        total_rot.RotateZ(rotPhi);
        TLorentzVector pair13r = p1r + p3r;
        TLorentzVector pair24r = p2r + p4r;
        TLorentzVector pair14r = p1r + p4r;
        TLorentzVector pair23r = p2r + p3r;
        TVector3 boostVec13 = -pair13r.BoostVector();
        TVector3 boostVec24 = -pair24r.BoostVector();
        TVector3 boostVec14 = -pair14r.BoostVector();
        TVector3 boostVec23 = -pair23r.BoostVector();
        TVector3 boostVec = -total_rot.BoostVector();
        pair13r.Boost(boostVec);
        pair24r.Boost(boostVec);
        pair14r.Boost(boostVec);
        pair23r.Boost(boostVec);
        pair_theta[0] = pair13r.Theta();
        pair_theta[1] = pair24r.Theta();
        pair_theta[2] = pair14r.Theta();
        pair_theta[3] = pair23r.Theta();
        pair_phi[0] = pair13r.Phi();
        pair_phi[1] = pair24r.Phi();
        pair_phi[2] = pair14r.Phi();
        pair_phi[3] = pair23r.Phi();
        total_pT = total.Pt();
        total_mass = total.M();

        p1r=p3r;p2r=p4r;
        p1r.Boost(boostVec13);
        p2r.Boost(boostVec24);
        p4r.Boost(boostVec14);
        p3r.Boost(boostVec23);
        pi_theta[0] = p1r.Theta();
        pi_theta[1] = p2r.Theta();
        pi_theta[2] = p4r.Theta();
        pi_theta[3] = p3r.Theta();
        pi_phi[0] = p1r.Phi();
        pi_phi[1] = p2r.Phi();
        pi_phi[2] = p4r.Phi();
        pi_phi[3] = p3r.Phi();

        TLorentzVector pair13x = p1x + p3x;
        TLorentzVector pair24x = p2x + p4x;
        TLorentzVector pair14x = p1x + p4x;
        TLorentzVector pair23x = p2x + p3x;
        double rotPhix1 = -pair13x.Phi();
        double rotPhix2 = -pair24x.Phi();
        double rotPhix3 = -pair14x.Phi();
        double rotPhix4 = -pair23x.Phi();
        double rotPhix = -total.Phi();
        p1x=p3x;p2x=p4x;// keep only positive pions
        p1x.RotateZ(rotPhix1); p2x.RotateZ(rotPhix2); p4x.RotateZ(rotPhix3); p3x.RotateZ(rotPhix4);
        pair13x.RotateZ(rotPhix1); pair24x.RotateZ(rotPhix2); pair14x.RotateZ(rotPhix3); pair23x.RotateZ(rotPhix4);
        TVector3 boostVec13x = -pair13x.BoostVector();
        TVector3 boostVec24x = -pair24x.BoostVector();
        TVector3 boostVec14x = -pair14x.BoostVector();
        TVector3 boostVec23x = -pair23x.BoostVector();
        p1x.Boost(boostVec13x);
        p2x.Boost(boostVec24x);
        p4x.Boost(boostVec14x);
        p3x.Boost(boostVec23x);
        pi_theta_x[0] = p1x.Theta();
        pi_theta_x[1] = p2x.Theta();
        pi_theta_x[2] = p4x.Theta();
        pi_theta_x[3] = p3x.Theta();
        pi_phi_x[0] = p1x.Phi();
        pi_phi_x[1] = p2x.Phi();
        pi_phi_x[2] = p4x.Phi();
        pi_phi_x[3] = p3x.Phi();

        TLorentzVector pair13a = p1a + p3a;
        TLorentzVector pair24a = p2a + p4a;
        TLorentzVector pair14a = p1a + p4a;
        TLorentzVector pair23a = p2a + p3a;
        double rotPhia1 = -pair13a.Phi();
        double rotPhia2 = -pair24a.Phi();
        double rotPhia3 = -pair14a.Phi();
        double rotPhia4 = -pair23a.Phi();
        p1a=p3a;p2a=p4a;// keep only positive pions
        p1a.RotateZ(rotPhia1); p2a.RotateZ(rotPhia2); p4a.RotateZ(rotPhia3); p3a.RotateZ(rotPhia4);
        pair13a.RotateZ(rotPhia1); pair24a.RotateZ(rotPhia2); pair14a.RotateZ(rotPhia3); pair23a.RotateZ(rotPhia4);
        TVector3 boostVec13a = -pair13a.BoostVector();
        TVector3 boostVec24a = -pair24a.BoostVector();
        TVector3 boostVec14a = -pair14a.BoostVector();
        TVector3 boostVec23a = -pair23a.BoostVector();
        p1a.Boost(boostVec13a);
        p2a.Boost(boostVec24a);
        p4a.Boost(boostVec14a);
        p3a.Boost(boostVec23a);
        pi_theta_a[0] = p1a.Theta();
        pi_theta_a[1] = p2a.Theta();
        pi_theta_a[2] = p4a.Theta();
        pi_theta_a[3] = p3a.Theta();
        pi_phi_a[0] = p1a.Phi();
        pi_phi_a[1] = p2a.Phi();
        pi_phi_a[2] = p4a.Phi();
        pi_phi_a[3] = p3a.Phi();
        
        tout->Fill();
    }

    if (export_pT && (generate_png || generate_pdf)) {
        std::filesystem::path total_pT_dir = "fourpair_png";
        if (generate_png) {
            try {
                std::filesystem::create_directories(total_pT_dir);
            } catch (const std::exception& ex) {
                std::cerr << "Warning: failed to create total pT output directory: " << ex.what() << std::endl;
            }
        }

        TCanvas *c_total = new TCanvas("c_total_pT", "Total pT: all events", 900, 700);
        h_total_pT->SetLineColor(kBlack);
        h_total_pT->SetLineWidth(2);
        h_total_pT->SetMarkerStyle(20);
        h_total_pT->SetMarkerSize(0.8);
        h_total_pT->SetMarkerColor(kBlack);
        float ymax_total = 1.2 * h_total_pT->GetMaximum();
        h_total_pT->SetMaximum(ymax_total);
        h_total_pT->Draw("E1P");
        auto leg_total = new TLegend(0.64, 0.72, 0.84, 0.88);
        leg_total->SetTextSize(0.03);
        leg_total->SetFillStyle(0);
        leg_total->SetBorderSize(0);
        leg_total->AddEntry(h_total_pT, "All events", "lep");
        leg_total->Draw();
        if (generate_png) c_total->SaveAs((total_pT_dir / "total_pT_all_events.png").string().c_str());
        if (generate_pdf) c_total->SaveAs((total_pT_dir / "total_pT_all_events.pdf").string().c_str());
    }

    std::cout << "Finished processing entries. Writing output..." << std::endl;
    fout->Write();
    fout->Close();
    fin->Close();
    std::cout << "Output file written successfully." << std::endl;

    // Plot inv_masses as four overlaid histograms and save image
    std::cout << "Creating mass histograms and overlay plot..." << std::endl;
    TFile *fplot = TFile::Open("new_tree.root", "UPDATE");
    if (!fplot || fplot->IsZombie()) {
        std::cerr << "Error: cannot open new_tree.root for plotting" << std::endl;
        return;
    }
    TTree *tplot = (TTree*)fplot->Get("NewTree");
    if (!tplot) {
        std::cerr << "Error: TTree 'NewTree' not found in new_tree.root" << std::endl;
        fplot->Close();
        return;
    }
    // Prepare histograms for the four fixed pairs
    TH1F *h_m13 = new TH1F("h_m13", "Invariant Mass;M (GeV);Counts", massBins, massMin, massMax);
    TH1F *h_m24 = new TH1F("h_m24", "Invariant Mass;M (GeV);Counts", massBins, massMin, massMax);
    TH1F *h_m14 = new TH1F("h_m14", "Invariant Mass;M (GeV);Counts", massBins, massMin, massMax);
    TH1F *h_m23 = new TH1F("h_m23", "Invariant Mass;M (GeV);Counts", massBins, massMin, massMax);
    TH1F *h_total_mass = new TH1F("h_total_mass", "Four-pion invariant mass;M (GeV);Counts", massBins, massMin, massMax);
    h_m13->SetStats(0);
    h_m24->SetStats(0);
    h_m14->SetStats(0);
    h_m23->SetStats(0);
    h_total_mass->SetStats(0);

    TH2F *h_pairMass_A1B1_vs_A2B2 = new TH2F(
        "h_pairMass_A1B1_vs_A2B2",
        "Pair1 mass vs Pair2 mass;M_{A1/B1} (GeV);M_{A2/B2} (GeV)",
        massBins, massMin, massMax,
        massBins, massMin, massMax
    );
    h_pairMass_A1B1_vs_A2B2->SetStats(0);

    const int nPairMassCats = 4;
    std::array<const char*, nPairMassCats> pairMassCatTags = {"A", "B", "C", "D"};
    std::array<const char*, nPairMassCats> pairMassCatTitles = {
        "A: 0.6<M_{pair1}<0.9 && 0.6<M_{pair2}<0.9",
        "B: M_{pair1}<0.6 && 0.6<M_{pair2}<0.9",
        "C: M_{pair1}<0.6 && M_{pair2}<0.6",
        "D: 0.6<M_{pair1}<0.9 && M_{pair2}<0.6"
    };
    std::array<TH1F*, nPairMassCats> h_pi_phi_pair1_cat;
    std::array<TH1F*, nPairMassCats> h_pi_phi_pair2_cat;
    std::array<TH1F*, nPairMassCats> h_pi_phi_pair1_cat_ptcut;
    std::array<TH1F*, nPairMassCats> h_pi_phi_pair2_cat_ptcut;
    TH1F *h_pi_phi_pair_cat_A_combined = nullptr;
    TH1F *h_pi_phi_pair_cat_C_combined = nullptr;
    TH1F *h_pi_phi_pair_cat_B1D2_combined = nullptr;
    TH1F *h_pi_phi_pair_cat_B2D1_combined = nullptr;
    TH1F *h_pi_a_phi_pair_cat_A_combined = nullptr;
    TH1F *h_pi_a_phi_pair_cat_C_combined = nullptr;
    TH1F *h_pi_a_phi_pair_cat_B1D2_combined = nullptr;
    TH1F *h_pi_a_phi_pair_cat_B2D1_combined = nullptr;
    std::array<TH1F*, nPairMassCats> h_pT_pair1_cat;
    std::array<TH1F*, nPairMassCats> h_pT_pair2_cat;
    for (int i = 0; i < nPairMassCats; ++i) {
        h_pi_phi_pair1_cat[i] = new TH1F(
            Form("h_pi_phi_pair1_cat_%s", pairMassCatTags[i]),
            ";#Delta#varphi;Counts",
            piPhiBins, phiMin, phiMax
        );
        h_pi_phi_pair2_cat[i] = new TH1F(
            Form("h_pi_phi_pair2_cat_%s", pairMassCatTags[i]),
            ";#Delta#varphi;Counts",
            piPhiBins, phiMin, phiMax
        );
        h_pi_phi_pair1_cat_ptcut[i] = new TH1F(
            Form("h_pi_phi_pair1_cat_%s_ptcut", pairMassCatTags[i]),
            ";#Delta#varphi;Counts",
            piPhiBins, phiMin, phiMax
        );
        h_pi_phi_pair2_cat_ptcut[i] = new TH1F(
            Form("h_pi_phi_pair2_cat_%s_ptcut", pairMassCatTags[i]),
            ";#Delta#varphi;Counts",
            piPhiBins, phiMin, phiMax
        );
        h_pT_pair1_cat[i] = new TH1F(
            Form("h_pT_pair1_cat_%s", pairMassCatTags[i]),
            Form("Pair p_{T} %s;p_{T} (GeV);Counts", pairMassCatTitles[i]),
            pTBins, pTMin, pTMax
        );
        h_pT_pair2_cat[i] = new TH1F(
            Form("h_pT_pair2_cat_%s", pairMassCatTags[i]),
            Form("Pair p_{T} (pair2) %s;p_{T} (GeV);Counts", pairMassCatTitles[i]),
            pTBins, pTMin, pTMax
        );
        h_pi_phi_pair1_cat[i]->SetStats(0);
        h_pi_phi_pair2_cat[i]->SetStats(0);
        h_pi_phi_pair1_cat_ptcut[i]->SetStats(0);
        h_pi_phi_pair2_cat_ptcut[i]->SetStats(0);
        h_pT_pair1_cat[i]->SetStats(0);
        h_pT_pair2_cat[i]->SetStats(0);
    }
    if (export_pairmass_phi_combined) {
        h_pi_phi_pair_cat_A_combined = new TH1F(
            "h_pi_phi_pair_cat_A_combined",
            ";#Delta#varphi;Counts",
            piPhiBins, phiMin, phiMax
        );
        h_pi_phi_pair_cat_C_combined = new TH1F(
            "h_pi_phi_pair_cat_C_combined",
            ";#Delta#varphi;Counts",
            piPhiBins, phiMin, phiMax
        );
        h_pi_phi_pair_cat_B1D2_combined = new TH1F(
            "h_pi_phi_pair_cat_B1D2_combined",
            ";#Delta#varphi;Counts",
            piPhiBins, phiMin, phiMax
        );
        h_pi_phi_pair_cat_B2D1_combined = new TH1F(
            "h_pi_phi_pair_cat_B2D1_combined",
            ";#Delta#varphi;Counts",
            piPhiBins, phiMin, phiMax
        );
        h_pi_phi_pair_cat_A_combined->SetStats(0);
        h_pi_phi_pair_cat_C_combined->SetStats(0);
        h_pi_phi_pair_cat_B1D2_combined->SetStats(0);
        h_pi_phi_pair_cat_B2D1_combined->SetStats(0);
    }
    if (export_pairmass_phi_a_combined) {
        h_pi_a_phi_pair_cat_A_combined = new TH1F(
            "h_pi_a_phi_pair_cat_A_combined",
            ";#Delta#varphi;Counts",
            piPhiBins, phiMin, phiMax
        );
        h_pi_a_phi_pair_cat_C_combined = new TH1F(
            "h_pi_a_phi_pair_cat_C_combined",
            ";#Delta#varphi;Counts",
            piPhiBins, phiMin, phiMax
        );
        h_pi_a_phi_pair_cat_B1D2_combined = new TH1F(
            "h_pi_a_phi_pair_cat_B1D2_combined",
            ";#Delta#varphi;Counts",
            piPhiBins, phiMin, phiMax
        );
        h_pi_a_phi_pair_cat_B2D1_combined = new TH1F(
            "h_pi_a_phi_pair_cat_B2D1_combined",
            ";#Delta#varphi;Counts",
            piPhiBins, phiMin, phiMax
        );
        h_pi_a_phi_pair_cat_A_combined->SetStats(0);
        h_pi_a_phi_pair_cat_C_combined->SetStats(0);
        h_pi_a_phi_pair_cat_B1D2_combined->SetStats(0);
        h_pi_a_phi_pair_cat_B2D1_combined->SetStats(0);
    }

    TH1F *h_pT13 = new TH1F("h_pT13", "Pair p_{T};p_{T} (GeV);Counts", pTBins, pTMin, pTMax);
    TH1F *h_pT24 = new TH1F("h_pT24", "Pair p_{T};p_{T} (GeV);Counts", pTBins, pTMin, pTMax);
    TH1F *h_pT14 = new TH1F("h_pT14", "Pair p_{T};p_{T} (GeV);Counts", pTBins, pTMin, pTMax);
    TH1F *h_pT23 = new TH1F("h_pT23", "Pair p_{T};p_{T} (GeV);Counts", pTBins, pTMin, pTMax);
    h_pT13->SetStats(0);
    h_pT24->SetStats(0);
    h_pT14->SetStats(0);
    h_pT23->SetStats(0);
    TH1F *h_pT_sum = new TH1F("h_pT_sum", "Sum of pair p_{T} hist counts; p_{T} (GeV);Counts", pTBins, pTMin, pTMax);
    h_pT_sum->SetStats(0);
    TH1F *h_pT_sum_totalpt = new TH1F("h_pT_sum_totalpt", "pair p_{T} (total p_{T}<0.1); p_{T} (GeV);Counts", pTBins, pTMin, pTMax);
    TH1F *h_pT_sum_totalpt_mid = new TH1F("h_pT_sum_totalpt_mid", "pair p_{T} (0.1<=total p_{T}<0.2); p_{T} (GeV);Counts", pTBins, pTMin, pTMax);
    TH1F *h_pT_sum_totalpt_high = new TH1F("h_pT_sum_totalpt_high", "pair p_{T} (total p_{T}>=0.2); p_{T} (GeV);Counts", pTBins, pTMin, pTMax);
    h_pT_sum_totalpt->SetStats(0);
    h_pT_sum_totalpt_mid->SetStats(0);
    h_pT_sum_totalpt_high->SetStats(0);
    TH1F *h_pT_sum_filtered = new TH1F("h_pT_sum_filtered", "Filtered sum of pair p_{T} (per A/B rule); p_{T} (GeV);Counts", pTBins, pTMin, pTMax);
    h_pT_sum_filtered->SetStats(0);
    TH1F *h_pT_sum_minpair = new TH1F("h_pT_sum_minpair", "Min-pair selection (|p_{T1}|+|p_{T2}|); p_{T} (GeV);Counts", pTBins, pTMin, pTMax);
    h_pT_sum_minpair->SetStats(0);
    TH1F *h_pT_sum_win_totalpt = new TH1F("h_pT_sum_win_totalpt", "Sum of pair p_{T} (0.6<M<0.9, total p_{T}<0.1); p_{T} (GeV);Counts", pTBins, pTMin, pTMax);
    TH1F *h_pT_sum_win_totalpt_mid = new TH1F("h_pT_sum_win_totalpt_mid", "Sum of pair p_{T} (0.6<M<0.9, 0.1<=total p_{T}<0.2); p_{T} (GeV);Counts", pTBins, pTMin, pTMax);
    TH1F *h_pT_sum_win_totalpt_high = new TH1F("h_pT_sum_win_totalpt_high", "Sum of pair p_{T} (0.6<M<0.9, total p_{T}>=0.2); p_{T} (GeV);Counts", pTBins, pTMin, pTMax);
    h_pT_sum_win_totalpt->SetStats(0);
    h_pT_sum_win_totalpt_mid->SetStats(0);
    h_pT_sum_win_totalpt_high->SetStats(0);

    // Angle histograms (0<theta<pi, -pi<phi<pi)
    TH1F *h_theta_A1 = new TH1F("h_theta_A1", "Theta Pair A;#theta;Counts", thetaBins, thetaMin, thetaMax);
    TH1F *h_theta_A2 = new TH1F("h_theta_A2", "Theta Pair A2;#theta;Counts", thetaBins, thetaMin, thetaMax);
    TH1F *h_theta_B1 = new TH1F("h_theta_B1", "Theta Pair B;#theta;Counts", thetaBins, thetaMin, thetaMax);
    TH1F *h_theta_B2 = new TH1F("h_theta_B2", "Theta Pair B2;#theta;Counts", thetaBins, thetaMin, thetaMax);
    TH1F *h_phi_A1 = new TH1F("h_phi_A1", ";#Delta#varphi;Counts", phiBins, phiMin, phiMax);
    TH1F *h_phi_A2 = new TH1F("h_phi_A2", ";#Delta#varphi;Counts", phiBins, phiMin, phiMax);
    TH1F *h_phi_B1 = new TH1F("h_phi_B1", ";#Delta#varphi;Counts", phiBins, phiMin, phiMax);
    TH1F *h_phi_B2 = new TH1F("h_phi_B2", ";#Delta#varphi;Counts", phiBins, phiMin, phiMax);
    TH1F *h_pair_phi_sel = new TH1F("h_pair_phi_sel", ";#Delta#varphi;Counts", phiBins, phiMin, phiMax);
    TH1F *h_pair_xphi_sel = new TH1F("h_pair_xphi_sel", ";#Delta#varphi;Counts", phiBins, phiMin, phiMax);
    h_theta_A1->SetStats(0);
    h_theta_A2->SetStats(0);
    h_theta_B1->SetStats(0);
    h_theta_B2->SetStats(0);
    h_phi_A1->SetStats(0);
    h_phi_A2->SetStats(0);
    h_phi_B1->SetStats(0);
    h_phi_B2->SetStats(0);
    h_pair_phi_sel->SetStats(0);
    h_pair_xphi_sel->SetStats(0);

    // Per-pion angle histograms (filled only if pi_theta/pi_phi branches exist)
    TH1F *h_pi_theta_1 = new TH1F("h_pi_theta_1", "Pion #theta;#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_theta_2 = new TH1F("h_pi_theta_2", "Pion #theta;#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_theta_3 = new TH1F("h_pi_theta_3", "Pion #theta;#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_theta_4 = new TH1F("h_pi_theta_4", "Pion #theta;#theta;Counts", piThetaBins, thetaMin, thetaMax);
    h_pi_theta_1->SetStats(0);
    h_pi_theta_2->SetStats(0);
    h_pi_theta_3->SetStats(0);
    h_pi_theta_4->SetStats(0);
    TH1F *h_pi_theta_1_win = new TH1F("h_pi_theta_1_win", "Pion #theta (0.6<M<0.9);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_theta_2_win = new TH1F("h_pi_theta_2_win", "Pion #theta (0.6<M<0.9);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_theta_3_win = new TH1F("h_pi_theta_3_win", "Pion #theta (0.6<M<0.9);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_theta_4_win = new TH1F("h_pi_theta_4_win", "Pion #theta (0.6<M<0.9);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    h_pi_theta_1_win->SetStats(0);
    h_pi_theta_2_win->SetStats(0);
    h_pi_theta_3_win->SetStats(0);
    h_pi_theta_4_win->SetStats(0);
    TH1F *h_pi_phi_1 = new TH1F("h_pi_phi_1", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_2 = new TH1F("h_pi_phi_2", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_3 = new TH1F("h_pi_phi_3", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_4 = new TH1F("h_pi_phi_4", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    h_pi_phi_1->SetStats(0);
    h_pi_phi_2->SetStats(0);
    h_pi_phi_3->SetStats(0);
    h_pi_phi_4->SetStats(0);

    // Per-pion phi histograms with mass window (0.6<M<0.9) applied per corresponding pair
    TH1F *h_pi_phi_1_win = new TH1F("h_pi_phi_1_win", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_2_win = new TH1F("h_pi_phi_2_win", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_3_win = new TH1F("h_pi_phi_3_win", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_4_win = new TH1F("h_pi_phi_4_win", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    h_pi_phi_1_win->SetStats(0);
    h_pi_phi_2_win->SetStats(0);
    h_pi_phi_3_win->SetStats(0);
    h_pi_phi_4_win->SetStats(0);

    // Summed filtered pion angle histograms (A/B rule and min-pair selection)
    TH1F *h_pi_theta_ab_sum = new TH1F("h_pi_theta_ab_sum", "Sum #pi^{+} #theta (A/B rule);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_phi_ab_sum   = new TH1F("h_pi_phi_ab_sum",   ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    h_pi_theta_ab_sum->SetStats(0);
    h_pi_phi_ab_sum->SetStats(0);
    TH1F *h_pi_theta_min_sum = new TH1F("h_pi_theta_min_sum", "Sum #pi^{+} #theta (min-pair);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_phi_min_sum   = new TH1F("h_pi_phi_min_sum",   ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    h_pi_theta_min_sum->SetStats(0);
    h_pi_phi_min_sum->SetStats(0);

    // Rho-closest selection histograms (best, counterpart, and remaining two pions)
    TH1F *h_pi_theta_rho_best = new TH1F("h_pi_theta_rho_best", "#pi^{+} #theta (closest to #rho);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_theta_rho_other = new TH1F("h_pi_theta_rho_other", "#pi^{+} #theta (counterpart of #rho-closest);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_theta_rho_rest = new TH1F("h_pi_theta_rho_rest", "#pi^{+} #theta (remaining two pions);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_phi_rho_best = new TH1F("h_pi_phi_rho_best", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_rho_other = new TH1F("h_pi_phi_rho_other", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_rho_rest = new TH1F("h_pi_phi_rho_rest", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_theta_x_rho_best = new TH1F("h_pi_theta_x_rho_best", "#pi^{+} #theta (closest to #rho);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_theta_x_rho_other = new TH1F("h_pi_theta_x_rho_other", "#pi^{+} #theta (counterpart of #rho-closest);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_theta_x_rho_rest = new TH1F("h_pi_theta_x_rho_rest", "#pi^{+} #theta (remaining two pions);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_phi_x_rho_best = new TH1F("h_pi_phi_x_rho_best", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_x_rho_other = new TH1F("h_pi_phi_x_rho_other", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_x_rho_rest = new TH1F("h_pi_phi_x_rho_rest", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_theta_x_rho_best_pt = new TH1F("h_pi_theta_x_rho_best_pt", "#pi^{+} #theta (closest #rho, p_{T}<0.1);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_theta_x_rho_best_pt_mid = new TH1F("h_pi_theta_x_rho_best_pt_mid", "#pi^{+} #theta (closest #rho, 0.1<p_{T}<0.2);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_theta_x_rho_best_pt_high = new TH1F("h_pi_theta_x_rho_best_pt_high", "#pi^{+} #theta (closest #rho, p_{T}>0.2);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_theta_x_rho_other_pt = new TH1F("h_pi_theta_x_rho_other_pt", "#pi^{+} #theta (counterpart, p_{T}<0.1);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_theta_x_rho_other_pt_mid = new TH1F("h_pi_theta_x_rho_other_pt_mid", "#pi^{+} #theta (counterpart, 0.1<p_{T}<0.2);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_theta_x_rho_other_pt_high = new TH1F("h_pi_theta_x_rho_other_pt_high", "#pi^{+} #theta (counterpart, p_{T}>0.2);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_theta_x_rho_rest_pt = new TH1F("h_pi_theta_x_rho_rest_pt", "#pi^{+} #theta (remaining, p_{T}<0.1);#theta;Counts", piThetaBins, thetaMin, thetaMax);
    TH1F *h_pi_phi_x_rho_best_pt = new TH1F("h_pi_phi_x_rho_best_pt", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_x_rho_other_pt = new TH1F("h_pi_phi_x_rho_other_pt", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_x_rho_rest_pt = new TH1F("h_pi_phi_x_rho_rest_pt", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_x_rho_best_pt_mid = new TH1F("h_pi_phi_x_rho_best_pt_mid", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_x_rho_best_pt_high = new TH1F("h_pi_phi_x_rho_best_pt_high", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_x_rho_other_pt_mid = new TH1F("h_pi_phi_x_rho_other_pt_mid", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    TH1F *h_pi_phi_x_rho_other_pt_high = new TH1F("h_pi_phi_x_rho_other_pt_high", ";#Delta#varphi;Counts", piPhiBins, phiMin, phiMax);
    h_pi_theta_rho_best->SetStats(0);
    h_pi_theta_rho_other->SetStats(0);
    h_pi_theta_rho_rest->SetStats(0);
    h_pi_phi_rho_best->SetStats(0);
    h_pi_phi_rho_other->SetStats(0);
    h_pi_phi_rho_rest->SetStats(0);
    h_pi_theta_x_rho_best->SetStats(0);
    h_pi_theta_x_rho_other->SetStats(0);
    h_pi_theta_x_rho_rest->SetStats(0);
    h_pi_phi_x_rho_best->SetStats(0);
    h_pi_phi_x_rho_other->SetStats(0);
    h_pi_phi_x_rho_rest->SetStats(0);
    h_pi_theta_x_rho_best_pt->SetStats(0);
    h_pi_theta_x_rho_best_pt_mid->SetStats(0);
    h_pi_theta_x_rho_best_pt_high->SetStats(0);
    h_pi_theta_x_rho_other_pt->SetStats(0);
    h_pi_theta_x_rho_other_pt_mid->SetStats(0);
    h_pi_theta_x_rho_other_pt_high->SetStats(0);
    h_pi_theta_x_rho_rest_pt->SetStats(0);
    h_pi_phi_x_rho_best_pt->SetStats(0);
    h_pi_phi_x_rho_other_pt->SetStats(0);
    h_pi_phi_x_rho_rest_pt->SetStats(0);
    h_pi_phi_x_rho_best_pt_mid->SetStats(0);
    h_pi_phi_x_rho_best_pt_high->SetStats(0);
    h_pi_phi_x_rho_other_pt_mid->SetStats(0);
    h_pi_phi_x_rho_other_pt_high->SetStats(0);

    TH1F *h_pT13_win = new TH1F("h_pT13_win", "Pair p_{T} (0.6<M<0.9);p_{T} (GeV);Counts", pTBins, pTMin, pTMax);
    TH1F *h_pT24_win = new TH1F("h_pT24_win", "Pair p_{T} (0.6<M<0.9);p_{T} (GeV);Counts", pTBins, pTMin, pTMax);
    TH1F *h_pT14_win = new TH1F("h_pT14_win", "Pair p_{T} (0.6<M<0.9);p_{T} (GeV);Counts", pTBins, pTMin, pTMax);
    TH1F *h_pT23_win = new TH1F("h_pT23_win", "Pair p_{T} (0.6<M<0.9);p_{T} (GeV);Counts", pTBins, pTMin, pTMax);
    h_pT13_win->SetStats(0);
    h_pT24_win->SetStats(0);
    h_pT14_win->SetStats(0);
    h_pT23_win->SetStats(0);
    TH1F *h_pT_sum_win = new TH1F("h_pT_sum_win", "Sum of pair p_{T} hist counts (0.6<M<0.9); p_{T} (GeV);Counts", pTBins, pTMin, pTMax);
    h_pT_sum_win->SetStats(0);

    // Total pT histogram from events before the total pT cut
    TH1F *h_total_pT_win = new TH1F("h_total_pT_win", "Total p_{T} (0.6<M<0.9);p_{T} (GeV);Counts", totalPtBins, totalPtPlotMin, totalPtPlotMax);
    h_total_pT_win->SetStats(0);

    // Read inv_masses branch
    std::vector<float> *inv_masses_read = nullptr;
    std::vector<float> *pair_pT_read = nullptr;
    std::vector<float> *pair_pT_scalar_sum_read = nullptr;
    std::vector<float> *pair_theta_read = nullptr;
    std::vector<float> *pair_phi_read = nullptr;
    std::vector<float> *pi_theta_read = nullptr;
    std::vector<float> *pi_phi_read = nullptr;
    std::vector<float> *pi_theta_x_read = nullptr;
    std::vector<float> *pi_phi_x_read = nullptr;
    std::vector<float> *pi_theta_a_read = nullptr;
    std::vector<float> *pi_phi_a_read = nullptr;
    float total_pT_read = 0.f;
    float total_mass_read = 0.f;
    static TRandom3 rng_sel(12345);
    tplot->SetBranchAddress("inv_masses", &inv_masses_read);
    tplot->SetBranchAddress("pair_pT", &pair_pT_read);
    if (tplot->GetBranch("pair_pT_scalar_sum")) {
        tplot->SetBranchAddress("pair_pT_scalar_sum", &pair_pT_scalar_sum_read);
    }
    if (tplot->GetBranch("pair_theta")) {
        tplot->SetBranchAddress("pair_theta", &pair_theta_read);
    }
    if (tplot->GetBranch("pair_phi")) {
        tplot->SetBranchAddress("pair_phi", &pair_phi_read);
    }
    if (tplot->GetBranch("pi_theta")) {
        tplot->SetBranchAddress("pi_theta", &pi_theta_read);
    }
    if (tplot->GetBranch("pi_phi")) {
        tplot->SetBranchAddress("pi_phi", &pi_phi_read);
    }
    if (tplot->GetBranch("pi_theta_x")) {
        tplot->SetBranchAddress("pi_theta_x", &pi_theta_x_read);
    }
    if (tplot->GetBranch("pi_phi_x")) {
        tplot->SetBranchAddress("pi_phi_x", &pi_phi_x_read);
    }
    if (tplot->GetBranch("pi_theta_a")) {
        tplot->SetBranchAddress("pi_theta_a", &pi_theta_a_read);
    }
    if (tplot->GetBranch("pi_phi_a")) {
        tplot->SetBranchAddress("pi_phi_a", &pi_phi_a_read);
    }
    if (tplot->GetBranch("total_pT")) {
        tplot->SetBranchAddress("total_pT", &total_pT_read);
    }
    if (tplot->GetBranch("total_mass")) {
        tplot->SetBranchAddress("total_mass", &total_mass_read);
    }



    Long64_t nplot = tplot->GetEntries();
    for (Long64_t i = 0; i < nplot; ++i) {
        tplot->GetEntry(i);
        bool passTotalPtCutPairMass = (!tplot->GetBranch("total_pT")) || (total_pT_read < 0.1f);
        if (inv_masses_read && inv_masses_read->size() == 4 && passTotalPtCutPairMass) {
        h_m13->Fill((*inv_masses_read)[0]);
        h_m24->Fill((*inv_masses_read)[1]);
        h_m14->Fill((*inv_masses_read)[2]);
        h_m23->Fill((*inv_masses_read)[3]);
        }
        if (tplot->GetBranch("total_mass") && passTotalPtCutPairMass) {
            h_total_mass->Fill(total_mass_read);
        }
        
        double pair1Mass_A1 = (*inv_masses_read)[0];
        double pair2Mass_A2 = (*inv_masses_read)[1];
        double pair1Mass_B1 = (*inv_masses_read)[2];
        double pair2Mass_B2 = (*inv_masses_read)[3];
        
        if (passTotalPtCutPairMass) {
            h_pairMass_A1B1_vs_A2B2->Fill(pair1Mass_A1, pair2Mass_A2);
            h_pairMass_A1B1_vs_A2B2->Fill(pair1Mass_B1, pair2Mass_B2);
        }
        bool passTotalPtCutPairPhi = (!tplot->GetBranch("total_pT")) || (total_pT_read < 0.1f);
        if (export_pairmass_phi && passTotalPtCutPairPhi && pair_pT_read && pair_pT_read->size() == 4 &&
            pi_phi_x_read && pi_phi_x_read->size() >= 4) {
            auto catIndex = [&](double m1, double m2) {
                bool m1win = (m1 > massWinMin && m1 < massWinMax);
                bool m2win = (m2 > massWinMin && m2 < massWinMax);
                if (m1win && m2win) return 0;               // A
                if (m1 < massWinMin && m2win) return 1;     // B
                if (m1 < massWinMin && m2 < massWinMin) return 2; // C
                if (m1win && m2 < massWinMin) return 3;     // D
                return -1;
            };
            int catA = catIndex(pair1Mass_A1, pair2Mass_A2);
            if (catA >= 0) {
                if (export_pairmass_phi_angles) {
                    h_pi_phi_pair1_cat[catA]->Fill((*pi_phi_x_read)[0]);
                    h_pi_phi_pair2_cat[catA]->Fill((*pi_phi_x_read)[1]);
                    if (pair1Mass_A1 > massWinMin && pair1Mass_A1 < massWinMax && (*pair_pT_read)[0] < 0.1f) {
                        h_pi_phi_pair1_cat_ptcut[catA]->Fill((*pi_phi_x_read)[0]);
                    }
                    if (pair2Mass_A2 > massWinMin && pair2Mass_A2 < massWinMax && (*pair_pT_read)[1] < 0.1f) {
                        h_pi_phi_pair2_cat_ptcut[catA]->Fill((*pi_phi_x_read)[1]);
                    }
                }
                if (export_pairmass_phi_combined && catA == 0) {
                    //if ((*pair_pT_read)[0] < 0.1f) 
                        h_pi_phi_pair_cat_A_combined->Fill((*pi_phi_read)[0]);
                    //if ((*pair_pT_read)[1] < 0.1f) 
                        h_pi_phi_pair_cat_A_combined->Fill((*pi_phi_read)[1]);
                }
                if (export_pairmass_phi_combined && catA == 2) {
                    h_pi_phi_pair_cat_C_combined->Fill((*pi_phi_read)[0]);
                    h_pi_phi_pair_cat_C_combined->Fill((*pi_phi_read)[1]);
                }
                if (export_pairmass_phi_combined && catA == 1) {
                    h_pi_phi_pair_cat_B1D2_combined->Fill((*pi_phi_read)[0]);
                    //if ((*pair_pT_read)[1] < 0.1f) 
                        h_pi_phi_pair_cat_B2D1_combined->Fill((*pi_phi_read)[1]);
                }
                if (export_pairmass_phi_combined && catA == 3) {
                    //if ((*pair_pT_read)[0] < 0.1f) 
                        h_pi_phi_pair_cat_B2D1_combined->Fill((*pi_phi_read)[0]);
                    h_pi_phi_pair_cat_B1D2_combined->Fill((*pi_phi_read)[1]);
                }
                if (export_pairmass_phi_a_combined && pi_phi_a_read && pi_phi_a_read->size() >= 4) {
                    if (catA == 0) {
                        if ((*pair_pT_read)[0] < 0.1f) h_pi_a_phi_pair_cat_A_combined->Fill((*pi_phi_a_read)[0]);
                        if ((*pair_pT_read)[1] < 0.1f) h_pi_a_phi_pair_cat_A_combined->Fill((*pi_phi_a_read)[1]);
                    }
                    if (catA == 2) {
                        h_pi_a_phi_pair_cat_C_combined->Fill((*pi_phi_a_read)[0]);
                        h_pi_a_phi_pair_cat_C_combined->Fill((*pi_phi_a_read)[1]);
                    }
                    if (catA == 1) {
                        h_pi_a_phi_pair_cat_B1D2_combined->Fill((*pi_phi_a_read)[0]);
                        if ((*pair_pT_read)[1] < 0.1f) h_pi_a_phi_pair_cat_B2D1_combined->Fill((*pi_phi_a_read)[1]);
                    }
                    if (catA == 3) {
                        if ((*pair_pT_read)[0] < 0.1f) h_pi_a_phi_pair_cat_B2D1_combined->Fill((*pi_phi_a_read)[0]);
                        h_pi_a_phi_pair_cat_B1D2_combined->Fill((*pi_phi_a_read)[1]);
                    }
                }
                if (export_pairmass_pair_pt) {
                    h_pT_pair1_cat[catA]->Fill((*pair_pT_read)[0]);
                    h_pT_pair2_cat[catA]->Fill((*pair_pT_read)[1]);
                }
            }
            int catB = catIndex(pair1Mass_B1, pair2Mass_B2);
            if (catB >= 0) {
                if (export_pairmass_phi_angles) {
                    h_pi_phi_pair1_cat[catB]->Fill((*pi_phi_read)[2]);
                    h_pi_phi_pair2_cat[catB]->Fill((*pi_phi_read)[3]);
                    if (pair1Mass_B1 > massWinMin && pair1Mass_B1 < massWinMax && (*pair_pT_read)[2] < 0.1f) {
                        h_pi_phi_pair1_cat_ptcut[catB]->Fill((*pi_phi_read)[2]);
                    }
                    if (pair2Mass_B2 > massWinMin && pair2Mass_B2 < massWinMax && (*pair_pT_read)[3] < 0.1f) {
                        h_pi_phi_pair2_cat_ptcut[catB]->Fill((*pi_phi_read)[3]);
                    }
                }
                if (export_pairmass_phi_combined && catB == 1) {
                    h_pi_phi_pair_cat_B1D2_combined->Fill((*pi_phi_read)[2]);
                    //if ((*pair_pT_read)[3] < 0.1f) 
                        h_pi_phi_pair_cat_B2D1_combined->Fill((*pi_phi_read)[3]);
                }
                if (export_pairmass_phi_combined && catB == 3) {
                    //if ((*pair_pT_read)[2] < 0.1f) 
                        h_pi_phi_pair_cat_B2D1_combined->Fill((*pi_phi_read)[2]);
                    h_pi_phi_pair_cat_B1D2_combined->Fill((*pi_phi_read)[3]);
                }
                if (export_pairmass_phi_combined && catB == 2) {
                    h_pi_phi_pair_cat_C_combined->Fill((*pi_phi_read)[2]);
                    h_pi_phi_pair_cat_C_combined->Fill((*pi_phi_read)[3]);
                }
                if (export_pairmass_phi_combined && catB == 0) {
                    //if ((*pair_pT_read)[2] < 0.1f) 
                        h_pi_phi_pair_cat_A_combined->Fill((*pi_phi_read)[2]);
                    //if ((*pair_pT_read)[3] < 0.1f) 
                        h_pi_phi_pair_cat_A_combined->Fill((*pi_phi_read)[3]);
                }
                if (export_pairmass_phi_a_combined && pi_phi_a_read && pi_phi_a_read->size() >= 4) {
                    if (catB == 1) {
                        h_pi_a_phi_pair_cat_B1D2_combined->Fill((*pi_phi_a_read)[2]);
                        if ((*pair_pT_read)[3] < 0.1f) h_pi_a_phi_pair_cat_B2D1_combined->Fill((*pi_phi_a_read)[3]);
                    }
                    if (catB == 3) {
                        if ((*pair_pT_read)[2] < 0.1f) h_pi_a_phi_pair_cat_B2D1_combined->Fill((*pi_phi_a_read)[2]);
                        h_pi_a_phi_pair_cat_B1D2_combined->Fill((*pi_phi_a_read)[3]);
                    }
                    if (catB == 2) {
                        h_pi_a_phi_pair_cat_C_combined->Fill((*pi_phi_a_read)[2]);
                        h_pi_a_phi_pair_cat_C_combined->Fill((*pi_phi_a_read)[3]);
                    }
                    if (catB == 0) {
                        if ((*pair_pT_read)[2] < 0.1f) h_pi_a_phi_pair_cat_A_combined->Fill((*pi_phi_a_read)[2]);
                        if ((*pair_pT_read)[3] < 0.1f) h_pi_a_phi_pair_cat_A_combined->Fill((*pi_phi_a_read)[3]);
                    }
                }
                if (export_pairmass_pair_pt) {
                    h_pT_pair1_cat[catB]->Fill((*pair_pT_read)[2]);
                    h_pT_pair2_cat[catB]->Fill((*pair_pT_read)[3]);
                }
            }
        }
        if (pair_pT_read && pair_pT_read->size() == 4) {
            h_pT13->Fill((*pair_pT_read)[0]);
            h_pT24->Fill((*pair_pT_read)[1]);
            h_pT14->Fill((*pair_pT_read)[2]);
            h_pT23->Fill((*pair_pT_read)[3]);
            bool hasTotalPt = tplot->GetBranch("total_pT");
            bool passTotalPtCutForSum = (!hasTotalPt) || (total_pT_read < 0.1f);
            if (passTotalPtCutForSum) {
                h_pT_sum_totalpt->Fill((*pair_pT_read)[0]);
                h_pT_sum_totalpt->Fill((*pair_pT_read)[1]);
                h_pT_sum_totalpt->Fill((*pair_pT_read)[2]);
                h_pT_sum_totalpt->Fill((*pair_pT_read)[3]);
                if ((*inv_masses_read)[0] > massWinMin && (*inv_masses_read)[0] < massWinMax) h_pT_sum_win_totalpt->Fill((*pair_pT_read)[0]);
                if ((*inv_masses_read)[1] > massWinMin && (*inv_masses_read)[1] < massWinMax) h_pT_sum_win_totalpt->Fill((*pair_pT_read)[1]);
                if ((*inv_masses_read)[2] > massWinMin && (*inv_masses_read)[2] < massWinMax) h_pT_sum_win_totalpt->Fill((*pair_pT_read)[2]);
                if ((*inv_masses_read)[3] > massWinMin && (*inv_masses_read)[3] < massWinMax) h_pT_sum_win_totalpt->Fill((*pair_pT_read)[3]);
            } else if (hasTotalPt && total_pT_read >= 0.1f && total_pT_read < 0.2f) {
                h_pT_sum_totalpt_mid->Fill((*pair_pT_read)[0]);
                h_pT_sum_totalpt_mid->Fill((*pair_pT_read)[1]);
                h_pT_sum_totalpt_mid->Fill((*pair_pT_read)[2]);
                h_pT_sum_totalpt_mid->Fill((*pair_pT_read)[3]);
                if ((*inv_masses_read)[0] > massWinMin && (*inv_masses_read)[0] < massWinMax) h_pT_sum_win_totalpt_mid->Fill((*pair_pT_read)[0]);
                if ((*inv_masses_read)[1] > massWinMin && (*inv_masses_read)[1] < massWinMax) h_pT_sum_win_totalpt_mid->Fill((*pair_pT_read)[1]);
                if ((*inv_masses_read)[2] > massWinMin && (*inv_masses_read)[2] < massWinMax) h_pT_sum_win_totalpt_mid->Fill((*pair_pT_read)[2]);
                if ((*inv_masses_read)[3] > massWinMin && (*inv_masses_read)[3] < massWinMax) h_pT_sum_win_totalpt_mid->Fill((*pair_pT_read)[3]);
            } else if (hasTotalPt && total_pT_read >= 0.2f) {
                h_pT_sum_totalpt_high->Fill((*pair_pT_read)[0]);
                h_pT_sum_totalpt_high->Fill((*pair_pT_read)[1]);
                h_pT_sum_totalpt_high->Fill((*pair_pT_read)[2]);
                h_pT_sum_totalpt_high->Fill((*pair_pT_read)[3]);
                if ((*inv_masses_read)[0] > massWinMin && (*inv_masses_read)[0] < massWinMax) h_pT_sum_win_totalpt_high->Fill((*pair_pT_read)[0]);
                if ((*inv_masses_read)[1] > massWinMin && (*inv_masses_read)[1] < massWinMax) h_pT_sum_win_totalpt_high->Fill((*pair_pT_read)[1]);
                if ((*inv_masses_read)[2] > massWinMin && (*inv_masses_read)[2] < massWinMax) h_pT_sum_win_totalpt_high->Fill((*pair_pT_read)[2]);
                if ((*inv_masses_read)[3] > massWinMin && (*inv_masses_read)[3] < massWinMax) h_pT_sum_win_totalpt_high->Fill((*pair_pT_read)[3]);
            }
            // Mass-filtered pT (0.6 < M < 0.9)
            if ((*inv_masses_read)[0] > massWinMin && (*inv_masses_read)[0] < massWinMax) h_pT13_win->Fill((*pair_pT_read)[0]);
            if ((*inv_masses_read)[1] > massWinMin && (*inv_masses_read)[1] < massWinMax) h_pT24_win->Fill((*pair_pT_read)[1]);
            if ((*inv_masses_read)[2] > massWinMin && (*inv_masses_read)[2] < massWinMax) h_pT14_win->Fill((*pair_pT_read)[2]);
            if ((*inv_masses_read)[3] > massWinMin && (*inv_masses_read)[3] < massWinMax) h_pT23_win->Fill((*pair_pT_read)[3]);

            // A/B selection rule for filtered sum
            auto in_cond = [&](int idx){
                return ((*inv_masses_read)[idx] > massWinMin && (*inv_masses_read)[idx] < massWinMax && (*pair_pT_read)[idx] < 0.15);
            };
            bool a1 = in_cond(0);
            bool a2 = in_cond(1);
            bool b1 = in_cond(2);
            bool b2 = in_cond(3);
            // Decide which to keep for A group
            int keepA = -1;
            if (a1 && !a2) keepA = 1;       // keep A2 only
            else if (a2 && !a1) keepA = 0;  // keep A1 only
            else if (a1 && a2) keepA = (rng_sel.Integer(2) == 0 ? 0 : 1); // both: random one
            // Decide which to keep for B group
            int keepB = -1;
            if (b1 && !b2) keepB = 3;       // keep B2 only
            else if (b2 && !b1) keepB = 2;  // keep B1 only
            else if (b1 && b2) keepB = (rng_sel.Integer(2) == 0 ? 2 : 3); // both: random one

            if (keepA >= 0) h_pT_sum_filtered->Fill((*pair_pT_read)[keepA]);
            if (keepB >= 0) h_pT_sum_filtered->Fill((*pair_pT_read)[keepB]);

            if (export_pion_angles_ab) {
                if (keepA >= 0) {
                    h_pi_theta_ab_sum->Fill((*pi_theta_x_read)[keepA]);
                    h_pi_phi_ab_sum->Fill((*pi_phi_x_read)[keepA]);
                }
                if (keepB >= 0) {
                    h_pi_theta_ab_sum->Fill((*pi_theta_x_read)[keepB]);
                    h_pi_phi_ab_sum->Fill((*pi_phi_x_read)[keepB]);
                }
            }

            // Min-pair method: compare total |pT| sums of (A1+A2) vs (B1+B2); keep the smaller group
            if (pair_pT_scalar_sum_read && pair_pT_scalar_sum_read->size() == 4) {
                double sumA = (*pair_pT_scalar_sum_read)[0] + (*pair_pT_scalar_sum_read)[1];
                double sumB = (*pair_pT_scalar_sum_read)[2] + (*pair_pT_scalar_sum_read)[3];
                bool keepAgroup = (sumA <= sumB); // keep A on tie
                if (keepAgroup) {
                    h_pT_sum_minpair->Fill((*pair_pT_read)[0]);
                    h_pT_sum_minpair->Fill((*pair_pT_read)[1]);
                    if (export_pion_angles_min) {
                        if (pi_theta_x_read && pi_theta_x_read->size() >= 4 && pi_phi_x_read && pi_phi_x_read->size() >= 4) {
                            h_pi_theta_min_sum->Fill((*pi_theta_x_read)[0]);
                            h_pi_theta_min_sum->Fill((*pi_theta_x_read)[1]);
                            h_pi_phi_min_sum->Fill((*pi_phi_x_read)[0]);
                            h_pi_phi_min_sum->Fill((*pi_phi_x_read)[1]);
                        }
                    }
                } else {
                    h_pT_sum_minpair->Fill((*pair_pT_read)[2]);
                    h_pT_sum_minpair->Fill((*pair_pT_read)[3]);
                    if (export_pion_angles_min) {
                        if (pi_theta_x_read && pi_theta_x_read->size() >= 4 && pi_phi_x_read && pi_phi_x_read->size() >= 4) {
                            h_pi_theta_min_sum->Fill((*pi_theta_x_read)[2]);
                            h_pi_theta_min_sum->Fill((*pi_theta_x_read)[3]);
                            h_pi_phi_min_sum->Fill((*pi_phi_x_read)[2]);
                            h_pi_phi_min_sum->Fill((*pi_phi_x_read)[3]);
                        }
                    }
                }
            }

            // Rho-closest selection using pion angles (pi1 for pairs 13/14, pi2 for pairs 24/23)
            if (export_pion_angles_rho && inv_masses_read && inv_masses_read->size() >= 4) {
                constexpr double mrho = 0.770; // GeV
                double d0 = std::abs((*inv_masses_read)[0] - mrho);
                double d1 = std::abs((*inv_masses_read)[1] - mrho);
                double d2 = std::abs((*inv_masses_read)[2] - mrho);
                double d3 = std::abs((*inv_masses_read)[3] - mrho);
                double deltas[4] = {d0, d1, d2, d3};
                int bestIdx = 0;
                for (int k = 1; k < 4; ++k) {
                    if (deltas[k] < deltas[bestIdx]) bestIdx = k;
                }
                bool used[4] = {false, false, false, false};
                used[bestIdx] = true;
                auto pionIndexForPair = [](int pairIdx) {
                    if (pairIdx == 0 || pairIdx == 3) return 2; // pairs 13,23 use pion 3
                    return 3; // pairs 24,14 use pion 4
                };
                auto counterpartPair = [](int pairIdx) {
                    if (pairIdx == 0) return 1; // A1 -> A2
                    if (pairIdx == 1) return 0; // A2 -> A1
                    if (pairIdx == 2) return 3; // B1 -> B2
                    return 2; // B2 -> B1
                };
                //int pionBest = pionIndexForPair(bestIdx);
                int otherPair = counterpartPair(bestIdx);
                used[otherPair] = true;
                //int pionOther = pionIndexForPair(otherPair);
                double bestPairPt = (pair_pT_read && pair_pT_read->size() == 4) ? (*pair_pT_read)[bestIdx] : -1.0;
                bool passTotalPtCut = (!tplot->GetBranch("total_pT")) || (total_pT_read < 0.1f);
                bool passPtCut = export_pion_angles_rho_x_ptcut && passTotalPtCut && bestPairPt >= 0.0 && bestPairPt < 0.1;
                bool passPtCutMid = export_pion_angles_rho_x_ptcut && export_pion_angles_rho_x_ptcut_mid_high && passTotalPtCut && bestPairPt >= 0.1 && bestPairPt < 0.2;
                bool passPtCutHigh = export_pion_angles_rho_x_ptcut && export_pion_angles_rho_x_ptcut_mid_high && passTotalPtCut && bestPairPt >= 0.2;
                if (pi_theta_read && pi_theta_read->size() >= 4) {
                    h_pi_theta_rho_best->Fill((*pi_theta_read)[bestIdx]);
                    h_pi_theta_rho_other->Fill((*pi_theta_read)[otherPair]);
                    if (export_pion_angles_rho_x && pi_theta_x_read && pi_theta_x_read->size() >= 4) {
                        h_pi_theta_x_rho_best->Fill((*pi_theta_x_read)[bestIdx]);
                        h_pi_theta_x_rho_other->Fill((*pi_theta_x_read)[otherPair]);
                        if (passPtCut) {
                            h_pi_theta_x_rho_best_pt->Fill((*pi_theta_x_read)[bestIdx]);
                            h_pi_theta_x_rho_other_pt->Fill((*pi_theta_x_read)[otherPair]);
                        }
                        if (passPtCutMid) {
                            h_pi_theta_x_rho_best_pt_mid->Fill((*pi_theta_x_read)[bestIdx]);
                            h_pi_theta_x_rho_other_pt_mid->Fill((*pi_theta_x_read)[otherPair]);
                        }
                        if (passPtCutHigh) {
                            h_pi_theta_x_rho_best_pt_high->Fill((*pi_theta_x_read)[bestIdx]);
                            h_pi_theta_x_rho_other_pt_high->Fill((*pi_theta_x_read)[otherPair]);
                        }
                    }
                    bool fill_theta_rest = export_pion_angles_rho_rest || (export_pion_angles_rho_x && export_pion_angles_rho_x_rest);
                    if (fill_theta_rest) {
                        // Remaining two pions (those not in best or counterpart)
                        for (int idx = 0; idx < 4; ++idx) {
                            if (used[idx]) continue;
                            //int pionIdx = pionIndexForPair(idx);
                            if (export_pion_angles_rho_rest) h_pi_theta_rho_rest->Fill((*pi_theta_read)[idx]);
                            if (export_pion_angles_rho_x && export_pion_angles_rho_x_rest && pi_theta_x_read && pi_theta_x_read->size() >= 4) {
                                h_pi_theta_x_rho_rest->Fill((*pi_theta_x_read)[idx]);
                                if (passPtCut) h_pi_theta_x_rho_rest_pt->Fill((*pi_theta_x_read)[idx]);
                            }
                        }
                    }
                }
                if (pi_phi_read && pi_phi_read->size() >= 4) {
                    h_pi_phi_rho_best->Fill((*pi_phi_read)[bestIdx]);
                    h_pi_phi_rho_other->Fill((*pi_phi_read)[otherPair]);
                    if (export_pion_angles_rho_x && pi_phi_x_read && pi_phi_x_read->size() >= 4) {
                        h_pi_phi_x_rho_best->Fill((*pi_phi_x_read)[bestIdx]);
                        h_pi_phi_x_rho_other->Fill((*pi_phi_x_read)[otherPair]);
                        if (passPtCut) {
                            h_pi_phi_x_rho_best_pt->Fill((*pi_phi_x_read)[bestIdx]);
                            h_pi_phi_x_rho_other_pt->Fill((*pi_phi_x_read)[otherPair]);
                        }
                        if (passPtCutMid) {
                            h_pi_phi_x_rho_best_pt_mid->Fill((*pi_phi_x_read)[bestIdx]);
                            h_pi_phi_x_rho_other_pt_mid->Fill((*pi_phi_x_read)[otherPair]);
                        }
                        if (passPtCutHigh) {
                            h_pi_phi_x_rho_best_pt_high->Fill((*pi_phi_x_read)[bestIdx]);
                            h_pi_phi_x_rho_other_pt_high->Fill((*pi_phi_x_read)[otherPair]);
                        }
                    }
                    bool fill_phi_rest = export_pion_angles_rho_rest || (export_pion_angles_rho_x && export_pion_angles_rho_x_rest);
                    if (fill_phi_rest) {
                        for (int idx = 0; idx < 4; ++idx) {
                            if (used[idx]) continue;
                            //int pionIdx = pionIndexForPair(idx);
                            if (export_pion_angles_rho_rest) h_pi_phi_rho_rest->Fill((*pi_phi_read)[idx]);
                            if (export_pion_angles_rho_x && export_pion_angles_rho_x_rest && pi_phi_x_read && pi_phi_x_read->size() >= 4) {
                                h_pi_phi_x_rho_rest->Fill((*pi_phi_x_read)[idx]);
                                if (passPtCut) h_pi_phi_x_rho_rest_pt->Fill((*pi_phi_x_read)[idx]);
                            }
                        }
                    }
                }
            }
        }
        if (pair_theta_read && pair_phi_read && pair_theta_read->size() == 4 && pair_phi_read->size() == 4) {
            h_theta_A1->Fill((*pair_theta_read)[0]);
            h_theta_A2->Fill((*pair_theta_read)[1]);
            h_theta_B1->Fill((*pair_theta_read)[2]);
            h_theta_B2->Fill((*pair_theta_read)[3]);
            h_phi_A1->Fill((*pair_phi_read)[0]);
            h_phi_A2->Fill((*pair_phi_read)[1]);
            h_phi_B1->Fill((*pair_phi_read)[2]);
            h_phi_B2->Fill((*pair_phi_read)[3]);
        }
        if (export_pair_phi_xphi_selected && pair_phi_read && pair_phi_read->size() == 4 &&
            inv_masses_read && inv_masses_read->size() == 4 && pair_pT_read && pair_pT_read->size() == 4) {
            bool passTotalPtCutPairPhi = (!tplot->GetBranch("total_pT")) || (total_pT_read < 0.1f);
            if (passTotalPtCutPairPhi) {
                bool sel13 = ((*inv_masses_read)[0] > massWinMin && (*inv_masses_read)[0] < massWinMax && (*pair_pT_read)[0] < 0.1f);
                bool sel24 = ((*inv_masses_read)[1] > massWinMin && (*inv_masses_read)[1] < massWinMax && (*pair_pT_read)[1] < 0.1f);
                bool anyPrimary = sel13 || sel24;
                auto fillSel = [&](int idx) {
                    h_pair_phi_sel->Fill((*pair_phi_read)[idx]);
                    if (pi_phi_x_read && pi_phi_x_read->size() == 4) {
                        h_pair_xphi_sel->Fill((*pi_phi_x_read)[idx]);
                    }
                };
                if (anyPrimary) {
                    if (sel13) fillSel(0);
                    if (sel24) fillSel(1);
                } else {
                    bool sel14 = ((*inv_masses_read)[2] > massWinMin && (*inv_masses_read)[2] < massWinMax && (*pair_pT_read)[2] < 0.1f);
                    bool sel23 = ((*inv_masses_read)[3] > massWinMin && (*inv_masses_read)[3] < massWinMax && (*pair_pT_read)[3] < 0.1f);
                    if (sel14) fillSel(2);
                    if (sel23) fillSel(3);
                }
            }
        }
        if (pi_theta_read && pi_theta_read->size() == 4) {
            h_pi_theta_1->Fill((*pi_theta_read)[0]);
            h_pi_theta_2->Fill((*pi_theta_read)[1]);
            h_pi_theta_3->Fill((*pi_theta_read)[2]);
            h_pi_theta_4->Fill((*pi_theta_read)[3]);
            if (inv_masses_read && inv_masses_read->size() >= 4) {
                if ((*inv_masses_read)[0] > massWinMin && (*inv_masses_read)[0] < massWinMax) h_pi_theta_1_win->Fill((*pi_theta_read)[0]);
                if ((*inv_masses_read)[1] > massWinMin && (*inv_masses_read)[1] < massWinMax) h_pi_theta_2_win->Fill((*pi_theta_read)[1]);
                if ((*inv_masses_read)[2] > massWinMin && (*inv_masses_read)[2] < massWinMax) h_pi_theta_3_win->Fill((*pi_theta_read)[2]);
                if ((*inv_masses_read)[3] > massWinMin && (*inv_masses_read)[3] < massWinMax) h_pi_theta_4_win->Fill((*pi_theta_read)[3]);
            }
        }
        if (pi_phi_read && pi_phi_read->size() == 4) {
            // Unfiltered per-pion phi
            h_pi_phi_1->Fill((*pi_phi_read)[0]);
            h_pi_phi_2->Fill((*pi_phi_read)[1]);
            h_pi_phi_3->Fill((*pi_phi_read)[2]);
            h_pi_phi_4->Fill((*pi_phi_read)[3]);
            // Mass-window filtered per-pion phi using the corresponding pair mass
            if (inv_masses_read && inv_masses_read->size() >= 4) {
                if ((*inv_masses_read)[0] > massWinMin && (*inv_masses_read)[0] < massWinMax) h_pi_phi_1_win->Fill((*pi_phi_read)[0]);
                if ((*inv_masses_read)[1] > massWinMin && (*inv_masses_read)[1] < massWinMax) h_pi_phi_2_win->Fill((*pi_phi_read)[1]);
                if ((*inv_masses_read)[2] > massWinMin && (*inv_masses_read)[2] < massWinMax) h_pi_phi_3_win->Fill((*pi_phi_read)[2]);
                if ((*inv_masses_read)[3] > massWinMin && (*inv_masses_read)[3] < massWinMax) h_pi_phi_4_win->Fill((*pi_phi_read)[3]);
            }
        }
        if (tplot->GetBranch("total_pT")) {
            // Use "any pair in mass window" as the filter condition for total pT
            bool any_in_window = (((*inv_masses_read)[0] > massWinMin && (*inv_masses_read)[0] < massWinMax) &&
                                 ((*inv_masses_read)[1] > massWinMin && (*inv_masses_read)[1] < massWinMax)) ||
                                 (((*inv_masses_read)[2] > massWinMin && (*inv_masses_read)[2] < massWinMax) &&
                                 ((*inv_masses_read)[3] > massWinMin && (*inv_masses_read)[3] < massWinMax));
            if (any_in_window) h_total_pT_win->Fill(total_pT_read);
        }
    }

    // Build summed histograms from individual pair histograms
    h_pT_sum->Add(h_pT13); h_pT_sum->Add(h_pT24); h_pT_sum->Add(h_pT14); h_pT_sum->Add(h_pT23);
    h_pT_sum_win->Add(h_pT13_win); h_pT_sum_win->Add(h_pT24_win); h_pT_sum_win->Add(h_pT14_win); h_pT_sum_win->Add(h_pT23_win);
    // h_pT_sum_filtered and h_pT_sum_minpair already filled during loop

    // Save PNG outputs into a dedicated local folder.
    std::filesystem::path original_workdir = std::filesystem::current_path();
    bool switched_to_png_dir = false;
    if (generate_png) {
        try {
            std::filesystem::path png_output_dir = original_workdir / "fourpair_png";
            std::filesystem::create_directories(png_output_dir);
            std::filesystem::current_path(png_output_dir);
            switched_to_png_dir = true;
            std::cout << "PNG output directory: " << png_output_dir << std::endl;
        } catch (const std::exception& ex) {
            std::cerr << "Warning: failed to switch to PNG output directory: " << ex.what() << std::endl;
        }
    }

    // Create single overlay canvas
    if (export_mass && (generate_png || generate_pdf)) {
        TCanvas *c = new TCanvas("c_pairs_overlay", "Pairs 13,24,14,23 Mass Overlay", 900, 700);
        h_m13->SetLineColor(kBlue);
        h_m24->SetLineColor(kRed);
        h_m14->SetLineColor(kGreen+2);
        h_m23->SetLineColor(kMagenta);
        h_m13->SetLineWidth(2);
        h_m24->SetLineWidth(2);
        h_m14->SetLineWidth(2);
        h_m23->SetLineWidth(2);

        float ymax = 1.2 * std::max(std::max(h_m13->GetMaximum(), h_m24->GetMaximum()), std::max(h_m14->GetMaximum(), h_m23->GetMaximum()));
        h_m13->SetMaximum(ymax);

        h_m13->Draw();
        h_m24->Draw("SAME");
        h_m14->Draw("SAME");
        h_m23->Draw("SAME");

        auto leg = new TLegend(0.64, 0.72, 0.84, 0.88);
        leg->SetTextSize(0.03);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->AddEntry(h_m13, "A1", "l");
        leg->AddEntry(h_m24, "A2", "l");
        leg->AddEntry(h_m14, "B1", "l");
        leg->AddEntry(h_m23, "B2", "l");
        leg->Draw();

        if (generate_png) c->SaveAs("pairs_13_24_14_23_overlay.png");
        if (generate_pdf) c->SaveAs("pairs_13_24_14_23_overlay.pdf");

        TCanvas *c_total_mass = new TCanvas("c_pairs_total_mass", "Four-pion invariant mass", 900, 700);
        c_total_mass->SetLogy();
        h_total_mass->SetLineColor(kBlue+2);
        h_total_mass->SetLineWidth(2);
        h_total_mass->SetMarkerStyle(20);
        h_total_mass->SetMarkerSize(0.9);
        h_total_mass->SetMarkerColor(kBlue+2);
        h_total_mass->Draw("E1P");
        if (generate_png) c_total_mass->SaveAs("pairs_total_mass.png");
        if (generate_pdf) c_total_mass->SaveAs("pairs_total_mass.pdf");
    }

    if (export_mass_2d && (generate_png || generate_pdf)) {
        TCanvas *c_pairMass = new TCanvas("c_pairMass_A1B1_vs_A2B2", "Pair1 vs Pair2 mass (A1+B1 vs A2+B2)", 900, 750);
        h_pairMass_A1B1_vs_A2B2->Draw("COLZ");
        TLine *line_x_rho = new TLine(0.77, massMin, 0.77, massMax);
        TLine *line_y_rho = new TLine(massMin, 0.77, massMax, 0.77);
        line_x_rho->SetLineStyle(2);
        line_x_rho->SetLineWidth(2);
        line_y_rho->SetLineStyle(2);
        line_y_rho->SetLineWidth(2);
        line_x_rho->Draw("SAME");
        line_y_rho->Draw("SAME");
        if (generate_png) c_pairMass->SaveAs("pairMass_A1B1_vs_A2B2_colz.png");
        if (generate_pdf) c_pairMass->SaveAs("pairMass_A1B1_vs_A2B2_colz.pdf");

        TCanvas *c_pairMass_cats = new TCanvas("c_pairMass_A1B1_vs_A2B2_cats", "Pair1 vs Pair2 mass categories", 900, 750);
        h_pairMass_A1B1_vs_A2B2->Draw("COLZ");
        TLine *line_x_06 = new TLine(massWinMin, massMin, massWinMin, massMax);
        TLine *line_y_06 = new TLine(massMin, massWinMin, massMax, massWinMin);
        TLine *line_x_09 = new TLine(massWinMax, massMin, massWinMax, massMax);
        TLine *line_y_09 = new TLine(massMin, massWinMax, massMax, massWinMax);
        line_x_06->SetLineStyle(2);
        line_x_06->SetLineWidth(2);
        line_y_06->SetLineStyle(2);
        line_y_06->SetLineWidth(2);
        line_x_09->SetLineStyle(2);
        line_x_09->SetLineWidth(2);
        line_y_09->SetLineStyle(2);
        line_y_09->SetLineWidth(2);
        line_x_06->Draw("SAME");
        line_y_06->Draw("SAME");
        line_x_09->Draw("SAME");
        line_y_09->Draw("SAME");
        TLatex catLabel;
        catLabel.SetTextSize(0.04);
        catLabel.DrawLatex(0.72, 0.78, "A");
        catLabel.DrawLatex(0.30, 0.78, "B");
        catLabel.DrawLatex(0.30, 0.30, "C");
        catLabel.DrawLatex(0.72, 0.30, "D");
        if (generate_png) c_pairMass_cats->SaveAs("pairMass_A1B1_vs_A2B2_colz_cats.png");
        if (generate_pdf) c_pairMass_cats->SaveAs("pairMass_A1B1_vs_A2B2_colz_cats.pdf");

    }

    // pT overlay
    if (export_pT && (generate_png || generate_pdf)) {
        TCanvas *c_pT = new TCanvas("c_pairs_pT_overlay", "Pairs 13,24,14,23 pT Overlay", 900, 700);
        h_pT13->SetLineColor(kBlue);
        h_pT24->SetLineColor(kRed);
        h_pT14->SetLineColor(kGreen+2);
        h_pT23->SetLineColor(kMagenta);
        h_pT13->SetLineWidth(2);
        h_pT24->SetLineWidth(2);
        h_pT14->SetLineWidth(2);
        h_pT23->SetLineWidth(2);
        float ymax_pT = 1.2 * std::max(std::max(h_pT13->GetMaximum(), h_pT24->GetMaximum()), std::max(h_pT14->GetMaximum(), h_pT23->GetMaximum()));
        h_pT13->SetMaximum(ymax_pT);
        h_pT13->Draw();
        h_pT24->Draw("SAME");
        h_pT14->Draw("SAME");
        h_pT23->Draw("SAME");
        auto leg_pT = new TLegend(0.64, 0.72, 0.84, 0.88);
        leg_pT->SetTextSize(0.03);
        leg_pT->SetFillStyle(0);
        leg_pT->SetBorderSize(0);
        leg_pT->AddEntry(h_pT13, "A1", "l");
        leg_pT->AddEntry(h_pT24, "A2", "l");
        leg_pT->AddEntry(h_pT14, "B1", "l");
        leg_pT->AddEntry(h_pT23, "B2", "l");
        leg_pT->Draw();
        if (generate_png) c_pT->SaveAs("pairs_pT_13_24_14_23_overlay.png");
        if (generate_pdf) c_pT->SaveAs("pairs_pT_13_24_14_23_overlay.pdf");
    }

    // Per-pair pT overlays (unfiltered vs filtered) in one multipage PDF
    std::array<TH1F*,4> pT_all   = {h_pT13, h_pT24, h_pT14, h_pT23};
    std::array<TH1F*,4> pT_win   = {h_pT13_win, h_pT24_win, h_pT14_win, h_pT23_win};
    std::array<const char*,4> pair_tags   = {"13", "24", "14", "23"};
    std::array<const char*,4> pair_labels = {"A1", "A2", "B1", "B2"};
    std::string per_pair_pdf = "pairs_pT_per_pair.pdf";

    if (export_pT && (generate_png || generate_pdf)) {
        for (size_t idx = 0; idx < 4; ++idx) {
            TCanvas *c_pair = new TCanvas(Form("c_pT_pair_%s", pair_tags[idx]), Form("pT %s: unfiltered vs filtered", pair_tags[idx]), 800, 650);
            pT_all[idx]->SetLineColor(kBlack);
            pT_all[idx]->SetLineWidth(2);
            pT_win[idx]->SetLineColor(kRed);
            pT_win[idx]->SetLineStyle(2);
            pT_win[idx]->SetLineWidth(2);
            float ymax_pair = 1.2 * std::max(pT_all[idx]->GetMaximum(), pT_win[idx]->GetMaximum());
            pT_all[idx]->SetMaximum(ymax_pair);
            pT_all[idx]->Draw();
            pT_win[idx]->Draw("SAME");
            auto leg_pair = new TLegend(0.64, 0.74, 0.86, 0.88);
            leg_pair->SetTextSize(0.03);
            leg_pair->SetFillStyle(0);
            leg_pair->SetBorderSize(0);
            leg_pair->AddEntry(pT_all[idx], pair_labels[idx], "l");
            leg_pair->AddEntry(pT_win[idx], "Filtered 0.6<M<0.9", "l");
            leg_pair->Draw();
            std::string pdfPage = per_pair_pdf;
            if (idx == 0) pdfPage += "(";           // open multipage
            else if (idx == 3) pdfPage += ")";      // close multipage
            if (generate_pdf) c_pair->SaveAs(pdfPage.c_str());
            if (generate_png) c_pair->SaveAs(Form("pairs_pT_%s_unfiltered_vs_filtered.png", pair_tags[idx]));
            delete c_pair;
        }
    }

    // pT overlay with mass window 0.6-0.9
    if (export_pT && (generate_png || generate_pdf)) {
        TCanvas *c_pT_win = new TCanvas("c_pairs_pT_overlay_win", "Pairs 13,24,14,23 pT Overlay (0.6<M<0.9)", 900, 700);
        h_pT13_win->SetLineColor(kBlue);
        h_pT24_win->SetLineColor(kRed);
        h_pT14_win->SetLineColor(kGreen+2);
        h_pT23_win->SetLineColor(kMagenta);
        h_pT13_win->SetLineStyle(1);
        h_pT24_win->SetLineStyle(1);
        h_pT14_win->SetLineStyle(1);
        h_pT23_win->SetLineStyle(1);
        h_pT13_win->SetLineWidth(2);
        h_pT24_win->SetLineWidth(2);
        h_pT14_win->SetLineWidth(2);
        h_pT23_win->SetLineWidth(2);
        float ymax_pT_win = 1.2 * std::max(std::max(h_pT13_win->GetMaximum(), h_pT24_win->GetMaximum()), std::max(h_pT14_win->GetMaximum(), h_pT23_win->GetMaximum()));
        h_pT13_win->SetMaximum(ymax_pT_win);
        h_pT13_win->Draw();
        h_pT24_win->Draw("SAME");
        h_pT14_win->Draw("SAME");
        h_pT23_win->Draw("SAME");
        auto leg_pT_win = new TLegend(0.64, 0.72, 0.84, 0.88);
        leg_pT_win->SetTextSize(0.03);
        leg_pT_win->SetFillStyle(0);
        leg_pT_win->SetBorderSize(0);
        leg_pT_win->AddEntry(h_pT13_win, "A1", "l");
        leg_pT_win->AddEntry(h_pT24_win, "A2", "l");
        leg_pT_win->AddEntry(h_pT14_win, "B1", "l");
        leg_pT_win->AddEntry(h_pT23_win, "B2", "l");
        leg_pT_win->Draw();
        if (generate_png) c_pT_win->SaveAs("pairs_pT_13_24_14_23_overlay_M0p6_0p9.png");
        if (generate_pdf) c_pT_win->SaveAs("pairs_pT_13_24_14_23_overlay_M0p6_0p9.pdf");
    }

    // Summed pair pT overlay (sum of counts of four pair histos)
    if (export_pT_sum && (generate_png || generate_pdf)) {
        auto drawSumPt = [&](TH1F *histAll, TH1F *histWin, const char *tag, const char *title) {
            TCanvas *c_pT_sum = new TCanvas(Form("c_pairs_pT_sum_%s", tag), title, 900, 700);
            histAll->SetLineColor(kBlack);
            histAll->SetLineWidth(2);
            histWin->SetLineColor(kRed);
            histWin->SetLineStyle(2);
            histWin->SetLineWidth(2);
            float ymax_pT_sum = 1.2f * std::max(static_cast<float>(histAll->GetMaximum()), static_cast<float>(histWin->GetMaximum()));
            histAll->SetMaximum(ymax_pT_sum);
            histAll->Draw();
            histWin->Draw("SAME");
            auto leg_pT_sum = new TLegend(0.60, 0.72, 0.86, 0.88);
            leg_pT_sum->SetTextSize(0.03);
            leg_pT_sum->SetFillStyle(0);
            leg_pT_sum->SetBorderSize(0);
            leg_pT_sum->AddEntry(histAll, "All", "l");
            leg_pT_sum->AddEntry(histWin, "0.6<M<0.9", "l");
            leg_pT_sum->Draw();
            if (generate_png) c_pT_sum->SaveAs(Form("pairs_pT_sum_overlay_%s.png", tag));
            if (generate_pdf) c_pT_sum->SaveAs(Form("pairs_pT_sum_overlay_%s.pdf", tag));
        };

        drawSumPt(h_pT_sum, h_pT_sum_win, "all", "Sum of pair pT (no total pT cut)");
        drawSumPt(h_pT_sum_totalpt, h_pT_sum_win_totalpt, "totalpt_lt0p1", "Sum of pair pT (total pT<0.1)");
        drawSumPt(h_pT_sum_totalpt_mid, h_pT_sum_win_totalpt_mid, "totalpt_0p1_0p2", "Sum of pair pT (0.1<=total pT<0.2)");
        drawSumPt(h_pT_sum_totalpt_high, h_pT_sum_win_totalpt_high, "totalpt_ge0p2", "Sum of pair pT (total pT>=0.2)");
    }

    // Summed pair pT drawn as histograms
    if (export_pT_sum && (generate_png || generate_pdf)) {
        TCanvas *c_pT_sum_pts = new TCanvas("c_pairs_pT_sum_points", "Sum pair pT (histogram)", 900, 700);
        h_pT_sum_totalpt->SetLineColor(kBlack);
        h_pT_sum_totalpt->SetLineWidth(2);
        h_pT_sum_win_totalpt->SetLineColor(kRed);
        h_pT_sum_win_totalpt->SetLineStyle(2);
        h_pT_sum_win_totalpt->SetLineWidth(2);
        float ymax_pT_sum_pts = 1.2 * std::max(h_pT_sum_totalpt->GetMaximum(), h_pT_sum_win_totalpt->GetMaximum());
        h_pT_sum_totalpt->SetMaximum(ymax_pT_sum_pts*0.8);
        h_pT_sum_totalpt->Draw("HIST");
        h_pT_sum_win_totalpt->Draw("HIST SAME");
        auto leg_pT_sum_pts = new TLegend(0.54, 0.72, 0.76, 0.88);
        leg_pT_sum_pts->SetTextSize(0.03);
        leg_pT_sum_pts->SetFillStyle(0);
        leg_pT_sum_pts->SetBorderSize(0);
        leg_pT_sum_pts->AddEntry(h_pT_sum_totalpt, "Total p_{T}<0.1 (sum)", "l");
        leg_pT_sum_pts->AddEntry(h_pT_sum_win_totalpt, "0.6<M<0.9, total p_{T}<0.1 (sum)", "l");
        leg_pT_sum_pts->Draw();
        {
            TLatex note;
            note.SetNDC();
            note.SetTextSize(0.03);
            note.DrawLatex(0.60, 0.66, "total p_{T}<0.1");
        }
        if (generate_png) c_pT_sum_pts->SaveAs("pairs_pT_sum_overlay_points.png");
        if (generate_pdf) c_pT_sum_pts->SaveAs("pairs_pT_sum_overlay_points.pdf");
    }

    // Filtered sum (per A/B drop rule) with error bars
    if (export_sum_filtered && (generate_png || generate_pdf)) {
        TCanvas *c_pT_sum_filt_pts = new TCanvas("c_pairs_pT_sum_filtered_points", "Filtered sum pair pT (A/B rule)", 900, 700);
        int nb_f = h_pT_sum_filtered->GetNbinsX();
        std::vector<double> xf(nb_f), yf(nb_f), exf(nb_f, 0.0), eyf(nb_f);
        int firstNonZero_f = -1, lastNonZero_f = -1;
        for (int i=1; i<=nb_f; ++i) {
            double xc = h_pT_sum_filtered->GetBinCenter(i);
            double c = h_pT_sum_filtered->GetBinContent(i);
            xf[i-1] = xc; yf[i-1] = c; eyf[i-1] = (c > 0.0) ? std::sqrt(c) : 0.0;
            if (c > 0.0) {
                if (firstNonZero_f < 0) firstNonZero_f = i;
                lastNonZero_f = i;
            }
        }
        TGraphErrors *g_sum_filt = new TGraphErrors(nb_f, xf.data(), yf.data(), exf.data(), eyf.data());
        g_sum_filt->SetTitle("Filtered sum (A/B rule);p_{T} (GeV);Counts");
        g_sum_filt->SetLineColor(kBlack);
        g_sum_filt->SetMarkerColor(kBlack);
        g_sum_filt->SetMarkerStyle(21);
        g_sum_filt->SetMarkerSize(1.1);
        g_sum_filt->SetLineWidth(2);
        if (firstNonZero_f > 0 && lastNonZero_f > 0 && lastNonZero_f >= firstNonZero_f) {
            double xmin = h_pT_sum_filtered->GetBinLowEdge(firstNonZero_f);
            double xmax = h_pT_sum_filtered->GetBinLowEdge(lastNonZero_f) + h_pT_sum_filtered->GetBinWidth(lastNonZero_f);
            double marginLeft = 0.02 * (xmax - xmin);
            xmin = std::max(0.0, xmin - marginLeft);
            g_sum_filt->GetXaxis()->SetLimits(xmin, xmax);
        }
        g_sum_filt->Draw("AP");
        auto leg_pT_sum_filt = new TLegend(0.66, 0.72, 0.86, 0.84);
        leg_pT_sum_filt->SetTextSize(0.03);
        leg_pT_sum_filt->SetFillStyle(0);
        leg_pT_sum_filt->SetBorderSize(0);
        leg_pT_sum_filt->AddEntry(g_sum_filt, "A/B rule selection", "p");
        leg_pT_sum_filt->Draw();
        if (generate_png) c_pT_sum_filt_pts->SaveAs("pairs_pT_sum_filtered_points.png");
        if (generate_pdf) c_pT_sum_filt_pts->SaveAs("pairs_pT_sum_filtered_points.pdf");
    }

    // Min-pair selection (|pT1|+|pT2|) points with errors
    if (export_sum_minpair && (generate_png || generate_pdf)) {
        TCanvas *c_pT_sum_min_pts = new TCanvas("c_pairs_pT_sum_minpair_points", "Min-pair sum selection", 900, 700);
        int nb_m = h_pT_sum_minpair->GetNbinsX();
        std::vector<double> xm(nb_m), ym(nb_m), exm(nb_m, 0.0), eym(nb_m);
        int firstNonZero_m = -1, lastNonZero_m = -1;
        for (int i=1; i<=nb_m; ++i) {
            double xc = h_pT_sum_minpair->GetBinCenter(i);
            double c = h_pT_sum_minpair->GetBinContent(i);
            xm[i-1] = xc; ym[i-1] = c; eym[i-1] = (c > 0.0) ? std::sqrt(c) : 0.0;
            if (c > 0.0) {
                if (firstNonZero_m < 0) firstNonZero_m = i;
                lastNonZero_m = i;
            }
        }
        TGraphErrors *g_sum_min = new TGraphErrors(nb_m, xm.data(), ym.data(), exm.data(), eym.data());
        g_sum_min->SetTitle("Min-pair (|p_{T1}|+|p_{T2}|);p_{T} (GeV);Counts");
        g_sum_min->SetLineColor(kBlue+2);
        g_sum_min->SetMarkerColor(kBlue+2);
        g_sum_min->SetMarkerStyle(22);
        g_sum_min->SetMarkerSize(1.1);
        g_sum_min->SetLineWidth(2);
        if (firstNonZero_m > 0 && lastNonZero_m > 0 && lastNonZero_m >= firstNonZero_m) {
            double xmin = h_pT_sum_minpair->GetBinLowEdge(firstNonZero_m);
            double xmax = h_pT_sum_minpair->GetBinLowEdge(lastNonZero_m) + h_pT_sum_minpair->GetBinWidth(lastNonZero_m);
            double marginLeft = 0.02 * (xmax - xmin);
            xmin = std::max(0.0, xmin - marginLeft);
            g_sum_min->GetXaxis()->SetLimits(xmin, xmax);
        }
        g_sum_min->Draw("AP");
        auto leg_pT_sum_min = new TLegend(0.16, 0.72, 0.36, 0.84);
        leg_pT_sum_min->SetTextSize(0.03);
        leg_pT_sum_min->SetFillStyle(0);
        leg_pT_sum_min->SetBorderSize(0);
        leg_pT_sum_min->AddEntry(g_sum_min, "Min(|pT1|+|pT2|) selection", "p");
        leg_pT_sum_min->Draw();
        if (generate_png) c_pT_sum_min_pts->SaveAs("pairs_pT_sum_minpair_points.png");
        if (generate_pdf) c_pT_sum_min_pts->SaveAs("pairs_pT_sum_minpair_points.pdf");
    }

    
        // Theta overlay
       /* TCanvas *c_theta = new TCanvas("c_theta_overlay", "Pair angles: theta", 900, 700);
        h_theta_A1->SetLineColor(kBlue);
        h_theta_A2->SetLineColor(kRed);
        h_theta_B1->SetLineColor(kGreen+2);
        h_theta_B2->SetLineColor(kMagenta);
        h_theta_A1->SetLineWidth(2);
        h_theta_A2->SetLineWidth(2);
        h_theta_B1->SetLineWidth(2);
        h_theta_B2->SetLineWidth(2);
        double ymax_theta = 1.2 * std::max(std::max(h_theta_A1->GetMaximum(), h_theta_A2->GetMaximum()), std::max(h_theta_B1->GetMaximum(), h_theta_B2->GetMaximum()));
        h_theta_A1->SetMaximum(ymax_theta);
        h_theta_A1->Draw();
        h_theta_A2->Draw("SAME");
        h_theta_B1->Draw("SAME");
        h_theta_B2->Draw("SAME");
        auto leg_theta = new TLegend(0.64, 0.72, 0.84, 0.88);
        leg_theta->SetTextSize(0.03);
        leg_theta->SetFillStyle(0);
        leg_theta->SetBorderSize(0);
        leg_theta->AddEntry(h_theta_A1, "Pair A1", "l");
        leg_theta->AddEntry(h_theta_A2, "Pair A2", "l");
        leg_theta->AddEntry(h_theta_B1, "Pair B1", "l");
        leg_theta->AddEntry(h_theta_B2, "Pair B2", "l");
        leg_theta->Draw();
        if (generate_png) c_theta->SaveAs("pairs_theta_overlay.png");
        if (generate_pdf) c_theta->SaveAs("pairs_theta_overlay.pdf");

        // Phi overlay
        TCanvas *c_phi = new TCanvas("c_phi_overlay", "Pair angles: phi", 900, 700);
        h_phi_A1->SetLineColor(kBlue);
        h_phi_A2->SetLineColor(kRed);
        h_phi_B1->SetLineColor(kGreen+2);
        h_phi_B2->SetLineColor(kMagenta);
        h_phi_A1->SetLineWidth(2);
        h_phi_A2->SetLineWidth(2);
        h_phi_B1->SetLineWidth(2);
        h_phi_B2->SetLineWidth(2);
        double ymax_phi = 1.2 * std::max(std::max(h_phi_A1->GetMaximum(), h_phi_A2->GetMaximum()), std::max(h_phi_B1->GetMaximum(), h_phi_B2->GetMaximum()));
        h_phi_A1->SetMaximum(ymax_phi);
        h_phi_A1->Draw();
        h_phi_A2->Draw("SAME");
        h_phi_B1->Draw("SAME");
        h_phi_B2->Draw("SAME");
        auto leg_phi = new TLegend(0.64, 0.72, 0.84, 0.88);
        leg_phi->SetTextSize(0.03);
        leg_phi->SetFillStyle(0);
        leg_phi->SetBorderSize(0);
        leg_phi->AddEntry(h_phi_A1, "Pair A1", "l");
        leg_phi->AddEntry(h_phi_A2, "Pair A2", "l");
        leg_phi->AddEntry(h_phi_B1, "Pair B1", "l");
        leg_phi->AddEntry(h_phi_B2, "Pair B2", "l");
        leg_phi->Draw();
        if (generate_png) c_phi->SaveAs("pairs_phi_overlay.png");
        if (generate_pdf) c_phi->SaveAs("pairs_phi_overlay.pdf");*/

        // Helper lambda to convert TH1F to TGraphErrors with sqrt(N) errors
        auto histToGraph = [](TH1F *h, int markerStyle, Color_t color) {
            int n = h->GetNbinsX();
            std::vector<double> x(n), y(n), ex(n, 0.0), ey(n);
            double ymax = 0.0;
            for (int i = 1; i <= n; ++i) {
                double xc = h->GetBinCenter(i);
                double c = h->GetBinContent(i);
                x[i-1] = xc; y[i-1] = c; ey[i-1] = (c > 0.0) ? std::sqrt(c) : 0.0;
                ymax = std::max(ymax, c);
            }
            auto g = new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data());
            g->SetTitle(h->GetTitle());
            g->SetMarkerStyle(markerStyle);
            g->SetMarkerColor(color);
            g->SetMarkerSize(1.1);
            g->SetLineColor(color);
            g->SetLineWidth(2);
            g->GetXaxis()->SetLimits(h->GetBinLowEdge(1), h->GetBinLowEdge(n) + h->GetBinWidth(n));
            g->SetMinimum(0.0);
            g->SetMaximum(1.2 * (ymax > 0.0 ? ymax : 1.0));
            return g;
        };

        
    // Angle overlays (theta and phi)
    if (export_pair_angles && export_theta_output && (generate_png || generate_pdf)) {
        // Theta: A-only and B-only, points with error bars
        TCanvas *c_theta_A_pts = new TCanvas("c_theta_A_points", "Theta A (points)", 900, 700);
        auto g_theta_A1 = histToGraph(h_theta_A1, 20, kBlue);
        g_theta_A1->Draw("AP");
        auto leg_theta_A = new TLegend(0.64, 0.72, 0.84, 0.82);
        leg_theta_A->SetTextSize(0.03);
        leg_theta_A->SetFillStyle(0);
        leg_theta_A->SetBorderSize(0);
        leg_theta_A->AddEntry(g_theta_A1, "A", "p");
        leg_theta_A->Draw();
        if (generate_png) c_theta_A_pts->SaveAs("pairs_theta_A_points.png");
        if (generate_pdf) c_theta_A_pts->SaveAs("pairs_theta_A_points.pdf");

        TCanvas *c_theta_B_pts = new TCanvas("c_theta_B_points", "Theta B (points)", 900, 700);
        auto g_theta_B1 = histToGraph(h_theta_B1, 21, kGreen+2);
        g_theta_B1->Draw("AP");
        auto leg_theta_B = new TLegend(0.64, 0.72, 0.84, 0.82);
        leg_theta_B->SetTextSize(0.03);
        leg_theta_B->SetFillStyle(0);
        leg_theta_B->SetBorderSize(0);
        leg_theta_B->AddEntry(g_theta_B1, "B", "p");
        leg_theta_B->Draw();
        if (generate_png) c_theta_B_pts->SaveAs("pairs_theta_B_points.png");
        if (generate_pdf) c_theta_B_pts->SaveAs("pairs_theta_B_points.pdf");
    }

    if (export_pair_angles && export_phi_output && (generate_png || generate_pdf)) {
        // Phi: A-only and B-only, points with error bars
        TCanvas *c_phi_A_pts = new TCanvas("c_phi_A_points", "", 900, 700);
        auto g_phi_A1 = histToGraph(h_phi_A1, 20, kBlue);
        g_phi_A1->SetTitle(";#Delta#varphi;Counts");
        g_phi_A1->Draw("AP");
        if (generate_png) c_phi_A_pts->SaveAs("pairs_phi_A_points.png");
        if (generate_pdf) c_phi_A_pts->SaveAs("pairs_phi_A_points.pdf");

        TCanvas *c_phi_B_pts = new TCanvas("c_phi_B_points", "", 900, 700);
        auto g_phi_B1 = histToGraph(h_phi_B1, 21, kGreen+2);
        g_phi_B1->SetTitle(";#Delta#varphi;Counts");
        g_phi_B1->Draw("AP");
        if (generate_png) c_phi_B_pts->SaveAs("pairs_phi_B_points.png");
        if (generate_pdf) c_phi_B_pts->SaveAs("pairs_phi_B_points.pdf");
    }

        // Per-pion angles: draw as points with error bars if available
        if (export_pion_angles && export_theta_output && ((pi_theta_read && pi_theta_read->size() == 4) || (h_pi_theta_1->GetEntries() > 0))) {
            TCanvas *c_pi_theta = new TCanvas("c_pi_theta_points", "Pion theta (points)", 900, 700);
            auto histToGraph = [](TH1F *h, int markerStyle, Color_t color) {
                int n = h->GetNbinsX();
                std::vector<double> x(n), y(n), ex(n, 0.0), ey(n);
                double ymax = 0.0;
                for (int i = 1; i <= n; ++i) {
                    double xc = h->GetBinCenter(i);
                    double c = h->GetBinContent(i);
                    x[i-1] = xc; y[i-1] = c; ey[i-1] = (c > 0.0) ? std::sqrt(c) : 0.0;
                    ymax = std::max(ymax, c);
                }
                auto g = new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data());
                g->SetTitle(h->GetTitle());
                g->SetMarkerStyle(markerStyle);
                g->SetMarkerColor(color);
                g->SetMarkerSize(1.1);
                g->SetLineColor(color);
                g->SetLineWidth(2);
                g->GetXaxis()->SetLimits(h->GetBinLowEdge(1), h->GetBinLowEdge(n) + h->GetBinWidth(n));
                g->SetMinimum(0.0);
                g->SetMaximum(1.2 * (ymax > 0.0 ? ymax : 1.0));
                return g;
            };
            auto g_pi_th1 = histToGraph(h_pi_theta_1, 20, kBlack);
            auto g_pi_th2 = histToGraph(h_pi_theta_2, 24, kRed);
            auto g_pi_th3 = histToGraph(h_pi_theta_3, 21, kBlue+2);
            auto g_pi_th4 = histToGraph(h_pi_theta_4, 22, kGreen+2);
            g_pi_th1->Draw("AP");
            g_pi_th2->Draw("P SAME");
            g_pi_th3->Draw("P SAME");
            g_pi_th4->Draw("P SAME");
            auto leg_pi_theta = new TLegend(0.64, 0.72, 0.86, 0.88);
            leg_pi_theta->SetTextSize(0.03);
            leg_pi_theta->SetFillStyle(0);
            leg_pi_theta->SetBorderSize(0);
            leg_pi_theta->AddEntry(g_pi_th1, "pi 1", "p");
            leg_pi_theta->AddEntry(g_pi_th2, "pi 2", "p");
            leg_pi_theta->AddEntry(g_pi_th3, "pi 3", "p");
            leg_pi_theta->AddEntry(g_pi_th4, "pi 4", "p");
            leg_pi_theta->Draw();
            if (generate_png) c_pi_theta->SaveAs("pi_theta_points.png");
            if (generate_pdf) c_pi_theta->SaveAs("pi_theta_points.pdf");
        }

        if (export_pion_angles && export_phi_output && ((pi_phi_read && pi_phi_read->size() == 4) || (h_pi_phi_1->GetEntries() > 0))) {
            TCanvas *c_pi_phi = new TCanvas("c_pi_phi_points", "", 900, 700);
            auto histToGraph2 = [](TH1F *h, int markerStyle, Color_t color) {
                int n = h->GetNbinsX();
                std::vector<double> x(n), y(n), ex(n, 0.0), ey(n);
                double ymax = 0.0;
                for (int i = 1; i <= n; ++i) {
                    double xc = h->GetBinCenter(i);
                    double c = h->GetBinContent(i);
                    x[i-1] = xc; y[i-1] = c; ey[i-1] = (c > 0.0) ? std::sqrt(c) : 0.0;
                    ymax = std::max(ymax, c);
                }
                auto g = new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data());
                g->SetTitle(h->GetTitle());
                g->SetMarkerStyle(markerStyle);
                g->SetMarkerColor(color);
                g->SetMarkerSize(1.1);
                g->SetLineColor(color);
                g->SetLineWidth(2);
                g->GetXaxis()->SetLimits(h->GetBinLowEdge(1), h->GetBinLowEdge(n) + h->GetBinWidth(n));
                g->SetMinimum(0.0);
                g->SetMaximum(1.2 * (ymax > 0.0 ? ymax : 1.0));
                return g;
            };
            auto g_pi_ph1 = histToGraph2(h_pi_phi_1, 20, kBlack);
            auto g_pi_ph2 = histToGraph2(h_pi_phi_2, 24, kRed);
            auto g_pi_ph3 = histToGraph2(h_pi_phi_3, 21, kBlue+2);
            auto g_pi_ph4 = histToGraph2(h_pi_phi_4, 22, kGreen+2);
            g_pi_ph1->SetTitle(";#Delta#varphi;Counts");
            g_pi_ph1->Draw("AP");
            g_pi_ph2->Draw("P SAME");
            g_pi_ph3->Draw("P SAME");
            g_pi_ph4->Draw("P SAME");
            if (generate_png) c_pi_phi->SaveAs("pi_phi_points.png");
            if (generate_pdf) c_pi_phi->SaveAs("pi_phi_points.pdf");
        }

        // NEW: Summed bin-count graphs for pair angles (A1+A2+B1+B2)
        // Build summed theta and phi histograms from the four pair histograms
        TH1F *h_theta_sum = (TH1F*)h_theta_A1->Clone("h_theta_sum");
        h_theta_sum->Reset();
        h_theta_sum->SetTitle("Sum of pair #theta (A1+A2+B1+B2);#theta;Counts");
        h_theta_sum->Add(h_theta_A1);
        h_theta_sum->Add(h_theta_A2);
        h_theta_sum->Add(h_theta_B1);
        h_theta_sum->Add(h_theta_B2);

        TH1F *h_phi_sum = (TH1F*)h_phi_A1->Clone("h_phi_sum");
        h_phi_sum->Reset();
        h_phi_sum->SetTitle(";#Delta#varphi;Counts");
        h_phi_sum->Add(h_phi_A1);
        h_phi_sum->Add(h_phi_A2);
        h_phi_sum->Add(h_phi_B1);
        h_phi_sum->Add(h_phi_B2);

        // Convert to TGraphErrors with sqrt(N) error bars and save
        if (export_theta_output) {
            TCanvas *c_theta_sum_pts = new TCanvas("c_theta_sum_points", "Theta sum (A1+A2+B1+B2)", 900, 700);
            auto g_theta_sum = histToGraph(h_theta_sum, 20, kBlack);
            g_theta_sum->Draw("AP");
            if (generate_png) c_theta_sum_pts->SaveAs("pairs_theta_sum_points.png");
            if (generate_pdf) c_theta_sum_pts->SaveAs("pairs_theta_sum_points.pdf");
        }

        if (export_phi_output) {
            TCanvas *c_phi_sum_pts = new TCanvas("c_phi_sum_points", "", 900, 700);
            auto g_phi_sum = histToGraph(h_phi_sum, 20, kBlack);
            g_phi_sum->SetTitle(";#Delta#varphi;Counts");
            g_phi_sum->Draw("AP");
            if (generate_png) c_phi_sum_pts->SaveAs("pairs_phi_sum_points.png");
            if (generate_pdf) c_phi_sum_pts->SaveAs("pairs_phi_sum_points.pdf");
        }

        // Persist summed histograms into the ROOT file as well
        fplot->cd();
        if (export_theta_output) h_theta_sum->Write();
        if (export_phi_output) h_phi_sum->Write();

        // NEW: Summed bin-count graphs for per-pion angles (A1+A2+B1+B2)
        TH1F *h_pi_theta_sum = (TH1F*)h_pi_theta_1->Clone("h_pi_theta_sum");
        h_pi_theta_sum->Reset();
        h_pi_theta_sum->SetTitle("Sum of #pi^{+} #theta (A1+A2+B1+B2);#theta;Counts");
        h_pi_theta_sum->Add(h_pi_theta_1);
        h_pi_theta_sum->Add(h_pi_theta_2);
        h_pi_theta_sum->Add(h_pi_theta_3);
        h_pi_theta_sum->Add(h_pi_theta_4);

        TH1F *h_pi_theta_sum_win = (TH1F*)h_pi_theta_1_win->Clone("h_pi_theta_sum_win");
        h_pi_theta_sum_win->Reset();
        h_pi_theta_sum_win->SetTitle("Sum of #pi^{+} #theta (0.6<M<0.9);#theta;Counts");
        h_pi_theta_sum_win->Add(h_pi_theta_1_win);
        h_pi_theta_sum_win->Add(h_pi_theta_2_win);
        h_pi_theta_sum_win->Add(h_pi_theta_3_win);
        h_pi_theta_sum_win->Add(h_pi_theta_4_win);

        TH1F *h_pi_phi_sum = (TH1F*)h_pi_phi_1->Clone("h_pi_phi_sum");
        h_pi_phi_sum->Reset();
        h_pi_phi_sum->SetTitle(";#Delta#varphi;Counts");
        h_pi_phi_sum->Add(h_pi_phi_1);
        h_pi_phi_sum->Add(h_pi_phi_2);
        h_pi_phi_sum->Add(h_pi_phi_3);
        h_pi_phi_sum->Add(h_pi_phi_4);

        TH1F *h_pi_phi_sum_win = (TH1F*)h_pi_phi_1_win->Clone("h_pi_phi_sum_win");
        h_pi_phi_sum_win->Reset();
        h_pi_phi_sum_win->SetTitle(";#Delta#varphi;Counts");
        h_pi_phi_sum_win->Add(h_pi_phi_1_win);
        h_pi_phi_sum_win->Add(h_pi_phi_2_win);
        h_pi_phi_sum_win->Add(h_pi_phi_3_win);
        h_pi_phi_sum_win->Add(h_pi_phi_4_win);

        // Flexible cos(n*phi) modulation fit: f(phi) = a0 + sum_{n=1..N} a_n cos(n*phi)
        auto buildCosSeries = [&](int order) {
            std::string expr = "[0]*(1";
            for (int n = 1; n <= order; ++n) {
                if (!only_even_cosine || (n % 2 == 0))
                expr += Form("+[%d]*cos(%d*x)", n, n);
            }
            expr += ")";
            return expr;
        };
        int fitOrder = std::max(1, pionPhiFitOrder);
        std::string cosExpr = buildCosSeries(fitOrder);
        if (export_pair_phi_xphi_selected && (generate_png || generate_pdf)) {
            auto drawPairPhiFit = [&](TH1F *hist, const char *tag) {
                if (!hist || hist->GetEntries() <= 0) return;
                TCanvas *c_pair = new TCanvas(Form("c_%s", tag), "", 900, 700);
                auto g_pair = histToGraph(hist, 24, kBlue+2);
                g_pair->SetTitle(";#Delta#varphi;Counts");
                g_pair->Draw("AP");
                TF1 *f_pair = new TF1(Form("f_%s", tag), cosExpr.c_str(), phiMin, phiMax);
                f_pair->SetNpx(500);
                f_pair->SetParName(0, "A0");
                f_pair->SetParameter(0, hist->GetMaximum());
                for (int n = 1; n <= fitOrder; ++n) {
                    if (!only_even_cosine || (n % 2 == 0)) {
                        f_pair->SetParName(n, Form("A%d", n));
                        f_pair->SetParameter(n, 0.05);
                    }
                }
                hist->Fit(f_pair, "Q0");
                f_pair->SetLineColor(kRed);
                f_pair->SetLineWidth(2);
                f_pair->Draw("SAME");
                double chi2_pair = f_pair->GetChisquare();
                double ndf_pair = static_cast<double>(f_pair->GetNDF());
                auto leg_pair = new TLegend(0.13, 0.14, 0.43, 0.40);
                leg_pair->SetTextSize(0.03);
                leg_pair->SetFillStyle(0);
                leg_pair->SetBorderSize(0);
                for (int n = 0; n <= fitOrder; ++n) {
                    if (!only_even_cosine || (n % 2 == 0)) {
                        leg_pair->AddEntry((TObject*)nullptr,
                                           Form("A%d = %.3f #pm %.3f", n, f_pair->GetParameter(n), f_pair->GetParError(n)),
                                           "");
                    }
                }
                leg_pair->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_pair, (int)ndf_pair), "");
                leg_pair->Draw();
                {
                    TLatex note;
                    note.SetNDC();
                    note.SetTextSize(0.03);
                    note.DrawLatex(0.58, 0.86, "0.6<M<0.9, p_{T}<0.1");
                }
                if (generate_png) c_pair->SaveAs(Form("%s.png", tag));
                if (generate_pdf) c_pair->SaveAs(Form("%s.pdf", tag));
            };

            drawPairPhiFit(h_pair_phi_sel, "pair_phi_masswin_ptcut");
            drawPairPhiFit(h_pair_xphi_sel, "pair_xphi_masswin_ptcut");
        }
        if (export_theta_output) {
            TCanvas *c_pi_theta_sum_pts = new TCanvas("c_pi_theta_sum_points", "Pion theta sum (A1+A2+B1+B2)", 900, 700);
            auto g_pi_theta_sum = histToGraph(h_pi_theta_sum, 20, kBlack);
            g_pi_theta_sum->Draw("AP");

            // Old cosine-series fit kept for reference (not executed)
            TF1 *f_pi_theta_mod_series = new TF1("f_pi_theta_mod_series", cosExpr.c_str(), thetaMin, thetaMax);
            f_pi_theta_mod_series->SetNpx(500);
            f_pi_theta_mod_series->SetParName(0, "A0");
            f_pi_theta_mod_series->SetParameter(0, h_pi_theta_sum->GetMaximum());
            for (int n = 1; n <= fitOrder; ++n) {
                if (!only_even_cosine || (n % 2 == 0)) {
                    f_pi_theta_mod_series->SetParName(n, Form("A%d", n));
                    f_pi_theta_mod_series->SetParameter(n, 0.05);
                }
            }
            // h_pi_theta_sum->Fit(f_pi_theta_mod_series, "Q0"); // not used

            // Active theta fit: A*(1 + B*cos^2(theta))
            TF1 *f_pi_theta_mod = new TF1("f_pi_theta_mod", "[0]*(1+[1]*pow(cos(x),2))", thetaMin, thetaMax);
            f_pi_theta_mod->SetNpx(500);
            f_pi_theta_mod->SetParName(0, "A");
            f_pi_theta_mod->SetParName(1, "B");
            f_pi_theta_mod->SetParameter(0, h_pi_theta_sum->GetMaximum());
            f_pi_theta_mod->SetParameter(1, 0.0);
            h_pi_theta_sum->Fit(f_pi_theta_mod, "Q0");
            f_pi_theta_mod->SetLineColor(kRed);
            f_pi_theta_mod->SetLineWidth(2);
            f_pi_theta_mod->Draw("SAME");

            // Print fit results
            std::cout << "pi_theta_sum fit: A(1+B cos^{2} theta)" << std::endl;
            std::cout << "  A = " << f_pi_theta_mod->GetParameter(0)
                      << " +/- " << f_pi_theta_mod->GetParError(0) << std::endl;
            std::cout << "  B = " << f_pi_theta_mod->GetParameter(1)
                      << " +/- " << f_pi_theta_mod->GetParError(1) << std::endl;

            // Fit quality
            double chi2_th = f_pi_theta_mod->GetChisquare();
            double ndf_th = static_cast<double>(f_pi_theta_mod->GetNDF());

            auto leg_fit_th = new TLegend(0.12, 0.60, 0.58, 0.86);
            leg_fit_th->SetTextSize(0.02);
            leg_fit_th->SetFillStyle(0);
            leg_fit_th->SetBorderSize(0);
            leg_fit_th->AddEntry(g_pi_theta_sum, "#pi^{+} #theta sum", "p");
            leg_fit_th->AddEntry(f_pi_theta_mod, "Fit: A(1+B cos^{2}#theta)", "l");
            leg_fit_th->AddEntry((TObject*)nullptr,
                                 Form("A = %.3f #pm %.3f", f_pi_theta_mod->GetParameter(0), f_pi_theta_mod->GetParError(0)),
                                 "");
            leg_fit_th->AddEntry((TObject*)nullptr,
                                 Form("B = %.3f #pm %.3f", f_pi_theta_mod->GetParameter(1), f_pi_theta_mod->GetParError(1)),
                                 "");
            leg_fit_th->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_th, (int)ndf_th), "");
            leg_fit_th->Draw();

            if (generate_png) c_pi_theta_sum_pts->SaveAs("pi_theta_sum_points.png");
            if (generate_pdf) c_pi_theta_sum_pts->SaveAs("pi_theta_sum_points.pdf");
        }

        if (export_phi_output) {
            TF1 *f_pi_phi_mod = new TF1("f_pi_phi_mod", cosExpr.c_str(), phiMin, phiMax);
            f_pi_phi_mod->SetNpx(500);
            f_pi_phi_mod->SetParName(0, "A0");
            for (int n = 1; n <= fitOrder; ++n) {
                if (!only_even_cosine || (n % 2 == 0)){
                    f_pi_phi_mod->SetParName(n, Form("A%d", n));
                    f_pi_phi_mod->SetParameter(n, 0.05); // mild modulation seed
                }
            }
            f_pi_phi_mod->SetParameter(0, h_pi_phi_sum->GetMaximum());
            h_pi_phi_sum->Fit(f_pi_phi_mod, "Q0"); // quiet, do not draw on hist

            TCanvas *c_pi_phi_sum_pts = new TCanvas("c_pi_phi_sum_points", "", 900, 700);
            auto g_pi_phi_sum = histToGraph(h_pi_phi_sum, 20, kBlack);
            g_pi_phi_sum->SetTitle(";#Delta#varphi;Counts");
            g_pi_phi_sum->Draw("AP");
            TGraphErrors *g_pi_phi_sum_win = nullptr;
            if (export_pion_angles_win) {
                g_pi_phi_sum_win = histToGraph(h_pi_phi_sum_win, 24, kBlue+2);
                g_pi_phi_sum_win->SetTitle(";#Delta#varphi;Counts");
                g_pi_phi_sum_win->Draw("P SAME");
            }
            f_pi_phi_mod->SetLineColor(kRed);
            f_pi_phi_mod->SetLineWidth(2);
            f_pi_phi_mod->Draw("SAME");
            std::cout << "pi_phi_sum cos-series fit (order " << fitOrder << ")" << std::endl;
            for (int n = 0; n <= fitOrder; ++n) {
                if (!only_even_cosine || (n % 2 == 0)) {
                std::cout << "  A" << n << " = " << f_pi_phi_mod->GetParameter(n)
                          << " +/- " << f_pi_phi_mod->GetParError(n) << std::endl;
                }
            }

            // Fit quality
            double chi2_phi = f_pi_phi_mod->GetChisquare();
            double ndf_phi = static_cast<double>(f_pi_phi_mod->GetNDF());

            // Legend with fit expression, parameters, and correlation coefficient
            auto leg_fit = new TLegend(0.18, 0.14, 0.48, 0.40);
            leg_fit->SetTextSize(0.03);
            leg_fit->SetFillStyle(0);
            leg_fit->SetBorderSize(0);
            for (int n = 0; n <= fitOrder; ++n) {
                if (!only_even_cosine || (n % 2 == 0)) {
                leg_fit->AddEntry((TObject*)nullptr,
                                  Form("A%d = %.3f #pm %.3f", n, f_pi_phi_mod->GetParameter(n), f_pi_phi_mod->GetParError(n)),
                                  "");
                }
            }
            leg_fit->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi, (int)ndf_phi), "");
            leg_fit->Draw();
            if (generate_png) c_pi_phi_sum_pts->SaveAs("pi_phi_sum_points.png");
            if (generate_pdf) c_pi_phi_sum_pts->SaveAs("pi_phi_sum_points.pdf");
        }

        if (export_pion_angles_win && export_theta_output && (generate_png || generate_pdf)) {
            // Filtered-only pi theta with fit
            TCanvas *c_pi_theta_sum_win_pts = new TCanvas("c_pi_theta_sum_win_points", "Pion theta sum (0.6<M<0.9)", 900, 700);
            auto g_pi_theta_sum_win = histToGraph(h_pi_theta_sum_win, 24, kBlue+2);
            g_pi_theta_sum_win->Draw("AP");
            TF1 *f_pi_theta_mod_win = new TF1("f_pi_theta_mod_win", "[0]*(1+[1]*pow(cos(x),2))", thetaMin, thetaMax);
            f_pi_theta_mod_win->SetNpx(500);
            f_pi_theta_mod_win->SetParName(0, "A");
            f_pi_theta_mod_win->SetParName(1, "B");
            f_pi_theta_mod_win->SetParameter(0, h_pi_theta_sum_win->GetMaximum());
            f_pi_theta_mod_win->SetParameter(1, 0.0);
            h_pi_theta_sum_win->Fit(f_pi_theta_mod_win, "Q0");
            f_pi_theta_mod_win->SetLineColor(kRed);
            f_pi_theta_mod_win->SetLineWidth(2);
            f_pi_theta_mod_win->Draw("SAME");
            double chi2_th_win = f_pi_theta_mod_win->GetChisquare();
            double ndf_th_win = static_cast<double>(f_pi_theta_mod_win->GetNDF());
            auto leg_fit_th_win = new TLegend(0.12, 0.60, 0.58, 0.86);
            leg_fit_th_win->SetTextSize(0.02);
            leg_fit_th_win->SetFillStyle(0);
            leg_fit_th_win->SetBorderSize(0);
            leg_fit_th_win->AddEntry(g_pi_theta_sum_win, "#pi^{+} #theta sum (0.6<M<0.9)", "p");
            leg_fit_th_win->AddEntry(f_pi_theta_mod_win, "Fit: A(1+B cos^{2}#theta)", "l");
            leg_fit_th_win->AddEntry((TObject*)nullptr,
                                     Form("A = %.3f #pm %.3f", f_pi_theta_mod_win->GetParameter(0), f_pi_theta_mod_win->GetParError(0)),
                                     "");
            leg_fit_th_win->AddEntry((TObject*)nullptr,
                                     Form("B = %.3f #pm %.3f", f_pi_theta_mod_win->GetParameter(1), f_pi_theta_mod_win->GetParError(1)),
                                     "");
            leg_fit_th_win->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_th_win, (int)ndf_th_win), "");
            leg_fit_th_win->Draw();
            if (generate_png) c_pi_theta_sum_win_pts->SaveAs("pi_theta_sum_win_points.png");
            if (generate_pdf) c_pi_theta_sum_win_pts->SaveAs("pi_theta_sum_win_points.pdf");
        }

        if (export_pion_angles_win && export_phi_output && (generate_png || generate_pdf)) {
            // Filtered-only pi phi with fit
            TCanvas *c_pi_phi_sum_win_pts = new TCanvas("c_pi_phi_sum_win_points", "", 900, 700);
            auto g_pi_phi_sum_win_only = histToGraph(h_pi_phi_sum_win, 24, kBlue+2);
            g_pi_phi_sum_win_only->SetTitle(";#Delta#varphi;Counts");
            g_pi_phi_sum_win_only->Draw("AP");
            TF1 *f_pi_phi_mod_win = new TF1("f_pi_phi_mod_win", cosExpr.c_str(), phiMin, phiMax);
            f_pi_phi_mod_win->SetNpx(500);
            f_pi_phi_mod_win->SetParName(0, "A0");
            f_pi_phi_mod_win->SetParameter(0, h_pi_phi_sum_win->GetMaximum());
            for (int n = 1; n <= fitOrder; ++n) {
                if (!only_even_cosine || (n % 2 == 0)) {
                    f_pi_phi_mod_win->SetParName(n, Form("A%d", n));
                    f_pi_phi_mod_win->SetParameter(n, 0.05);
                }
            }
            h_pi_phi_sum_win->Fit(f_pi_phi_mod_win, "Q0");
            f_pi_phi_mod_win->SetLineColor(kRed);
            f_pi_phi_mod_win->SetLineWidth(2);
            f_pi_phi_mod_win->Draw("SAME");
            double chi2_phi_win = f_pi_phi_mod_win->GetChisquare();
            double ndf_phi_win = static_cast<double>(f_pi_phi_mod_win->GetNDF());
            auto leg_fit_phi_win = new TLegend(0.18, 0.14, 0.48, 0.40);
            leg_fit_phi_win->SetTextSize(0.03);
            leg_fit_phi_win->SetFillStyle(0);
            leg_fit_phi_win->SetBorderSize(0);
            for (int n = 0; n <= fitOrder; ++n) {
                if (!only_even_cosine || (n % 2 == 0)) {
                leg_fit_phi_win->AddEntry((TObject*)nullptr,
                                          Form("A%d = %.3f #pm %.3f", n, f_pi_phi_mod_win->GetParameter(n), f_pi_phi_mod_win->GetParError(n)),
                                          "");
                }
            }
            leg_fit_phi_win->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_win, (int)ndf_phi_win), "");
            leg_fit_phi_win->Draw();
            if (generate_png) c_pi_phi_sum_win_pts->SaveAs("pi_phi_sum_win_points.png");
            if (generate_pdf) c_pi_phi_sum_win_pts->SaveAs("pi_phi_sum_win_points.pdf");
        }

        if (export_pion_angles_rho && export_theta_output && (generate_png || generate_pdf)) {
            // Rho-closest pi theta (best)
            TCanvas *c_pi_theta_rho_best = new TCanvas("c_pi_theta_rho_best", "Pion #theta (closest to #rho)", 900, 700);
            auto g_pi_theta_rho_best = histToGraph(h_pi_theta_rho_best, 24, kBlue+2);
            g_pi_theta_rho_best->Draw("AP");
            TF1 *f_pi_theta_rho_best = new TF1("f_pi_theta_rho_best", "[0]*(1+[1]*pow(cos(x),2))", thetaMin, thetaMax);
            f_pi_theta_rho_best->SetNpx(500);
            f_pi_theta_rho_best->SetParName(0, "A");
            f_pi_theta_rho_best->SetParName(1, "B");
            f_pi_theta_rho_best->SetParameter(0, h_pi_theta_rho_best->GetMaximum());
            f_pi_theta_rho_best->SetParameter(1, 0.0);
            h_pi_theta_rho_best->Fit(f_pi_theta_rho_best, "Q0");
            f_pi_theta_rho_best->SetLineColor(kRed);
            f_pi_theta_rho_best->SetLineWidth(2);
            f_pi_theta_rho_best->Draw("SAME");
            double chi2_theta_rho_best = f_pi_theta_rho_best->GetChisquare();
            double ndf_theta_rho_best = static_cast<double>(f_pi_theta_rho_best->GetNDF());
            auto leg_theta_rho_best = new TLegend(0.12, 0.60, 0.58, 0.86);
            leg_theta_rho_best->SetTextSize(0.02);
            leg_theta_rho_best->SetFillStyle(0);
            leg_theta_rho_best->SetBorderSize(0);
            leg_theta_rho_best->AddEntry(g_pi_theta_rho_best, "#pi^{+} #theta (closest #rho)", "p");
            leg_theta_rho_best->AddEntry(f_pi_theta_rho_best, "Fit: A(1+B cos^{2}#theta)", "l");
            leg_theta_rho_best->AddEntry((TObject*)nullptr,
                                         Form("A = %.3f #pm %.3f", f_pi_theta_rho_best->GetParameter(0), f_pi_theta_rho_best->GetParError(0)),
                                         "");
            leg_theta_rho_best->AddEntry((TObject*)nullptr,
                                         Form("B = %.3f #pm %.3f", f_pi_theta_rho_best->GetParameter(1), f_pi_theta_rho_best->GetParError(1)),
                                         "");
            leg_theta_rho_best->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_theta_rho_best, (int)ndf_theta_rho_best), "");
            leg_theta_rho_best->Draw();
            if (generate_png) c_pi_theta_rho_best->SaveAs("pi_theta_rho_best.png");
            if (generate_pdf) c_pi_theta_rho_best->SaveAs("pi_theta_rho_best.pdf");

            // Rho-closest counterpart pi theta
            TCanvas *c_pi_theta_rho_other = new TCanvas("c_pi_theta_rho_other", "Pion #theta (counterpart of #rho-closest)", 900, 700);
            auto g_pi_theta_rho_other = histToGraph(h_pi_theta_rho_other, 24, kMagenta+2);
            g_pi_theta_rho_other->Draw("AP");
            TF1 *f_pi_theta_rho_other = new TF1("f_pi_theta_rho_other", "[0]*(1+[1]*pow(cos(x),2))", thetaMin, thetaMax);
            f_pi_theta_rho_other->SetNpx(500);
            f_pi_theta_rho_other->SetParName(0, "A");
            f_pi_theta_rho_other->SetParName(1, "B");
            f_pi_theta_rho_other->SetParameter(0, h_pi_theta_rho_other->GetMaximum());
            f_pi_theta_rho_other->SetParameter(1, 0.0);
            h_pi_theta_rho_other->Fit(f_pi_theta_rho_other, "Q0");
            f_pi_theta_rho_other->SetLineColor(kRed);
            f_pi_theta_rho_other->SetLineWidth(2);
            f_pi_theta_rho_other->Draw("SAME");
            double chi2_theta_rho_other = f_pi_theta_rho_other->GetChisquare();
            double ndf_theta_rho_other = static_cast<double>(f_pi_theta_rho_other->GetNDF());
            auto leg_theta_rho_other = new TLegend(0.12, 0.60, 0.58, 0.86);
            leg_theta_rho_other->SetTextSize(0.02);
            leg_theta_rho_other->SetFillStyle(0);
            leg_theta_rho_other->SetBorderSize(0);
            leg_theta_rho_other->AddEntry(g_pi_theta_rho_other, "#pi^{+} #theta (counterpart)", "p");
            leg_theta_rho_other->AddEntry(f_pi_theta_rho_other, "Fit: A(1+B cos^{2}#theta)", "l");
            leg_theta_rho_other->AddEntry((TObject*)nullptr,
                                          Form("A = %.3f #pm %.3f", f_pi_theta_rho_other->GetParameter(0), f_pi_theta_rho_other->GetParError(0)),
                                          "");
            leg_theta_rho_other->AddEntry((TObject*)nullptr,
                                          Form("B = %.3f #pm %.3f", f_pi_theta_rho_other->GetParameter(1), f_pi_theta_rho_other->GetParError(1)),
                                          "");
            leg_theta_rho_other->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_theta_rho_other, (int)ndf_theta_rho_other), "");
            leg_theta_rho_other->Draw();
            if (generate_png) c_pi_theta_rho_other->SaveAs("pi_theta_rho_other.png");
            if (generate_pdf) c_pi_theta_rho_other->SaveAs("pi_theta_rho_other.pdf");

            if (export_pion_angles_rho_rest) {
                // Rho-closest remaining-two pi theta
                TCanvas *c_pi_theta_rho_rest = new TCanvas("c_pi_theta_rho_rest", "Pion #theta (remaining two, #rho selection)", 900, 700);
                auto g_pi_theta_rho_rest = histToGraph(h_pi_theta_rho_rest, 24, kOrange+7);
                g_pi_theta_rho_rest->Draw("AP");
                TF1 *f_pi_theta_rho_rest = new TF1("f_pi_theta_rho_rest", "[0]*(1+[1]*pow(cos(x),2))", thetaMin, thetaMax);
                f_pi_theta_rho_rest->SetNpx(500);
                f_pi_theta_rho_rest->SetParName(0, "A");
                f_pi_theta_rho_rest->SetParName(1, "B");
                f_pi_theta_rho_rest->SetParameter(0, h_pi_theta_rho_rest->GetMaximum());
                f_pi_theta_rho_rest->SetParameter(1, 0.0);
                h_pi_theta_rho_rest->Fit(f_pi_theta_rho_rest, "Q0");
                f_pi_theta_rho_rest->SetLineColor(kRed);
                f_pi_theta_rho_rest->SetLineWidth(2);
                f_pi_theta_rho_rest->Draw("SAME");
                double chi2_theta_rho_rest = f_pi_theta_rho_rest->GetChisquare();
                double ndf_theta_rho_rest = static_cast<double>(f_pi_theta_rho_rest->GetNDF());
                auto leg_theta_rho_rest = new TLegend(0.12, 0.60, 0.58, 0.86);
                leg_theta_rho_rest->SetTextSize(0.02);
                leg_theta_rho_rest->SetFillStyle(0);
                leg_theta_rho_rest->SetBorderSize(0);
                leg_theta_rho_rest->AddEntry(g_pi_theta_rho_rest, "#pi^{+} #theta (remaining)", "p");
                leg_theta_rho_rest->AddEntry(f_pi_theta_rho_rest, "Fit: A(1+B cos^{2}#theta)", "l");
                leg_theta_rho_rest->AddEntry((TObject*)nullptr,
                                              Form("A = %.3f #pm %.3f", f_pi_theta_rho_rest->GetParameter(0), f_pi_theta_rho_rest->GetParError(0)),
                                              "");
                leg_theta_rho_rest->AddEntry((TObject*)nullptr,
                                              Form("B = %.3f #pm %.3f", f_pi_theta_rho_rest->GetParameter(1), f_pi_theta_rho_rest->GetParError(1)),
                                              "");
                leg_theta_rho_rest->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_theta_rho_rest, (int)ndf_theta_rho_rest), "");
                leg_theta_rho_rest->Draw();
                if (generate_png) c_pi_theta_rho_rest->SaveAs("pi_theta_rho_rest.png");
                if (generate_pdf) c_pi_theta_rho_rest->SaveAs("pi_theta_rho_rest.pdf");
            }

            if (export_pion_angles_rho_x && (generate_png || generate_pdf)) {
                // Rho-closest pi theta_x (best)
                TCanvas *c_pi_theta_x_rho_best = new TCanvas("c_pi_theta_x_rho_best", "Pion #theta (closest to #rho)", 900, 700);
                auto g_pi_theta_x_rho_best = histToGraph(h_pi_theta_x_rho_best, 24, kBlue+2);
                g_pi_theta_x_rho_best->SetMinimum(0.0);
                g_pi_theta_x_rho_best->SetMaximum(h_pi_theta_x_rho_best->GetMaximum() * 2.0);
                g_pi_theta_x_rho_best->Draw("AP");
                TF1 *f_pi_theta_x_rho_best = new TF1("f_pi_theta_x_rho_best", "[0]*(1+[1]*pow(cos(x),2))", thetaMin, thetaMax);
                f_pi_theta_x_rho_best->SetNpx(500);
                f_pi_theta_x_rho_best->SetParName(0, "A");
                f_pi_theta_x_rho_best->SetParName(1, "B");
                f_pi_theta_x_rho_best->SetParameter(0, h_pi_theta_x_rho_best->GetMaximum());
                f_pi_theta_x_rho_best->SetParameter(1, 0.0);
                h_pi_theta_x_rho_best->Fit(f_pi_theta_x_rho_best, "Q0");
                f_pi_theta_x_rho_best->SetLineColor(kRed);
                f_pi_theta_x_rho_best->SetLineWidth(2);
                f_pi_theta_x_rho_best->Draw("SAME");
                double chi2_theta_x_rho_best = f_pi_theta_x_rho_best->GetChisquare();
                double ndf_theta_x_rho_best = static_cast<double>(f_pi_theta_x_rho_best->GetNDF());
                auto leg_theta_x_rho_best = new TLegend(0.12, 0.72, 0.58, 0.90);
                leg_theta_x_rho_best->SetTextSize(0.02);
                leg_theta_x_rho_best->SetFillStyle(0);
                leg_theta_x_rho_best->SetBorderSize(0);
                leg_theta_x_rho_best->AddEntry(g_pi_theta_x_rho_best, "#pi^{+} #theta (closest #rho)", "p");
                leg_theta_x_rho_best->AddEntry(f_pi_theta_x_rho_best, "Fit: A(1+B cos^{2}#theta)", "l");
                leg_theta_x_rho_best->AddEntry((TObject*)nullptr,
                                             Form("A = %.3f #pm %.3f", f_pi_theta_x_rho_best->GetParameter(0), f_pi_theta_x_rho_best->GetParError(0)),
                                             "");
                leg_theta_x_rho_best->AddEntry((TObject*)nullptr,
                                             Form("B = %.3f #pm %.3f", f_pi_theta_x_rho_best->GetParameter(1), f_pi_theta_x_rho_best->GetParError(1)),
                                             "");
                leg_theta_x_rho_best->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_theta_x_rho_best, (int)ndf_theta_x_rho_best), "");
                leg_theta_x_rho_best->Draw();
                if (generate_png) c_pi_theta_x_rho_best->SaveAs("pi_theta_x_rho_best.png");
                if (generate_pdf) c_pi_theta_x_rho_best->SaveAs("pi_theta_x_rho_best.pdf");

                // Rho-closest counterpart pi theta_x
                TCanvas *c_pi_theta_x_rho_other = new TCanvas("c_pi_theta_x_rho_other", "Pion #theta (counterpart of #rho-closest)", 900, 700);
                auto g_pi_theta_x_rho_other = histToGraph(h_pi_theta_x_rho_other, 24, kMagenta+2);
                g_pi_theta_x_rho_other->SetMinimum(0.0);
                g_pi_theta_x_rho_other->SetMaximum(h_pi_theta_x_rho_other->GetMaximum() * 2.0);
                g_pi_theta_x_rho_other->Draw("AP");
                TF1 *f_pi_theta_x_rho_other = new TF1("f_pi_theta_x_rho_other", "[0]*(1+[1]*pow(cos(x),2))", thetaMin, thetaMax);
                f_pi_theta_x_rho_other->SetNpx(500);
                f_pi_theta_x_rho_other->SetParName(0, "A");
                f_pi_theta_x_rho_other->SetParName(1, "B");
                f_pi_theta_x_rho_other->SetParameter(0, h_pi_theta_x_rho_other->GetMaximum());
                f_pi_theta_x_rho_other->SetParameter(1, 0.0);
                h_pi_theta_x_rho_other->Fit(f_pi_theta_x_rho_other, "Q0");
                f_pi_theta_x_rho_other->SetLineColor(kRed);
                f_pi_theta_x_rho_other->SetLineWidth(2);
                f_pi_theta_x_rho_other->Draw("SAME");
                double chi2_theta_x_rho_other = f_pi_theta_x_rho_other->GetChisquare();
                double ndf_theta_x_rho_other = static_cast<double>(f_pi_theta_x_rho_other->GetNDF());
                auto leg_theta_x_rho_other = new TLegend(0.12, 0.72, 0.58, 0.90);
                leg_theta_x_rho_other->SetTextSize(0.02);
                leg_theta_x_rho_other->SetFillStyle(0);
                leg_theta_x_rho_other->SetBorderSize(0);
                leg_theta_x_rho_other->AddEntry(g_pi_theta_x_rho_other, "#pi^{+} #theta (counterpart)", "p");
                leg_theta_x_rho_other->AddEntry(f_pi_theta_x_rho_other, "Fit: A(1+B cos^{2}#theta)", "l");
                leg_theta_x_rho_other->AddEntry((TObject*)nullptr,
                                              Form("A = %.3f #pm %.3f", f_pi_theta_x_rho_other->GetParameter(0), f_pi_theta_x_rho_other->GetParError(0)),
                                              "");
                leg_theta_x_rho_other->AddEntry((TObject*)nullptr,
                                              Form("B = %.3f #pm %.3f", f_pi_theta_x_rho_other->GetParameter(1), f_pi_theta_x_rho_other->GetParError(1)),
                                              "");
                leg_theta_x_rho_other->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_theta_x_rho_other, (int)ndf_theta_x_rho_other), "");
                leg_theta_x_rho_other->Draw();
                if (generate_png) c_pi_theta_x_rho_other->SaveAs("pi_theta_x_rho_other.png");
                if (generate_pdf) c_pi_theta_x_rho_other->SaveAs("pi_theta_x_rho_other.pdf");

                if (export_pion_angles_rho_x_rest) {
                    // Rho-closest remaining-two pi theta_x
                    TCanvas *c_pi_theta_x_rho_rest = new TCanvas("c_pi_theta_x_rho_rest", "Pion #theta (remaining two, #rho selection)", 900, 700);
                    auto g_pi_theta_x_rho_rest = histToGraph(h_pi_theta_x_rho_rest, 24, kOrange+7);
                    g_pi_theta_x_rho_rest->SetMinimum(0.0);
                    g_pi_theta_x_rho_rest->SetMaximum(h_pi_theta_x_rho_rest->GetMaximum() * 2.0);
                    g_pi_theta_x_rho_rest->Draw("AP");
                    TF1 *f_pi_theta_x_rho_rest = new TF1("f_pi_theta_x_rho_rest", "[0]*(1+[1]*pow(cos(x),2))", thetaMin, thetaMax);
                    f_pi_theta_x_rho_rest->SetNpx(500);
                    f_pi_theta_x_rho_rest->SetParName(0, "A");
                    f_pi_theta_x_rho_rest->SetParName(1, "B");
                    f_pi_theta_x_rho_rest->SetParameter(0, h_pi_theta_x_rho_rest->GetMaximum());
                    f_pi_theta_x_rho_rest->SetParameter(1, 0.0);
                    h_pi_theta_x_rho_rest->Fit(f_pi_theta_x_rho_rest, "Q0");
                    f_pi_theta_x_rho_rest->SetLineColor(kRed);
                    f_pi_theta_x_rho_rest->SetLineWidth(2);
                    f_pi_theta_x_rho_rest->Draw("SAME");
                    double chi2_theta_x_rho_rest = f_pi_theta_x_rho_rest->GetChisquare();
                    double ndf_theta_x_rho_rest = static_cast<double>(f_pi_theta_x_rho_rest->GetNDF());
                    auto leg_theta_x_rho_rest = new TLegend(0.12, 0.72, 0.58, 0.90);
                    leg_theta_x_rho_rest->SetTextSize(0.02);
                    leg_theta_x_rho_rest->SetFillStyle(0);
                    leg_theta_x_rho_rest->SetBorderSize(0);
                    leg_theta_x_rho_rest->AddEntry(g_pi_theta_x_rho_rest, "#pi^{+} #theta (remaining)", "p");
                    leg_theta_x_rho_rest->AddEntry(f_pi_theta_x_rho_rest, "Fit: A(1+B cos^{2}#theta)", "l");
                    leg_theta_x_rho_rest->AddEntry((TObject*)nullptr,
                                                  Form("A = %.3f #pm %.3f", f_pi_theta_x_rho_rest->GetParameter(0), f_pi_theta_x_rho_rest->GetParError(0)),
                                                  "");
                    leg_theta_x_rho_rest->AddEntry((TObject*)nullptr,
                                                  Form("B = %.3f #pm %.3f", f_pi_theta_x_rho_rest->GetParameter(1), f_pi_theta_x_rho_rest->GetParError(1)),
                                                  "");
                    leg_theta_x_rho_rest->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_theta_x_rho_rest, (int)ndf_theta_x_rho_rest), "");
                    leg_theta_x_rho_rest->Draw();
                    if (generate_png) c_pi_theta_x_rho_rest->SaveAs("pi_theta_x_rho_rest.png");
                    if (generate_pdf) c_pi_theta_x_rho_rest->SaveAs("pi_theta_x_rho_rest.pdf");
                }
            }

        }

        if (export_pion_angles_rho && export_phi_output && (generate_png || generate_pdf)) {
            // Rho-closest pi phi (best)
            TCanvas *c_pi_phi_rho_best = new TCanvas("c_pi_phi_rho_best", "", 900, 700);
            auto g_pi_phi_rho_best = histToGraph(h_pi_phi_rho_best, 24, kBlue+2);
            g_pi_phi_rho_best->SetTitle(";#Delta#varphi;Counts");
            g_pi_phi_rho_best->Draw("AP");
            TF1 *f_pi_phi_rho_best = new TF1("f_pi_phi_rho_best", cosExpr.c_str(), phiMin, phiMax);
            f_pi_phi_rho_best->SetNpx(500);
            f_pi_phi_rho_best->SetParName(0, "A0");
            f_pi_phi_rho_best->SetParameter(0, h_pi_phi_rho_best->GetMaximum());
            for (int n = 1; n <= fitOrder; ++n) {
                if (!only_even_cosine || (n % 2 == 0)) {
                    f_pi_phi_rho_best->SetParName(n, Form("A%d", n));
                    f_pi_phi_rho_best->SetParameter(n, 0.05);
                }
            }
            h_pi_phi_rho_best->Fit(f_pi_phi_rho_best, "Q0");
            f_pi_phi_rho_best->SetLineColor(kRed);
            f_pi_phi_rho_best->SetLineWidth(2);
            f_pi_phi_rho_best->Draw("SAME");
            double chi2_phi_rho_best = f_pi_phi_rho_best->GetChisquare();
            double ndf_phi_rho_best = static_cast<double>(f_pi_phi_rho_best->GetNDF());
            auto leg_phi_rho_best = new TLegend(0.13, 0.14, 0.43, 0.40);
            leg_phi_rho_best->SetTextSize(0.03);
            leg_phi_rho_best->SetFillStyle(0);
            leg_phi_rho_best->SetBorderSize(0);
            for (int n = 0; n <= fitOrder; ++n) {
                if (!only_even_cosine || (n % 2 == 0)) {
                leg_phi_rho_best->AddEntry((TObject*)nullptr,
                                           Form("A%d = %.3f #pm %.3f", n, f_pi_phi_rho_best->GetParameter(n), f_pi_phi_rho_best->GetParError(n)),
                                           "");
                }
            }
            leg_phi_rho_best->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_rho_best, (int)ndf_phi_rho_best), "");
            leg_phi_rho_best->Draw();
            if (generate_png) c_pi_phi_rho_best->SaveAs("pi_phi_rho_best.png");
            if (generate_pdf) c_pi_phi_rho_best->SaveAs("pi_phi_rho_best.pdf");

            // Rho-closest counterpart pi phi
            TCanvas *c_pi_phi_rho_other = new TCanvas("c_pi_phi_rho_other", "", 900, 700);
            auto g_pi_phi_rho_other = histToGraph(h_pi_phi_rho_other, 24, kMagenta+2);
            g_pi_phi_rho_other->SetTitle(";#Delta#varphi;Counts");
            g_pi_phi_rho_other->Draw("AP");
            TF1 *f_pi_phi_rho_other = new TF1("f_pi_phi_rho_other", cosExpr.c_str(), phiMin, phiMax);
            f_pi_phi_rho_other->SetNpx(500);
            f_pi_phi_rho_other->SetParName(0, "A0");
            f_pi_phi_rho_other->SetParameter(0, h_pi_phi_rho_other->GetMaximum());
            for (int n = 1; n <= fitOrder; ++n) {
                if (!only_even_cosine || (n % 2 == 0)) {
                    f_pi_phi_rho_other->SetParName(n, Form("A%d", n));
                    f_pi_phi_rho_other->SetParameter(n, 0.05);
                }
            }
            h_pi_phi_rho_other->Fit(f_pi_phi_rho_other, "Q0");
            f_pi_phi_rho_other->SetLineColor(kRed);
            f_pi_phi_rho_other->SetLineWidth(2);
            f_pi_phi_rho_other->Draw("SAME");
            double chi2_phi_rho_other = f_pi_phi_rho_other->GetChisquare();
            double ndf_phi_rho_other = static_cast<double>(f_pi_phi_rho_other->GetNDF());
            auto leg_phi_rho_other = new TLegend(0.13, 0.14, 0.43, 0.40);
            leg_phi_rho_other->SetTextSize(0.03);
            leg_phi_rho_other->SetFillStyle(0);
            leg_phi_rho_other->SetBorderSize(0);
            for (int n = 0; n <= fitOrder; ++n) {
                if (!only_even_cosine || (n % 2 == 0)) {
                leg_phi_rho_other->AddEntry((TObject*)nullptr,
                                            Form("A%d = %.3f #pm %.3f", n, f_pi_phi_rho_other->GetParameter(n), f_pi_phi_rho_other->GetParError(n)),
                                            "");
                }
            }
            leg_phi_rho_other->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_rho_other, (int)ndf_phi_rho_other), "");
            leg_phi_rho_other->Draw();
            if (generate_png) c_pi_phi_rho_other->SaveAs("pi_phi_rho_other.png");
            if (generate_pdf) c_pi_phi_rho_other->SaveAs("pi_phi_rho_other.pdf");

            if (export_pion_angles_rho_rest) {
                // Rho-closest remaining-two pi phi
                TCanvas *c_pi_phi_rho_rest = new TCanvas("c_pi_phi_rho_rest", "", 900, 700);
                auto g_pi_phi_rho_rest = histToGraph(h_pi_phi_rho_rest, 24, kOrange+7);
                g_pi_phi_rho_rest->SetTitle(";#Delta#varphi;Counts");
                g_pi_phi_rho_rest->Draw("AP");
                TF1 *f_pi_phi_rho_rest = new TF1("f_pi_phi_rho_rest", cosExpr.c_str(), phiMin, phiMax);
                f_pi_phi_rho_rest->SetNpx(500);
                f_pi_phi_rho_rest->SetParName(0, "A0");
                f_pi_phi_rho_rest->SetParameter(0, h_pi_phi_rho_rest->GetMaximum());
                for (int n = 1; n <= fitOrder; ++n) {
                    if (!only_even_cosine || (n % 2 == 0)) {
                        f_pi_phi_rho_rest->SetParName(n, Form("A%d", n));
                        f_pi_phi_rho_rest->SetParameter(n, 0.05);
                    }
                }
                h_pi_phi_rho_rest->Fit(f_pi_phi_rho_rest, "Q0");
                f_pi_phi_rho_rest->SetLineColor(kRed);
                f_pi_phi_rho_rest->SetLineWidth(2);
                f_pi_phi_rho_rest->Draw("SAME");
                double chi2_phi_rho_rest = f_pi_phi_rho_rest->GetChisquare();
                double ndf_phi_rho_rest = static_cast<double>(f_pi_phi_rho_rest->GetNDF());
                auto leg_phi_rho_rest = new TLegend(0.13, 0.14, 0.43, 0.40);
                leg_phi_rho_rest->SetTextSize(0.03);
                leg_phi_rho_rest->SetFillStyle(0);
                leg_phi_rho_rest->SetBorderSize(0);
                for (int n = 0; n <= fitOrder; ++n) {
                    if (!only_even_cosine || (n % 2 == 0)) {
                    leg_phi_rho_rest->AddEntry((TObject*)nullptr,
                                               Form("A%d = %.3f #pm %.3f", n, f_pi_phi_rho_rest->GetParameter(n), f_pi_phi_rho_rest->GetParError(n)),
                                               "");
                    }
                }
                leg_phi_rho_rest->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_rho_rest, (int)ndf_phi_rho_rest), "");
                leg_phi_rho_rest->Draw();
                if (generate_png) c_pi_phi_rho_rest->SaveAs("pi_phi_rho_rest.png");
                if (generate_pdf) c_pi_phi_rho_rest->SaveAs("pi_phi_rho_rest.pdf");
            }

            // Rho-closest pi phi_x (best)
            TCanvas *c_pi_phi_x_rho_best = new TCanvas("c_pi_phi_x_rho_best", "", 900, 700);
            auto g_pi_phi_x_rho_best = histToGraph(h_pi_phi_x_rho_best, 24, kBlue+2);
            g_pi_phi_x_rho_best->SetMinimum(0.0);
            g_pi_phi_x_rho_best->SetMaximum(h_pi_phi_x_rho_best->GetMaximum() * 1.6);
            g_pi_phi_x_rho_best->SetTitle(";#Delta#varphi;Counts");
            g_pi_phi_x_rho_best->Draw("AP");
            TF1 *f_pi_phi_x_rho_best = new TF1("f_pi_phi_x_rho_best", cosExpr.c_str(), phiMin, phiMax);
            f_pi_phi_x_rho_best->SetNpx(500);
            f_pi_phi_x_rho_best->SetParName(0, "A0");
            f_pi_phi_x_rho_best->SetParameter(0, h_pi_phi_x_rho_best->GetMaximum());
            for (int n = 1; n <= fitOrder; ++n) {
                if (!only_even_cosine || (n % 2 == 0)) {
                    f_pi_phi_x_rho_best->SetParName(n, Form("A%d", n));
                    f_pi_phi_x_rho_best->SetParameter(n, 0.05);
                }
            }
            h_pi_phi_x_rho_best->Fit(f_pi_phi_x_rho_best, "Q0");
            f_pi_phi_x_rho_best->SetLineColor(kRed);
            f_pi_phi_x_rho_best->SetLineWidth(2);
            f_pi_phi_x_rho_best->Draw("SAME");
            double chi2_phi_x_rho_best = f_pi_phi_x_rho_best->GetChisquare();
            double ndf_phi_x_rho_best = static_cast<double>(f_pi_phi_x_rho_best->GetNDF());
            auto leg_phi_x_rho_best = new TLegend(0.13, 0.63, 0.43, 0.88);
            leg_phi_x_rho_best->SetTextSize(0.03);
            leg_phi_x_rho_best->SetFillStyle(0);
            leg_phi_x_rho_best->SetBorderSize(0);
            for (int n = 0; n <= fitOrder; ++n) {
                if (!only_even_cosine || (n % 2 == 0)) {
                leg_phi_x_rho_best->AddEntry((TObject*)nullptr,
                                           Form("A%d = %.3f #pm %.3f", n, f_pi_phi_x_rho_best->GetParameter(n), f_pi_phi_x_rho_best->GetParError(n)),
                                           "");
                }
            }
            leg_phi_x_rho_best->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_x_rho_best, (int)ndf_phi_x_rho_best), "");
            leg_phi_x_rho_best->Draw();
            if (generate_png) c_pi_phi_x_rho_best->SaveAs("pi_phi_x_rho_best.png");
            if (generate_pdf) c_pi_phi_x_rho_best->SaveAs("pi_phi_x_rho_best.pdf");

            // Rho-closest counterpart pi phi_x
            TCanvas *c_pi_phi_x_rho_other = new TCanvas("c_pi_phi_x_rho_other", "", 900, 700);
            auto g_pi_phi_x_rho_other = histToGraph(h_pi_phi_x_rho_other, 24, kMagenta+2);
            g_pi_phi_x_rho_other->SetMinimum(0.0);
            g_pi_phi_x_rho_other->SetMaximum(h_pi_phi_x_rho_other->GetMaximum() * 2.0);
            g_pi_phi_x_rho_other->SetTitle(";#Delta#varphi;Counts");
            g_pi_phi_x_rho_other->Draw("AP");
            TF1 *f_pi_phi_x_rho_other = new TF1("f_pi_phi_x_rho_other", cosExpr.c_str(), phiMin, phiMax);
            f_pi_phi_x_rho_other->SetNpx(500);
            f_pi_phi_x_rho_other->SetParName(0, "A0");
            f_pi_phi_x_rho_other->SetParameter(0, h_pi_phi_x_rho_other->GetMaximum());
            for (int n = 1; n <= fitOrder; ++n) {
                if (!only_even_cosine || (n % 2 == 0)) {
                    f_pi_phi_x_rho_other->SetParName(n, Form("A%d", n));
                    f_pi_phi_x_rho_other->SetParameter(n, 0.05);
                }
            }
            h_pi_phi_x_rho_other->Fit(f_pi_phi_x_rho_other, "Q0");
            f_pi_phi_x_rho_other->SetLineColor(kRed);
            f_pi_phi_x_rho_other->SetLineWidth(2);
            f_pi_phi_x_rho_other->Draw("SAME");
            double chi2_phi_x_rho_other = f_pi_phi_x_rho_other->GetChisquare();
            double ndf_phi_x_rho_other = static_cast<double>(f_pi_phi_x_rho_other->GetNDF());
            auto leg_phi_x_rho_other = new TLegend(0.13, 0.58, 0.43, 0.88);
            leg_phi_x_rho_other->SetTextSize(0.03);
            leg_phi_x_rho_other->SetFillStyle(0);
            leg_phi_x_rho_other->SetBorderSize(0);
            for (int n = 0; n <= fitOrder; ++n) {
                if (!only_even_cosine || (n % 2 == 0)) {
                leg_phi_x_rho_other->AddEntry((TObject*)nullptr,
                                               Form("A%d = %.3f #pm %.3f", n, f_pi_phi_x_rho_other->GetParameter(n), f_pi_phi_x_rho_other->GetParError(n)),
                                               "");
                }
            }
            leg_phi_x_rho_other->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_x_rho_other, (int)ndf_phi_x_rho_other), "");
            leg_phi_x_rho_other->Draw();
            if (generate_png) c_pi_phi_x_rho_other->SaveAs("pi_phi_x_rho_other.png");
            if (generate_pdf) c_pi_phi_x_rho_other->SaveAs("pi_phi_x_rho_other.pdf");

            if (export_pion_angles_rho_x_rest) {
                // Rho-closest remaining-two pi phi_x
                TCanvas *c_pi_phi_x_rho_rest = new TCanvas("c_pi_phi_x_rho_rest", "", 900, 700);
                auto g_pi_phi_x_rho_rest = histToGraph(h_pi_phi_x_rho_rest, 24, kOrange+7);
                g_pi_phi_x_rho_rest->SetMinimum(0.0);
                g_pi_phi_x_rho_rest->SetMaximum(h_pi_phi_x_rho_rest->GetMaximum() * 2.0);
                g_pi_phi_x_rho_rest->SetTitle(";#Delta#varphi;Counts");
                g_pi_phi_x_rho_rest->Draw("AP");
                TF1 *f_pi_phi_x_rho_rest = new TF1("f_pi_phi_x_rho_rest", cosExpr.c_str(), phiMin, phiMax);
                f_pi_phi_x_rho_rest->SetNpx(500);
                f_pi_phi_x_rho_rest->SetParName(0, "A0");
                f_pi_phi_x_rho_rest->SetParameter(0, h_pi_phi_x_rho_rest->GetMaximum());
                for (int n = 1; n <= fitOrder; ++n) {
                    if (!only_even_cosine || (n % 2 == 0)) {
                        f_pi_phi_x_rho_rest->SetParName(n, Form("A%d", n));
                        f_pi_phi_x_rho_rest->SetParameter(n, 0.05);
                    }
                }
                h_pi_phi_x_rho_rest->Fit(f_pi_phi_x_rho_rest, "Q0");
                f_pi_phi_x_rho_rest->SetLineColor(kRed);
                f_pi_phi_x_rho_rest->SetLineWidth(2);
                f_pi_phi_x_rho_rest->Draw("SAME");
                double chi2_phi_x_rho_rest = f_pi_phi_x_rho_rest->GetChisquare();
                double ndf_phi_x_rho_rest = static_cast<double>(f_pi_phi_x_rho_rest->GetNDF());
                auto leg_phi_x_rho_rest = new TLegend(0.13, 0.58, 0.43, 0.88);
                leg_phi_x_rho_rest->SetTextSize(0.03);
                leg_phi_x_rho_rest->SetFillStyle(0);
                leg_phi_x_rho_rest->SetBorderSize(0);
                for (int n = 0; n <= fitOrder; ++n) {
                    if (!only_even_cosine || (n % 2 == 0)) {
                    leg_phi_x_rho_rest->AddEntry((TObject*)nullptr,
                                               Form("A%d = %.3f #pm %.3f", n, f_pi_phi_x_rho_rest->GetParameter(n), f_pi_phi_x_rho_rest->GetParError(n)),
                                               "");
                    }
                }
                leg_phi_x_rho_rest->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_x_rho_rest, (int)ndf_phi_x_rho_rest), "");
                leg_phi_x_rho_rest->Draw();
                if (generate_png) c_pi_phi_x_rho_rest->SaveAs("pi_phi_x_rho_rest.png");
                if (generate_pdf) c_pi_phi_x_rho_rest->SaveAs("pi_phi_x_rho_rest.pdf");
            }

            if (export_pion_angles_rho_x_ptcut && (generate_png || generate_pdf)) {
                // pT-filtered x-frame theta (best/other)
                TCanvas *c_pi_theta_x_rho_best_pt = new TCanvas("c_pi_theta_x_rho_best_pt", "Pion #theta (closest #rho, pair p_{T}<0.1)", 900, 700);
                auto g_pi_theta_x_rho_best_pt = histToGraph(h_pi_theta_x_rho_best_pt, 24, kBlue+2);
                g_pi_theta_x_rho_best_pt->SetMinimum(0.0);
                g_pi_theta_x_rho_best_pt->SetMaximum(h_pi_theta_x_rho_best_pt->GetMaximum() * 2.0);
                g_pi_theta_x_rho_best_pt->Draw("AP");
                TF1 *f_pi_theta_x_rho_best_pt = new TF1("f_pi_theta_x_rho_best_pt", "[0]*(1+[1]*pow(cos(x),2))", thetaMin, thetaMax);
                f_pi_theta_x_rho_best_pt->SetNpx(500);
                f_pi_theta_x_rho_best_pt->SetParName(0, "A");
                f_pi_theta_x_rho_best_pt->SetParName(1, "B");
                f_pi_theta_x_rho_best_pt->SetParameter(0, h_pi_theta_x_rho_best_pt->GetMaximum());
                f_pi_theta_x_rho_best_pt->SetParameter(1, 0.0);
                h_pi_theta_x_rho_best_pt->Fit(f_pi_theta_x_rho_best_pt, "Q0");
                f_pi_theta_x_rho_best_pt->SetLineColor(kRed);
                f_pi_theta_x_rho_best_pt->SetLineWidth(2);
                f_pi_theta_x_rho_best_pt->Draw("SAME");
                double chi2_theta_x_rho_best_pt = f_pi_theta_x_rho_best_pt->GetChisquare();
                double ndf_theta_x_rho_best_pt = static_cast<double>(f_pi_theta_x_rho_best_pt->GetNDF());
                auto leg_theta_x_rho_best_pt = new TLegend(0.13, 0.58, 0.43, 0.88);
                leg_theta_x_rho_best_pt->SetTextSize(0.03);
                leg_theta_x_rho_best_pt->SetFillStyle(0);
                leg_theta_x_rho_best_pt->SetBorderSize(0);
                leg_theta_x_rho_best_pt->AddEntry(g_pi_theta_x_rho_best_pt, "#pi^{+} #theta (closest #rho, pair p_{T}<0.1)", "p");
                leg_theta_x_rho_best_pt->AddEntry(f_pi_theta_x_rho_best_pt, "Fit: A(1+B cos^{2}#theta)", "l");
                leg_theta_x_rho_best_pt->AddEntry((TObject*)nullptr,
                                               Form("A = %.3f #pm %.3f", f_pi_theta_x_rho_best_pt->GetParameter(0), f_pi_theta_x_rho_best_pt->GetParError(0)),
                                               "");
                leg_theta_x_rho_best_pt->AddEntry((TObject*)nullptr,
                                               Form("B = %.3f #pm %.3f", f_pi_theta_x_rho_best_pt->GetParameter(1), f_pi_theta_x_rho_best_pt->GetParError(1)),
                                               "");
                leg_theta_x_rho_best_pt->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_theta_x_rho_best_pt, (int)ndf_theta_x_rho_best_pt), "");
                leg_theta_x_rho_best_pt->Draw();
                {
                    TLatex note;
                    note.SetNDC();
                    note.SetTextSize(0.03);
                    note.DrawLatex(0.60, 0.86, "pair p_{T}<0.1");
                }
                if (generate_png) c_pi_theta_x_rho_best_pt->SaveAs("pi_theta_x_rho_best_pt.png");
                if (generate_pdf) c_pi_theta_x_rho_best_pt->SaveAs("pi_theta_x_rho_best_pt.pdf");

                if (export_pion_angles_rho_x_ptcut_mid_high) {
                    // pT-filtered x-frame theta (best, 0.1<pT<0.2)
                    TCanvas *c_pi_theta_x_rho_best_pt_mid = new TCanvas("c_pi_theta_x_rho_best_pt_mid", "Pion #theta (closest #rho, 0.1<pair p_{T}<0.2)", 900, 700);
                    auto g_pi_theta_x_rho_best_pt_mid = histToGraph(h_pi_theta_x_rho_best_pt_mid, 24, kBlue+2);
                    g_pi_theta_x_rho_best_pt_mid->SetMinimum(0.0);
                    g_pi_theta_x_rho_best_pt_mid->SetMaximum(h_pi_theta_x_rho_best_pt_mid->GetMaximum() * 2.0);
                    g_pi_theta_x_rho_best_pt_mid->Draw("AP");
                    TF1 *f_pi_theta_x_rho_best_pt_mid = new TF1("f_pi_theta_x_rho_best_pt_mid", "[0]*(1+[1]*pow(cos(x),2))", thetaMin, thetaMax);
                    f_pi_theta_x_rho_best_pt_mid->SetNpx(500);
                    f_pi_theta_x_rho_best_pt_mid->SetParName(0, "A");
                    f_pi_theta_x_rho_best_pt_mid->SetParName(1, "B");
                    f_pi_theta_x_rho_best_pt_mid->SetParameter(0, h_pi_theta_x_rho_best_pt_mid->GetMaximum());
                    f_pi_theta_x_rho_best_pt_mid->SetParameter(1, 0.0);
                    h_pi_theta_x_rho_best_pt_mid->Fit(f_pi_theta_x_rho_best_pt_mid, "Q0");
                    f_pi_theta_x_rho_best_pt_mid->SetLineColor(kRed);
                    f_pi_theta_x_rho_best_pt_mid->SetLineWidth(2);
                    f_pi_theta_x_rho_best_pt_mid->Draw("SAME");
                    double chi2_theta_x_rho_best_pt_mid = f_pi_theta_x_rho_best_pt_mid->GetChisquare();
                    double ndf_theta_x_rho_best_pt_mid = static_cast<double>(f_pi_theta_x_rho_best_pt_mid->GetNDF());
                    auto leg_theta_x_rho_best_pt_mid = new TLegend(0.13, 0.58, 0.43, 0.88);
                    leg_theta_x_rho_best_pt_mid->SetTextSize(0.03);
                    leg_theta_x_rho_best_pt_mid->SetFillStyle(0);
                    leg_theta_x_rho_best_pt_mid->SetBorderSize(0);
                    leg_theta_x_rho_best_pt_mid->AddEntry(g_pi_theta_x_rho_best_pt_mid, "#pi^{+} #theta (closest #rho, 0.1<pair p_{T}<0.2)", "p");
                    leg_theta_x_rho_best_pt_mid->AddEntry(f_pi_theta_x_rho_best_pt_mid, "Fit: A(1+B cos^{2}#theta)", "l");
                    leg_theta_x_rho_best_pt_mid->AddEntry((TObject*)nullptr,
                                                         Form("A = %.3f #pm %.3f", f_pi_theta_x_rho_best_pt_mid->GetParameter(0), f_pi_theta_x_rho_best_pt_mid->GetParError(0)),
                                                         "");
                    leg_theta_x_rho_best_pt_mid->AddEntry((TObject*)nullptr,
                                                         Form("B = %.3f #pm %.3f", f_pi_theta_x_rho_best_pt_mid->GetParameter(1), f_pi_theta_x_rho_best_pt_mid->GetParError(1)),
                                                         "");
                    leg_theta_x_rho_best_pt_mid->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_theta_x_rho_best_pt_mid, (int)ndf_theta_x_rho_best_pt_mid), "");
                    leg_theta_x_rho_best_pt_mid->Draw();
                    {
                        TLatex note;
                        note.SetNDC();
                        note.SetTextSize(0.03);
                        note.DrawLatex(0.60, 0.86, "0.1<pair p_{T}<0.2");
                    }
                    if (generate_png) c_pi_theta_x_rho_best_pt_mid->SaveAs("pi_theta_x_rho_best_pt_0p1_0p2.png");
                    if (generate_pdf) c_pi_theta_x_rho_best_pt_mid->SaveAs("pi_theta_x_rho_best_pt_0p1_0p2.pdf");

                    // pT-filtered x-frame theta (best, pT>0.2)
                    TCanvas *c_pi_theta_x_rho_best_pt_high = new TCanvas("c_pi_theta_x_rho_best_pt_high", "Pion #theta (closest #rho, pair p_{T}>0.2)", 900, 700);
                    auto g_pi_theta_x_rho_best_pt_high = histToGraph(h_pi_theta_x_rho_best_pt_high, 24, kBlue+2);
                    g_pi_theta_x_rho_best_pt_high->SetMinimum(0.0);
                    g_pi_theta_x_rho_best_pt_high->SetMaximum(h_pi_theta_x_rho_best_pt_high->GetMaximum() * 2.0);
                    g_pi_theta_x_rho_best_pt_high->Draw("AP");
                    TF1 *f_pi_theta_x_rho_best_pt_high = new TF1("f_pi_theta_x_rho_best_pt_high", "[0]*(1+[1]*pow(cos(x),2))", thetaMin, thetaMax);
                    f_pi_theta_x_rho_best_pt_high->SetNpx(500);
                    f_pi_theta_x_rho_best_pt_high->SetParName(0, "A");
                    f_pi_theta_x_rho_best_pt_high->SetParName(1, "B");
                    f_pi_theta_x_rho_best_pt_high->SetParameter(0, h_pi_theta_x_rho_best_pt_high->GetMaximum());
                    f_pi_theta_x_rho_best_pt_high->SetParameter(1, 0.0);
                    h_pi_theta_x_rho_best_pt_high->Fit(f_pi_theta_x_rho_best_pt_high, "Q0");
                    f_pi_theta_x_rho_best_pt_high->SetLineColor(kRed);
                    f_pi_theta_x_rho_best_pt_high->SetLineWidth(2);
                    f_pi_theta_x_rho_best_pt_high->Draw("SAME");
                    double chi2_theta_x_rho_best_pt_high = f_pi_theta_x_rho_best_pt_high->GetChisquare();
                    double ndf_theta_x_rho_best_pt_high = static_cast<double>(f_pi_theta_x_rho_best_pt_high->GetNDF());
                    auto leg_theta_x_rho_best_pt_high = new TLegend(0.13, 0.58, 0.43, 0.88);
                    leg_theta_x_rho_best_pt_high->SetTextSize(0.03);
                    leg_theta_x_rho_best_pt_high->SetFillStyle(0);
                    leg_theta_x_rho_best_pt_high->SetBorderSize(0);
                    leg_theta_x_rho_best_pt_high->AddEntry(g_pi_theta_x_rho_best_pt_high, "#pi^{+} #theta (closest #rho, pair p_{T}>0.2)", "p");
                    leg_theta_x_rho_best_pt_high->AddEntry(f_pi_theta_x_rho_best_pt_high, "Fit: A(1+B cos^{2}#theta)", "l");
                    leg_theta_x_rho_best_pt_high->AddEntry((TObject*)nullptr,
                                                          Form("A = %.3f #pm %.3f", f_pi_theta_x_rho_best_pt_high->GetParameter(0), f_pi_theta_x_rho_best_pt_high->GetParError(0)),
                                                          "");
                    leg_theta_x_rho_best_pt_high->AddEntry((TObject*)nullptr,
                                                          Form("B = %.3f #pm %.3f", f_pi_theta_x_rho_best_pt_high->GetParameter(1), f_pi_theta_x_rho_best_pt_high->GetParError(1)),
                                                          "");
                    leg_theta_x_rho_best_pt_high->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_theta_x_rho_best_pt_high, (int)ndf_theta_x_rho_best_pt_high), "");
                    leg_theta_x_rho_best_pt_high->Draw();
                    {
                        TLatex note;
                        note.SetNDC();
                        note.SetTextSize(0.03);
                        note.DrawLatex(0.60, 0.86, "pair p_{T}>0.2");
                    }
                    if (generate_png) c_pi_theta_x_rho_best_pt_high->SaveAs("pi_theta_x_rho_best_pt_gt_0p2.png");
                    if (generate_pdf) c_pi_theta_x_rho_best_pt_high->SaveAs("pi_theta_x_rho_best_pt_gt_0p2.pdf");
                }

                TCanvas *c_pi_theta_x_rho_other_pt = new TCanvas("c_pi_theta_x_rho_other_pt", "Pion #theta (counterpart, pair p_{T}<0.1)", 900, 700);
                auto g_pi_theta_x_rho_other_pt = histToGraph(h_pi_theta_x_rho_other_pt, 24, kMagenta+2);
                g_pi_theta_x_rho_other_pt->SetMinimum(0.0);
                g_pi_theta_x_rho_other_pt->SetMaximum(h_pi_theta_x_rho_other_pt->GetMaximum() * 2.0);
                g_pi_theta_x_rho_other_pt->Draw("AP");
                TF1 *f_pi_theta_x_rho_other_pt = new TF1("f_pi_theta_x_rho_other_pt", "[0]*(1+[1]*pow(cos(x),2))", thetaMin, thetaMax);
                f_pi_theta_x_rho_other_pt->SetNpx(500);
                f_pi_theta_x_rho_other_pt->SetParName(0, "A");
                f_pi_theta_x_rho_other_pt->SetParName(1, "B");
                f_pi_theta_x_rho_other_pt->SetParameter(0, h_pi_theta_x_rho_other_pt->GetMaximum());
                f_pi_theta_x_rho_other_pt->SetParameter(1, 0.0);
                h_pi_theta_x_rho_other_pt->Fit(f_pi_theta_x_rho_other_pt, "Q0");
                f_pi_theta_x_rho_other_pt->SetLineColor(kRed);
                f_pi_theta_x_rho_other_pt->SetLineWidth(2);
                f_pi_theta_x_rho_other_pt->Draw("SAME");
                double chi2_theta_x_rho_other_pt = f_pi_theta_x_rho_other_pt->GetChisquare();
                double ndf_theta_x_rho_other_pt = static_cast<double>(f_pi_theta_x_rho_other_pt->GetNDF());
                auto leg_theta_x_rho_other_pt = new TLegend(0.13, 0.58, 0.43, 0.88);
                leg_theta_x_rho_other_pt->SetTextSize(0.03);
                leg_theta_x_rho_other_pt->SetFillStyle(0);
                leg_theta_x_rho_other_pt->SetBorderSize(0);
                leg_theta_x_rho_other_pt->AddEntry(g_pi_theta_x_rho_other_pt, "#pi^{+} #theta (counterpart, pair p_{T}<0.1)", "p");
                leg_theta_x_rho_other_pt->AddEntry(f_pi_theta_x_rho_other_pt, "Fit: A(1+B cos^{2}#theta)", "l");
                leg_theta_x_rho_other_pt->AddEntry((TObject*)nullptr,
                                                  Form("A = %.3f #pm %.3f", f_pi_theta_x_rho_other_pt->GetParameter(0), f_pi_theta_x_rho_other_pt->GetParError(0)),
                                                  "");
                leg_theta_x_rho_other_pt->AddEntry((TObject*)nullptr,
                                                  Form("B = %.3f #pm %.3f", f_pi_theta_x_rho_other_pt->GetParameter(1), f_pi_theta_x_rho_other_pt->GetParError(1)),
                                                  "");
                leg_theta_x_rho_other_pt->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_theta_x_rho_other_pt, (int)ndf_theta_x_rho_other_pt), "");
                leg_theta_x_rho_other_pt->Draw();
                {
                    TLatex note;
                    note.SetNDC();
                    note.SetTextSize(0.03);
                    note.DrawLatex(0.60, 0.86, "pair p_{T}<0.1");
                }
                if (generate_png) c_pi_theta_x_rho_other_pt->SaveAs("pi_theta_x_rho_other_pt.png");
                if (generate_pdf) c_pi_theta_x_rho_other_pt->SaveAs("pi_theta_x_rho_other_pt.pdf");

                if (export_pion_angles_rho_x_ptcut_mid_high) {
                    // pT-filtered x-frame theta (counterpart, 0.1<pT<0.2)
                    TCanvas *c_pi_theta_x_rho_other_pt_mid = new TCanvas("c_pi_theta_x_rho_other_pt_mid", "Pion #theta (counterpart, 0.1<pair p_{T}<0.2)", 900, 700);
                    auto g_pi_theta_x_rho_other_pt_mid = histToGraph(h_pi_theta_x_rho_other_pt_mid, 24, kMagenta+2);
                    g_pi_theta_x_rho_other_pt_mid->SetMinimum(0.0);
                    g_pi_theta_x_rho_other_pt_mid->SetMaximum(h_pi_theta_x_rho_other_pt_mid->GetMaximum() * 2.0);
                    g_pi_theta_x_rho_other_pt_mid->Draw("AP");
                    TF1 *f_pi_theta_x_rho_other_pt_mid = new TF1("f_pi_theta_x_rho_other_pt_mid", "[0]*(1+[1]*pow(cos(x),2))", thetaMin, thetaMax);
                    f_pi_theta_x_rho_other_pt_mid->SetNpx(500);
                    f_pi_theta_x_rho_other_pt_mid->SetParName(0, "A");
                    f_pi_theta_x_rho_other_pt_mid->SetParName(1, "B");
                    f_pi_theta_x_rho_other_pt_mid->SetParameter(0, h_pi_theta_x_rho_other_pt_mid->GetMaximum());
                    f_pi_theta_x_rho_other_pt_mid->SetParameter(1, 0.0);
                    h_pi_theta_x_rho_other_pt_mid->Fit(f_pi_theta_x_rho_other_pt_mid, "Q0");
                    f_pi_theta_x_rho_other_pt_mid->SetLineColor(kRed);
                    f_pi_theta_x_rho_other_pt_mid->SetLineWidth(2);
                    f_pi_theta_x_rho_other_pt_mid->Draw("SAME");
                    double chi2_theta_x_rho_other_pt_mid = f_pi_theta_x_rho_other_pt_mid->GetChisquare();
                    double ndf_theta_x_rho_other_pt_mid = static_cast<double>(f_pi_theta_x_rho_other_pt_mid->GetNDF());
                    auto leg_theta_x_rho_other_pt_mid = new TLegend(0.13, 0.58, 0.43, 0.88);
                    leg_theta_x_rho_other_pt_mid->SetTextSize(0.03);
                    leg_theta_x_rho_other_pt_mid->SetFillStyle(0);
                    leg_theta_x_rho_other_pt_mid->SetBorderSize(0);
                    leg_theta_x_rho_other_pt_mid->AddEntry(g_pi_theta_x_rho_other_pt_mid, "#pi^{+} #theta (counterpart, 0.1<pair p_{T}<0.2)", "p");
                    leg_theta_x_rho_other_pt_mid->AddEntry(f_pi_theta_x_rho_other_pt_mid, "Fit: A(1+B cos^{2}#theta)", "l");
                    leg_theta_x_rho_other_pt_mid->AddEntry((TObject*)nullptr,
                                                          Form("A = %.3f #pm %.3f", f_pi_theta_x_rho_other_pt_mid->GetParameter(0), f_pi_theta_x_rho_other_pt_mid->GetParError(0)),
                                                          "");
                    leg_theta_x_rho_other_pt_mid->AddEntry((TObject*)nullptr,
                                                          Form("B = %.3f #pm %.3f", f_pi_theta_x_rho_other_pt_mid->GetParameter(1), f_pi_theta_x_rho_other_pt_mid->GetParError(1)),
                                                          "");
                    leg_theta_x_rho_other_pt_mid->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_theta_x_rho_other_pt_mid, (int)ndf_theta_x_rho_other_pt_mid), "");
                    leg_theta_x_rho_other_pt_mid->Draw();
                    {
                        TLatex note;
                        note.SetNDC();
                        note.SetTextSize(0.03);
                        note.DrawLatex(0.60, 0.86, "0.1<pair p_{T}<0.2");
                    }
                    if (generate_png) c_pi_theta_x_rho_other_pt_mid->SaveAs("pi_theta_x_rho_other_pt_0p1_0p2.png");
                    if (generate_pdf) c_pi_theta_x_rho_other_pt_mid->SaveAs("pi_theta_x_rho_other_pt_0p1_0p2.pdf");

                    // pT-filtered x-frame theta (counterpart, pT>0.2)
                    TCanvas *c_pi_theta_x_rho_other_pt_high = new TCanvas("c_pi_theta_x_rho_other_pt_high", "Pion #theta (counterpart, pair p_{T}>0.2)", 900, 700);
                    auto g_pi_theta_x_rho_other_pt_high = histToGraph(h_pi_theta_x_rho_other_pt_high, 24, kMagenta+2);
                    g_pi_theta_x_rho_other_pt_high->SetMinimum(0.0);
                    g_pi_theta_x_rho_other_pt_high->SetMaximum(h_pi_theta_x_rho_other_pt_high->GetMaximum() * 2.0);
                    g_pi_theta_x_rho_other_pt_high->Draw("AP");
                    TF1 *f_pi_theta_x_rho_other_pt_high = new TF1("f_pi_theta_x_rho_other_pt_high", "[0]*(1+[1]*pow(cos(x),2))", thetaMin, thetaMax);
                    f_pi_theta_x_rho_other_pt_high->SetNpx(500);
                    f_pi_theta_x_rho_other_pt_high->SetParName(0, "A");
                    f_pi_theta_x_rho_other_pt_high->SetParName(1, "B");
                    f_pi_theta_x_rho_other_pt_high->SetParameter(0, h_pi_theta_x_rho_other_pt_high->GetMaximum());
                    f_pi_theta_x_rho_other_pt_high->SetParameter(1, 0.0);
                    h_pi_theta_x_rho_other_pt_high->Fit(f_pi_theta_x_rho_other_pt_high, "Q0");
                    f_pi_theta_x_rho_other_pt_high->SetLineColor(kRed);
                    f_pi_theta_x_rho_other_pt_high->SetLineWidth(2);
                    f_pi_theta_x_rho_other_pt_high->Draw("SAME");
                    double chi2_theta_x_rho_other_pt_high = f_pi_theta_x_rho_other_pt_high->GetChisquare();
                    double ndf_theta_x_rho_other_pt_high = static_cast<double>(f_pi_theta_x_rho_other_pt_high->GetNDF());
                    auto leg_theta_x_rho_other_pt_high = new TLegend(0.13, 0.58, 0.43, 0.88);
                    leg_theta_x_rho_other_pt_high->SetTextSize(0.03);
                    leg_theta_x_rho_other_pt_high->SetFillStyle(0);
                    leg_theta_x_rho_other_pt_high->SetBorderSize(0);
                    leg_theta_x_rho_other_pt_high->AddEntry(g_pi_theta_x_rho_other_pt_high, "#pi^{+} #theta (counterpart, pair p_{T}>0.2)", "p");
                    leg_theta_x_rho_other_pt_high->AddEntry(f_pi_theta_x_rho_other_pt_high, "Fit: A(1+B cos^{2}#theta)", "l");
                    leg_theta_x_rho_other_pt_high->AddEntry((TObject*)nullptr,
                                                           Form("A = %.3f #pm %.3f", f_pi_theta_x_rho_other_pt_high->GetParameter(0), f_pi_theta_x_rho_other_pt_high->GetParError(0)),
                                                           "");
                    leg_theta_x_rho_other_pt_high->AddEntry((TObject*)nullptr,
                                                           Form("B = %.3f #pm %.3f", f_pi_theta_x_rho_other_pt_high->GetParameter(1), f_pi_theta_x_rho_other_pt_high->GetParError(1)),
                                                           "");
                    leg_theta_x_rho_other_pt_high->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_theta_x_rho_other_pt_high, (int)ndf_theta_x_rho_other_pt_high), "");
                    leg_theta_x_rho_other_pt_high->Draw();
                    {
                        TLatex note;
                        note.SetNDC();
                        note.SetTextSize(0.03);
                        note.DrawLatex(0.60, 0.86, "pair p_{T}>0.2");
                    }
                    if (generate_png) c_pi_theta_x_rho_other_pt_high->SaveAs("pi_theta_x_rho_other_pt_gt_0p2.png");
                    if (generate_pdf) c_pi_theta_x_rho_other_pt_high->SaveAs("pi_theta_x_rho_other_pt_gt_0p2.pdf");
                }

                if (export_pion_angles_rho_x_rest) {
                    TCanvas *c_pi_theta_x_rho_rest_pt = new TCanvas("c_pi_theta_x_rho_rest_pt", "Pion #theta (remaining, pair p_{T}<0.1)", 900, 700);
                    auto g_pi_theta_x_rho_rest_pt = histToGraph(h_pi_theta_x_rho_rest_pt, 24, kOrange+7);
                    g_pi_theta_x_rho_rest_pt->SetMinimum(0.0);
                    g_pi_theta_x_rho_rest_pt->SetMaximum(h_pi_theta_x_rho_rest_pt->GetMaximum() * 2.0);
                    g_pi_theta_x_rho_rest_pt->Draw("AP");
                    TF1 *f_pi_theta_x_rho_rest_pt = new TF1("f_pi_theta_x_rho_rest_pt", "[0]*(1+[1]*pow(cos(x),2))", thetaMin, thetaMax);
                    f_pi_theta_x_rho_rest_pt->SetNpx(500);
                    f_pi_theta_x_rho_rest_pt->SetParName(0, "A");
                    f_pi_theta_x_rho_rest_pt->SetParName(1, "B");
                    f_pi_theta_x_rho_rest_pt->SetParameter(0, h_pi_theta_x_rho_rest_pt->GetMaximum());
                    f_pi_theta_x_rho_rest_pt->SetParameter(1, 0.0);
                    h_pi_theta_x_rho_rest_pt->Fit(f_pi_theta_x_rho_rest_pt, "Q0");
                    f_pi_theta_x_rho_rest_pt->SetLineColor(kRed);
                    f_pi_theta_x_rho_rest_pt->SetLineWidth(2);
                    f_pi_theta_x_rho_rest_pt->Draw("SAME");
                    double chi2_theta_x_rho_rest_pt = f_pi_theta_x_rho_rest_pt->GetChisquare();
                    double ndf_theta_x_rho_rest_pt = static_cast<double>(f_pi_theta_x_rho_rest_pt->GetNDF());
                    auto leg_theta_x_rho_rest_pt = new TLegend(0.13, 0.58, 0.43, 0.88);
                    leg_theta_x_rho_rest_pt->SetTextSize(0.03);
                    leg_theta_x_rho_rest_pt->SetFillStyle(0);
                    leg_theta_x_rho_rest_pt->SetBorderSize(0);
                    leg_theta_x_rho_rest_pt->AddEntry(g_pi_theta_x_rho_rest_pt, "#pi^{+} #theta (remaining, pair p_{T}<0.1)", "p");
                    leg_theta_x_rho_rest_pt->AddEntry(f_pi_theta_x_rho_rest_pt, "Fit: A(1+B cos^{2}#theta)", "l");
                    leg_theta_x_rho_rest_pt->AddEntry((TObject*)nullptr,
                                                      Form("A = %.3f #pm %.3f", f_pi_theta_x_rho_rest_pt->GetParameter(0), f_pi_theta_x_rho_rest_pt->GetParError(0)),
                                                      "");
                    leg_theta_x_rho_rest_pt->AddEntry((TObject*)nullptr,
                                                      Form("B = %.3f #pm %.3f", f_pi_theta_x_rho_rest_pt->GetParameter(1), f_pi_theta_x_rho_rest_pt->GetParError(1)),
                                                      "");
                    leg_theta_x_rho_rest_pt->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_theta_x_rho_rest_pt, (int)ndf_theta_x_rho_rest_pt), "");
                    leg_theta_x_rho_rest_pt->Draw();
                    {
                        TLatex note;
                        note.SetNDC();
                        note.SetTextSize(0.03);
                        note.DrawLatex(0.60, 0.86, "pair p_{T}<0.1");
                    }
                    if (generate_png) c_pi_theta_x_rho_rest_pt->SaveAs("pi_theta_x_rho_rest_pt.png");
                    if (generate_pdf) c_pi_theta_x_rho_rest_pt->SaveAs("pi_theta_x_rho_rest_pt.pdf");
                }

                // pT-filtered x-frame phi (best/other)
                TCanvas *c_pi_phi_x_rho_best_pt = new TCanvas("c_pi_phi_x_rho_best_pt", "", 900, 700);
                auto g_pi_phi_x_rho_best_pt = histToGraph(h_pi_phi_x_rho_best_pt, 24, kBlue+2);
                g_pi_phi_x_rho_best_pt->SetMinimum(0.0);
                g_pi_phi_x_rho_best_pt->SetMaximum(h_pi_phi_x_rho_best_pt->GetMaximum() * 2.0);
                g_pi_phi_x_rho_best_pt->SetTitle(";#Delta#varphi;Counts");
                g_pi_phi_x_rho_best_pt->Draw("AP");
                TF1 *f_pi_phi_x_rho_best_pt = new TF1("f_pi_phi_x_rho_best_pt", cosExpr.c_str(), phiMin, phiMax);
                f_pi_phi_x_rho_best_pt->SetNpx(500);
                f_pi_phi_x_rho_best_pt->SetParName(0, "A0");
                f_pi_phi_x_rho_best_pt->SetParameter(0, h_pi_phi_x_rho_best_pt->GetMaximum());
                for (int n = 1; n <= fitOrder; ++n) {
                    if (!only_even_cosine || (n % 2 == 0)) {
                        f_pi_phi_x_rho_best_pt->SetParName(n, Form("A%d", n));
                        f_pi_phi_x_rho_best_pt->SetParameter(n, 0.05);
                    }
                }
                h_pi_phi_x_rho_best_pt->Fit(f_pi_phi_x_rho_best_pt, "Q0");
                f_pi_phi_x_rho_best_pt->SetLineColor(kRed);
                f_pi_phi_x_rho_best_pt->SetLineWidth(2);
                f_pi_phi_x_rho_best_pt->Draw("SAME");
                double chi2_phi_x_rho_best_pt = f_pi_phi_x_rho_best_pt->GetChisquare();
                double ndf_phi_x_rho_best_pt = static_cast<double>(f_pi_phi_x_rho_best_pt->GetNDF());
                auto leg_phi_x_rho_best_pt = new TLegend(0.13, 0.58, 0.43, 0.88);
                leg_phi_x_rho_best_pt->SetTextSize(0.03);
                leg_phi_x_rho_best_pt->SetFillStyle(0);
                leg_phi_x_rho_best_pt->SetBorderSize(0);
                for (int n = 0; n <= fitOrder; ++n) {
                    if (!only_even_cosine || (n % 2 == 0)) {
                        leg_phi_x_rho_best_pt->AddEntry((TObject*)nullptr,
                                                       Form("A%d = %.3f #pm %.3f", n, f_pi_phi_x_rho_best_pt->GetParameter(n), f_pi_phi_x_rho_best_pt->GetParError(n)),
                                                       "");
                    }
                }
                leg_phi_x_rho_best_pt->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_x_rho_best_pt, (int)ndf_phi_x_rho_best_pt), "");
                leg_phi_x_rho_best_pt->Draw();
                {
                    TLatex note;
                    note.SetNDC();
                    note.SetTextSize(0.03);
                    note.DrawLatex(0.60, 0.86, "pair p_{T}<0.1");
                }
                if (generate_png) c_pi_phi_x_rho_best_pt->SaveAs("pi_phi_x_rho_best_pt.png");
                if (generate_pdf) c_pi_phi_x_rho_best_pt->SaveAs("pi_phi_x_rho_best_pt.pdf");

                TCanvas *c_pi_phi_x_rho_other_pt = new TCanvas("c_pi_phi_x_rho_other_pt", "", 900, 700);
                auto g_pi_phi_x_rho_other_pt = histToGraph(h_pi_phi_x_rho_other_pt, 24, kMagenta+2);
                g_pi_phi_x_rho_other_pt->SetMinimum(0.0);
                g_pi_phi_x_rho_other_pt->SetMaximum(h_pi_phi_x_rho_other_pt->GetMaximum() * 2.0);
                g_pi_phi_x_rho_other_pt->SetTitle(";#Delta#varphi;Counts");
                g_pi_phi_x_rho_other_pt->Draw("AP");
                TF1 *f_pi_phi_x_rho_other_pt = new TF1("f_pi_phi_x_rho_other_pt", cosExpr.c_str(), phiMin, phiMax);
                f_pi_phi_x_rho_other_pt->SetNpx(500);
                f_pi_phi_x_rho_other_pt->SetParName(0, "A0");
                f_pi_phi_x_rho_other_pt->SetParameter(0, h_pi_phi_x_rho_other_pt->GetMaximum());
                for (int n = 1; n <= fitOrder; ++n) {
                    if (!only_even_cosine || (n % 2 == 0)) {
                        f_pi_phi_x_rho_other_pt->SetParName(n, Form("A%d", n));
                        f_pi_phi_x_rho_other_pt->SetParameter(n, 0.05);
                    }
                }
                h_pi_phi_x_rho_other_pt->Fit(f_pi_phi_x_rho_other_pt, "Q0");
                f_pi_phi_x_rho_other_pt->SetLineColor(kRed);
                f_pi_phi_x_rho_other_pt->SetLineWidth(2);
                f_pi_phi_x_rho_other_pt->Draw("SAME");
                double chi2_phi_x_rho_other_pt = f_pi_phi_x_rho_other_pt->GetChisquare();
                double ndf_phi_x_rho_other_pt = static_cast<double>(f_pi_phi_x_rho_other_pt->GetNDF());
                auto leg_phi_x_rho_other_pt = new TLegend(0.13, 0.58, 0.43, 0.88);
                leg_phi_x_rho_other_pt->SetTextSize(0.03);
                leg_phi_x_rho_other_pt->SetFillStyle(0);
                leg_phi_x_rho_other_pt->SetBorderSize(0);
                for (int n = 0; n <= fitOrder; ++n) {
                    if (!only_even_cosine || (n % 2 == 0)) {
                        leg_phi_x_rho_other_pt->AddEntry((TObject*)nullptr,
                                                        Form("A%d = %.3f #pm %.3f", n, f_pi_phi_x_rho_other_pt->GetParameter(n), f_pi_phi_x_rho_other_pt->GetParError(n)),
                                                        "");
                    }
                }
                leg_phi_x_rho_other_pt->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_x_rho_other_pt, (int)ndf_phi_x_rho_other_pt), "");
                leg_phi_x_rho_other_pt->Draw();
                {
                    TLatex note;
                    note.SetNDC();
                    note.SetTextSize(0.03);
                    note.DrawLatex(0.60, 0.86, "pair p_{T}<0.1");
                }
                if (generate_png) c_pi_phi_x_rho_other_pt->SaveAs("pi_phi_x_rho_other_pt.png");
                if (generate_pdf) c_pi_phi_x_rho_other_pt->SaveAs("pi_phi_x_rho_other_pt.pdf");

                if (export_pion_angles_rho_x_ptcut_mid_high) {
                    // pT-filtered x-frame phi (best, 0.1<pT<0.2)
                    TCanvas *c_pi_phi_x_rho_best_pt_mid = new TCanvas("c_pi_phi_x_rho_best_pt_mid", "", 900, 700);
                    auto g_pi_phi_x_rho_best_pt_mid = histToGraph(h_pi_phi_x_rho_best_pt_mid, 24, kBlue+2);
                    g_pi_phi_x_rho_best_pt_mid->SetMinimum(0.0);
                    g_pi_phi_x_rho_best_pt_mid->SetMaximum(h_pi_phi_x_rho_best_pt_mid->GetMaximum() * 2.0);
                    g_pi_phi_x_rho_best_pt_mid->SetTitle(";#Delta#varphi;Counts");
                    g_pi_phi_x_rho_best_pt_mid->Draw("AP");
                    TF1 *f_pi_phi_x_rho_best_pt_mid = new TF1("f_pi_phi_x_rho_best_pt_mid", cosExpr.c_str(), phiMin, phiMax);
                    f_pi_phi_x_rho_best_pt_mid->SetNpx(500);
                    f_pi_phi_x_rho_best_pt_mid->SetParName(0, "A0");
                    f_pi_phi_x_rho_best_pt_mid->SetParameter(0, h_pi_phi_x_rho_best_pt_mid->GetMaximum());
                    for (int n = 1; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            f_pi_phi_x_rho_best_pt_mid->SetParName(n, Form("A%d", n));
                            f_pi_phi_x_rho_best_pt_mid->SetParameter(n, 0.05);
                        }
                    }
                    h_pi_phi_x_rho_best_pt_mid->Fit(f_pi_phi_x_rho_best_pt_mid, "Q0");
                    f_pi_phi_x_rho_best_pt_mid->SetLineColor(kRed);
                    f_pi_phi_x_rho_best_pt_mid->SetLineWidth(2);
                    f_pi_phi_x_rho_best_pt_mid->Draw("SAME");
                    double chi2_phi_x_rho_best_pt_mid = f_pi_phi_x_rho_best_pt_mid->GetChisquare();
                    double ndf_phi_x_rho_best_pt_mid = static_cast<double>(f_pi_phi_x_rho_best_pt_mid->GetNDF());
                    auto leg_phi_x_rho_best_pt_mid = new TLegend(0.13, 0.58, 0.43, 0.88);
                    leg_phi_x_rho_best_pt_mid->SetTextSize(0.03);
                    leg_phi_x_rho_best_pt_mid->SetFillStyle(0);
                    leg_phi_x_rho_best_pt_mid->SetBorderSize(0);
                    for (int n = 0; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            leg_phi_x_rho_best_pt_mid->AddEntry((TObject*)nullptr,
                                                               Form("A%d = %.3f #pm %.3f", n, f_pi_phi_x_rho_best_pt_mid->GetParameter(n), f_pi_phi_x_rho_best_pt_mid->GetParError(n)),
                                                               "");
                        }
                    }
                    leg_phi_x_rho_best_pt_mid->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_x_rho_best_pt_mid, (int)ndf_phi_x_rho_best_pt_mid), "");
                    leg_phi_x_rho_best_pt_mid->Draw();
                    {
                        TLatex note;
                        note.SetNDC();
                        note.SetTextSize(0.03);
                        note.DrawLatex(0.60, 0.86, "0.1<pair p_{T}<0.2");
                    }
                    if (generate_png) c_pi_phi_x_rho_best_pt_mid->SaveAs("pi_phi_x_rho_best_pt_0p1_0p2.png");
                    if (generate_pdf) c_pi_phi_x_rho_best_pt_mid->SaveAs("pi_phi_x_rho_best_pt_0p1_0p2.pdf");

                    // pT-filtered x-frame phi (best, pT>0.2)
                    TCanvas *c_pi_phi_x_rho_best_pt_high = new TCanvas("c_pi_phi_x_rho_best_pt_high", "", 900, 700);
                    auto g_pi_phi_x_rho_best_pt_high = histToGraph(h_pi_phi_x_rho_best_pt_high, 24, kBlue+2);
                    g_pi_phi_x_rho_best_pt_high->SetMinimum(0.0);
                    g_pi_phi_x_rho_best_pt_high->SetMaximum(h_pi_phi_x_rho_best_pt_high->GetMaximum() * 2.0);
                    g_pi_phi_x_rho_best_pt_high->SetTitle(";#Delta#varphi;Counts");
                    g_pi_phi_x_rho_best_pt_high->Draw("AP");
                    TF1 *f_pi_phi_x_rho_best_pt_high = new TF1("f_pi_phi_x_rho_best_pt_high", cosExpr.c_str(), phiMin, phiMax);
                    f_pi_phi_x_rho_best_pt_high->SetNpx(500);
                    f_pi_phi_x_rho_best_pt_high->SetParName(0, "A0");
                    f_pi_phi_x_rho_best_pt_high->SetParameter(0, h_pi_phi_x_rho_best_pt_high->GetMaximum());
                    for (int n = 1; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            f_pi_phi_x_rho_best_pt_high->SetParName(n, Form("A%d", n));
                            f_pi_phi_x_rho_best_pt_high->SetParameter(n, 0.05);
                        }
                    }
                    h_pi_phi_x_rho_best_pt_high->Fit(f_pi_phi_x_rho_best_pt_high, "Q0");
                    f_pi_phi_x_rho_best_pt_high->SetLineColor(kRed);
                    f_pi_phi_x_rho_best_pt_high->SetLineWidth(2);
                    f_pi_phi_x_rho_best_pt_high->Draw("SAME");
                    double chi2_phi_x_rho_best_pt_high = f_pi_phi_x_rho_best_pt_high->GetChisquare();
                    double ndf_phi_x_rho_best_pt_high = static_cast<double>(f_pi_phi_x_rho_best_pt_high->GetNDF());
                    auto leg_phi_x_rho_best_pt_high = new TLegend(0.13, 0.58, 0.43, 0.88);
                    leg_phi_x_rho_best_pt_high->SetTextSize(0.03);
                    leg_phi_x_rho_best_pt_high->SetFillStyle(0);
                    leg_phi_x_rho_best_pt_high->SetBorderSize(0);
                    for (int n = 0; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            leg_phi_x_rho_best_pt_high->AddEntry((TObject*)nullptr,
                                                                Form("A%d = %.3f #pm %.3f", n, f_pi_phi_x_rho_best_pt_high->GetParameter(n), f_pi_phi_x_rho_best_pt_high->GetParError(n)),
                                                                "");
                        }
                    }
                    leg_phi_x_rho_best_pt_high->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_x_rho_best_pt_high, (int)ndf_phi_x_rho_best_pt_high), "");
                    leg_phi_x_rho_best_pt_high->Draw();
                    {
                        TLatex note;
                        note.SetNDC();
                        note.SetTextSize(0.03);
                        note.DrawLatex(0.60, 0.86, "pair p_{T}>0.2");
                    }
                    if (generate_png) c_pi_phi_x_rho_best_pt_high->SaveAs("pi_phi_x_rho_best_pt_gt_0p2.png");
                    if (generate_pdf) c_pi_phi_x_rho_best_pt_high->SaveAs("pi_phi_x_rho_best_pt_gt_0p2.pdf");

                    // pT-filtered x-frame phi (counterpart, 0.1<pT<0.2)
                    TCanvas *c_pi_phi_x_rho_other_pt_mid = new TCanvas("c_pi_phi_x_rho_other_pt_mid", "", 900, 700);
                    auto g_pi_phi_x_rho_other_pt_mid = histToGraph(h_pi_phi_x_rho_other_pt_mid, 24, kMagenta+2);
                    g_pi_phi_x_rho_other_pt_mid->SetMinimum(0.0);
                    g_pi_phi_x_rho_other_pt_mid->SetMaximum(h_pi_phi_x_rho_other_pt_mid->GetMaximum() * 2.0);
                    g_pi_phi_x_rho_other_pt_mid->SetTitle(";#Delta#varphi;Counts");
                    g_pi_phi_x_rho_other_pt_mid->Draw("AP");
                    TF1 *f_pi_phi_x_rho_other_pt_mid = new TF1("f_pi_phi_x_rho_other_pt_mid", cosExpr.c_str(), phiMin, phiMax);
                    f_pi_phi_x_rho_other_pt_mid->SetNpx(500);
                    f_pi_phi_x_rho_other_pt_mid->SetParName(0, "A0");
                    f_pi_phi_x_rho_other_pt_mid->SetParameter(0, h_pi_phi_x_rho_other_pt_mid->GetMaximum());
                    for (int n = 1; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            f_pi_phi_x_rho_other_pt_mid->SetParName(n, Form("A%d", n));
                            f_pi_phi_x_rho_other_pt_mid->SetParameter(n, 0.05);
                        }
                    }
                    h_pi_phi_x_rho_other_pt_mid->Fit(f_pi_phi_x_rho_other_pt_mid, "Q0");
                    f_pi_phi_x_rho_other_pt_mid->SetLineColor(kRed);
                    f_pi_phi_x_rho_other_pt_mid->SetLineWidth(2);
                    f_pi_phi_x_rho_other_pt_mid->Draw("SAME");
                    double chi2_phi_x_rho_other_pt_mid = f_pi_phi_x_rho_other_pt_mid->GetChisquare();
                    double ndf_phi_x_rho_other_pt_mid = static_cast<double>(f_pi_phi_x_rho_other_pt_mid->GetNDF());
                    auto leg_phi_x_rho_other_pt_mid = new TLegend(0.13, 0.58, 0.43, 0.88);
                    leg_phi_x_rho_other_pt_mid->SetTextSize(0.03);
                    leg_phi_x_rho_other_pt_mid->SetFillStyle(0);
                    leg_phi_x_rho_other_pt_mid->SetBorderSize(0);
                    for (int n = 0; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            leg_phi_x_rho_other_pt_mid->AddEntry((TObject*)nullptr,
                                                                Form("A%d = %.3f #pm %.3f", n, f_pi_phi_x_rho_other_pt_mid->GetParameter(n), f_pi_phi_x_rho_other_pt_mid->GetParError(n)),
                                                                "");
                        }
                    }
                    leg_phi_x_rho_other_pt_mid->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_x_rho_other_pt_mid, (int)ndf_phi_x_rho_other_pt_mid), "");
                    leg_phi_x_rho_other_pt_mid->Draw();
                    {
                        TLatex note;
                        note.SetNDC();
                        note.SetTextSize(0.03);
                        note.DrawLatex(0.60, 0.86, "0.1<pair p_{T}<0.2");
                    }
                    if (generate_png) c_pi_phi_x_rho_other_pt_mid->SaveAs("pi_phi_x_rho_other_pt_0p1_0p2.png");
                    if (generate_pdf) c_pi_phi_x_rho_other_pt_mid->SaveAs("pi_phi_x_rho_other_pt_0p1_0p2.pdf");

                    // pT-filtered x-frame phi (counterpart, pT>0.2)
                    TCanvas *c_pi_phi_x_rho_other_pt_high = new TCanvas("c_pi_phi_x_rho_other_pt_high", "", 900, 700);
                    auto g_pi_phi_x_rho_other_pt_high = histToGraph(h_pi_phi_x_rho_other_pt_high, 24, kMagenta+2);
                    g_pi_phi_x_rho_other_pt_high->SetMinimum(0.0);
                    g_pi_phi_x_rho_other_pt_high->SetMaximum(h_pi_phi_x_rho_other_pt_high->GetMaximum() * 2.0);
                    g_pi_phi_x_rho_other_pt_high->SetTitle(";#Delta#varphi;Counts");
                    g_pi_phi_x_rho_other_pt_high->Draw("AP");
                    TF1 *f_pi_phi_x_rho_other_pt_high = new TF1("f_pi_phi_x_rho_other_pt_high", cosExpr.c_str(), phiMin, phiMax);
                    f_pi_phi_x_rho_other_pt_high->SetNpx(500);
                    f_pi_phi_x_rho_other_pt_high->SetParName(0, "A0");
                    f_pi_phi_x_rho_other_pt_high->SetParameter(0, h_pi_phi_x_rho_other_pt_high->GetMaximum());
                    for (int n = 1; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            f_pi_phi_x_rho_other_pt_high->SetParName(n, Form("A%d", n));
                            f_pi_phi_x_rho_other_pt_high->SetParameter(n, 0.05);
                        }
                    }
                    h_pi_phi_x_rho_other_pt_high->Fit(f_pi_phi_x_rho_other_pt_high, "Q0");
                    f_pi_phi_x_rho_other_pt_high->SetLineColor(kRed);
                    f_pi_phi_x_rho_other_pt_high->SetLineWidth(2);
                    f_pi_phi_x_rho_other_pt_high->Draw("SAME");
                    double chi2_phi_x_rho_other_pt_high = f_pi_phi_x_rho_other_pt_high->GetChisquare();
                    double ndf_phi_x_rho_other_pt_high = static_cast<double>(f_pi_phi_x_rho_other_pt_high->GetNDF());
                    auto leg_phi_x_rho_other_pt_high = new TLegend(0.13, 0.58, 0.43, 0.88);
                    leg_phi_x_rho_other_pt_high->SetTextSize(0.03);
                    leg_phi_x_rho_other_pt_high->SetFillStyle(0);
                    leg_phi_x_rho_other_pt_high->SetBorderSize(0);
                    for (int n = 0; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            leg_phi_x_rho_other_pt_high->AddEntry((TObject*)nullptr,
                                                                 Form("A%d = %.3f #pm %.3f", n, f_pi_phi_x_rho_other_pt_high->GetParameter(n), f_pi_phi_x_rho_other_pt_high->GetParError(n)),
                                                                 "");
                        }
                    }
                    leg_phi_x_rho_other_pt_high->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_x_rho_other_pt_high, (int)ndf_phi_x_rho_other_pt_high), "");
                    leg_phi_x_rho_other_pt_high->Draw();
                    {
                        TLatex note;
                        note.SetNDC();
                        note.SetTextSize(0.03);
                        note.DrawLatex(0.60, 0.86, "pair p_{T}>0.2");
                    }
                    if (generate_png) c_pi_phi_x_rho_other_pt_high->SaveAs("pi_phi_x_rho_other_pt_gt_0p2.png");
                    if (generate_pdf) c_pi_phi_x_rho_other_pt_high->SaveAs("pi_phi_x_rho_other_pt_gt_0p2.pdf");
                }

                if (export_pion_angles_rho_x_rest) {
                    TCanvas *c_pi_phi_x_rho_rest_pt = new TCanvas("c_pi_phi_x_rho_rest_pt", "", 900, 700);
                    auto g_pi_phi_x_rho_rest_pt = histToGraph(h_pi_phi_x_rho_rest_pt, 24, kOrange+7);
                    g_pi_phi_x_rho_rest_pt->SetMinimum(0.0);
                    g_pi_phi_x_rho_rest_pt->SetMaximum(h_pi_phi_x_rho_rest_pt->GetMaximum() * 2.0);
                    g_pi_phi_x_rho_rest_pt->SetTitle(";#Delta#varphi;Counts");
                    g_pi_phi_x_rho_rest_pt->Draw("AP");
                    TF1 *f_pi_phi_x_rho_rest_pt = new TF1("f_pi_phi_x_rho_rest_pt", cosExpr.c_str(), phiMin, phiMax);
                    f_pi_phi_x_rho_rest_pt->SetNpx(500);
                    f_pi_phi_x_rho_rest_pt->SetParName(0, "A0");
                    f_pi_phi_x_rho_rest_pt->SetParameter(0, h_pi_phi_x_rho_rest_pt->GetMaximum());
                    for (int n = 1; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            f_pi_phi_x_rho_rest_pt->SetParName(n, Form("A%d", n));
                            f_pi_phi_x_rho_rest_pt->SetParameter(n, 0.05);
                        }
                    }
                    h_pi_phi_x_rho_rest_pt->Fit(f_pi_phi_x_rho_rest_pt, "Q0");
                    f_pi_phi_x_rho_rest_pt->SetLineColor(kRed);
                    f_pi_phi_x_rho_rest_pt->SetLineWidth(2);
                    f_pi_phi_x_rho_rest_pt->Draw("SAME");
                    double chi2_phi_x_rho_rest_pt = f_pi_phi_x_rho_rest_pt->GetChisquare();
                    double ndf_phi_x_rho_rest_pt = static_cast<double>(f_pi_phi_x_rho_rest_pt->GetNDF());
                    auto leg_phi_x_rho_rest_pt = new TLegend(0.13, 0.58, 0.43, 0.88);
                    leg_phi_x_rho_rest_pt->SetTextSize(0.03);
                    leg_phi_x_rho_rest_pt->SetFillStyle(0);
                    leg_phi_x_rho_rest_pt->SetBorderSize(0);
                    for (int n = 0; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            leg_phi_x_rho_rest_pt->AddEntry((TObject*)nullptr,
                                                           Form("A%d = %.3f #pm %.3f", n, f_pi_phi_x_rho_rest_pt->GetParameter(n), f_pi_phi_x_rho_rest_pt->GetParError(n)),
                                                           "");
                        }
                    }
                    leg_phi_x_rho_rest_pt->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_x_rho_rest_pt, (int)ndf_phi_x_rho_rest_pt), "");
                    leg_phi_x_rho_rest_pt->Draw();
                    {
                        TLatex note;
                        note.SetNDC();
                        note.SetTextSize(0.03);
                        note.DrawLatex(0.60, 0.86, "pair p_{T}<0.1");
                    }
                    if (generate_png) c_pi_phi_x_rho_rest_pt->SaveAs("pi_phi_x_rho_rest_pt.png");
                    if (generate_pdf) c_pi_phi_x_rho_rest_pt->SaveAs("pi_phi_x_rho_rest_pt.pdf");
                }
            }
        }

        if (export_pion_angles_ab && (generate_png || generate_pdf)) {
            // A/B rule filtered pi theta with fit
            TCanvas *c_pi_theta_ab_pts = new TCanvas("c_pi_theta_ab_points", "Pion theta sum (A/B rule)", 900, 700);
            auto g_pi_theta_ab = histToGraph(h_pi_theta_ab_sum, 24, kBlue+2);
            g_pi_theta_ab->Draw("AP");
            TF1 *f_pi_theta_ab = new TF1("f_pi_theta_ab", "[0]*(1+[1]*pow(cos(x),2))", thetaMin, thetaMax);
            f_pi_theta_ab->SetNpx(500);
            f_pi_theta_ab->SetParName(0, "A");
            f_pi_theta_ab->SetParName(1, "B");
            f_pi_theta_ab->SetParameter(0, h_pi_theta_ab_sum->GetMaximum());
            f_pi_theta_ab->SetParameter(1, 0.0);
            h_pi_theta_ab_sum->Fit(f_pi_theta_ab, "Q0");
            f_pi_theta_ab->SetLineColor(kRed);
            f_pi_theta_ab->SetLineWidth(2);
            f_pi_theta_ab->Draw("SAME");
            double chi2_theta_ab = f_pi_theta_ab->GetChisquare();
            double ndf_theta_ab = static_cast<double>(f_pi_theta_ab->GetNDF());
            auto leg_theta_ab = new TLegend(0.12, 0.60, 0.58, 0.86);
            leg_theta_ab->SetTextSize(0.02);
            leg_theta_ab->SetFillStyle(0);
            leg_theta_ab->SetBorderSize(0);
            leg_theta_ab->AddEntry(g_pi_theta_ab, "#pi^{+} #theta (A/B rule)", "p");
            leg_theta_ab->AddEntry(f_pi_theta_ab, "Fit: A(1+B cos^{2}#theta)", "l");
            leg_theta_ab->AddEntry((TObject*)nullptr,
                                   Form("A = %.3f #pm %.3f", f_pi_theta_ab->GetParameter(0), f_pi_theta_ab->GetParError(0)),
                                   "");
            leg_theta_ab->AddEntry((TObject*)nullptr,
                                   Form("B = %.3f #pm %.3f", f_pi_theta_ab->GetParameter(1), f_pi_theta_ab->GetParError(1)),
                                   "");
            leg_theta_ab->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_theta_ab, (int)ndf_theta_ab), "");
            leg_theta_ab->Draw();
            if (generate_png) c_pi_theta_ab_pts->SaveAs("pi_theta_ab_points.png");
            if (generate_pdf) c_pi_theta_ab_pts->SaveAs("pi_theta_ab_points.pdf");

            // A/B rule filtered pi phi with fit
            TCanvas *c_pi_phi_ab_pts = new TCanvas("c_pi_phi_ab_points", "", 900, 700);
            auto g_pi_phi_ab = histToGraph(h_pi_phi_ab_sum, 24, kBlue+2);
            g_pi_phi_ab->SetTitle(";#Delta#varphi;Counts");
            g_pi_phi_ab->Draw("AP");
            TF1 *f_pi_phi_ab = new TF1("f_pi_phi_ab", cosExpr.c_str(), phiMin, phiMax);
            f_pi_phi_ab->SetNpx(500);
            f_pi_phi_ab->SetParName(0, "A0");
            f_pi_phi_ab->SetParameter(0, h_pi_phi_ab_sum->GetMaximum());
            for (int n = 1; n <= fitOrder; ++n) {
                if (!only_even_cosine || (n % 2 == 0)) {
                    f_pi_phi_ab->SetParName(n, Form("A%d", n));
                    f_pi_phi_ab->SetParameter(n, 0.05);
                }
            }
            h_pi_phi_ab_sum->Fit(f_pi_phi_ab, "Q0");
            f_pi_phi_ab->SetLineColor(kRed);
            f_pi_phi_ab->SetLineWidth(2);
            f_pi_phi_ab->Draw("SAME");
            double chi2_phi_ab = f_pi_phi_ab->GetChisquare();
            double ndf_phi_ab = static_cast<double>(f_pi_phi_ab->GetNDF());
            auto leg_phi_ab = new TLegend(0.18, 0.14, 0.48, 0.40);
            leg_phi_ab->SetTextSize(0.03);
            leg_phi_ab->SetFillStyle(0);
            leg_phi_ab->SetBorderSize(0);
            for (int n = 0; n <= fitOrder; ++n) {
                if (!only_even_cosine || (n % 2 == 0)) {
                    leg_phi_ab->AddEntry((TObject*)nullptr,
                                          Form("A%d = %.3f #pm %.3f", n, f_pi_phi_ab->GetParameter(n), f_pi_phi_ab->GetParError(n)),
                                          "");
                }
            }
            leg_phi_ab->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_ab, (int)ndf_phi_ab), "");
            leg_phi_ab->Draw();
            if (generate_png) c_pi_phi_ab_pts->SaveAs("pi_phi_ab_points.png");
            if (generate_pdf) c_pi_phi_ab_pts->SaveAs("pi_phi_ab_points.pdf");
        }

        if (export_pion_angles_min && (generate_png || generate_pdf)) {
            // Min-pair filtered pi theta with fit
            TCanvas *c_pi_theta_min_pts = new TCanvas("c_pi_theta_min_points", "Pion theta sum (min-pair)", 900, 700);
            auto g_pi_theta_min = histToGraph(h_pi_theta_min_sum, 24, kBlue+2);
            g_pi_theta_min->Draw("AP");
            TF1 *f_pi_theta_min = new TF1("f_pi_theta_min", "[0]*(1+[1]*pow(cos(x),2))", thetaMin, thetaMax);
            f_pi_theta_min->SetNpx(500);
            f_pi_theta_min->SetParName(0, "A");
            f_pi_theta_min->SetParName(1, "B");
            f_pi_theta_min->SetParameter(0, h_pi_theta_min_sum->GetMaximum());
            f_pi_theta_min->SetParameter(1, 0.0);
            h_pi_theta_min_sum->Fit(f_pi_theta_min, "Q0");
            f_pi_theta_min->SetLineColor(kRed);
            f_pi_theta_min->SetLineWidth(2);
            f_pi_theta_min->Draw("SAME");
            double chi2_theta_min = f_pi_theta_min->GetChisquare();
            double ndf_theta_min = static_cast<double>(f_pi_theta_min->GetNDF());
            auto leg_theta_min = new TLegend(0.12, 0.60, 0.58, 0.86);
            leg_theta_min->SetTextSize(0.02);
            leg_theta_min->SetFillStyle(0);
            leg_theta_min->SetBorderSize(0);
            leg_theta_min->AddEntry(g_pi_theta_min, "#pi^{+} #theta (min-pair)", "p");
            leg_theta_min->AddEntry(f_pi_theta_min, "Fit: A(1+B cos^{2}#theta)", "l");
            leg_theta_min->AddEntry((TObject*)nullptr,
                                    Form("A = %.3f #pm %.3f", f_pi_theta_min->GetParameter(0), f_pi_theta_min->GetParError(0)),
                                    "");
            leg_theta_min->AddEntry((TObject*)nullptr,
                                    Form("B = %.3f #pm %.3f", f_pi_theta_min->GetParameter(1), f_pi_theta_min->GetParError(1)),
                                    "");
            leg_theta_min->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_theta_min, (int)ndf_theta_min), "");
            leg_theta_min->Draw();
            if (generate_png) c_pi_theta_min_pts->SaveAs("pi_theta_min_points.png");
            if (generate_pdf) c_pi_theta_min_pts->SaveAs("pi_theta_min_points.pdf");

            // Min-pair filtered pi phi with fit
            TCanvas *c_pi_phi_min_pts = new TCanvas("c_pi_phi_min_points", "", 900, 700);
            auto g_pi_phi_min = histToGraph(h_pi_phi_min_sum, 24, kBlue+2);
            g_pi_phi_min->SetTitle(";#Delta#varphi;Counts");
            g_pi_phi_min->Draw("AP");
            TF1 *f_pi_phi_min = new TF1("f_pi_phi_min", cosExpr.c_str(), phiMin, phiMax);
            f_pi_phi_min->SetNpx(500);
            f_pi_phi_min->SetParName(0, "A0");
            f_pi_phi_min->SetParameter(0, h_pi_phi_min_sum->GetMaximum());
            for (int n = 1; n <= fitOrder; ++n) {
                if (!only_even_cosine || (n % 2 == 0)) {
                    f_pi_phi_min->SetParName(n, Form("A%d", n));
                    f_pi_phi_min->SetParameter(n, 0.05);
                }
            }
            h_pi_phi_min_sum->Fit(f_pi_phi_min, "Q0");
            f_pi_phi_min->SetLineColor(kRed);
            f_pi_phi_min->SetLineWidth(2);
            f_pi_phi_min->Draw("SAME");
            double chi2_phi_min = f_pi_phi_min->GetChisquare();
            double ndf_phi_min = static_cast<double>(f_pi_phi_min->GetNDF());
            auto leg_phi_min = new TLegend(0.18, 0.14, 0.48, 0.40);
            leg_phi_min->SetTextSize(0.03);
            leg_phi_min->SetFillStyle(0);
            leg_phi_min->SetBorderSize(0);
            for (int n = 0; n <= fitOrder; ++n) {
                if (!only_even_cosine || (n % 2 == 0)) {
                leg_phi_min->AddEntry((TObject*)nullptr,
                                      Form("A%d = %.3f #pm %.3f", n, f_pi_phi_min->GetParameter(n), f_pi_phi_min->GetParError(n)),
                                      "");
                }
            }
            leg_phi_min->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_min, (int)ndf_phi_min), "");
            leg_phi_min->Draw();
            if (generate_png) c_pi_phi_min_pts->SaveAs("pi_phi_min_points.png");
            if (generate_pdf) c_pi_phi_min_pts->SaveAs("pi_phi_min_points.pdf");
        }

        if (export_pairmass_phi && (generate_png || generate_pdf)) {
            std::string cosExpr = buildCosSeries(fitOrder);
            for (int i = 0; i < nPairMassCats; ++i) {
                if (export_pairmass_phi_angles) {
                    // Pair1 pi-phi
                    TCanvas *c_pair1_phi = new TCanvas(Form("c_pi_phi_pair1_cat_%s", pairMassCatTags[i]),
                                                       "", 900, 700);
                    auto g_phi_pair1 = histToGraph(h_pi_phi_pair1_cat[i], 24, kBlue+2);
                    g_phi_pair1->SetTitle(";#Delta#varphi;Counts");
                    g_phi_pair1->Draw("AP");
                    TF1 *f_phi_pair1 = new TF1(Form("f_pi_phi_pair1_cat_%s", pairMassCatTags[i]), cosExpr.c_str(), phiMin, phiMax);
                    f_phi_pair1->SetNpx(500);
                    f_phi_pair1->SetParName(0, "A0");
                    f_phi_pair1->SetParameter(0, h_pi_phi_pair1_cat[i]->GetMaximum());
                    for (int n = 1; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            f_phi_pair1->SetParName(n, Form("A%d", n));
                            f_phi_pair1->SetParameter(n, 0.05);
                        }
                    }
                    h_pi_phi_pair1_cat[i]->Fit(f_phi_pair1, "Q0");
                    f_phi_pair1->SetLineColor(kRed);
                    f_phi_pair1->SetLineWidth(2);
                    f_phi_pair1->Draw("SAME");
                    double chi2_phi_pair1 = f_phi_pair1->GetChisquare();
                    double ndf_phi_pair1 = static_cast<double>(f_phi_pair1->GetNDF());
                    auto leg_phi_pair1 = new TLegend(0.13, 0.14, 0.43, 0.40);
                    leg_phi_pair1->SetTextSize(0.03);
                    leg_phi_pair1->SetFillStyle(0);
                    leg_phi_pair1->SetBorderSize(0);
                    for (int n = 0; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            leg_phi_pair1->AddEntry((TObject*)nullptr,
                                                   Form("A%d = %.3f #pm %.3f", n, f_phi_pair1->GetParameter(n), f_phi_pair1->GetParError(n)),
                                                   "");
                        }
                    }
                    leg_phi_pair1->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_pair1, (int)ndf_phi_pair1), "");
                    leg_phi_pair1->Draw();
                    if (generate_png) c_pair1_phi->SaveAs(Form("pi_phi_pair1_cat_%s.png", pairMassCatTags[i]));
                    if (generate_pdf) c_pair1_phi->SaveAs(Form("pi_phi_pair1_cat_%s.pdf", pairMassCatTags[i]));

                    // Pair2 pi-phi
                    TCanvas *c_pair2_phi = new TCanvas(Form("c_pi_phi_pair2_cat_%s", pairMassCatTags[i]),
                                                       "", 900, 700);
                    auto g_phi_pair2 = histToGraph(h_pi_phi_pair2_cat[i], 24, kMagenta+2);
                    g_phi_pair2->SetTitle(";#Delta#varphi;Counts");
                    g_phi_pair2->Draw("AP");
                    TF1 *f_phi_pair2 = new TF1(Form("f_pi_phi_pair2_cat_%s", pairMassCatTags[i]), cosExpr.c_str(), phiMin, phiMax);
                    f_phi_pair2->SetNpx(500);
                    f_phi_pair2->SetParName(0, "A0");
                    f_phi_pair2->SetParameter(0, h_pi_phi_pair2_cat[i]->GetMaximum());
                    for (int n = 1; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            f_phi_pair2->SetParName(n, Form("A%d", n));
                            f_phi_pair2->SetParameter(n, 0.05);
                        }
                    }
                    h_pi_phi_pair2_cat[i]->Fit(f_phi_pair2, "Q0");
                    f_phi_pair2->SetLineColor(kRed);
                    f_phi_pair2->SetLineWidth(2);
                    f_phi_pair2->Draw("SAME");
                    double chi2_phi_pair2 = f_phi_pair2->GetChisquare();
                    double ndf_phi_pair2 = static_cast<double>(f_phi_pair2->GetNDF());
                    auto leg_phi_pair2 = new TLegend(0.13, 0.14, 0.43, 0.40);
                    leg_phi_pair2->SetTextSize(0.03);
                    leg_phi_pair2->SetFillStyle(0);
                    leg_phi_pair2->SetBorderSize(0);
                    for (int n = 0; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            leg_phi_pair2->AddEntry((TObject*)nullptr,
                                                   Form("A%d = %.3f #pm %.3f", n, f_phi_pair2->GetParameter(n), f_phi_pair2->GetParError(n)),
                                                   "");
                        }
                    }
                    leg_phi_pair2->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_pair2, (int)ndf_phi_pair2), "");
                    leg_phi_pair2->Draw();
                    if (generate_png) c_pair2_phi->SaveAs(Form("pi_phi_pair2_cat_%s.png", pairMassCatTags[i]));
                    if (generate_pdf) c_pair2_phi->SaveAs(Form("pi_phi_pair2_cat_%s.pdf", pairMassCatTags[i]));

                    // Pair1 pi-phi with pair pT<0.1
                    TCanvas *c_pair1_phi_ptcut = new TCanvas(Form("c_pi_phi_pair1_cat_%s_ptcut", pairMassCatTags[i]),
                                                             "", 900, 700);
                    auto g_phi_pair1_ptcut = histToGraph(h_pi_phi_pair1_cat_ptcut[i], 24, kBlue+2);
                    g_phi_pair1_ptcut->SetTitle(";#Delta#varphi;Counts");
                    g_phi_pair1_ptcut->Draw("AP");
                    TF1 *f_phi_pair1_ptcut = new TF1(Form("f_pi_phi_pair1_cat_%s_ptcut", pairMassCatTags[i]), cosExpr.c_str(), phiMin, phiMax);
                    f_phi_pair1_ptcut->SetNpx(500);
                    f_phi_pair1_ptcut->SetParName(0, "A0");
                    f_phi_pair1_ptcut->SetParameter(0, h_pi_phi_pair1_cat_ptcut[i]->GetMaximum());
                    for (int n = 1; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            f_phi_pair1_ptcut->SetParName(n, Form("A%d", n));
                            f_phi_pair1_ptcut->SetParameter(n, 0.05);
                        }
                    }
                    h_pi_phi_pair1_cat_ptcut[i]->Fit(f_phi_pair1_ptcut, "Q0");
                    f_phi_pair1_ptcut->SetLineColor(kRed);
                    f_phi_pair1_ptcut->SetLineWidth(2);
                    f_phi_pair1_ptcut->Draw("SAME");
                    double chi2_phi_pair1_ptcut = f_phi_pair1_ptcut->GetChisquare();
                    double ndf_phi_pair1_ptcut = static_cast<double>(f_phi_pair1_ptcut->GetNDF());
                    auto leg_phi_pair1_ptcut = new TLegend(0.10, 0.12, 0.43, 0.35);
                    leg_phi_pair1_ptcut->SetTextSize(0.03);
                    leg_phi_pair1_ptcut->SetFillStyle(0);
                    leg_phi_pair1_ptcut->SetBorderSize(0);
                    for (int n = 0; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            leg_phi_pair1_ptcut->AddEntry((TObject*)nullptr,
                                                         Form("A%d = %.3f #pm %.3f", n, f_phi_pair1_ptcut->GetParameter(n), f_phi_pair1_ptcut->GetParError(n)),
                                                         "");
                        }
                    }
                    leg_phi_pair1_ptcut->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_pair1_ptcut, (int)ndf_phi_pair1_ptcut), "");
                    leg_phi_pair1_ptcut->Draw();
                    {
                        TLatex note;
                        note.SetNDC();
                        note.SetTextSize(0.03);
                        note.DrawLatex(0.60, 0.86, "pair p_{T}<0.1");
                    }
                    if (generate_png) c_pair1_phi_ptcut->SaveAs(Form("pi_phi_pair1_cat_%s_ptcut.png", pairMassCatTags[i]));
                    if (generate_pdf) c_pair1_phi_ptcut->SaveAs(Form("pi_phi_pair1_cat_%s_ptcut.pdf", pairMassCatTags[i]));

                    // Pair2 pi-phi with pair pT<0.1
                    TCanvas *c_pair2_phi_ptcut = new TCanvas(Form("c_pi_phi_pair2_cat_%s_ptcut", pairMassCatTags[i]),
                                                             "", 900, 700);
                    auto g_phi_pair2_ptcut = histToGraph(h_pi_phi_pair2_cat_ptcut[i], 24, kMagenta+2);
                    g_phi_pair2_ptcut->SetTitle(";#Delta#varphi;Counts");
                    g_phi_pair2_ptcut->Draw("AP");
                    TF1 *f_phi_pair2_ptcut = new TF1(Form("f_pi_phi_pair2_cat_%s_ptcut", pairMassCatTags[i]), cosExpr.c_str(), phiMin, phiMax);
                    f_phi_pair2_ptcut->SetNpx(500);
                    f_phi_pair2_ptcut->SetParName(0, "A0");
                    f_phi_pair2_ptcut->SetParameter(0, h_pi_phi_pair2_cat_ptcut[i]->GetMaximum());
                    for (int n = 1; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            f_phi_pair2_ptcut->SetParName(n, Form("A%d", n));
                            f_phi_pair2_ptcut->SetParameter(n, 0.05);
                        }
                    }
                    h_pi_phi_pair2_cat_ptcut[i]->Fit(f_phi_pair2_ptcut, "Q0");
                    f_phi_pair2_ptcut->SetLineColor(kRed);
                    f_phi_pair2_ptcut->SetLineWidth(2);
                    f_phi_pair2_ptcut->Draw("SAME");
                    double chi2_phi_pair2_ptcut = f_phi_pair2_ptcut->GetChisquare();
                    double ndf_phi_pair2_ptcut = static_cast<double>(f_phi_pair2_ptcut->GetNDF());
                    auto leg_phi_pair2_ptcut = new TLegend(0.13, 0.14, 0.43, 0.40);
                    leg_phi_pair2_ptcut->SetTextSize(0.03);
                    leg_phi_pair2_ptcut->SetFillStyle(0);
                    leg_phi_pair2_ptcut->SetBorderSize(0);
                    for (int n = 0; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            leg_phi_pair2_ptcut->AddEntry((TObject*)nullptr,
                                                         Form("A%d = %.3f #pm %.3f", n, f_phi_pair2_ptcut->GetParameter(n), f_phi_pair2_ptcut->GetParError(n)),
                                                         "");
                        }
                    }
                    leg_phi_pair2_ptcut->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_phi_pair2_ptcut, (int)ndf_phi_pair2_ptcut), "");
                    leg_phi_pair2_ptcut->Draw();
                    {
                        TLatex note;
                        note.SetNDC();
                        note.SetTextSize(0.03);
                        note.DrawLatex(0.60, 0.86, "pair p_{T}<0.1");
                    }
                    if (generate_png) c_pair2_phi_ptcut->SaveAs(Form("pi_phi_pair2_cat_%s_ptcut.png", pairMassCatTags[i]));
                    if (generate_pdf) c_pair2_phi_ptcut->SaveAs(Form("pi_phi_pair2_cat_%s_ptcut.pdf", pairMassCatTags[i]));
                }

                if (export_pairmass_pair_pt) {
                    // Pair1 pT
                    TCanvas *c_pair1_pT_cat = new TCanvas(Form("c_pair1_pT_cat_%s", pairMassCatTags[i]),
                                                          Form("Pair p_{T} %s", pairMassCatTitles[i]), 900, 700);
                    h_pT_pair1_cat[i]->SetLineColor(kBlue+2);
                    h_pT_pair1_cat[i]->SetLineWidth(2);
                    h_pT_pair1_cat[i]->Draw();
                    auto leg_pair1_pt = new TLegend(0.64, 0.72, 0.86, 0.88);
                    leg_pair1_pt->SetTextSize(0.03);
                    leg_pair1_pt->SetFillStyle(0);
                    leg_pair1_pt->SetBorderSize(0);
                    leg_pair1_pt->AddEntry(h_pT_pair1_cat[i], "Pair1", "l");
                    leg_pair1_pt->Draw();
                    if (generate_png) c_pair1_pT_cat->SaveAs(Form("pair1_pT_cat_%s.png", pairMassCatTags[i]));
                    if (generate_pdf) c_pair1_pT_cat->SaveAs(Form("pair1_pT_cat_%s.pdf", pairMassCatTags[i]));

                    // Pair2 pT
                    TCanvas *c_pair2_pT_cat = new TCanvas(Form("c_pair2_pT_cat_%s", pairMassCatTags[i]),
                                                          Form("Pair2 p_{T} %s", pairMassCatTitles[i]), 900, 700);
                    h_pT_pair2_cat[i]->SetLineColor(kMagenta+2);
                    h_pT_pair2_cat[i]->SetLineWidth(2);
                    h_pT_pair2_cat[i]->Draw();
                    auto leg_pair2_pt = new TLegend(0.64, 0.72, 0.86, 0.88);
                    leg_pair2_pt->SetTextSize(0.03);
                    leg_pair2_pt->SetFillStyle(0);
                    leg_pair2_pt->SetBorderSize(0);
                    leg_pair2_pt->AddEntry(h_pT_pair2_cat[i], "Pair2", "l");
                    leg_pair2_pt->Draw();
                    if (generate_png) c_pair2_pT_cat->SaveAs(Form("pair2_pT_cat_%s.png", pairMassCatTags[i]));
                    if (generate_pdf) c_pair2_pT_cat->SaveAs(Form("pair2_pT_cat_%s.pdf", pairMassCatTags[i]));

                    // Pair1+Pair2 overlay pT
                    TCanvas *c_pair_pT_cat = new TCanvas(Form("c_pair_pT_cat_%s", pairMassCatTags[i]),
                                                         Form("Pair p_{T} %s", pairMassCatTitles[i]), 900, 700);
                    h_pT_pair1_cat[i]->SetLineColor(kBlue+2);
                    h_pT_pair1_cat[i]->SetLineWidth(2);
                    h_pT_pair2_cat[i]->SetLineColor(kMagenta+2);
                    h_pT_pair2_cat[i]->SetLineWidth(2);
                    double max_pair_pT = std::max(h_pT_pair1_cat[i]->GetMaximum(), h_pT_pair2_cat[i]->GetMaximum());
                    h_pT_pair1_cat[i]->SetMaximum(max_pair_pT * 1.15);
                    h_pT_pair1_cat[i]->Draw();
                    h_pT_pair2_cat[i]->Draw("SAME");
                    auto leg_pair_pt = new TLegend(0.64, 0.72, 0.86, 0.88);
                    leg_pair_pt->SetTextSize(0.03);
                    leg_pair_pt->SetFillStyle(0);
                    leg_pair_pt->SetBorderSize(0);
                    leg_pair_pt->AddEntry(h_pT_pair1_cat[i], "Pair1", "l");
                    leg_pair_pt->AddEntry(h_pT_pair2_cat[i], "Pair2", "l");
                    leg_pair_pt->Draw();
                    if (generate_png) c_pair_pT_cat->SaveAs(Form("pair_pT_cat_%s.png", pairMassCatTags[i]));
                    if (generate_pdf) c_pair_pT_cat->SaveAs(Form("pair_pT_cat_%s.pdf", pairMassCatTags[i]));
                }
            }
            if (export_pairmass_phi_combined && export_pairmass_phi_angles) {
                auto drawCombinedPhi = [&](TH1F *hist, const char *tag, const char *title, bool showPtCut) {
                    if (!hist) return;
                    TCanvas *c_comb = new TCanvas(Form("c_%s", tag), title, 900, 700);
                    auto g_comb = histToGraph(hist, 24, kBlue+2);
                    g_comb->SetTitle(";#Delta#varphi;Counts");
                    g_comb->Draw("AP");
                    TF1 *f_comb = new TF1(Form("f_%s", tag), cosExpr.c_str(), phiMin, phiMax);
                    f_comb->SetNpx(500);
                    f_comb->SetParName(0, "A0");
                    f_comb->SetParameter(0, hist->GetMaximum());
                    for (int n = 1; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            f_comb->SetParName(n, Form("A%d", n));
                            f_comb->SetParameter(n, 0.05);
                        }
                    }
                    hist->Fit(f_comb, "Q0");
                    f_comb->SetLineColor(kRed);
                    f_comb->SetLineWidth(2);
                    f_comb->Draw("SAME");
                    double chi2_comb = f_comb->GetChisquare();
                    double ndf_comb = static_cast<double>(f_comb->GetNDF());
                    auto leg_comb = new TLegend(0.13, 0.14, 0.43, 0.40);
                    leg_comb->SetTextSize(0.03);
                    leg_comb->SetFillStyle(0);
                    leg_comb->SetBorderSize(0);
                    for (int n = 0; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            leg_comb->AddEntry((TObject*)nullptr,
                                               Form("A%d = %.3f #pm %.3f", n, f_comb->GetParameter(n), f_comb->GetParError(n)),
                                               "");
                        }
                    }
                    leg_comb->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_comb, (int)ndf_comb), "");
                    leg_comb->Draw();
                    if (showPtCut) {
                        TLatex note;
                        note.SetNDC();
                        note.SetTextSize(0.03);
                        note.DrawLatex(0.60, 0.86, "pair p_{T}<0.1");
                    }
                    if (generate_png) c_comb->SaveAs(Form("%s.png", tag));
                    if (generate_pdf) c_comb->SaveAs(Form("%s.pdf", tag));
                };

                drawCombinedPhi(h_pi_phi_pair_cat_A_combined,
                                "pi_phi_pair_cat_A_combined",
                                "#pi^{+} #varphi (A)", false);
                drawCombinedPhi(h_pi_phi_pair_cat_C_combined,
                                "pi_phi_pair_cat_C_combined",
                                "#pi^{+} #varphi (C)", false);
                drawCombinedPhi(h_pi_phi_pair_cat_B1D2_combined,
                                "pi_phi_pair_cat_B1D2_combined",
                                "#pi^{+} #varphi (B+D, M<0.6)", false);
                drawCombinedPhi(h_pi_phi_pair_cat_B2D1_combined,
                                "pi_phi_pair_cat_B2D1_combined",
                                "#pi^{+} #varphi (B+D, 0.6<M<0.9)", false);
            }
            if (export_pairmass_phi_a_combined && export_pairmass_phi_angles) {
                auto drawCombinedPhiA = [&](TH1F *hist, const char *tag, const char *title, bool showPtCut) {
                    if (!hist) return;
                    TCanvas *c_comb = new TCanvas(Form("c_%s", tag), title, 900, 700);
                    auto g_comb = histToGraph(hist, 24, kBlue+2);
                    g_comb->SetTitle(";#Delta#varphi;Counts");
                    g_comb->Draw("AP");
                    TF1 *f_comb = new TF1(Form("f_%s", tag), cosExpr.c_str(), phiMin, phiMax);
                    f_comb->SetNpx(500);
                    f_comb->SetParName(0, "A0");
                    f_comb->SetParameter(0, hist->GetMaximum());
                    for (int n = 1; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            f_comb->SetParName(n, Form("A%d", n));
                            f_comb->SetParameter(n, 0.05);
                        }
                    }
                    hist->Fit(f_comb, "Q0");
                    f_comb->SetLineColor(kRed);
                    f_comb->SetLineWidth(2);
                    f_comb->Draw("SAME");
                    double chi2_comb = f_comb->GetChisquare();
                    double ndf_comb = static_cast<double>(f_comb->GetNDF());
                    auto leg_comb = new TLegend(0.13, 0.14, 0.43, 0.40);
                    leg_comb->SetTextSize(0.03);
                    leg_comb->SetFillStyle(0);
                    leg_comb->SetBorderSize(0);
                    for (int n = 0; n <= fitOrder; ++n) {
                        if (!only_even_cosine || (n % 2 == 0)) {
                            leg_comb->AddEntry((TObject*)nullptr,
                                               Form("A%d = %.3f #pm %.3f", n, f_comb->GetParameter(n), f_comb->GetParError(n)),
                                               "");
                        }
                    }
                    leg_comb->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f/%d", chi2_comb, (int)ndf_comb), "");
                    leg_comb->Draw();
                    if (showPtCut) {
                        TLatex note;
                        note.SetNDC();
                        note.SetTextSize(0.03);
                        note.DrawLatex(0.60, 0.86, "pair p_{T}<0.1");
                    }
                    if (generate_png) c_comb->SaveAs(Form("%s.png", tag));
                    if (generate_pdf) c_comb->SaveAs(Form("%s.pdf", tag));
                };

                drawCombinedPhiA(h_pi_a_phi_pair_cat_A_combined,
                                 "pi_a_phi_pair_cat_A_combined",
                                 "#pi^{+} #varphi (A)", true);
                drawCombinedPhiA(h_pi_a_phi_pair_cat_C_combined,
                                 "pi_a_phi_pair_cat_C_combined",
                                 "#pi^{+} #varphi (C)", false);
                drawCombinedPhiA(h_pi_a_phi_pair_cat_B1D2_combined,
                                 "pi_a_phi_pair_cat_B1D2_combined",
                                 "#pi^{+} #varphi (B+D, M<0.6)", false);
                drawCombinedPhiA(h_pi_a_phi_pair_cat_B2D1_combined,
                                 "pi_a_phi_pair_cat_B2D1_combined",
                                 "#pi^{+} #varphi (B+D, 0.6<M<0.9)", true);
            }
        }

        h_pi_theta_sum->Write();
        h_pi_phi_sum->Write();
        h_pi_theta_sum_win->Write();
        h_pi_phi_sum_win->Write();
        h_pi_theta_ab_sum->Write();
        h_pi_phi_ab_sum->Write();
        h_pi_theta_min_sum->Write();
        h_pi_phi_min_sum->Write();
        if (export_pion_angles_rho) {
            h_pi_theta_rho_best->Write();
            h_pi_theta_rho_other->Write();
            h_pi_phi_rho_best->Write();
            h_pi_phi_rho_other->Write();
            h_pi_theta_x_rho_best->Write();
            h_pi_theta_x_rho_other->Write();
            h_pi_phi_x_rho_best->Write();
            h_pi_phi_x_rho_other->Write();
            if (export_pion_angles_rho_rest) {
                h_pi_theta_rho_rest->Write();
                h_pi_phi_rho_rest->Write();
            }
            if (export_pion_angles_rho_x_rest) {
                h_pi_theta_x_rho_rest->Write();
                h_pi_phi_x_rho_rest->Write();
            }
        }
    

    // Save histograms to ROOT file
    fplot->cd();
    h_m13->Write(); h_m24->Write(); h_m14->Write(); h_m23->Write();
    h_pairMass_A1B1_vs_A2B2->Write();
    h_pT13->Write(); h_pT24->Write(); h_pT14->Write(); h_pT23->Write();
    h_pT13_win->Write(); h_pT24_win->Write(); h_pT14_win->Write(); h_pT23_win->Write();
    h_pT_sum->Write(); h_pT_sum_win->Write();
    h_pT_sum_totalpt->Write(); h_pT_sum_win_totalpt->Write();
    h_pT_sum_totalpt_mid->Write(); h_pT_sum_win_totalpt_mid->Write();
    h_pT_sum_totalpt_high->Write(); h_pT_sum_win_totalpt_high->Write();
    h_pT_sum_filtered->Write();
    h_pT_sum_minpair->Write();
    h_theta_A1->Write(); h_theta_A2->Write(); h_theta_B1->Write(); h_theta_B2->Write();
    h_phi_A1->Write(); h_phi_A2->Write(); h_phi_B1->Write(); h_phi_B2->Write();
    h_pi_theta_1->Write(); h_pi_theta_2->Write(); h_pi_theta_3->Write(); h_pi_theta_4->Write();
    h_pi_phi_1->Write(); h_pi_phi_2->Write(); h_pi_phi_3->Write(); h_pi_phi_4->Write();
    if (export_pion_angles_rho) {
        h_pi_theta_rho_best->Write();
        h_pi_theta_rho_other->Write();
        h_pi_phi_rho_best->Write();
        h_pi_phi_rho_other->Write();
        h_pi_theta_x_rho_best->Write();
        h_pi_theta_x_rho_other->Write();
        h_pi_phi_x_rho_best->Write();
        h_pi_phi_x_rho_other->Write();
        if (export_pion_angles_rho_x_ptcut) {
            h_pi_theta_x_rho_best_pt->Write();
            h_pi_theta_x_rho_other_pt->Write();
            if (export_pion_angles_rho_x_ptcut_mid_high) {
                h_pi_theta_x_rho_best_pt_mid->Write();
                h_pi_theta_x_rho_best_pt_high->Write();
                h_pi_theta_x_rho_other_pt_mid->Write();
                h_pi_theta_x_rho_other_pt_high->Write();
                h_pi_phi_x_rho_best_pt_mid->Write();
                h_pi_phi_x_rho_best_pt_high->Write();
                h_pi_phi_x_rho_other_pt_mid->Write();
                h_pi_phi_x_rho_other_pt_high->Write();
            }
            h_pi_phi_x_rho_best_pt->Write();
            h_pi_phi_x_rho_other_pt->Write();
        }
        if (export_pion_angles_rho_rest) {
            h_pi_theta_rho_rest->Write();
            h_pi_phi_rho_rest->Write();
        }
        if (export_pion_angles_rho_x_rest) {
            h_pi_theta_x_rho_rest->Write();
            h_pi_phi_x_rho_rest->Write();
            if (export_pion_angles_rho_x_ptcut) {
                h_pi_theta_x_rho_rest_pt->Write();
                h_pi_phi_x_rho_rest_pt->Write();
            }
        }
    }
    if (export_pairmass_phi) {
        for (int i = 0; i < nPairMassCats; ++i) {
            if (export_pairmass_phi_angles) {
                h_pi_phi_pair1_cat[i]->Write();
                h_pi_phi_pair2_cat[i]->Write();
                h_pi_phi_pair1_cat_ptcut[i]->Write();
                h_pi_phi_pair2_cat_ptcut[i]->Write();
            }
            if (export_pairmass_pair_pt) {
                h_pT_pair1_cat[i]->Write();
                h_pT_pair2_cat[i]->Write();
            }
        }
        if (export_pairmass_phi_combined) {
            if (h_pi_phi_pair_cat_A_combined) h_pi_phi_pair_cat_A_combined->Write();
            if (h_pi_phi_pair_cat_C_combined) h_pi_phi_pair_cat_C_combined->Write();
            if (h_pi_phi_pair_cat_B1D2_combined) h_pi_phi_pair_cat_B1D2_combined->Write();
            if (h_pi_phi_pair_cat_B2D1_combined) h_pi_phi_pair_cat_B2D1_combined->Write();
        }
        if (export_pair_phi_xphi_selected) {
            if (h_pair_phi_sel) h_pair_phi_sel->Write();
            if (h_pair_xphi_sel) h_pair_xphi_sel->Write();
        }
        if (export_pairmass_phi_a_combined) {
            if (h_pi_a_phi_pair_cat_A_combined) h_pi_a_phi_pair_cat_A_combined->Write();
            if (h_pi_a_phi_pair_cat_C_combined) h_pi_a_phi_pair_cat_C_combined->Write();
            if (h_pi_a_phi_pair_cat_B1D2_combined) h_pi_a_phi_pair_cat_B1D2_combined->Write();
            if (h_pi_a_phi_pair_cat_B2D1_combined) h_pi_a_phi_pair_cat_B2D1_combined->Write();
        }
    }
    h_total_pT->Write(); h_total_pT_win->Write();
    fplot->Close();
    if (switched_to_png_dir) {
        std::filesystem::current_path(original_workdir);
    }
    std::cout << "Overlay plot created: pairs_13_24_14_23_overlay" << (generate_png?".png":"") << std::endl;
}

// Auto-run when interpreted as a ROOT macro
//fourpair();
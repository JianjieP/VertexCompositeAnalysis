// Enhanced Four-Pion Pairing Analysis with Vector Angular Distributions
// This macro analyzes four-pion events using multiple pairing methods:
// 1. Mass difference minimization
// 2. pT difference minimization  
// 3. pT sum maximization
// 4. Random pairing
//
// Additionally computes phi and theta distributions for relative momentum
// vectors and their angular differences.
//
// NAMING CONVENTION:
// - pairAB: Difference between Pair A and Pair B momentum vectors
// - pairAT: Difference between Pair A and Total 4π momentum vectors  
// - pairBT: Difference between Pair B and Total 4π momentum vectors
//
// PHYSICS QUANTITIES:
// - cos_theta1/2: Cosine of angle between pair momentum and z-axis in 4π rest frame
// - phi_diff: Azimuthal angle difference between the two pairs
// - phi_prime: Angle between transverse momentum vectors of the two pairs

// Required includes
#include <vector>
#include <cmath>
#include <string>
#include <functional>
#include <random>
#include <filesystem>
#include <iostream>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"

// Structure to hold results for all pairing methods
struct PairingMethod {
    std::vector<float> masses;
    std::vector<float> pair_pT;
    std::vector<float> phi_diff; // Angle difference within each pair
    float phi_prime;  // Angle between the two pT pairs
    float cos_theta1;  // cos(theta1): pair polarization - angle between pair and z-axis in 4pi rest frame (only one value since cos_theta1[0] = -cos_theta1[1])
    std::vector<float> cos_theta2;  // cos(theta2): single pion polarization - angle between pion and z-axis in pair rest frame
    std::vector<float> vec_phi; // phi of the relative momentum vector for each pair (pairA, pairB)
    std::vector<float> vec_theta; // theta of the relative momentum vector for each pair (pairA, pairB)
    std::string name;
    PairingMethod(const std::string& n) : masses(2, 0), pair_pT(2, 0), phi_diff(2, 0), phi_prime(0), cos_theta1(0), cos_theta2(2, 0), vec_phi(2,0), vec_theta(2,0), name(n) {}
};

void twopairs() {
    // Control flags for output generation
    bool generate_pdf = false;   // Set to false to skip PDF generation
    bool generate_png = true;   // Set to false to skip PNG generation
    
    std::cout << "=== Output Generation Settings ===" << std::endl;
    std::cout << "PDF generation: " << (generate_pdf ? "ENABLED" : "DISABLED") << std::endl;
    std::cout << "PNG generation: " << (generate_png ? "ENABLED" : "DISABLED") << std::endl;
    std::cout << "==================================" << std::endl;

    // ----- User-selectable flags -----
    // Control which pairing methods to adopt (order: mass, pTdiff, pTsum, random)
    std::array<bool,4> method_enabled = { true, false, false, false }; // Only mass method enabled by default

    // Control which results to compute/store/plot
    bool save_phi_diff = true;    // phi differences within each pair
    bool save_vec_angles = true;  // relative 3-vector phi/theta for each pair
    bool save_cos_theta1 = true;  // cos(theta1) pair polarization
    bool save_cos_theta2 = true;  // cos(theta2) single pion polarization
    bool save_phi_prime = true;   // phi prime (angle between the two pT pairs)
    bool save_four_vectors = true; // four_pi phi/theta
    bool save_delta_distributions = true; // delta phi and delta theta distributions

    // You can toggle these flags at the top of the macro or expose them as function args.
    // If true, delete existing .png files in the working directory before generating new ones
    bool clean_png_before = true;
    std::cout << "Cleanup existing PNGs before generating: " << (clean_png_before ? "YES" : "NO") << std::endl;
    // Clarify delta histogram index mapping used later in the code
    std::cout << "Delta histogram index mapping: 0 = pair A vs pair B, 1 = pair A vs Total 4π, 2 = pair B vs Total 4π" << std::endl;

    // Save PNG files under a dedicated local folder
    std::filesystem::path png_output_dir = std::filesystem::path("twopairs_png");
    if (generate_png) {
        try {
            std::filesystem::create_directories(png_output_dir);
            std::cout << "PNG output directory: " << png_output_dir << std::endl;
        } catch (const std::filesystem::filesystem_error& ex) {
            std::cout << "Error: Could not create PNG output directory: " << ex.what() << std::endl;
            return;
        }
    }
    
    TFile *fin = TFile::Open("/eos/user/j/jianjie/fourpi/HIRun2025A_FourPi_0410/filtered_chained_ParticleTree.root");
    if (!fin || fin->IsZombie()) {
        std::cout << "Error: Could not open input file /eos/user/j/jianjie/fourpi/HIRun2025A_FourPi_0410/filtered_chained_ParticleTree.root" << std::endl;
        return;
    }
    TTree *tin = (TTree*)fin->Get("ParticleTree");
    if (!tin) {
        std::cout << "Error: Could not find ParticleTree in input file" << std::endl;
        fin->Close();
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
    tout->SetAutoFlush(10000);   // Reduced from 30000 for better memory usage
    tout->SetAutoSave(100000);   // Auto-save every 100k entries
    fout->SetCompressionLevel(1);  // Fast compression

    std::vector<PairingMethod> methods = {
        PairingMethod("mass"),
        PairingMethod("pTdiff"),
        PairingMethod("pTsum"),
        PairingMethod("random")
    };
    
    std::vector<float> inv_masses(4, 0);
    tout->Branch("inv_masses", &inv_masses);
    
    // Add branches for four-pion kinematics
    float four_pi_mass = 0;
    float four_pi_pT = 0;
    float four_phi = 0;
    float four_theta = 0;
    std::vector<float> individual_pT(4, 0);
    tout->Branch("four_pi_mass", &four_pi_mass);
    tout->Branch("four_pi_pT", &four_pi_pT);
    tout->Branch("four_phi", &four_phi);
    tout->Branch("four_theta", &four_theta);
    tout->Branch("individual_pT", &individual_pT);
    
    // Create branches for all methods
    for (auto& method : methods) {
        tout->Branch(("inv_mass_pair_" + method.name).c_str(), &method.masses);
        tout->Branch(("pair_pT_" + method.name).c_str(), &method.pair_pT);
        tout->Branch(("phi_diff_" + method.name).c_str(), &method.phi_diff);
        tout->Branch(("vec_phi_" + method.name).c_str(), &method.vec_phi);
        tout->Branch(("vec_theta_" + method.name).c_str(), &method.vec_theta);
        tout->Branch(("phi_prime_" + method.name).c_str(), &method.phi_prime);
        tout->Branch(("cos_theta1_" + method.name).c_str(), &method.cos_theta1);
        tout->Branch(("cos_theta2_" + method.name).c_str(), &method.cos_theta2);
    }

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

    // Random number generator for random pairing method (initialize once)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

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

        TLorentzVector p1 = get_pion(1);
        TLorentzVector p2 = get_pion(2);
        TLorentzVector p3 = get_pion(3);
        TLorentzVector p4 = get_pion(4);

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

        // Calculate four-pion kinematics
        TLorentzVector total = p1 + p2 + p3 + p4;
        four_pi_mass = total.M();
        four_pi_pT = total.Pt();
    TVector3 totalVec = total.Vect();
    four_phi = totalVec.Phi();
    four_theta = totalVec.Theta();
        individual_pT[0] = p1.Pt();
        individual_pT[1] = p2.Pt();
        individual_pT[2] = p3.Pt();
        individual_pT[3] = p4.Pt();

        // Pre-calculate pT values once
        float pT_pair13 = pair13.Pt();
        float pT_pair24 = pair24.Pt();
        float pT_pair14 = pair14.Pt();
        float pT_pair23 = pair23.Pt();

        // Function to select pairing and calculate additional properties
        auto processPairing = [&](bool use1324, TLorentzVector p1, TLorentzVector p2, 
                                 TLorentzVector p3, TLorentzVector p4) -> PairingMethod {
            PairingMethod result("");
            TLorentzVector pairA, pairB;
            
            if (use1324) {
                // Pairing 1&3, 2&4
                pairA = p1 + p3;
                pairB = p2 + p4;
                result.masses[0] = pairA.M();
                result.masses[1] = pairB.M();
                result.pair_pT[0] = pairA.Pt();
                result.pair_pT[1] = pairB.Pt();
                // Phi difference within each pair
                result.phi_diff[0] = fabs(p1.Phi() - p3.Phi());
                result.phi_diff[1] = fabs(p2.Phi() - p4.Phi());
                // relative momentum vectors: p1 - p3, p2 - p4
                TVector3 vecA = p1.Vect() - p3.Vect();
                TVector3 vecB = p2.Vect() - p4.Vect();
                if (vecA.Mag() > 1e-6) {
                    result.vec_phi[0] = vecA.Phi();
                    result.vec_theta[0] = vecA.Theta();
                } else { result.vec_phi[0] = 0.; result.vec_theta[0] = 0.; }
                if (vecB.Mag() > 1e-6) {
                    result.vec_phi[1] = vecB.Phi();
                    result.vec_theta[1] = vecB.Theta();
                } else { result.vec_phi[1] = 0.; result.vec_theta[1] = 0.; }
            } else {
                // Pairing 1&4, 2&3
                pairA = p1 + p4;
                pairB = p2 + p3;
                result.masses[0] = pairA.M();
                result.masses[1] = pairB.M();
                result.pair_pT[0] = pairA.Pt();
                result.pair_pT[1] = pairB.Pt();
                // Phi difference within each pair
                result.phi_diff[0] = fabs(p1.Phi() - p4.Phi());
                result.phi_diff[1] = fabs(p2.Phi() - p3.Phi());
                // relative momentum vectors: p1 - p4, p2 - p3
                TVector3 vecA = p1.Vect() - p4.Vect();
                TVector3 vecB = p2.Vect() - p3.Vect();
                if (vecA.Mag() > 1e-6) {
                    result.vec_phi[0] = vecA.Phi();
                    result.vec_theta[0] = vecA.Theta();
                } else { result.vec_phi[0] = 0.; result.vec_theta[0] = 0.; }
                if (vecB.Mag() > 1e-6) {
                    result.vec_phi[1] = vecB.Phi();
                    result.vec_theta[1] = vecB.Theta();
                } else { result.vec_phi[1] = 0.; result.vec_theta[1] = 0.; }
            }
            
            // Ensure phi differences are in [0, pi] range
            for (int i = 0; i < 2; i++) {
                if (result.phi_diff[i] > M_PI) {
                    result.phi_diff[i] = 2 * M_PI - result.phi_diff[i];
                }
            }
            
            // Calculate phi prime: angle between the two pT pairs
            result.phi_prime = fabs(pairA.Phi() - pairB.Phi());
            if (result.phi_prime > M_PI) {
                result.phi_prime = 2 * M_PI - result.phi_prime;
            }
            
            // Calculate cos(theta1) and cos(theta2) for pair and single pion polarizations
            TLorentzVector fourPion = p1 + p2 + p3 + p4;  // Total four-pion system
            
            if (use1324) {
                // For pairing 1&3, 2&4
                TLorentzVector pairA_rest4pi = pairA;
                TLorentzVector pairB_rest4pi = pairB;
                pairA_rest4pi.Boost(-fourPion.BoostVector());  // Boost pairA to four-pion rest frame
                pairB_rest4pi.Boost(-fourPion.BoostVector());  // Boost pairB to four-pion rest frame
                
                // theta1: pair polarization - angle between pair momentum and z-axis in 4pi rest frame
                // Only store one value since cos_theta1[0] = -cos_theta1[1] in 4pi rest frame
                result.cos_theta1 = pairA_rest4pi.Pz() / pairA_rest4pi.P();  // cos(theta1) for pair A (1&3)
                
                // theta2: single pion polarization - angle between pion and z-axis in pair rest frame
                TLorentzVector p1_restA = p1;
                TLorentzVector p2_restB = p2;
                p1_restA.Boost(-pairA.BoostVector());  // Boost p1 to pairA rest frame
                p2_restB.Boost(-pairB.BoostVector());  // Boost p2 to pairB rest frame
                result.cos_theta2[0] = p1_restA.Pz() / p1_restA.P();  // cos(theta2) for p1 in pairA rest frame
                result.cos_theta2[1] = p2_restB.Pz() / p2_restB.P();  // cos(theta2) for p2 in pairB rest frame
            } else {
                // For pairing 1&4, 2&3
                TLorentzVector pairA_rest4pi = pairA;
                TLorentzVector pairB_rest4pi = pairB;
                pairA_rest4pi.Boost(-fourPion.BoostVector());  // Boost pairA to four-pion rest frame
                pairB_rest4pi.Boost(-fourPion.BoostVector());  // Boost pairB to four-pion rest frame
                
                // theta1: pair polarization - angle between pair momentum and z-axis in 4pi rest frame
                // Only store one value since cos_theta1[0] = -cos_theta1[1] in 4pi rest frame
                result.cos_theta1 = pairA_rest4pi.Pz() / pairA_rest4pi.P();  // cos(theta1) for pair A (1&4)
                
                // theta2: single pion polarization - angle between pion and z-axis in pair rest frame
                TLorentzVector p1_restA = p1;
                TLorentzVector p2_restB = p2;
                p1_restA.Boost(-pairA.BoostVector());  // Boost p1 to pairA rest frame
                p2_restB.Boost(-pairB.BoostVector());  // Boost p2 to pairB rest frame
                result.cos_theta2[0] = p1_restA.Pz() / p1_restA.P();  // cos(theta2) for p1 in pairA rest frame
                result.cos_theta2[1] = p2_restB.Pz() / p2_restB.P();  // cos(theta2) for p2 in pairB rest frame
            }
            
            return result;
        };

        // Method 1: Compare mass differences
        if (method_enabled[0]) {
            float diff_1324 = fabs(m13 - m24);
            float diff_1423 = fabs(m14 - m23);
            auto result_mass = processPairing(diff_1324 > diff_1423, p1, p2, p3, p4);
            methods[0].masses = result_mass.masses;
            methods[0].pair_pT = result_mass.pair_pT;
            methods[0].phi_diff = result_mass.phi_diff;
            methods[0].phi_prime = result_mass.phi_prime;
            methods[0].cos_theta1 = result_mass.cos_theta1;
            methods[0].cos_theta2 = result_mass.cos_theta2;
            methods[0].vec_phi = result_mass.vec_phi;
            methods[0].vec_theta = result_mass.vec_theta;
        }

        // Method 2: Compare pT magnitude differences
        if (method_enabled[1]) {
            float diff_pT_1324 = fabs(pT_pair13 - pT_pair24);
            float diff_pT_1423 = fabs(pT_pair14 - pT_pair23);
            auto result_pT = processPairing(diff_pT_1324 < diff_pT_1423, p1, p2, p3, p4);
            methods[1].masses = result_pT.masses;
            methods[1].pair_pT = result_pT.pair_pT;
            methods[1].phi_diff = result_pT.phi_diff;
            methods[1].phi_prime = result_pT.phi_prime;
            methods[1].cos_theta1 = result_pT.cos_theta1;
            methods[1].cos_theta2 = result_pT.cos_theta2;
            methods[1].vec_phi = result_pT.vec_phi;
            methods[1].vec_theta = result_pT.vec_theta;
        }

        // Method 3: Compare pT magnitude sums
        if (method_enabled[2]) {
            float pTsum_1324 = pT_pair13 + pT_pair24;
            float pTsum_1423 = pT_pair14 + pT_pair23;
            auto result_pTsum = processPairing(pTsum_1324 < pTsum_1423, p1, p2, p3, p4);
            methods[2].masses = result_pTsum.masses;
            methods[2].pair_pT = result_pTsum.pair_pT;
            methods[2].phi_diff = result_pTsum.phi_diff;
            methods[2].phi_prime = result_pTsum.phi_prime;
            methods[2].cos_theta1 = result_pTsum.cos_theta1;
            methods[2].cos_theta2 = result_pTsum.cos_theta2;
            methods[2].vec_phi = result_pTsum.vec_phi;
            methods[2].vec_theta = result_pTsum.vec_theta;
        }

        // Method 4: Random selection between neutral charge combinations
        if (method_enabled[3]) {
            bool random_choice = dis(gen) < 0.5;  // 50% probability for each choice
            auto result_random = processPairing(random_choice, p1, p2, p3, p4);
            methods[3].masses = result_random.masses;
            methods[3].pair_pT = result_random.pair_pT;
            methods[3].phi_diff = result_random.phi_diff;
            methods[3].phi_prime = result_random.phi_prime;
            methods[3].cos_theta1 = result_random.cos_theta1;
            methods[3].cos_theta2 = result_random.cos_theta2;
            methods[3].vec_phi = result_random.vec_phi;
            methods[3].vec_theta = result_random.vec_theta;
        }

        tout->Fill();
    }

    std::cout << "Finished processing entries. Writing output..." << std::endl;
    fout->Write();
    fout->Close();
    fin->Close();
    std::cout << "Output file written successfully." << std::endl;

    // Plot inv_mass_pair as histogram and save image
    std::cout << "Creating histograms and plots..." << std::endl;
    TFile *fplot = TFile::Open("new_tree.root", "UPDATE");
    TTree *tplot = (TTree*)fplot->Get("NewTree");
    
    // Enable faster tree reading
    tplot->SetCacheSize(50000000);  // 50MB cache
    tplot->AddBranchToCache("*", kTRUE);  // Cache all branches
    
    // Create arrays to hold pointers to the data
    std::vector<std::vector<float>*> mass_data(4, nullptr);
    std::vector<std::vector<float>*> pT_data(4, nullptr);
    std::vector<std::vector<float>*> phi_data(4, nullptr);
    std::vector<std::vector<float>*> vec_phi_data(4, nullptr);
    std::vector<std::vector<float>*> vec_theta_data(4, nullptr);
    std::vector<float> cos_theta1_data(4, 0);  // Pair polarization (single value per method)
    std::vector<std::vector<float>*> cos_theta2_data(4, nullptr);  // Single pion polarization
    std::vector<float> phi_prime_data(4, 0);  // Direct values, not pointers
    
    // Four-pion kinematic data - scalar branches, not pointers
    float four_pi_mass_data = 0;
    float four_pi_pT_data = 0;
    float four_phi_data = 0;
    float four_theta_data = 0;
    std::vector<float> *individual_pT_data = nullptr;
    
    // Set branch addresses for all methods
    for (size_t i = 0; i < methods.size(); ++i) {
        std::string mass_branch = "inv_mass_pair_" + methods[i].name;
        std::string pT_branch = "pair_pT_" + methods[i].name;
        std::string phi_branch = "phi_diff_" + methods[i].name;
        std::string vec_phi_branch = "vec_phi_" + methods[i].name;
        std::string vec_theta_branch = "vec_theta_" + methods[i].name;
        std::string phi_prime_branch = "phi_prime_" + methods[i].name;
        std::string cos_theta1_branch = "cos_theta1_" + methods[i].name;
        std::string cos_theta2_branch = "cos_theta2_" + methods[i].name;
        tplot->SetBranchAddress(mass_branch.c_str(), &mass_data[i]);
        tplot->SetBranchAddress(pT_branch.c_str(), &pT_data[i]);
        tplot->SetBranchAddress(phi_branch.c_str(), &phi_data[i]);
        tplot->SetBranchAddress(vec_phi_branch.c_str(), &vec_phi_data[i]);
        tplot->SetBranchAddress(vec_theta_branch.c_str(), &vec_theta_data[i]);
        tplot->SetBranchAddress(phi_prime_branch.c_str(), &phi_prime_data[i]);
        tplot->SetBranchAddress(cos_theta1_branch.c_str(), &cos_theta1_data[i]);
        tplot->SetBranchAddress(cos_theta2_branch.c_str(), &cos_theta2_data[i]);
    }
    
    // Set branch addresses for four-pion kinematics  
    tplot->SetBranchAddress("four_pi_mass", &four_pi_mass_data);
    tplot->SetBranchAddress("four_pi_pT", &four_pi_pT_data);
    tplot->SetBranchAddress("four_phi", &four_phi_data);
    tplot->SetBranchAddress("four_theta", &four_theta_data);
    tplot->SetBranchAddress("individual_pT", &individual_pT_data);

    // Create histograms for all methods
    std::vector<std::vector<TH1F*>> mass_hists(4);
    std::vector<std::vector<TH1F*>> pT_hists(4);
    std::vector<std::vector<TH1F*>> phi_hists(4);
    std::vector<std::vector<TH1F*>> vec_phi_hists(4);
    std::vector<std::vector<TH1F*>> vec_theta_hists(4);
    // delta histograms between the three vectors: (A,B), (A,total), (B,total)
    std::vector<std::vector<TH1F*>> delta_phi_hists(4);
    std::vector<std::vector<TH1F*>> delta_theta_hists(4);
    std::vector<TH1F*> cos_theta1_hists(4);  // Pair polarization (single histogram per method)
    std::vector<std::vector<TH1F*>> cos_theta2_hists(4);  // Single pion polarization
    std::vector<TH1F*> phi_prime_hists(4);
    
    std::vector<std::string> method_titles = {"Mass Difference", "pT Difference", "pT Sum", "Random"};
    
    // Create four-pion kinematics histograms
    TH1F *h_four_pi_mass = new TH1F("h_four_pi_mass", "Four-Pion Mass Distribution;M_{4#pi} (GeV);Events", 200, 0, 10);
    h_four_pi_mass->SetStats(0);
    
    TH1F *h_four_pi_pT = new TH1F("h_four_pi_pT", "Four-Pion pT Distribution;pT_{4#pi} (GeV);Events", 200, 0, 1);
    h_four_pi_pT->SetStats(0);
    
    // Individual pion pT histograms (all four pions in one histogram)
    TH1F *h_individual_pT = new TH1F("h_individual_pT", "Individual Pion pT Distribution;pT_{#pi} (GeV);Events", 100, 0, 5);
    h_individual_pT->SetStats(0);
    
    for (size_t i = 0; i < methods.size(); ++i) {
        // Mass histograms (separate for plotting)
        for (int j = 0; j < 2; ++j) {
            std::string name = "h_mass_" + methods[i].name + "_" + std::to_string(j);
            std::string title = "Invariant Mass Pair (" + method_titles[i] + " Selection);Mass (GeV);Events";
            mass_hists[i].push_back(new TH1F(name.c_str(), title.c_str(), 200, 0, 3));
            mass_hists[i][j]->SetStats(0);
        }
        
        // pT histograms (separate for plotting)
        for (int j = 0; j < 2; ++j) {
            std::string name = "h_pT_" + methods[i].name + "_" + std::to_string(j);
            std::string title = "Pair pT (" + method_titles[i] + " Selection);pT (GeV);Events";
            pT_hists[i].push_back(new TH1F(name.c_str(), title.c_str(), 100, 0, 5));
            pT_hists[i][j]->SetStats(0);
        }
        
        // Phi difference histograms (separate for plotting)
        for (int j = 0; j < 2; ++j) {
            std::string name = "h_phi_" + methods[i].name + "_" + std::to_string(j);
            std::string title = "Phi Difference Within Each Pair (" + method_titles[i] + " Selection);#Delta#phi within pair (rad);Events";
            phi_hists[i].push_back(new TH1F(name.c_str(), title.c_str(), 50, 0, M_PI));
            phi_hists[i][j]->SetStats(0);
        }

        // Relative momentum vector phi/theta histograms for pairA/pairB
        for (int j = 0; j < 2; ++j) {
            std::string name_phi = "h_vec_phi_" + methods[i].name + "_" + std::to_string(j);
            std::string title_phi = "Relative Momentum Phi (" + method_titles[i] + " Selection);#phi (rad);Events";
            vec_phi_hists[i].push_back(new TH1F(name_phi.c_str(), title_phi.c_str(), 64, -M_PI, M_PI));
            vec_phi_hists[i][j]->SetStats(0);

            std::string name_theta = "h_vec_theta_" + methods[i].name + "_" + std::to_string(j);
            std::string title_theta = "Relative Momentum Theta (" + method_titles[i] + " Selection);#theta (rad);Events";
            vec_theta_hists[i].push_back(new TH1F(name_theta.c_str(), title_theta.c_str(), 64, 0, M_PI));
            vec_theta_hists[i][j]->SetStats(0);
        }

        // Delta histograms (index mapping: 0 = pair A vs pair B, 1 = pair A vs total, 2 = pair B vs total)
        {
            const std::array<std::string,3> suffix = {"pairAB","pairAT","pairBT"};
            const std::array<std::string,3> label = {"Pair A vs Pair B","Pair A vs Total 4π","Pair B vs Total 4π"};
            for (int j = 0; j < 3; ++j) {
                std::string name_dphi = "h_delta_phi_" + methods[i].name + "_" + suffix[j];
                std::string title_dphi = std::string("Delta Phi (") + label[j] + ") (" + method_titles[i] + " Selection);#Delta#phi (rad);Events";
                delta_phi_hists[i].push_back(new TH1F(name_dphi.c_str(), title_dphi.c_str(), 64, 0, M_PI));
                delta_phi_hists[i][j]->SetStats(0);

                std::string name_dtheta = "h_delta_theta_" + methods[i].name + "_" + suffix[j];
                std::string title_dtheta = std::string("Delta Theta (") + label[j] + ") (" + method_titles[i] + " Selection);#Delta#theta (rad);Events";
                delta_theta_hists[i].push_back(new TH1F(name_dtheta.c_str(), title_dtheta.c_str(), 64, 0, M_PI));
                delta_theta_hists[i][j]->SetStats(0);
            }
        }
        
        // cos(theta1) histogram - pair polarization (single histogram per method)
        {
            std::string name = "h_cos_theta1_" + methods[i].name;
            std::string title = "cos(#theta_{1}) Pair Polarization - Pair Momentum vs Z-axis in 4#pi Rest Frame (" + method_titles[i] + " Selection);cos(#theta_{1});Events";
            cos_theta1_hists[i] = new TH1F(name.c_str(), title.c_str(), 100, -1, 1);
            cos_theta1_hists[i]->SetStats(0);
        }
        
        // cos(theta2) histograms - single pion polarization (separate for plotting)
        for (int j = 0; j < 2; ++j) {
            std::string name = "h_cos_theta2_" + methods[i].name + "_" + std::to_string(j);
            std::string title = "cos(#theta_{2}) Single Pion Polarization - Pion Momentum vs Z-axis in Pair Rest Frame (" + method_titles[i] + " Selection);cos(#theta_{2});Events";
            cos_theta2_hists[i].push_back(new TH1F(name.c_str(), title.c_str(), 100, -1, 1));
            cos_theta2_hists[i][j]->SetStats(0);
        }
        
        // Phi prime histograms (angle between pT pairs)
        {
            std::string name = "h_phi_prime_" + methods[i].name;
            std::string title = "Phi Prime - Angle Between Two Pair pT Vectors (" + method_titles[i] + " Selection);#phi' (rad);Events";
            phi_prime_hists[i] = new TH1F(name.c_str(), title.c_str(), 50, 0, M_PI);
            phi_prime_hists[i]->SetStats(0);
        }
    }

    Long64_t nplot = tplot->GetEntries();
    std::cout << "Filling histograms from " << nplot << " entries..." << std::endl;
    
    // Progress reporting for histogram filling
    Long64_t plot_progress = nplot / 10;  // Report every 10%
    if (plot_progress == 0) plot_progress = 1000;
    
    for (Long64_t i=0; i<nplot; ++i) {
        if (i % plot_progress == 0) {
            std::cout << "Filling histograms: " << (100.0 * i / nplot) << "% complete" << std::endl;
        }
        
        tplot->GetEntry(i);
        
        // Fill four-pion kinematics histograms - access as scalar values now
        if (save_four_vectors) {
            h_four_pi_mass->Fill(four_pi_mass_data);
            h_four_pi_pT->Fill(four_pi_pT_data);
        }
        
        if (individual_pT_data && individual_pT_data->size() == 4) {
            for (int j = 0; j < 4; ++j) {
                h_individual_pT->Fill((*individual_pT_data)[j]);
            }
        }
        
        // Fill histograms for all methods
        for (size_t m = 0; m < methods.size(); ++m) {
            if (mass_data[m] && mass_data[m]->size() == 2) {
                // Fill separate histograms (for plotting)
                mass_hists[m][0]->Fill((*mass_data[m])[0]);
                mass_hists[m][1]->Fill((*mass_data[m])[1]);
            }
            if (pT_data[m] && pT_data[m]->size() == 2) {
                // Fill separate histograms (for plotting)
                pT_hists[m][0]->Fill((*pT_data[m])[0]);
                pT_hists[m][1]->Fill((*pT_data[m])[1]);
            }
            if (save_phi_diff && phi_data[m] && phi_data[m]->size() == 2 && method_enabled[m]) {
                // Fill separate histograms (for plotting)
                phi_hists[m][0]->Fill((*phi_data[m])[0]);
                phi_hists[m][1]->Fill((*phi_data[m])[1]);
            }
            // Fill cos(theta1) histogram - pair polarization (single value)
            if (save_cos_theta1 && method_enabled[m]) cos_theta1_hists[m]->Fill(cos_theta1_data[m]);
            if (save_cos_theta2 && cos_theta2_data[m] && cos_theta2_data[m]->size() == 2 && method_enabled[m]) {
                // Fill cos(theta2) histograms - single pion polarization
                cos_theta2_hists[m][0]->Fill((*cos_theta2_data[m])[0]);
                cos_theta2_hists[m][1]->Fill((*cos_theta2_data[m])[1]);
            }
            // Fill phi_prime histogram (direct value, no pointer dereferencing)
            if (save_phi_prime && method_enabled[m]) phi_prime_hists[m]->Fill(phi_prime_data[m]);
            // Fill relative momentum phi/theta histograms for pairA/pairB if available
            if (vec_phi_data[m] && vec_phi_data[m]->size() == 2 && vec_theta_data[m] && vec_theta_data[m]->size() == 2) {
                float phiA = (*vec_phi_data[m])[0];
                float phiB = (*vec_phi_data[m])[1];
                float thetaA = (*vec_theta_data[m])[0];
                float thetaB = (*vec_theta_data[m])[1];
                if (save_vec_angles) {
                    vec_phi_hists[m][0]->Fill(phiA);
                    vec_phi_hists[m][1]->Fill(phiB);
                    vec_theta_hists[m][0]->Fill(thetaA);
                    vec_theta_hists[m][1]->Fill(thetaB);
                }

                // Delta phi: A vs B
                float dphiAB = fabs(phiA - phiB);
                if (dphiAB > M_PI) dphiAB = 2*M_PI - dphiAB;
                if (save_vec_angles) delta_phi_hists[m][0]->Fill(dphiAB);

                // Delta theta: A vs B
                float dthetaAB = fabs(thetaA - thetaB);
                if (save_vec_angles) delta_theta_hists[m][0]->Fill(dthetaAB);

                // A vs total and B vs total (use four_phi_data/four_theta_data)
                float phiTot = four_phi_data;
                float thetaTot = four_theta_data;
                float dphiAT = fabs(phiA - phiTot);
                if (dphiAT > M_PI) dphiAT = 2*M_PI - dphiAT;
                if (save_vec_angles) delta_phi_hists[m][1]->Fill(dphiAT);
                float dthetaAT = fabs(thetaA - thetaTot);
                if (save_vec_angles) delta_theta_hists[m][1]->Fill(dthetaAT);

                float dphiBT = fabs(phiB - phiTot);
                if (dphiBT > M_PI) dphiBT = 2*M_PI - dphiBT;
                if (save_vec_angles) delta_phi_hists[m][2]->Fill(dphiBT);
                float dthetaBT = fabs(thetaB - thetaTot);
                if (save_vec_angles) delta_theta_hists[m][2]->Fill(dthetaBT);
            }
        }
    }
    
    std::cout << "Histogram filling complete." << std::endl;
    
    // Clean up existing PNG files in the local output folder before generating new ones
    if (generate_png && clean_png_before) {
        std::cout << "Cleaning up existing PNG files in " << png_output_dir << " ..." << std::endl;
        try {
            for (const auto& entry : std::filesystem::directory_iterator(png_output_dir)) {
                if (entry.is_regular_file() && entry.path().extension() == ".png") {
                    std::filesystem::remove(entry.path());
                    std::cout << "Deleted: " << entry.path().filename() << std::endl;
                }
            }
        } catch (const std::filesystem::filesystem_error& ex) {
            std::cout << "Warning: Could not clean up PNG files: " << ex.what() << std::endl;
        }
    }
    
    if (generate_pdf || generate_png) {
        std::cout << "Creating plots..." << std::endl;
    }
    
    // Helper function to create and save plots for pair distributions
    auto createPlot = [&](const std::vector<TH1F*>& hists, const std::string& title, 
                         const std::string& filename) {
        TCanvas *c = new TCanvas(("c_" + filename).c_str(), title.c_str(), 800, 600);
        
        hists[0]->SetLineColor(kBlue);
        hists[1]->SetLineColor(kRed);
        
        float max0 = hists[0]->GetMaximum();
        float max1 = hists[1]->GetMaximum();
        float ymax = 1.2 * std::max(max0, max1);
        hists[0]->SetMaximum(ymax);
        hists[1]->SetMaximum(ymax);
        
        hists[0]->Draw();
        hists[1]->Draw("SAME");
        
        auto legend = new TLegend(0.7,0.7,0.9,0.9);
        legend->AddEntry(hists[0], "Pair 1", "l");
        legend->AddEntry(hists[1], "Pair 2", "l");
        legend->Draw();
        
        if (generate_png) {
            std::filesystem::path png_path = png_output_dir / (filename + ".png");
            c->SaveAs(png_path.string().c_str());
        }
        if (generate_pdf) {
            c->SaveAs((filename + ".pdf").c_str());
        }
        delete c;
    };

    // Helper function for single histogram plots (phi prime)
    auto createSinglePlot = [&](TH1F* hist, const std::string& title, 
                               const std::string& filename) {
        TCanvas *c = new TCanvas(("c_" + filename).c_str(), title.c_str(), 800, 600);
        
        hist->SetLineColor(kBlue);
        hist->SetLineWidth(2);
        hist->Draw();
        
        if (generate_png) {
            std::filesystem::path png_path = png_output_dir / (filename + ".png");
            c->SaveAs(png_path.string().c_str());
        }
        if (generate_pdf) {
            c->SaveAs((filename + ".pdf").c_str());
        }
        delete c;
    };

    // Create plots for all methods and all distributions
    std::vector<std::string> mass_filenames = {"inv_mass_pair_mass_selection", 
                                              "inv_mass_pair_pTdiff_selection", 
                                              "inv_mass_pair_pTsum_selection",
                                              "inv_mass_pair_random_selection"};
    std::vector<std::string> pT_filenames = {"pair_pT_mass_selection", 
                                            "pair_pT_pTdiff_selection", 
                                            "pair_pT_pTsum_selection",
                                            "pair_pT_random_selection"};
    std::vector<std::string> phi_filenames = {"phi_diff_mass_selection", 
                                             "phi_diff_pTdiff_selection", 
                                             "phi_diff_pTsum_selection",
                                             "phi_diff_random_selection"};
    std::vector<std::string> cos_theta1_filenames = {"cos_theta1_mass_selection", 
                                                    "cos_theta1_pTdiff_selection", 
                                                    "cos_theta1_pTsum_selection",
                                                    "cos_theta1_random_selection"};
    std::vector<std::string> cos_theta2_filenames = {"cos_theta2_mass_selection", 
                                                    "cos_theta2_pTdiff_selection", 
                                                    "cos_theta2_pTsum_selection",
                                                    "cos_theta2_random_selection"};
    std::vector<std::string> phi_prime_filenames = {"phi_prime_mass_selection", 
                                                   "phi_prime_pTdiff_selection", 
                                                   "phi_prime_pTsum_selection",
                                                   "phi_prime_random_selection"};
    
    // Helper function to create combined PDF for each method
    auto createCombinedPDF = [&](size_t method_idx) {
        std::string pdf_filename = "combined_plots_" + methods[method_idx].name + ".pdf";
        std::string pdf_open = pdf_filename + "(";
        std::string pdf_close = pdf_filename + ")";
        
        // Create combined canvas for method overview
        TCanvas *combined_canvas = new TCanvas(("combined_" + methods[method_idx].name).c_str(), 
                                             ("Combined Plots - " + method_titles[method_idx] + " Selection").c_str(), 
                                             1200, 800);
        combined_canvas->Divide(3, 2); // 3x2 grid for 6 plots
        
        // Page 1: Overview with all plots on one page
        combined_canvas->cd(1);
        mass_hists[method_idx][0]->SetLineColor(kBlue);
        mass_hists[method_idx][1]->SetLineColor(kRed);
        float max_mass = std::max(mass_hists[method_idx][0]->GetMaximum(), mass_hists[method_idx][1]->GetMaximum());
        mass_hists[method_idx][0]->SetMaximum(1.2 * max_mass);
        mass_hists[method_idx][0]->SetTitle("Invariant Mass");
        mass_hists[method_idx][0]->Draw();
        mass_hists[method_idx][1]->Draw("SAME");
        
        combined_canvas->cd(2);
        pT_hists[method_idx][0]->SetLineColor(kBlue);
        pT_hists[method_idx][1]->SetLineColor(kRed);
        float max_pT = std::max(pT_hists[method_idx][0]->GetMaximum(), pT_hists[method_idx][1]->GetMaximum());
        pT_hists[method_idx][0]->SetMaximum(1.2 * max_pT);
        pT_hists[method_idx][0]->SetTitle("Pair pT");
        pT_hists[method_idx][0]->Draw();
        pT_hists[method_idx][1]->Draw("SAME");
        
        combined_canvas->cd(3);
        phi_hists[method_idx][0]->SetLineColor(kBlue);
        phi_hists[method_idx][1]->SetLineColor(kRed);
        float max_phi = std::max(phi_hists[method_idx][0]->GetMaximum(), phi_hists[method_idx][1]->GetMaximum());
        phi_hists[method_idx][0]->SetMaximum(1.2 * max_phi);
        phi_hists[method_idx][0]->SetTitle("Phi Diff Within Pairs");
        phi_hists[method_idx][0]->Draw();
        phi_hists[method_idx][1]->Draw("SAME");
        
        combined_canvas->cd(4);
        cos_theta1_hists[method_idx]->SetLineColor(kBlue);
        cos_theta1_hists[method_idx]->SetTitle("cos(#theta_{1}) Pair vs Z in 4#pi Rest");
        cos_theta1_hists[method_idx]->Draw();
        
        combined_canvas->cd(5);
        cos_theta2_hists[method_idx][0]->SetLineColor(kBlue);
        cos_theta2_hists[method_idx][1]->SetLineColor(kRed);
        float max_theta2 = std::max(cos_theta2_hists[method_idx][0]->GetMaximum(), cos_theta2_hists[method_idx][1]->GetMaximum());
        cos_theta2_hists[method_idx][0]->SetMaximum(1.2 * max_theta2);
        cos_theta2_hists[method_idx][0]->SetTitle("cos(#theta_{2}) Pion vs Z in Pair Rest");
        cos_theta2_hists[method_idx][0]->Draw();
        cos_theta2_hists[method_idx][1]->Draw("SAME");
        
        combined_canvas->cd(6);
        phi_prime_hists[method_idx]->SetLineColor(kBlue);
        phi_prime_hists[method_idx]->SetTitle("Phi Prime - Angle Between Pair pT");
        phi_prime_hists[method_idx]->Draw();
        
        // Add legend to the combined plot
        combined_canvas->cd();
        auto legend = new TLegend(0.02, 0.02, 0.18, 0.12);
        legend->AddEntry(mass_hists[method_idx][0], "Pair A", "l");
        legend->AddEntry(mass_hists[method_idx][1], "Pair B", "l");
        legend->Draw();
        
        // Save overview page
        combined_canvas->SaveAs(pdf_open.c_str());
        
        // Page 2-6: Individual detailed plots
        // Mass plots
        TCanvas *mass_canvas = new TCanvas("mass_detailed", "Mass Distribution", 800, 600);
        mass_hists[method_idx][0]->SetLineColor(kBlue);
        mass_hists[method_idx][1]->SetLineColor(kRed);
        mass_hists[method_idx][0]->SetTitle(("Invariant Mass - " + method_titles[method_idx] + " Selection").c_str());
        mass_hists[method_idx][0]->Draw();
        mass_hists[method_idx][1]->Draw("SAME");
        auto leg1 = new TLegend(0.7,0.7,0.9,0.9);
        leg1->AddEntry(mass_hists[method_idx][0], "Pair A", "l");
        leg1->AddEntry(mass_hists[method_idx][1], "Pair B", "l");
        leg1->Draw();
        mass_canvas->SaveAs(pdf_filename.c_str());
        
        // pT plots
        TCanvas *pT_canvas = new TCanvas("pT_detailed", "pT Distribution", 800, 600);
        pT_hists[method_idx][0]->SetLineColor(kBlue);
        pT_hists[method_idx][1]->SetLineColor(kRed);
        pT_hists[method_idx][0]->SetTitle(("Pair pT - " + method_titles[method_idx] + " Selection").c_str());
        pT_hists[method_idx][0]->Draw();
        pT_hists[method_idx][1]->Draw("SAME");
        auto leg2 = new TLegend(0.7,0.7,0.9,0.9);
        leg2->AddEntry(pT_hists[method_idx][0], "Pair A", "l");
        leg2->AddEntry(pT_hists[method_idx][1], "Pair B", "l");
        leg2->Draw();
        pT_canvas->SaveAs(pdf_filename.c_str());
        
        // Phi difference plots
        TCanvas *phi_canvas = new TCanvas("phi_detailed", "Phi Difference", 800, 600);
        phi_hists[method_idx][0]->SetLineColor(kBlue);
        phi_hists[method_idx][1]->SetLineColor(kRed);
        phi_hists[method_idx][0]->SetTitle(("Phi Difference Within Pairs - " + method_titles[method_idx] + " Selection").c_str());
        phi_hists[method_idx][0]->Draw();
        phi_hists[method_idx][1]->Draw("SAME");
        auto leg3 = new TLegend(0.7,0.7,0.9,0.9);
        leg3->AddEntry(phi_hists[method_idx][0], "Pair A", "l");
        leg3->AddEntry(phi_hists[method_idx][1], "Pair B", "l");
        leg3->Draw();
        phi_canvas->SaveAs(pdf_filename.c_str());
        
        // cos(theta1) plot
        TCanvas *theta1_canvas = new TCanvas("theta1_detailed", "Theta1 Polarization", 800, 600);
        cos_theta1_hists[method_idx]->SetLineColor(kBlue);
        cos_theta1_hists[method_idx]->SetTitle(("cos(#theta_{1}) Pair Momentum vs Z-axis in 4#pi Rest Frame - " + method_titles[method_idx] + " Selection").c_str());
        cos_theta1_hists[method_idx]->Draw();
        theta1_canvas->SaveAs(pdf_filename.c_str());
        
        // cos(theta2) plot
        TCanvas *theta2_canvas = new TCanvas("theta2_detailed", "Theta2 Polarization", 800, 600);
        cos_theta2_hists[method_idx][0]->SetLineColor(kBlue);
        cos_theta2_hists[method_idx][1]->SetLineColor(kRed);
        cos_theta2_hists[method_idx][0]->SetTitle(("cos(#theta_{2}) Pion Momentum vs Z-axis in Pair Rest Frame - " + method_titles[method_idx] + " Selection").c_str());
        cos_theta2_hists[method_idx][0]->Draw();
        cos_theta2_hists[method_idx][1]->Draw("SAME");
        auto leg4 = new TLegend(0.7,0.7,0.9,0.9);
        leg4->AddEntry(cos_theta2_hists[method_idx][0], "Pair A", "l");
        leg4->AddEntry(cos_theta2_hists[method_idx][1], "Pair B", "l");
        leg4->Draw();
        theta2_canvas->SaveAs(pdf_filename.c_str());
        
        // Phi prime plot
        TCanvas *phi_prime_canvas = new TCanvas("phi_prime_detailed", "Phi Prime", 800, 600);
        phi_prime_hists[method_idx]->SetLineColor(kBlue);
        phi_prime_hists[method_idx]->SetTitle(("Phi Prime - Angle Between Two Pair pT Vectors - " + method_titles[method_idx] + " Selection").c_str());
        phi_prime_hists[method_idx]->Draw();
        phi_prime_canvas->SaveAs(pdf_filename.c_str());
        
        // Vector phi angles plot
        if (save_vec_angles && vec_phi_hists[method_idx].size() >= 2) {
            TCanvas *vec_phi_canvas = new TCanvas("vec_phi_detailed", "Vector Phi Angles", 800, 600);
            vec_phi_hists[method_idx][0]->SetLineColor(kBlue);
            vec_phi_hists[method_idx][1]->SetLineColor(kRed);
            float max_vec_phi = std::max(vec_phi_hists[method_idx][0]->GetMaximum(), vec_phi_hists[method_idx][1]->GetMaximum());
            vec_phi_hists[method_idx][0]->SetMaximum(1.2 * max_vec_phi);
            vec_phi_hists[method_idx][0]->SetTitle(("Phi Angles of Relative Momentum Vectors - " + method_titles[method_idx] + " Selection").c_str());
            vec_phi_hists[method_idx][0]->Draw();
            vec_phi_hists[method_idx][1]->Draw("SAME");
            auto leg_vphi = new TLegend(0.7,0.7,0.9,0.9);
            leg_vphi->AddEntry(vec_phi_hists[method_idx][0], "Pair A", "l");
            leg_vphi->AddEntry(vec_phi_hists[method_idx][1], "Pair B", "l");
            leg_vphi->Draw();
            vec_phi_canvas->SaveAs(pdf_filename.c_str());
            delete vec_phi_canvas;
        }
        
        // Vector theta angles plot
        if (save_vec_angles && vec_theta_hists[method_idx].size() >= 2) {
            TCanvas *vec_theta_canvas = new TCanvas("vec_theta_detailed", "Vector Theta Angles", 800, 600);
            vec_theta_hists[method_idx][0]->SetLineColor(kBlue);
            vec_theta_hists[method_idx][1]->SetLineColor(kRed);
            float max_vec_theta = std::max(vec_theta_hists[method_idx][0]->GetMaximum(), vec_theta_hists[method_idx][1]->GetMaximum());
            vec_theta_hists[method_idx][0]->SetMaximum(1.2 * max_vec_theta);
            vec_theta_hists[method_idx][0]->SetTitle(("Theta Angles of Relative Momentum Vectors - " + method_titles[method_idx] + " Selection").c_str());
            vec_theta_hists[method_idx][0]->Draw();
            vec_theta_hists[method_idx][1]->Draw("SAME");
            auto leg_vtheta = new TLegend(0.7,0.7,0.9,0.9);
            leg_vtheta->AddEntry(vec_theta_hists[method_idx][0], "Pair A", "l");
            leg_vtheta->AddEntry(vec_theta_hists[method_idx][1], "Pair B", "l");
            leg_vtheta->Draw();
            vec_theta_canvas->SaveAs(pdf_filename.c_str());
            delete vec_theta_canvas;
        }
        
        // Delta phi distributions plot
        if (save_delta_distributions && delta_phi_hists[method_idx].size() >= 3) {
            TCanvas *delta_phi_canvas = new TCanvas("delta_phi_detailed", "Delta Phi Distributions", 800, 600);
            delta_phi_hists[method_idx][0]->SetLineColor(kBlue);   // pairAB
            delta_phi_hists[method_idx][1]->SetLineColor(kRed);    // pairAT
            delta_phi_hists[method_idx][2]->SetLineColor(kGreen);  // pairBT
            float max_dphi = std::max({delta_phi_hists[method_idx][0]->GetMaximum(), 
                                      delta_phi_hists[method_idx][1]->GetMaximum(),
                                      delta_phi_hists[method_idx][2]->GetMaximum()});
            delta_phi_hists[method_idx][0]->SetMaximum(1.2 * max_dphi);
            delta_phi_hists[method_idx][0]->SetTitle(("Delta Phi Between Momentum Vectors - " + method_titles[method_idx] + " Selection").c_str());
            delta_phi_hists[method_idx][0]->Draw();
            delta_phi_hists[method_idx][1]->Draw("SAME");
            delta_phi_hists[method_idx][2]->Draw("SAME");
            auto leg_dphi = new TLegend(0.7,0.6,0.9,0.9);
            leg_dphi->AddEntry(delta_phi_hists[method_idx][0], "pairAB", "l");
            leg_dphi->AddEntry(delta_phi_hists[method_idx][1], "pairAT", "l");
            leg_dphi->AddEntry(delta_phi_hists[method_idx][2], "pairBT", "l");
            leg_dphi->Draw();
            delta_phi_canvas->SaveAs(pdf_filename.c_str());
            delete delta_phi_canvas;
        }
        
        // Delta theta distributions plot (final page)
        if (save_delta_distributions && delta_theta_hists[method_idx].size() >= 3) {
            TCanvas *delta_theta_canvas = new TCanvas("delta_theta_detailed", "Delta Theta Distributions", 800, 600);
            delta_theta_hists[method_idx][0]->SetLineColor(kBlue);   // pairAB
            delta_theta_hists[method_idx][1]->SetLineColor(kRed);    // pairAT
            delta_theta_hists[method_idx][2]->SetLineColor(kGreen);  // pairBT
            float max_dtheta = std::max({delta_theta_hists[method_idx][0]->GetMaximum(), 
                                        delta_theta_hists[method_idx][1]->GetMaximum(),
                                        delta_theta_hists[method_idx][2]->GetMaximum()});
            delta_theta_hists[method_idx][0]->SetMaximum(1.2 * max_dtheta);
            delta_theta_hists[method_idx][0]->SetTitle(("Delta Theta Between Momentum Vectors - " + method_titles[method_idx] + " Selection").c_str());
            delta_theta_hists[method_idx][0]->Draw();
            delta_theta_hists[method_idx][1]->Draw("SAME");
            delta_theta_hists[method_idx][2]->Draw("SAME");
            auto leg_dtheta = new TLegend(0.7,0.6,0.9,0.9);
            leg_dtheta->AddEntry(delta_theta_hists[method_idx][0], "pairAB", "l");
            leg_dtheta->AddEntry(delta_theta_hists[method_idx][1], "pairAT", "l");
            leg_dtheta->AddEntry(delta_theta_hists[method_idx][2], "pairBT", "l");
            leg_dtheta->Draw();
            delta_theta_canvas->SaveAs(pdf_close.c_str());  // Close PDF with final page
            delete delta_theta_canvas;
        } else {
            // If no delta distributions, close PDF here
            TCanvas *dummy_canvas = new TCanvas("dummy", "dummy", 800, 600);
            dummy_canvas->SaveAs(pdf_close.c_str());
            delete dummy_canvas;
        }
        
        // Clean up canvases
        delete combined_canvas;
        delete mass_canvas;
        delete pT_canvas;
        delete phi_canvas;
        delete theta1_canvas;
        delete theta2_canvas;
        delete phi_prime_canvas;
        
        std::cout << "Created combined PDF: " << pdf_filename << std::endl;
    };

    // Create combined PDFs and PNG plots for each method
    // Helper to check if histogram contains data
    auto histNonEmpty = [&](TH1F* h) {
        return (h && h->GetEntries() > 0);
    };

    for (size_t i = 0; i < methods.size(); ++i) {
        // Determine whether there is anything to plot for this method, respecting method_enabled and save flags
        bool anyToPlot = false;
        if (method_enabled[i]) {
            // mass / pT always present when method enabled
            if (histNonEmpty(mass_hists[i][0]) || histNonEmpty(mass_hists[i][1]) ||
                histNonEmpty(pT_hists[i][0]) || histNonEmpty(pT_hists[i][1])) anyToPlot = true;
            if (save_phi_diff && (histNonEmpty(phi_hists[i][0]) || histNonEmpty(phi_hists[i][1]))) anyToPlot = true;
            if (save_vec_angles) {
                if (histNonEmpty(vec_phi_hists[i][0]) || histNonEmpty(vec_phi_hists[i][1]) ||
                    histNonEmpty(delta_phi_hists[i][0]) || histNonEmpty(delta_phi_hists[i][1]) || histNonEmpty(delta_phi_hists[i][2])) anyToPlot = true;
            }
            if (save_cos_theta1 && histNonEmpty(cos_theta1_hists[i])) anyToPlot = true;
            if (save_cos_theta2 && (histNonEmpty(cos_theta2_hists[i][0]) || histNonEmpty(cos_theta2_hists[i][1]))) anyToPlot = true;
            if (save_phi_prime && histNonEmpty(phi_prime_hists[i])) anyToPlot = true;
        }

        if (generate_pdf && anyToPlot) {
            createCombinedPDF(i);
        }

        if ((generate_png || generate_pdf) && anyToPlot) {
            // Create individual plots (PNG and PDF) only if corresponding histograms contain data
            if (histNonEmpty(mass_hists[i][0]) || histNonEmpty(mass_hists[i][1]))
                createPlot(mass_hists[i], "Invariant Mass Pair Histograms (" + method_titles[i] + " Selection)", mass_filenames[i]);

            if (histNonEmpty(pT_hists[i][0]) || histNonEmpty(pT_hists[i][1]))
                createPlot(pT_hists[i], "Pair pT Histograms (" + method_titles[i] + " Selection)", pT_filenames[i]);

            if (save_phi_diff && (histNonEmpty(phi_hists[i][0]) || histNonEmpty(phi_hists[i][1])))
                createPlot(phi_hists[i], "Phi Difference Within Pairs (" + method_titles[i] + " Selection)", phi_filenames[i]);

            if (save_cos_theta1 && histNonEmpty(cos_theta1_hists[i]))
                createSinglePlot(cos_theta1_hists[i], "cos(#theta_{1}) Pair Momentum vs Z-axis in 4#pi Rest Frame (" + method_titles[i] + " Selection)", cos_theta1_filenames[i]);

            if (save_cos_theta2 && (histNonEmpty(cos_theta2_hists[i][0]) || histNonEmpty(cos_theta2_hists[i][1])))
                createPlot(cos_theta2_hists[i], "cos(#theta_{2}) Pion Momentum vs Z-axis in Pair Rest Frame (" + method_titles[i] + " Selection)", cos_theta2_filenames[i]);

            if (save_phi_prime && histNonEmpty(phi_prime_hists[i]))
                createSinglePlot(phi_prime_hists[i], "Phi Prime - Angle Between Two Pair pT Vectors (" + method_titles[i] + " Selection)", phi_prime_filenames[i]);
        }
    }

    // Create individual plots for four-pion kinematics (PNG and PDF)
    if ((generate_png || generate_pdf) && save_four_vectors) {
        auto histNonEmpty = [&](TH1F* h) { return (h && h->GetEntries() > 0); };
        if (histNonEmpty(h_four_pi_mass) || histNonEmpty(h_four_pi_pT) || histNonEmpty(h_individual_pT)) {
            if (histNonEmpty(h_four_pi_mass))
                createSinglePlot(h_four_pi_mass, "Four-Pion Mass Distribution", "four_pi_mass_distribution");
            if (histNonEmpty(h_four_pi_pT))
                createSinglePlot(h_four_pi_pT, "Four-Pion pT Distribution", "four_pi_pT_distribution");
            if (histNonEmpty(h_individual_pT))
                createSinglePlot(h_individual_pT, "Individual Pion pT Distribution", "individual_pion_pT_distribution");
        }
    }

    // Create comprehensive overall PDF with four-pion kinematics
    if (generate_pdf && save_four_vectors) {
        auto histNonEmpty = [&](TH1F* h) { return (h && h->GetEntries() > 0); };
        if (histNonEmpty(h_four_pi_mass) || histNonEmpty(h_four_pi_pT) || histNonEmpty(h_individual_pT)) {
            std::string overall_pdf = "overall_four_pion_analysis.pdf";
            std::string pdf_open = overall_pdf + "(";
            std::string pdf_close = overall_pdf + ")";
            
            // Create combined overview of four-pion analysis
            TCanvas *overall_canvas = new TCanvas("overall_analysis", "Four-Pion Analysis Overview", 1200, 800);
            overall_canvas->Divide(2, 2); // 2x2 grid
            
            if (histNonEmpty(h_four_pi_mass)) {
                overall_canvas->cd(1);
                h_four_pi_mass->SetLineColor(kBlue);
                h_four_pi_mass->Draw();
            }
            
            if (histNonEmpty(h_four_pi_pT)) {
                overall_canvas->cd(2);
                h_four_pi_pT->SetLineColor(kRed);
                h_four_pi_pT->Draw();
            }
            
            if (histNonEmpty(h_individual_pT)) {
                overall_canvas->cd(3);
                h_individual_pT->SetLineColor(kGreen);
                h_individual_pT->Draw();
            }
            
            overall_canvas->SaveAs(pdf_open.c_str());
            
            // Individual detailed plots
            if (histNonEmpty(h_four_pi_mass)) {
                TCanvas *mass_canvas = new TCanvas("four_pi_mass", "Four-Pion Mass", 800, 600);
                h_four_pi_mass->SetLineColor(kBlue);
                h_four_pi_mass->Draw();
                mass_canvas->SaveAs(overall_pdf.c_str());
                delete mass_canvas;
            }
            
            if (histNonEmpty(h_four_pi_pT)) {
                TCanvas *pT_canvas = new TCanvas("four_pi_pT", "Four-Pion pT", 800, 600);
                h_four_pi_pT->SetLineColor(kRed);
                h_four_pi_pT->Draw();
                pT_canvas->SaveAs(overall_pdf.c_str());
                delete pT_canvas;
            }
            
            if (histNonEmpty(h_individual_pT)) {
                TCanvas *indiv_canvas = new TCanvas("individual_pT", "Individual Pion pT", 800, 600);
                h_individual_pT->SetLineColor(kGreen);
                h_individual_pT->Draw();
                indiv_canvas->SaveAs(pdf_close.c_str());
                delete indiv_canvas;
            } else {
                // Close PDF if individual pT not available
                TCanvas *dummy_canvas = new TCanvas("dummy", "dummy", 800, 600);
                dummy_canvas->SaveAs(pdf_close.c_str());
                delete dummy_canvas;
            }
            
            delete overall_canvas;
            std::cout << "Created overall analysis PDF: " << overall_pdf << std::endl;
        }
    }

    // Save histograms to ROOT file
    std::cout << "Saving histograms to ROOT file..." << std::endl;
    fplot->cd();
    
    // Save four-pion kinematics histograms
    h_four_pi_mass->Write();
    h_four_pi_pT->Write();
    h_individual_pT->Write();
    
    // Save separate histograms (for plotting and analysis)
    for (size_t i = 0; i < methods.size(); ++i) {
        // Write relative-vector histograms and delta histograms (only if requested)
        if (save_vec_angles) {
            for (int j = 0; j < 2; ++j) {
                vec_phi_hists[i][j]->Write();
                vec_theta_hists[i][j]->Write();
            }
            for (int j = 0; j < 3; ++j) {
                delta_phi_hists[i][j]->Write();
                delta_theta_hists[i][j]->Write();
            }
        }
    if (save_phi_prime) phi_prime_hists[i]->Write();
        if (save_cos_theta1) cos_theta1_hists[i]->Write();  // Pair polarization (single histogram)
        for (int j = 0; j < 2; ++j) {
            mass_hists[i][j]->Write();
            pT_hists[i][j]->Write();
            if (save_phi_diff) phi_hists[i][j]->Write();
            if (save_cos_theta2) cos_theta2_hists[i][j]->Write();  // Single pion polarization
        }
    }

    fplot->Close();
    std::cout << "Analysis complete! All histograms saved to ROOT file." << std::endl;
    
    // Report what was generated based on flags
    std::string output_summary = "Created: ";
    if (generate_pdf && generate_png) {
        output_summary += "4 combined PDF files (one per method) + 20 individual PNG plots + 3 four-pion kinematic plots";
    } else if (generate_pdf && !generate_png) {
        output_summary += "4 combined PDF files (one per method) only";
    } else if (!generate_pdf && generate_png) {
        output_summary += "20 individual PNG plots + 3 four-pion kinematic plots only";
    } else {
        output_summary += "ROOT histograms only (no plots generated)";
    }
    std::cout << output_summary << std::endl;
    
    std::cout << "theta1: Pair polarization (pair momentum vs z-axis in 4π rest frame)" << std::endl;
    std::cout << "theta2: Single pion polarization (pion momentum vs z-axis in pair rest frame)" << std::endl;
}
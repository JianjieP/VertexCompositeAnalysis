#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TText.h"
#include "TF1.h"
#include "TMath.h"
#include "TString.h"
#include "TROOT.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

void samplaAna_general(const char* inputFile = "/eos/user/j/jianjie/fourpi/HIRun2025A_FourPi_0410/chained_ParticleTree_merged.root") {
    // Set ROOT to batch mode to avoid GUI issues
    gROOT->SetBatch(kTRUE);

    // Feature flags: keep filter logic on, turn off the non-filter outputs by default.
    const bool kMakeFilteredTree = true;
    const bool kMakeGraphs = false;
    
    // Define histogram bin numbers
    const Int_t nBins_mass = 300;
    const Int_t nBins_pt = 150;
    const Int_t nBins_pt2 = 150;
    const Int_t nBins_y = 150;
    const Int_t nBins_eta =150;
    const Int_t nBins_phi = 150;
    const Int_t nBins_cuts = 9;  // Increased to 9 for new filter
    
    // Define histogram ranges
    float pt_min = 0.0;
    float pt_max = 1.0;
    float pt2_min = 0.0;
    float pt2_max = 1.0;  // pt^2 range (GeV/c)^2
    float mass_min = 0.8;
    float mass_max = 8.0;
    cout << "Opening input file: " << inputFile << endl;
    
    // Open the chained ROOT file
    TFile* file = TFile::Open(inputFile, "READ");
    if (!file || file->IsZombie()) {
        cout << "Error: Could not open file " << inputFile << endl;
        return;
    }
    
    // Get the tree from the file
    TTree* tree = (TTree*)file->Get("ParticleTree");
    if (!tree) {
        cout << "Error: Could not find ParticleTree in file" << endl;
        file->Close();
        return;
    }
    
    cout << "Loaded tree with " << tree->GetEntries() << " entries" << endl;
    
    // Apply filtering conditions and fill histogram manually
    cout << "Applying filtering conditions and creating histograms for data extraction..." << endl;
    
    // Set up branches
    unsigned int RunNb;
    float TWHFmaxEMinus, TWHFmaxEPlus;
    std::vector<bool>* passHLT = nullptr;
    std::vector<bool>* evtSel = nullptr;
    std::vector<bool>* trk_isHP = nullptr;
    std::vector<float>* trk_xyDCASignificance = nullptr;
    std::vector<float>* trk_zDCASignificance = nullptr;
    std::vector<float>* cand_mass = nullptr;
    std::vector<float>* cand_pT = nullptr;
    std::vector<float>* cand_y = nullptr;
    std::vector<float>* cand_eta = nullptr;
    std::vector<float>* cand_phi = nullptr;
    
    tree->SetBranchAddress("RunNb", &RunNb);
    tree->SetBranchAddress("TWHFmaxEMinus", &TWHFmaxEMinus);
    tree->SetBranchAddress("TWHFmaxEPlus", &TWHFmaxEPlus);
    tree->SetBranchAddress("passHLT", &passHLT);
    tree->SetBranchAddress("evtSel", &evtSel);
    tree->SetBranchAddress("trk_isHP", &trk_isHP);
    tree->SetBranchAddress("trk_xyDCASignificance", &trk_xyDCASignificance);
    tree->SetBranchAddress("trk_zDCASignificance", &trk_zDCASignificance);
    tree->SetBranchAddress("cand_mass", &cand_mass);
    tree->SetBranchAddress("cand_pT", &cand_pT);
    tree->SetBranchAddress("cand_y", &cand_y);
    tree->SetBranchAddress("cand_eta", &cand_eta);
    tree->SetBranchAddress("cand_phi", &cand_phi);

    // Define filtering parameters
    unsigned int runNb_min = 399465;
    unsigned int runNb_max = 400426;
    float twhfMinus_max = 8.6;
    float twhfPlus_max = 9.2;
    float dca_significance_max = 3.0;
    float candidate_pt_min = 0.2;
    

    // Create histograms
    TH1F* hist_mass = new TH1F("hist_mass", "Mass (p_{T} < 0.2 GeV/c);Mass (GeV/c^{2});Raw Counts", nBins_mass, mass_min, mass_max);
    hist_mass->SetDirectory(0);  // Prevent ROOT from auto-managing this histogram
    
    // Create additional kinematic histograms
    TH1F* hist_pt = new TH1F("hist_pt", "p_{T};p_{T} (GeV/c);Raw Counts", nBins_pt, pt_min, pt_max);
    hist_pt->SetDirectory(0);
    TH1F* hist_pt_normalized = new TH1F("hist_pt_normalized", "N/p_{T} vs p_{T};p_{T} (GeV/c);N/p_{T} (GeV/c)^{-1}", nBins_pt, pt_min, pt_max);
    hist_pt_normalized->SetDirectory(0);
    TH1F* hist_pt2 = new TH1F("hist_pt2", "p_{T}^{2};p_{T}^{2} ((GeV/c)^{2});Raw Counts", nBins_pt2, pt2_min, pt2_max);
    hist_pt2->SetDirectory(0);
    TH1F* hist_y = new TH1F("hist_y", "Rapidity;y;Raw Counts", nBins_y, -3, 3);
    hist_y->SetDirectory(0);
    TH1F* hist_eta = new TH1F("hist_eta", "Pseudorapidity;#eta;Raw Counts", nBins_eta, -3, 3);
    hist_eta->SetDirectory(0);
    TH1F* hist_phi = new TH1F("hist_phi", "#phi;#phi;Raw Counts", nBins_phi, -TMath::Pi(), TMath::Pi());
    hist_phi->SetDirectory(0);
    
    // Create filtering steps histogram
    TH1F* hist_cuts = new TH1F("hist_cuts", "Event Selection Steps;Selection Step;Number of Events", nBins_cuts, 0, nBins_cuts);
    hist_cuts->SetDirectory(0);  // Prevent ROOT from auto-managing this histogram
    hist_cuts->GetXaxis()->SetBinLabel(1, "Total Events");
    hist_cuts->GetXaxis()->SetBinLabel(2, "Run Number");
    hist_cuts->GetXaxis()->SetBinLabel(3, "TWHFmaxEMinus");
    hist_cuts->GetXaxis()->SetBinLabel(4, "TWHFmaxEPlus");
    hist_cuts->GetXaxis()->SetBinLabel(5, "passHLT");
    hist_cuts->GetXaxis()->SetBinLabel(6, "Event Selection");
    hist_cuts->GetXaxis()->SetBinLabel(7, "Track Quality");
    hist_cuts->GetXaxis()->SetBinLabel(8, "Track DCA");
    hist_cuts->GetXaxis()->SetBinLabel(9, "Track p_{T} > 0.2");
    
    // Counters for filtering steps
    Long64_t count_total = 0;
    Long64_t count_runNb = 0;
    Long64_t count_twhfMinus = 0;
    Long64_t count_twhfPlus = 0;
    Long64_t count_passHLT = 0;
    Long64_t count_evtSel = 0;
    Long64_t count_trackQuality = 0;
    Long64_t count_trackDCA = 0;
    Long64_t count_final = 0;
    
    Long64_t nEntries = tree->GetEntries();
    cout << "Processing " << nEntries << " entries..." << endl;
    
    // Prepare output file and filtered tree only when requested
    TFile* filteredFile = nullptr;
    TTree* filteredTree = nullptr;
    if (kMakeFilteredTree) {
        TString filteredOutput = "/eos/user/j/jianjie/fourpi/HIRun2025A_FourPi_0410/filtered_chained_ParticleTree.root";
        filteredFile = new TFile(filteredOutput, "RECREATE");
        if (!filteredFile || filteredFile->IsZombie()) {
            cout << "Error: Could not create filtered output file " << filteredOutput << endl;
            return;
        }
        filteredTree = tree->CloneTree(0); // Empty tree with same structure
    }

    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % 100000 == 0) {  // Progress indicator
            cout << "Processing entry " << i << "/" << nEntries << " (" << (100.0*i/nEntries) << "%)" << endl;
        }
        
        tree->GetEntry(i);
        count_total++;
        
        // Event-level cuts
        if (!(RunNb >= runNb_min && RunNb <= runNb_max)) continue;
        count_runNb++;
        
        if (!(TWHFmaxEMinus <= twhfMinus_max)) continue;
        count_twhfMinus++;
        
        if (!(TWHFmaxEPlus <= twhfPlus_max)) continue;
        count_twhfPlus++;
        
        // passHLT is a vector<bool>, require all true and at least one element
        if (!passHLT || passHLT->empty()) continue;
        bool allHLTPass = true;
        for (size_t hlt_i = 0; hlt_i < passHLT->size(); ++hlt_i) {
            if (!(*passHLT)[hlt_i]) {
                allHLTPass = false;
                break;
            }
        }
        if (!allHLTPass) continue;
        count_passHLT++;
        
        // evtSel is a vector<bool> for primaryVertex&clusterCompatibility, require all true and at least one element
        if (!evtSel || evtSel->empty()) continue;
        bool allEvtSelPass = true;
        for (size_t evt_i = 0; evt_i < evtSel->size(); ++evt_i) {
            if (!(*evtSel)[evt_i]) {
                allEvtSelPass = false;
                break;
            }
        }
        if (!allEvtSelPass) continue;
        count_evtSel++;
        
        // Track-level cuts: all tracks must pass
        if (!trk_isHP || !trk_xyDCASignificance || !trk_zDCASignificance) continue;
        if (trk_isHP->size() != trk_xyDCASignificance->size() || trk_isHP->size() != trk_zDCASignificance->size()) continue;
        count_trackQuality++;
        
        bool allTracksPass = true;
        for (size_t j = 0; j < trk_isHP->size(); ++j) {
            if (!(*trk_isHP)[j]) { allTracksPass = false; break; }
            if (std::abs((*trk_xyDCASignificance)[j]) >= dca_significance_max) { allTracksPass = false; break; }
            if (std::abs((*trk_zDCASignificance)[j]) >= dca_significance_max) { allTracksPass = false; break; }
        }
        if (!allTracksPass) continue;
        count_trackDCA++;
        
        // cand_pT filter: check that the last four members have pT > candidate_pt_min
        if (!cand_pT || cand_pT->size() < 4) continue;
        bool lastFourPtPass = true;
        for (size_t k = cand_pT->size() - 4; k < cand_pT->size(); ++k) {
            if ((*cand_pT)[k] <= candidate_pt_min) {
                lastFourPtPass = false;
                break;
            }
        }
        if (!lastFourPtPass) continue;
        count_final++;

        // Fill the filtered tree with this passing event
        if (kMakeFilteredTree && filteredTree) {
            filteredTree->Fill();
        }

        // Fill kinematic histograms using first candidate
        // Apply pT < 0.2 restriction for mass histogram only
        if (cand_mass && !cand_mass->empty() && cand_pT && !cand_pT->empty() && (*cand_pT)[0] < 0.2) {
            hist_mass->Fill((*cand_mass)[0]);
        }
        if (cand_pT && !cand_pT->empty()) {
            hist_pt->Fill((*cand_pT)[0]);
            // Fill N/pT histogram - avoid division by zero
            if ((*cand_pT)[0] > 0) {
                hist_pt_normalized->Fill((*cand_pT)[0], 1.0/(*cand_pT)[0]);
            }
            hist_pt2->Fill((*cand_pT)[0] * (*cand_pT)[0]);  // Fill pT^2
        }
        if (cand_y && !cand_y->empty()) hist_y->Fill((*cand_y)[0]);
        if (cand_eta && !cand_eta->empty()) hist_eta->Fill((*cand_eta)[0]);
        if (cand_phi && !cand_phi->empty()) hist_phi->Fill((*cand_phi)[0]);
    }
    
    // Fill the filtering steps histogram
    hist_cuts->SetBinContent(1, count_total);
    hist_cuts->SetBinContent(2, count_runNb);
    hist_cuts->SetBinContent(3, count_twhfMinus);
    hist_cuts->SetBinContent(4, count_twhfPlus);
    hist_cuts->SetBinContent(5, count_passHLT);
    hist_cuts->SetBinContent(6, count_evtSel);
    hist_cuts->SetBinContent(7, count_trackQuality);
    hist_cuts->SetBinContent(8, count_trackDCA);
    hist_cuts->SetBinContent(9, count_final);
    
    cout << "Analysis complete!" << endl;
    cout << "Filtering summary:" << endl;
    cout << "  Total events: " << count_total << endl;
    cout << "  After Run Number cut: " << count_runNb << " (" << (100.0*count_runNb/count_total) << "%)" << endl;
    cout << "  After TWHFmaxEMinus cut: " << count_twhfMinus << " (" << (100.0*count_twhfMinus/count_total) << "%)" << endl;
    cout << "  After TWHFmaxEPlus cut: " << count_twhfPlus << " (" << (100.0*count_twhfPlus/count_total) << "%)" << endl;
    cout << "  After passHLT cut: " << count_passHLT << " (" << (100.0*count_passHLT/count_total) << "%)" << endl;
    cout << "  After Event Selection cut: " << count_evtSel << " (" << (100.0*count_evtSel/count_total) << "%)" << endl;
    cout << "  After Track Quality cut: " << count_trackQuality << " (" << (100.0*count_trackQuality/count_total) << "%)" << endl;
    cout << "  After Track DCA cut: " << count_trackDCA << " (" << (100.0*count_trackDCA/count_total) << "%)" << endl;
    cout << "  After Track p_{T} > 0.2 cut: " << count_final << " (" << (100.0*count_final/count_total) << "%)" << endl;

    if (kMakeGraphs && hist_mass->GetEntries() > 0) {
        cout << "Histograms created:" << endl;
        cout << "  Mass histogram entries: " << hist_mass->GetEntries() << endl;
        cout << "  p_{T} histogram entries: " << hist_pt->GetEntries() << endl;
        cout << "  N/p_{T} histogram entries: " << hist_pt_normalized->GetEntries() << endl;
        cout << "  p_{T}^{2} histogram entries: " << hist_pt2->GetEntries() << endl;
        cout << "  Rapidity histogram entries: " << hist_y->GetEntries() << endl;
        cout << "  Pseudorapidity histogram entries: " << hist_eta->GetEntries() << endl;
        cout << "  Phi histogram entries: " << hist_phi->GetEntries() << endl;
        
        // Convert histogram to TGraphErrors (dots with error bars)
        Int_t nBins = hist_mass->GetNbinsX();
        vector<Double_t> x, y, ex, ey;
        
        for (Int_t i = 1; i <= nBins; i++) {
            Double_t binContent = hist_mass->GetBinContent(i);
            Double_t binWidth = hist_mass->GetBinWidth(i);
            if (binContent > 0) {  // Only include bins with data
                Double_t binCenter = hist_mass->GetBinCenter(i);
                Double_t binError = hist_mass->GetBinError(i);
                x.push_back(binCenter);
                y.push_back(binContent); // Raw count
                // y.push_back(binContent / binWidth); // Differential count
                ex.push_back(binWidth/2.0);  // Horizontal error = half bin width
                ey.push_back(binError); // Raw error
                // ey.push_back(binError / binWidth); // Error also divided by bin width
            }
        }
        
        cout << "Creating graph with " << x.size() << " data points" << endl;
        
        // Create TGraphErrors
        TGraphErrors* graph = new TGraphErrors(x.size(), &x[0], &y[0], &ex[0], &ey[0]);
        graph->SetTitle("Mass Distribution (p_{T} < 0.2 GeV/c);Mass (GeV/c^{2});Raw Counts");
        graph->SetMarkerStyle(20);  // Filled circles
        graph->SetMarkerSize(0.8);
        graph->SetMarkerColor(kBlue);
        graph->SetLineColor(kBlue);
        
        // Create TGraphErrors for pT distribution
        Int_t nBins_pt = hist_pt->GetNbinsX();
        vector<Double_t> x_pt, y_pt, ex_pt, ey_pt;
        
        for (Int_t i = 1; i <= nBins_pt; i++) {
            Double_t binContent = hist_pt->GetBinContent(i);
            Double_t binWidth = hist_pt->GetBinWidth(i);
            if (binContent > 0) {  // Only include bins with data
                Double_t binCenter = hist_pt->GetBinCenter(i);
                Double_t binError = hist_pt->GetBinError(i);
                x_pt.push_back(binCenter);
                y_pt.push_back(binContent); // Raw count
                // y_pt.push_back(binContent / binWidth); // Differential count
                ex_pt.push_back(binWidth/2.0);
                ey_pt.push_back(binError); // Raw error
                // ey_pt.push_back(binError / binWidth); // Error also divided by bin width
            }
        }
        
        TGraphErrors* graph_pt = new TGraphErrors(x_pt.size(), &x_pt[0], &y_pt[0], &ex_pt[0], &ey_pt[0]);
        graph_pt->SetTitle("p_{T} Distribution;p_{T} (GeV/c);Raw Counts");
        graph_pt->SetMarkerStyle(20);
        graph_pt->SetMarkerSize(0.8);
        graph_pt->SetMarkerColor(kBlue);
        graph_pt->SetLineColor(kBlue);
        
        // Create TGraphErrors for N/pT distribution
        Int_t nBins_pt_normalized = hist_pt_normalized->GetNbinsX();
        vector<Double_t> x_pt_normalized, y_pt_normalized, ex_pt_normalized, ey_pt_normalized;
        
        for (Int_t i = 1; i <= nBins_pt_normalized; i++) {
            Double_t binContent = hist_pt_normalized->GetBinContent(i);
            Double_t binWidth = hist_pt_normalized->GetBinWidth(i);
            if (binContent > 0) {  // Only include bins with data
                Double_t binCenter = hist_pt_normalized->GetBinCenter(i);
                Double_t binError = hist_pt_normalized->GetBinError(i);
                x_pt_normalized.push_back(binCenter);
                y_pt_normalized.push_back(binContent);
                ex_pt_normalized.push_back(binWidth/2.0);
                ey_pt_normalized.push_back(binError);
            }
        }
        
        TGraphErrors* graph_pt_normalized = new TGraphErrors(x_pt_normalized.size(), &x_pt_normalized[0], &y_pt_normalized[0], &ex_pt_normalized[0], &ey_pt_normalized[0]);
        graph_pt_normalized->SetTitle("N/p_{T} Distribution;p_{T} (GeV/c);N/p_{T} (GeV/c)^{-1}");
        graph_pt_normalized->SetMarkerStyle(21);
        graph_pt_normalized->SetMarkerSize(0.8);
        graph_pt_normalized->SetMarkerColor(kRed);
        graph_pt_normalized->SetLineColor(kRed);
        
        // Create TGraphErrors for pT^2 distribution
        Int_t nBins_pt2 = hist_pt2->GetNbinsX();
        vector<Double_t> x_pt2, y_pt2, ex_pt2, ey_pt2;
        
        for (Int_t i = 1; i <= nBins_pt2; i++) {
            Double_t binContent = hist_pt2->GetBinContent(i);
            Double_t binWidth = hist_pt2->GetBinWidth(i);
            if (binContent > 0) {  // Only include bins with data
                Double_t binCenter = hist_pt2->GetBinCenter(i);
                Double_t binError = hist_pt2->GetBinError(i);
                x_pt2.push_back(binCenter);
                y_pt2.push_back(binContent); // Raw count
                // y_pt2.push_back(binContent / binWidth); // Differential count
                ex_pt2.push_back(binWidth/2.0);
                ey_pt2.push_back(binError); // Raw error
                // ey_pt2.push_back(binError / binWidth); // Error also divided by bin width
            }
        }
        
        TGraphErrors* graph_pt2 = new TGraphErrors(x_pt2.size(), &x_pt2[0], &y_pt2[0], &ex_pt2[0], &ey_pt2[0]);
        graph_pt2->SetTitle("p_{T}^{2} Distribution;p_{T}^{2} ((GeV/c)^{2});Raw Counts");
        graph_pt2->SetMarkerStyle(20);
        graph_pt2->SetMarkerSize(0.8);
        graph_pt2->SetMarkerColor(kCyan);
        graph_pt2->SetLineColor(kCyan);
        
        // Create TGraphErrors for rapidity distribution
        Int_t nBins_y = hist_y->GetNbinsX();
        vector<Double_t> x_y, y_y, ex_y, ey_y;
        
        for (Int_t i = 1; i <= nBins_y; i++) {
            Double_t binContent = hist_y->GetBinContent(i);
            Double_t binWidth = hist_y->GetBinWidth(i);
            if (binContent > 0) {  // Only include bins with data
                Double_t binCenter = hist_y->GetBinCenter(i);
                Double_t binError = hist_y->GetBinError(i);
                x_y.push_back(binCenter);
                y_y.push_back(binContent); // Raw count
                // y_y.push_back(binContent / binWidth); // Differential count
                ex_y.push_back(binWidth/2.0);
                ey_y.push_back(binError); // Raw error
                // ey_y.push_back(binError / binWidth); // Error also divided by bin width
            }
        }
        
        TGraphErrors* graph_y = new TGraphErrors(x_y.size(), &x_y[0], &y_y[0], &ex_y[0], &ey_y[0]);
        graph_y->SetTitle("Rapidity Distribution;y;Raw Counts");
        graph_y->SetMarkerStyle(20);
        graph_y->SetMarkerSize(0.8);
        graph_y->SetMarkerColor(kGreen);
        graph_y->SetLineColor(kGreen);
        
        // Create TGraphErrors for pseudorapidity distribution
        Int_t nBins_eta = hist_eta->GetNbinsX();
        vector<Double_t> x_eta, y_eta, ex_eta, ey_eta;
        
        for (Int_t i = 1; i <= nBins_eta; i++) {
            Double_t binContent = hist_eta->GetBinContent(i);
            Double_t binWidth = hist_eta->GetBinWidth(i);
            if (binContent > 0) {  // Only include bins with data
                Double_t binCenter = hist_eta->GetBinCenter(i);
                Double_t binError = hist_eta->GetBinError(i);
                x_eta.push_back(binCenter);
                y_eta.push_back(binContent); // Raw count
                // y_eta.push_back(binContent / binWidth); // Differential count
                ex_eta.push_back(binWidth/2.0);
                ey_eta.push_back(binError); // Raw error
                // ey_eta.push_back(binError / binWidth); // Error also divided by bin width
            }
        }
        
        TGraphErrors* graph_eta = new TGraphErrors(x_eta.size(), &x_eta[0], &y_eta[0], &ex_eta[0], &ey_eta[0]);
        graph_eta->SetTitle("Pseudorapidity Distribution;#eta;Raw Counts");
        graph_eta->SetMarkerStyle(20);
        graph_eta->SetMarkerSize(0.8);
        graph_eta->SetMarkerColor(kRed);
        graph_eta->SetLineColor(kRed);
        
        // Create TGraphErrors for phi distribution
        Int_t nBins_phi = hist_phi->GetNbinsX();
        vector<Double_t> x_phi, y_phi, ex_phi, ey_phi;
        
        for (Int_t i = 1; i <= nBins_phi; i++) {
            Double_t binContent = hist_phi->GetBinContent(i);
            Double_t binWidth = hist_phi->GetBinWidth(i);
            if (binContent > 0) {  // Only include bins with data
                Double_t binCenter = hist_phi->GetBinCenter(i);
                Double_t binError = hist_phi->GetBinError(i);
                x_phi.push_back(binCenter);
                y_phi.push_back(binContent); // Raw count
                // y_phi.push_back(binContent / binWidth); // Differential count
                ex_phi.push_back(binWidth/2.0);
                ey_phi.push_back(binError); // Raw error
                // ey_phi.push_back(binError / binWidth); // Error also divided by bin width
            }
        }
        
        TGraphErrors* graph_phi = new TGraphErrors(x_phi.size(), &x_phi[0], &y_phi[0], &ex_phi[0], &ey_phi[0]);
        graph_phi->SetTitle("#phi Distribution;#phi;Raw Counts");
        graph_phi->SetMarkerStyle(20);
        graph_phi->SetMarkerSize(0.8);
        graph_phi->SetMarkerColor(kMagenta);
        graph_phi->SetLineColor(kMagenta);
        
        // Save all histograms and graphs to file
        TFile* outFile = new TFile("cand_general.root", "RECREATE");
        if (outFile && !outFile->IsZombie()) {
            // Clone objects before writing to avoid ownership issues
            TH1F* hist_mass_copy = (TH1F*)hist_mass->Clone("hist_mass");
            TH1F* hist_cuts_copy = (TH1F*)hist_cuts->Clone("hist_cuts");
            TH1F* hist_pt_copy = (TH1F*)hist_pt->Clone("hist_pt");
            TH1F* hist_pt_normalized_copy = (TH1F*)hist_pt_normalized->Clone("hist_pt_normalized");
            TH1F* hist_pt2_copy = (TH1F*)hist_pt2->Clone("hist_pt2");
            TH1F* hist_y_copy = (TH1F*)hist_y->Clone("hist_y");
            TH1F* hist_eta_copy = (TH1F*)hist_eta->Clone("hist_eta");
            TH1F* hist_phi_copy = (TH1F*)hist_phi->Clone("hist_phi");
            TGraphErrors* graph_copy = (TGraphErrors*)graph->Clone("graph_mass");
            TGraphErrors* graph_pt_copy = (TGraphErrors*)graph_pt->Clone("graph_pt");
            TGraphErrors* graph_pt_normalized_copy = (TGraphErrors*)graph_pt_normalized->Clone("graph_pt_normalized");
            TGraphErrors* graph_pt2_copy = (TGraphErrors*)graph_pt2->Clone("graph_pt2");
            TGraphErrors* graph_y_copy = (TGraphErrors*)graph_y->Clone("graph_y");
            TGraphErrors* graph_eta_copy = (TGraphErrors*)graph_eta->Clone("graph_eta");
            TGraphErrors* graph_phi_copy = (TGraphErrors*)graph_phi->Clone("graph_phi");
            
            hist_mass_copy->Write();
            hist_cuts_copy->Write();
            hist_pt_copy->Write();
            hist_pt_normalized_copy->Write();
            hist_pt2_copy->Write();
            hist_y_copy->Write();
            hist_eta_copy->Write();
            hist_phi_copy->Write();
            graph_copy->Write();
            graph_pt_copy->Write();
            graph_pt_normalized_copy->Write();
            graph_pt2_copy->Write();
            graph_y_copy->Write();
            graph_eta_copy->Write();
            graph_phi_copy->Write();
            
            // Create TTree to store filtering parameters
            TTree* paramTree = new TTree("parameters", "Analysis Parameters");
            
            // Create branches using the original parameter variables directly
            paramTree->Branch("runNb_min", &runNb_min, "runNb_min/i");
            paramTree->Branch("runNb_max", &runNb_max, "runNb_max/i");
            paramTree->Branch("twhfMinus_max", &twhfMinus_max, "twhfMinus_max/F");
            paramTree->Branch("twhfPlus_max", &twhfPlus_max, "twhfPlus_max/F");
            paramTree->Branch("dca_significance_max", &dca_significance_max, "dca_significance_max/F");
            paramTree->Branch("candidate_pt_min", &candidate_pt_min, "candidate_pt_min/F");
            paramTree->Branch("mass_min", &mass_min, "mass_min/F");
            paramTree->Branch("mass_max", &mass_max, "mass_max/F");
            
            // Fill the tree with one entry containing all parameters
            paramTree->Fill();
            paramTree->Write();
            
            outFile->Close();
            delete outFile;
            cout << "All histograms, graphs, and parameters saved to cand_general.root" << endl;
        } else {
            cout << "Error creating output file!" << endl;
        }
        
        // Create canvas for mass distribution graph
        TCanvas* c1 = new TCanvas("c1", "Mass Distribution", 800, 600);
        c1->SetGrid();
        c1->SetLogy();  // Use log scale for y-axis
        graph->Draw("AP");  // "A" = axis, "P" = points with error bars
        graph->GetXaxis()->SetRangeUser(mass_min, mass_max);  // Set x-axis range to match histogram
        c1->SaveAs("cand_mass_graph.png");
        //c1->SaveAs("cand_mass_graph.pdf");
        cout << "Mass distribution graph saved as cand_mass_graph.png and cand_mass_graph.pdf" << endl;
        
        // Create canvas for p_{T} distribution histogram
        TCanvas* c2 = new TCanvas("c2", "p_{T} Distribution", 800, 600);
        c2->SetGrid();
        hist_pt->SetFillColor(kBlue);
        hist_pt->SetFillStyle(3001);
        hist_pt->Draw("HIST");
        c2->SaveAs("cand_pt_distribution.png");
        //c2->SaveAs("cand_pt_distribution.pdf");
        cout << "p_{T} distribution saved as cand_pt_distribution.png and cand_pt_distribution.pdf" << endl;
        
        // Create canvas for p_{T} distribution graph
        TCanvas* c2_graph = new TCanvas("c2_graph", "p_{T} Distribution Graph", 800, 600);
        c2_graph->SetGrid();
        graph_pt->Draw("AP");  // "A" = axis, "P" = points with error bars
        graph_pt->GetXaxis()->SetRangeUser(pt_min, pt_max);  // Set x-axis range to match histogram
        c2_graph->SaveAs("cand_pt_graph.png");
        //c2_graph->SaveAs("cand_pt_graph.pdf");
        cout << "p_{T} distribution graph saved as cand_pt_graph.png and cand_pt_graph.pdf" << endl;
        
        // Create canvas for N/pT distribution graph
        TCanvas* c2_normalized = new TCanvas("c2_normalized", "N/p_{T} Distribution Graph", 800, 600);
        c2_normalized->SetGrid();
        graph_pt_normalized->Draw("AP");  // "A" = axis, "P" = points with error bars
        graph_pt_normalized->GetXaxis()->SetRangeUser(pt_min, pt_max);  // Set x-axis range to match histogram
        graph_pt_normalized->GetYaxis()->SetMaxDigits(3);  // Use scientific notation for y-axis
        c2_normalized->SaveAs("cand_pt_normalized_graph.png");
        //c2_normalized->SaveAs("cand_pt_normalized_graph.pdf");
        cout << "N/p_{T} distribution graph saved as cand_pt_normalized_graph.png and cand_pt_normalized_graph.pdf" << endl;
        
        // Define fit function: a*exp(-b*pt^2) + c/(1+pt^2/d)^2
        TF1* fitFunc = new TF1("fitFunc", "[0]*exp(-[1]*x*x) + [2]/((1+x*x/[3])*(1+x*x/[3]))", pt_min, pt_max);
        fitFunc->SetParameters(800, 400, 100, 0.01);  // Initial parameter estimates based on data description
        fitFunc->SetParNames("a", "b", "c", "d");
        
        // Set reasonable parameter limits to help convergence
        fitFunc->SetParLimits(0, 100, 2000);    // a: amplitude of exponential (100-2000)
        fitFunc->SetParLimits(1, 10, 1000);     // b: exponential decay rate (10-1000)  
        fitFunc->SetParLimits(2, 10, 500);      // c: amplitude of power law (10-500)
        fitFunc->SetParLimits(3, 0.001, 0.1);   // d: power law scale parameter (0.001-0.1)
        
        // Perform the fit
        cout << "Fitting p_{T} distribution with function: a*exp(-b*pt^2) + c/(1+pt^2/d)^2" << endl;
        TFitResultPtr fitResult = graph_pt->Fit(fitFunc, "RS");
        
        // Create individual component functions
        TF1* expComponent = new TF1("expComponent", "[0]*exp(-[1]*x*x)", pt_min, pt_max);
        expComponent->SetParameters(fitFunc->GetParameter(0), fitFunc->GetParameter(1));
        expComponent->SetLineColor(kRed);
        expComponent->SetLineStyle(2);  // Dashed line
        expComponent->SetLineWidth(2);
        
        TF1* powerComponent = new TF1("powerComponent", "[0]/((1+x*x/[1])*(1+x*x/[1]))", pt_min, pt_max);
        powerComponent->SetParameters(fitFunc->GetParameter(2), fitFunc->GetParameter(3));
        powerComponent->SetLineColor(kGreen);
        powerComponent->SetLineStyle(3);  // Dotted line
        powerComponent->SetLineWidth(2);
        
        // Set fit function style
        fitFunc->SetLineColor(kMagenta);
        fitFunc->SetLineWidth(3);
        
        // Create canvas for p_{T} distribution with fit
        TCanvas* c2_fit = new TCanvas("c2_fit", "p_{T} Distribution with Fit", 800, 600);
        c2_fit->SetGrid();
        c2_fit->SetLogy();  // Use log scale for y-axis to better see the fit
        graph_pt->Draw("AP");  // "A" = axis, "P" = points with error bars
        graph_pt->GetXaxis()->SetRangeUser(pt_min, pt_max);  // Set x-axis range to match histogram
        fitFunc->Draw("same");
        expComponent->Draw("same");
        powerComponent->Draw("same");
        
        // Add legend
        TLegend* legend = new TLegend(0.5, 0.7, 0.9, 0.9);
        legend->AddEntry(graph_pt, "Data", "p");
        legend->AddEntry(fitFunc, "Total Fit: a*exp(-b*p_{T}^{2}) + c/(1+p_{T}^{2}/d)^{2}", "l");
        legend->AddEntry(expComponent, "Exponential: a*exp(-b*p_{T}^{2})", "l");
        legend->AddEntry(powerComponent, "Power Law: c/(1+p_{T}^{2}/d)^{2}", "l");
        legend->SetTextSize(0.03);
        legend->Draw();
        
        c2_fit->SaveAs("cand_pt_fit.png");
        c2_fit->SaveAs("cand_pt_fit.pdf");
        cout << "p_{T} distribution fit saved as cand_pt_fit.png and cand_pt_fit.pdf" << endl;
        
        // Print fit results
        cout << "Fit Results:" << endl;
        cout << "  a = " << fitFunc->GetParameter(0) << " ± " << fitFunc->GetParError(0) << endl;
        cout << "  b = " << fitFunc->GetParameter(1) << " ± " << fitFunc->GetParError(1) << endl;
        cout << "  c = " << fitFunc->GetParameter(2) << " ± " << fitFunc->GetParError(2) << endl;
        cout << "  d = " << fitFunc->GetParameter(3) << " ± " << fitFunc->GetParError(3) << endl;
        cout << "  Chi2/NDF = " << fitFunc->GetChisquare() << "/" << fitFunc->GetNDF() << " = " << fitFunc->GetChisquare()/fitFunc->GetNDF() << endl;
        
        // Save fit parameters to text file
        ofstream fitFile("pt_fit_parameters.txt");
        if (fitFile.is_open()) {
            fitFile << "p_T Distribution Fit Results" << endl;
            fitFile << "============================" << endl;
            fitFile << "Fit Function: a*exp(-b*p_T^2) + c/(1+p_T^2/d)^2" << endl;
            fitFile << endl;
            fitFile << "Parameters:" << endl;
            fitFile << "a = " << fitFunc->GetParameter(0) << " ± " << fitFunc->GetParError(0) << endl;
            fitFile << "b = " << fitFunc->GetParameter(1) << " ± " << fitFunc->GetParError(1) << endl;
            fitFile << "c = " << fitFunc->GetParameter(2) << " ± " << fitFunc->GetParError(2) << endl;
            fitFile << "d = " << fitFunc->GetParameter(3) << " ± " << fitFunc->GetParError(3) << endl;
            fitFile << endl;
            fitFile << "Goodness of Fit:" << endl;
            fitFile << "Chi2 = " << fitFunc->GetChisquare() << endl;
            fitFile << "NDF = " << fitFunc->GetNDF() << endl;
            fitFile << "Chi2/NDF = " << fitFunc->GetChisquare()/fitFunc->GetNDF() << endl;
            fitFile << "Probability = " << TMath::Prob(fitFunc->GetChisquare(), fitFunc->GetNDF()) << endl;
            fitFile << endl;
            fitFile << "Fit Range: " << pt_min << " to " << pt_max << " GeV/c" << endl;
            fitFile << "Number of data points: " << graph_pt->GetN() << endl;
            fitFile.close();
            cout << "Fit parameters saved to pt_fit_parameters.txt" << endl;
        } else {
            cout << "Error: Could not create pt_fit_parameters.txt file" << endl;
        }
        
        // Create canvas for p_{T}^{2} distribution histogram
        TCanvas* c2_pt2 = new TCanvas("c2_pt2", "p_{T}^{2} Distribution", 800, 600);
        c2_pt2->SetGrid();
        c2_pt2->SetLogx();  // Use log scale for x-axis
        c2_pt2->SetLogy();  // Use log scale for y-axis
        hist_pt2->SetFillColor(kCyan);
        hist_pt2->SetFillStyle(3001);
        hist_pt2->Draw("HIST");
        c2_pt2->SaveAs("cand_pt2_distribution.png");
        //c2_pt2->SaveAs("cand_pt2_distribution.pdf");
        cout << "p_{T}^{2} distribution saved as cand_pt2_distribution.png and cand_pt2_distribution.pdf" << endl;
        
        // Create canvas for p_{T}^{2} distribution graph
        TCanvas* c2_pt2_graph = new TCanvas("c2_pt2_graph", "p_{T}^{2} Distribution Graph", 800, 600);
        c2_pt2_graph->SetGrid();
        c2_pt2_graph->SetLogx();  // Use log scale for x-axis
        c2_pt2_graph->SetLogy();  // Use log scale for y-axis
        graph_pt2->Draw("AP");  // "A" = axis, "P" = points with error bars
        graph_pt2->GetXaxis()->SetRangeUser(pt2_min, pt2_max);  // Set x-axis range to match histogram
        c2_pt2_graph->SaveAs("cand_pt2_graph.png");
        //c2_pt2_graph->SaveAs("cand_pt2_graph.pdf");
        cout << "p_{T}^{2} distribution graph saved as cand_pt2_graph.png and cand_pt2_graph.pdf" << endl;
        
        // Create canvas for rapidity distribution histogram
        TCanvas* c3 = new TCanvas("c3", "Rapidity Distribution", 800, 600);
        c3->SetGrid();
        c3->SetLogy();  // Use log scale for y-axis
        hist_y->SetFillColor(kGreen);
        hist_y->SetFillStyle(3001);
        hist_y->Draw("HIST");
        c3->SaveAs("cand_y_distribution.png");
        //c3->SaveAs("cand_y_distribution.pdf");
        cout << "Rapidity distribution saved as cand_y_distribution.png and cand_y_distribution.pdf" << endl;
        
        // Create canvas for rapidity distribution graph
        TCanvas* c3_graph = new TCanvas("c3_graph", "Rapidity Distribution Graph", 800, 600);
        c3_graph->SetGrid();
        c3_graph->SetLogy();  // Use log scale for y-axis
        graph_y->Draw("AP");  // "A" = axis, "P" = points with error bars
        graph_y->GetXaxis()->SetRangeUser(-3, 3);  // Set x-axis range to match histogram
        c3_graph->SaveAs("cand_y_graph.png");
        //c3_graph->SaveAs("cand_y_graph.pdf");
        cout << "Rapidity distribution graph saved as cand_y_graph.png and cand_y_graph.pdf" << endl;
        
        // Create canvas for pseudorapidity distribution histogram
        TCanvas* c4 = new TCanvas("c4", "Pseudorapidity Distribution", 800, 600);
        c4->SetGrid();
        hist_eta->SetFillColor(kRed);
        hist_eta->SetFillStyle(3001);
        hist_eta->Draw("HIST");
        c4->SaveAs("cand_eta_distribution.png");
        //c4->SaveAs("cand_eta_distribution.pdf");
        cout << "Pseudorapidity distribution saved as cand_eta_distribution.png and cand_eta_distribution.pdf" << endl;
        
        // Create canvas for pseudorapidity distribution graph
        TCanvas* c4_graph = new TCanvas("c4_graph", "Pseudorapidity Distribution Graph", 800, 600);
        c4_graph->SetGrid();
        graph_eta->Draw("AP");  // "A" = axis, "P" = points with error bars
        graph_eta->GetXaxis()->SetRangeUser(-3, 3);  // Set x-axis range to match histogram
        c4_graph->SaveAs("cand_eta_graph.png");
        //c4_graph->SaveAs("cand_eta_graph.pdf");
        cout << "Pseudorapidity distribution graph saved as cand_eta_graph.png and cand_eta_graph.pdf" << endl;
        
        // Create canvas for phi distribution histogram
        TCanvas* c6 = new TCanvas("c6", "#phi Distribution", 800, 600);
        c6->SetGrid();
        hist_phi->SetFillColor(kMagenta);
        hist_phi->SetFillStyle(3001);
        hist_phi->Draw("HIST");
        // Set x-axis labels in relation to pi
        hist_phi->GetXaxis()->SetNdivisions(505);
        hist_phi->GetXaxis()->SetLabelSize(0.03);
        c6->SaveAs("cand_phi_distribution.png");
        //c6->SaveAs("cand_phi_distribution.pdf");
        cout << "#phi distribution saved as cand_phi_distribution.png and cand_phi_distribution.pdf" << endl;
        
        // Create canvas for phi distribution graph
        TCanvas* c6_graph = new TCanvas("c6_graph", "#phi Distribution Graph", 800, 600);
        c6_graph->SetGrid();
        graph_phi->Draw("AP");  // "A" = axis, "P" = points with error bars
        graph_phi->GetXaxis()->SetRangeUser(-TMath::Pi(), TMath::Pi());  // Set x-axis range to match histogram
        // Set x-axis labels in relation to pi for the graph
        graph_phi->GetXaxis()->SetNdivisions(505);
        graph_phi->GetXaxis()->SetLabelSize(0.03);
        c6_graph->SaveAs("cand_phi_graph.png");
        //c6_graph->SaveAs("cand_phi_graph.pdf");
        cout << "#phi distribution graph saved as cand_phi_graph.png and cand_phi_graph.pdf" << endl;
        
        // Create canvas for filtering steps histogram
        TCanvas* c5 = new TCanvas("c5", "Event Selection Steps", 1000, 600);
        //c5->SetLogy();  // Use log scale for y-axis due to large range
        hist_cuts->SetFillColor(kBlue);
        hist_cuts->SetFillStyle(3001);
        hist_cuts->SetStats(0);  // Disable statistics box
        
        // Set y-axis maximum to accommodate text labels
        Double_t maxBinContent = hist_cuts->GetMaximum();
        hist_cuts->SetMaximum(maxBinContent * 1.2);  // 60% extra space for text labels
        
        hist_cuts->Draw("HIST");
        
        // Add text annotations with percentages - use array instead of vector to avoid conflicts
        TText* textLabels[11];  // Updated size for 9 bins + final count + total count
        int nLabels = 0;
        
        // Add the total number of events on top of the first bin
        Double_t firstBinContent = hist_cuts->GetBinContent(1);
        TString totalLabel = Form("%lld", count_total);
        textLabels[nLabels] = new TText(0.5, firstBinContent * 1.12, totalLabel);
        textLabels[nLabels]->SetTextAlign(22);  // Center alignment
        textLabels[nLabels]->SetTextSize(0.04);
        textLabels[nLabels]->SetTextColor(kBlack);
        textLabels[nLabels]->Draw();
        nLabels++;
        
        for (int i = 1; i <= hist_cuts->GetNbinsX(); i++) {
            Double_t binContent = hist_cuts->GetBinContent(i);
            if (binContent > 0 && nLabels < 10) {
                Double_t percentage = (100.0 * binContent / count_total);
                TString label = Form("%.1f%%", percentage);
                textLabels[nLabels] = new TText(i-0.5, binContent +firstBinContent* 0.05, label);
                textLabels[nLabels]->SetTextAlign(22);  // Center alignment
                textLabels[nLabels]->SetTextSize(0.03);
                textLabels[nLabels]->Draw();
                nLabels++;
            }
        }
        
        // Add the number of remaining events on top of the last bin
        if (nLabels < 11) {
            Double_t lastBinContent = hist_cuts->GetBinContent(hist_cuts->GetNbinsX());
            TString finalLabel = Form("%lld", count_final);
            textLabels[nLabels] = new TText(hist_cuts->GetNbinsX()-0.5, lastBinContent + firstBinContent* 0.12, finalLabel);
            textLabels[nLabels]->SetTextAlign(22);  // Center alignment
            textLabels[nLabels]->SetTextSize(0.04);
            textLabels[nLabels]->SetTextColor(kRed);
            textLabels[nLabels]->Draw();
            nLabels++;
        }
        
        c5->SaveAs("filtering_steps.png");
        //c5->SaveAs("filtering_steps.pdf");
        cout << "Filtering steps histogram saved as filtering_steps.png and filtering_steps.pdf" << endl;
        
        // Print statistics before cleanup
        cout << "Graph Statistics:" << endl;
        cout << "  Mass graph data points: " << graph->GetN() << endl;
        cout << "  p_{T} graph data points: " << graph_pt->GetN() << endl;
        cout << "  p_{T}^{2} graph data points: " << graph_pt2->GetN() << endl;
        cout << "  Rapidity graph data points: " << graph_y->GetN() << endl;
        cout << "  Pseudorapidity graph data points: " << graph_eta->GetN() << endl;
        cout << "  Phi graph data points: " << graph_phi->GetN() << endl;
        Double_t xmin, ymin, xmax, ymax;
        graph->ComputeRange(xmin, ymin, xmax, ymax);
        cout << "  Mass X range: " << xmin << " to " << xmax << endl;
        cout << "  Mass Y range: " << ymin << " to " << ymax << endl;
        
        // Proper cleanup - delete text labels
        for (int i = 0; i < nLabels; i++) {
            delete textLabels[i];
        }
        
        // Don't delete canvases and objects that ROOT manages
        // Just clear the canvases and let ROOT handle cleanup
        c6_graph->Clear();
        c6->Clear();
        c5->Clear();
        c4_graph->Clear();
        c4->Clear();
        c3_graph->Clear();
        c3->Clear();
        c2_pt2_graph->Clear();
        c2_pt2->Clear();
        c2_fit->Clear();
        c2_graph->Clear();
        c2->Clear();
        c1->Clear();
        
        cout << "Analysis finished!" << endl;
    } 
    else {
        cout << "Failed to create histogram - no entries passed the cuts" << endl;
    }
    
    // Write and close filtered tree and file
    if (kMakeFilteredTree && filteredFile && filteredTree) {
        filteredFile->cd();
        filteredTree->Write();
        filteredFile->Close();
        delete filteredFile;
        cout << "Filtered tree saved to filtered_chained_ParticleTree.root" << endl;
    }

    // Clean up
    file->Close();
    delete file;
    cout << "Analysis finished!" << endl;
}

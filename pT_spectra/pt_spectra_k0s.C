// ============================================================================
// K0s_pT_Spectra_Fresh.C - Fresh start for K0s pT spectra comparison
// ============================================================================
// Compares PYTHIA (Monash & Junctions) vs ALICE for K0s
// Three multiplicity selections: 0-1%, 0-20%, 80-100%
//
// ALICE mapping (from HEPData Table 52):
//   0-1%:     Table 1 (V0M Class I: 0.0-0.95%)
//   0-20%:    Tables 1-5 averaged (V0M Classes I-V: 0.0-19.0%)
//   80-100%:  Table 10 (V0M Class X: 68.0-100.0%)
//
// Usage:
// root -l 'K0s_pT_Spectra_Fresh.C("spectra_monash_100M.root", "spectra_junctions_100M.root", "HEPData-ins1471838-v1-root.root")'
//

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

// ============================================================================
// Calculate percentile bins from multiplicity histogram
// ============================================================================
std::vector<int> CalculatePercentileBins(TH1D* hMult) {
    const int nbMult = hMult->GetNbinsX();
    const double total = hMult->Integral(1, nbMult);

    std::vector<int> binEdges;

    // We need edges for: 0-1%, 0-20%, 80-100%
    // Which gives us bins: [low_mult, 80%, 20%, 1%, high_mult]

    double acc = 0.0;
    int bin_80 = 1;   // 80th percentile (from high mult side)
    int bin_20 = 1;   // 20th percentile
    int bin_1 = 1;    // 1st percentile

    // Start from highest multiplicity, accumulate events
    for (int b = nbMult; b >= 1; --b) {
        acc += hMult->GetBinContent(b);
        double frac = acc / total;

        if (frac >= 0.01 && bin_1 == 1) bin_1 = b;
        if (frac >= 0.20 && bin_20 == 1) bin_20 = b;
        if (frac >= 0.80 && bin_80 == 1) bin_80 = b;
    }

    // Build bin edges: [lowest, 80%, 20%, 1%, highest+1]
    binEdges.push_back(1);           // Lowest mult (bin 1)
    binEdges.push_back(bin_80);      // 80th percentile
    binEdges.push_back(bin_20);      // 20th percentile
    binEdges.push_back(bin_1);       // 1st percentile
    binEdges.push_back(nbMult + 1);  // Highest mult + 1

    return binEdges;
}

// ============================================================================
// Get ALICE spectrum from single table or average of tables
// ============================================================================
TGraphErrors* GetAliceSpectrum(TFile* aliceFile, std::vector<int> tables) {
    std::vector<TH1*> histograms;

    for (int tableNum : tables) {
        TDirectory* tableDir = (TDirectory*)aliceFile->Get(Form("Table %d", tableNum));
        if (!tableDir) {
            std::cout << "  Warning: Table " << tableNum << " not found\n";
            continue;
        }

        auto* h = dynamic_cast<TH1*>(tableDir->Get("Hist1D_y1"));
        if (h) {
            histograms.push_back((TH1*)h->Clone(Form("hAlice_%d", tableNum)));
            std::cout << "  Loaded Table " << tableNum << "\n";
        }
    }

    if (histograms.empty()) return nullptr;

    // Average if multiple tables
    TH1* hAvg = (TH1*)histograms[0]->Clone("hAvg");
    for (size_t i = 1; i < histograms.size(); ++i) {
        hAvg->Add(histograms[i]);
    }
    if (histograms.size() > 1) {
        hAvg->Scale(1.0 / histograms.size());
        std::cout << "  Averaged " << histograms.size() << " tables\n";
    }

    // Convert to TGraph (no error bars for ALICE)
    std::vector<double> vx, vy, vex, vey;
    for (int b = 1; b <= hAvg->GetNbinsX(); ++b) {
        double pt = hAvg->GetBinCenter(b);
        double value = hAvg->GetBinContent(b);

        if (value > 0 && pt > 0) {
            vx.push_back(pt);
            vy.push_back(value);
            vex.push_back(0.0);
            vey.push_back(0.0);  // No error bars
        }
    }

    delete hAvg;
    for (auto* h : histograms) delete h;

    if (vx.empty()) return nullptr;

    return new TGraphErrors(vx.size(), vx.data(), vy.data(), vex.data(), vey.data());
}

// ============================================================================
// Get PYTHIA spectrum with simple normalization
// ============================================================================
TGraphErrors* GetPythiaSpectrum(TFile* pythiaFile, TH1D* hMult, int binLow, int binHigh, const char* name) {
    auto* h2D = dynamic_cast<TH2D*>(pythiaFile->Get("fHistPtK0s"));
    if (!h2D) return nullptr;

    // Project pT spectrum
    TH1D* hPt = h2D->ProjectionX(name, binLow, binHigh);
    double nEvents = hMult->Integral(binLow, binHigh);

    if (nEvents <= 0) {
        delete hPt;
        return nullptr;
    }

    std::cout << "  Events in range: " << nEvents << "\n";

    // Normalize: (1/Nev) × (d²N/dy dpT)
    // PYTHIA uses |η| < 4.0 (Δη = 8.0), ALICE uses |y| < 0.5 (Δy = 1.0)
    // Divide by 8.0 to convert from integrated counts to differential per unit rapidity
    const double acceptanceCorrection = 8.0;

    std::vector<double> vx, vy, vex, vey;

    for (int b = 1; b <= hPt->GetNbinsX(); ++b) {
        double pt = hPt->GetBinCenter(b);
        double binWidth = hPt->GetBinWidth(b);
        double counts = hPt->GetBinContent(b);

        if (counts > 0 && binWidth > 0) {
            // Divide by acceptance to match ALICE format: d²N/(dy dpT)
            double d2NdydpT = (counts / nEvents) / (binWidth * acceptanceCorrection);
            double d2NdydpT_err = (TMath::Sqrt(counts) / nEvents) / (binWidth * acceptanceCorrection);

            vx.push_back(pt);
            vy.push_back(d2NdydpT);
            vex.push_back(0.0);
            vey.push_back(d2NdydpT_err);
        }
    }

    delete hPt;

    if (vx.empty()) return nullptr;

    return new TGraphErrors(vx.size(), vx.data(), vy.data(), vex.data(), vey.data());
}

// ============================================================================
// Main function
// ============================================================================
void K0s_pT_Spectra_Fresh(
    const char* monashFile = "spectra_monash_100M.root",
    const char* junctionsFile = "spectra_junctions_100M.root",
    const char* aliceFile = "HEPData-ins1471838-v1-root.root")
{
    std::cout << "\n╔══════════════════════════════════════════════════════╗\n";
    std::cout << "║  K0s pT Spectra: 0-1%, 0-20%, 80-100% Multiplicity  ║\n";
    std::cout << "╚══════════════════════════════════════════════════════╝\n\n";

    gStyle->SetOptStat(0);
    gStyle->SetGridColor(kGray);
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(1);

    // Open files
    TFile* fMon = TFile::Open(monashFile);
    TFile* fJunc = TFile::Open(junctionsFile);
    TFile* fAlice = TFile::Open(aliceFile);

    if (!fMon || !fJunc || !fAlice) {
        std::cerr << "ERROR: Cannot open files!\n";
        return;
    }

    auto* hMultMon = dynamic_cast<TH1D*>(fMon->Get("fHistMultiplicity"));
    auto* hMultJunc = dynamic_cast<TH1D*>(fJunc->Get("fHistMultiplicity"));

    if (!hMultMon || !hMultJunc) {
        std::cerr << "ERROR: Missing multiplicity histograms!\n";
        return;
    }

    // Calculate percentile bins
    // binEdges: [lowest_mult, 80%, 20%, 1%, highest_mult+1]
    std::vector<int> binEdges = CalculatePercentileBins(hMultMon);

    std::cout << "Multiplicity bin edges:\n";
    std::cout << "  80-100%:  bins " << binEdges[0] << " to " << binEdges[1] << "\n";
    std::cout << "  0-20%:    bins " << binEdges[2] << " to " << binEdges[4] << "\n";
    std::cout << "  0-1%:     bins " << binEdges[3] << " to " << binEdges[4] << "\n\n";

    // Define multiplicity selections
    struct MultSelection {
        const char* name;
        const char* label;
        int binLow;
        int binHigh;
        std::vector<int> aliceTables;
    };

    std::vector<MultSelection> selections = {
        {"0_1",   "0-1%",    binEdges[3], binEdges[4], {1}},           // Class I (highest mult)
        {"0_20",  "0-20%",   binEdges[2], binEdges[4], {1,2,3,4,5}},   // Classes I-V
        {"80_100","80-100%", binEdges[0], binEdges[1], {10}}           // Class X (lowest mult)
    };

    // Create separate canvas for each multiplicity selection
    for (size_t i = 0; i < selections.size(); ++i) {
        auto& sel = selections[i];

        std::cout << "═══════════════════════════════════════════════════\n";
        std::cout << "Processing " << sel.label << "\n";
        std::cout << "  PYTHIA bins: " << sel.binLow << " to " << sel.binHigh << "\n";
        std::cout << "  ALICE tables: ";
        for (auto t : sel.aliceTables) std::cout << t << " ";
        std::cout << "\n";

        // Create individual canvas for this selection
        TCanvas* c = new TCanvas(Form("c_%s", sel.name), Form("K0s %s", sel.label), 800, 700);
        c->SetLogy();
        c->SetGrid();
        c->SetLeftMargin(0.14);
        c->SetRightMargin(0.05);
        c->SetBottomMargin(0.12);
        c->SetTopMargin(0.10);

        // Get PYTHIA spectra
        std::cout << "  Getting Monash spectrum...\n";
        auto* gMonash = GetPythiaSpectrum(fMon, hMultMon, sel.binLow, sel.binHigh, Form("gMon_%s", sel.name));

        std::cout << "  Getting Junctions spectrum...\n";
        auto* gJunctions = GetPythiaSpectrum(fJunc, hMultJunc, sel.binLow, sel.binHigh, Form("gJunc_%s", sel.name));

        // Get ALICE spectrum
        std::cout << "  Getting ALICE spectrum...\n";
        auto* gAlice = GetAliceSpectrum(fAlice, sel.aliceTables);

        if (!gMonash || !gJunctions) {
            std::cerr << "  ERROR: Missing PYTHIA data for " << sel.label << "\n";
            continue;
        }

        // Style Monash
        gMonash->SetLineColor(kBlue+1);
        gMonash->SetLineWidth(2);
        gMonash->SetMarkerColor(kBlue+1);
        gMonash->SetMarkerStyle(20);
        gMonash->SetMarkerSize(0.8);

        // Style Junctions
        gJunctions->SetLineColor(kRed+1);
        gJunctions->SetLineWidth(2);
        gJunctions->SetLineStyle(2);
        gJunctions->SetMarkerColor(kRed+1);
        gJunctions->SetMarkerStyle(21);
        gJunctions->SetMarkerSize(0.8);

        // Style ALICE
        if (gAlice) {
            gAlice->SetMarkerStyle(25);
            gAlice->SetMarkerColor(kBlack);
            gAlice->SetLineColor(kBlack);
            gAlice->SetMarkerSize(1.2);
        }

        // Draw
        gMonash->Draw("APL");

        // Set title with K0s and multiplicity label
        gMonash->SetTitle(Form("K^{0}_{s} Transverse Momentum Spectrum (%s);p_{T} (GeV/c);(1/N_{ev}) d^{2}N/(dy dp_{T}) (GeV/c)^{-1}", sel.label));
        gMonash->GetXaxis()->SetRangeUser(0, 10);
        gMonash->GetYaxis()->SetRangeUser(1e-6, 10);  // Extended from 1e-5 to 1e-6
        gMonash->GetXaxis()->SetTitleSize(0.045);
        gMonash->GetYaxis()->SetTitleSize(0.038);
        gMonash->GetXaxis()->SetLabelSize(0.04);
        gMonash->GetYaxis()->SetLabelSize(0.04);
        gMonash->GetYaxis()->SetTitleOffset(1.5);

        gJunctions->Draw("PL SAME");
        if (gAlice) gAlice->Draw("P SAME");

        // Legend
        TLegend* leg = new TLegend(0.60, 0.73, 0.88, 0.88);
        leg->SetBorderSize(1);
        leg->SetFillColor(kWhite);
        leg->SetTextSize(0.035);
        leg->AddEntry(gMonash, "Monash", "lp");
        leg->AddEntry(gJunctions, "Junctions", "lp");
        if (gAlice) leg->AddEntry(gAlice, "ALICE", "p");
        leg->Draw();

        // Save individual plot
        TString outName = Form("K0s_pT_Spectrum_%s.png", sel.name);
        c->SaveAs(outName);
        std::cout << "  ✓ Saved: " << outName << "\n\n";

        delete c;
    }

    std::cout << "═══════════════════════════════════════════════════\n";
    std::cout << "✓ All plots saved!\n\n";

    fMon->Close();
    fJunc->Close();
    fAlice->Close();
}

// ============================================================================
//  plot_pT_spectra_focused.C - Focused pT spectra for specific percentile bins
// ============================================================================
// Compares PYTHIA (Monash & Junctions) vs ALICE for ONE multiplicity bin at a time
// Uses YOUR existing binning and normalization scheme
//
// Usage: root -l -q 'plot_pT_spectra_focused.C("spectra_monash_100M.root", "spectra_junctions_100M.root", "HEPData-ins1471838-v1-root.root", "K0s", "0_20")'

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <iostream>
#include <vector>

// Map particle name to histogram and ALICE tables
struct ParticleConfig {
    const char* histName;      // e.g., "fHistPtK0s"
    const char* displayName;   // e.g., "K^{0}_{s}"
    int aliceTableStart;       // First table for this particle
    int color;
};

std::map<std::string, ParticleConfig> GetParticleConfigs() {
    std::map<std::string, ParticleConfig> configs;
    configs["K0s"]    = {"fHistPtK0s",     "K^{0}_{s}",    1,  kGreen+2};
    configs["Lambda"] = {"fHistPtLambdas", "#Lambda",      11, kRed+1};
    configs["Xi"]     = {"fHistPtXis",     "#Xi^{#pm}",   21, kViolet};
    configs["Omega"]  = {"fHistPtOmegas",  "#Omega^{#pm}",31, kOrange+4};
    return configs;
}

// Calculate percentile bin edges (same as drawMultiplicity.C)
std::vector<int> CalculatePercentileBins(TH1D* hMult) {
    const int nbMult = hMult->GetNbinsX();
    const double total = hMult->Integral(1, nbMult);

    std::vector<int> binEdges;
    binEdges.push_back(nbMult + 1);

    double acc = 0.0;
    int targetIdx = 1;
    const std::array<double,5> targets{0.20, 0.40, 0.60, 0.80, 1.00};  // 5 bins: 0-20, 20-40, etc.

    for (int b = nbMult; b >= 1 && targetIdx <= 5; --b) {
        acc += hMult->GetBinContent(b);
        double frac = acc / total;
        while (targetIdx <= 5 && frac >= targets[targetIdx-1]) {
            binEdges.push_back(b);
            ++targetIdx;
        }
    }
    while (binEdges.size() < 6) binEdges.push_back(1);
    std::reverse(binEdges.begin(), binEdges.end());

    return binEdges;
}

// Extract and sum ALICE spectra from multiple tables
// Converts ALICE format to YOUR simple dN/dpT format
TGraphErrors* GetAliceCombinedSpectrum(TFile* aliceFile, int tableStart, std::vector<int> tableOffsets) {
    std::vector<TH1*> histograms;

    // Load all requested tables
    for (int offset : tableOffsets) {
        int tableNum = tableStart + offset;
        TDirectory* tableDir = (TDirectory*)aliceFile->Get(Form("Table %d", tableNum));
        if (!tableDir) continue;

        auto* h = dynamic_cast<TH1*>(tableDir->Get("Hist1D_y1"));
        if (h) histograms.push_back((TH1*)h->Clone(Form("hAlice_table%d", tableNum)));
    }

    if (histograms.empty()) return nullptr;

    // Sum histograms (they have same binning)
    TH1* hSum = (TH1*)histograms[0]->Clone("hSum_alice");
    for (size_t i = 1; i < histograms.size(); ++i) {
        hSum->Add(histograms[i]);
    }

    // Convert ALICE format [1/(Nev·2πpT·ΔpT) d²N/dydpT] to per-event format [1/(Nev·ΔpT) dN/dpT]
    // ALICE has: (1/Nev) × (1/(2π·pT·ΔpT)) × d²N/dydpT
    // We want: (1/Nev) × (1/ΔpT) × dN/dpT
    // So multiply by (2π × pT) to remove those factors
    std::vector<double> vx, vy, vex, vey;
    for (int b = 1; b <= hSum->GetNbinsX(); ++b) {
        double pt = hSum->GetBinCenter(b);
        double binWidth = hSum->GetBinWidth(b);
        double aliceValue = hSum->GetBinContent(b);
        double aliceError = hSum->GetBinError(b);

        if (aliceValue > 0 && binWidth > 0 && pt > 0) {
            // ALICE Table 1 format (from HEPData): 1/(Nev·2πpT) d²N/(dydpT)
            // Just use raw ALICE value - no corrections applied
            double perEventYield = aliceValue;
            double perEventError = aliceError;

            vx.push_back(pt);
            vy.push_back(perEventYield);
            vex.push_back(0.0);
            vey.push_back(perEventError);
        }
    }

    delete hSum;

    if (vx.empty()) return nullptr;

    auto* g = new TGraphErrors((int)vx.size(), vx.data(), vy.data(), vex.data(), vey.data());
    return g;
}

// Get PYTHIA spectrum for a percentile bin using YOUR simple normalization
// Uses subsampling method to calculate errors
TH1D* GetPythiaSpectrum(TFile* pythiaFile, const char* histName, TH1D* hMult,
                        int binLow, int binHigh, const char* name) {
    auto* h2D = dynamic_cast<TH2D*>(pythiaFile->Get(histName));
    if (!h2D) return nullptr;

    // Get number of events in this multiplicity range
    double nEvents = hMult->Integral(binLow, binHigh);
    if (nEvents <= 0) return nullptr;

    // Subsampling parameters
    const int nSubsamples = 10;  // Number of subsamples
    std::vector<TH1D*> subsamples;

    // Calculate events per subsample
    int eventsPerSample = (binHigh - binLow + 1) / nSubsamples;
    if (eventsPerSample < 1) eventsPerSample = 1;

    // Create subsamples by dividing the multiplicity range
    for (int i = 0; i < nSubsamples; ++i) {
        int subBinLow = binLow + i * eventsPerSample;
        int subBinHigh = (i == nSubsamples - 1) ? binHigh : (binLow + (i + 1) * eventsPerSample - 1);

        if (subBinLow > binHigh) break;
        if (subBinHigh > binHigh) subBinHigh = binHigh;

        // Project subsample
        TH1D* hSub = h2D->ProjectionX(Form("%s_sub%d", name, i), subBinLow, subBinHigh);
        double nEvSub = hMult->Integral(subBinLow, subBinHigh);

        if (nEvSub > 0) {
            // Normalize: divide by events AND bin width to match ALICE per-event format
            for (int b = 1; b <= hSub->GetNbinsX(); ++b) {
                double binWidth = hSub->GetBinWidth(b);
                double content = hSub->GetBinContent(b);

                if (binWidth > 0) {
                    hSub->SetBinContent(b, content / (nEvSub * binWidth));
                }
            }
            subsamples.push_back(hSub);
        } else {
            delete hSub;
        }
    }

    if (subsamples.empty()) return nullptr;

    // Create final histogram (use first subsample as template)
    TH1D* hPt = (TH1D*)subsamples[0]->Clone(name);
    hPt->Reset();

    // Calculate mean and standard deviation across subsamples
    for (int b = 1; b <= hPt->GetNbinsX(); ++b) {
        std::vector<double> values;

        for (auto* hSub : subsamples) {
            double val = hSub->GetBinContent(b);
            if (val > 0) values.push_back(val);
        }

        if (values.empty()) continue;

        // Calculate mean
        double mean = 0;
        for (double v : values) mean += v;
        mean /= values.size();

        // Calculate standard deviation
        double variance = 0;
        for (double v : values) {
            variance += (v - mean) * (v - mean);
        }
        variance /= (values.size() > 1) ? (values.size() - 1) : 1;
        double stddev = TMath::Sqrt(variance);

        // Standard error of the mean (accounts for number of subsamples)
        double stderr = stddev / TMath::Sqrt(values.size());

        hPt->SetBinContent(b, mean);
        hPt->SetBinError(b, stderr);
    }

    // Clean up subsamples
    for (auto* hSub : subsamples) delete hSub;

    return hPt;
}

// ============================================================================
// Main function
// ============================================================================
void plot_pT_spectra_focused_v10(
    const char* monashFile = "spectra_monash_100M.root",
    const char* junctionsFile = "spectra_junctions_100M.root",
    const char* aliceFile = "HEPData-ins1471838-v1-root.root",
    const char* particle = "K0s",      // K0s, Lambda, Xi, Omega
    const char* multBin = "0_20")      // 0_20, 20_40, 40_60, 60_80, 80_100
{
    gStyle->SetOptStat(0);

    auto particleConfigs = GetParticleConfigs();
    if (particleConfigs.find(particle) == particleConfigs.end()) {
        std::cerr << "Unknown particle: " << particle << std::endl;
        return;
    }

    ParticleConfig pc = particleConfigs[particle];

    // Open files
    TFile* fMon = TFile::Open(monashFile);
    TFile* fJunc = TFile::Open(junctionsFile);
    TFile* fAlice = TFile::Open(aliceFile);

    if (!fMon || !fJunc || !fAlice) {
        std::cerr << "Error opening files!\n";
        return;
    }

    auto* hMultMon = dynamic_cast<TH1D*>(fMon->Get("fHistMultiplicity"));
    auto* hMultJunc = dynamic_cast<TH1D*>(fJunc->Get("fHistMultiplicity"));

    if (!hMultMon || !hMultJunc) {
        std::cerr << "Error: Missing fHistMultiplicity\n";
        return;
    }

    // Calculate percentile bins (YOUR binning scheme)
    std::vector<int> binEdges = CalculatePercentileBins(hMultMon);

    // Map multiplicity bin name to index and ALICE tables
    std::map<std::string, int> binIndex;
    std::map<std::string, std::vector<int>> aliceTables;
    std::map<std::string, const char*> binLabels;

    binIndex["0_20"]    = 0;
    binIndex["20_40"]   = 1;
    binIndex["40_60"]   = 2;
    binIndex["60_80"]   = 3;
    binIndex["80_100"]  = 4;

    binLabels["0_20"]    = "0-20%";
    binLabels["20_40"]   = "20-40%";
    binLabels["40_60"]   = "40-60%";
    binLabels["60_80"]   = "60-80%";
    binLabels["80_100"]  = "80-100%";

    // Map YOUR bins to ALICE classes (from Table 52)
    // Class I-V covers ~0-19% (matches your 0-20%)
    aliceTables["0_20"]   = {0, 1, 2, 3, 4};  // Classes I, II, III, IV, V
    // Class VI-VII covers ~19-38% (closest to your 20-40%)
    aliceTables["20_40"]  = {5, 6};           // Classes VI, VII
    // Class VIII covers ~38-48% (closest to your 40-60%)
    aliceTables["40_60"]  = {7};              // Class VIII
    // Class IX covers ~48-68% (overlaps your 60-80%)
    aliceTables["60_80"]  = {8};              // Class IX
    // Class X covers ~68-100% (matches your 80-100%)
    aliceTables["80_100"] = {9};              // Class X

    if (binIndex.find(multBin) == binIndex.end()) {
        std::cerr << "Unknown multiplicity bin: " << multBin << std::endl;
        return;
    }

    int idx = binIndex[multBin];
    int binLow = binEdges[idx];
    int binHigh = binEdges[idx + 1];

    std::cout << "\n=== Plotting " << pc.displayName << " for " << binLabels[multBin] << " ===\n";
    std::cout << "PYTHIA bins: " << binLow << " to " << binHigh << "\n";
    std::cout << "ALICE tables: ";
    for (int t : aliceTables[multBin]) std::cout << (pc.aliceTableStart + t) << " ";
    std::cout << "\n\n";

    // Get PYTHIA spectra (with proper normalization)
    TH1D* hMon = GetPythiaSpectrum(fMon, pc.histName, hMultMon, binLow, binHigh, "hMon");
    TH1D* hJunc = GetPythiaSpectrum(fJunc, pc.histName, hMultJunc, binLow, binHigh, "hJunc");

    // Get ALICE spectrum (combined from multiple classes if needed)
    TGraphErrors* gAlice = GetAliceCombinedSpectrum(fAlice, pc.aliceTableStart, aliceTables[multBin]);

    if (!hMon || !hJunc) {
        std::cerr << "Error loading PYTHIA data\n";
        return;
    }

    // Create canvas
    TCanvas* c = new TCanvas("c", Form("%s pT spectrum", pc.displayName), 900, 700);
    c->SetLogy();
    c->SetLeftMargin(0.14);
    c->SetBottomMargin(0.12);
    c->SetGrid();

    // Lighter gridlines
    gStyle->SetGridColor(kGray);
    gStyle->SetGridStyle(3);  // dotted
    gStyle->SetGridWidth(1);

    // Styling
    hMon->SetLineColor(kBlue+1);
    hMon->SetLineWidth(2);
    hMon->SetMarkerColor(kBlue+1);
    hMon->SetMarkerStyle(20);
    hMon->SetMarkerSize(0.8);

    hJunc->SetLineColor(kRed+1);
    hJunc->SetLineWidth(2);
    hJunc->SetLineStyle(2);
    hJunc->SetMarkerColor(kRed+1);
    hJunc->SetMarkerStyle(21);
    hJunc->SetMarkerSize(0.8);

    if (gAlice) {
        gAlice->SetMarkerStyle(25);
        gAlice->SetMarkerColor(kBlack);
        gAlice->SetLineColor(kBlack);
        gAlice->SetMarkerSize(1.0);
    }

    // Draw
    hMon->SetTitle(Form("%s Transverse Momentum Spectrum (%s)", pc.displayName, binLabels[multBin]));
    hMon->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hMon->GetYaxis()->SetTitle("dN/dp_{T}");
    hMon->GetXaxis()->SetTitleSize(0.045);
    hMon->GetYaxis()->SetTitleSize(0.038);
    hMon->GetXaxis()->SetLabelSize(0.04);
    hMon->GetYaxis()->SetLabelSize(0.04);
    hMon->GetYaxis()->SetTitleOffset(1.5);
    hMon->GetXaxis()->SetRangeUser(0, 10);
    hMon->GetYaxis()->SetRangeUser(1e-9, 1000);

    hMon->Draw("E");
    hJunc->Draw("E SAME");
    if (gAlice) gAlice->Draw("P SAME");

    // Legend
    TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->SetBorderSize(1);
    leg->SetFillColor(kWhite);
    leg->SetTextSize(0.028);
    leg->AddEntry(hMon, "Monash", "l");
    leg->AddEntry(hJunc, "Junctions", "l");
    if (gAlice) leg->AddEntry(gAlice, "ALICE", "p");
    leg->Draw();

    // Label (removed particle name and multiplicity bin labels)

    // Save
    TString outName = Form("pT_spectrum_%s_%s.png", particle, multBin);
    c->SaveAs(outName);
    std::cout << "Created: " << outName << "\n";

    delete c;
    fMon->Close();
    fJunc->Close();
    fAlice->Close();
}

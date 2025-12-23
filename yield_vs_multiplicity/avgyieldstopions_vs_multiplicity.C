// File: avgYields_vs_Multiplicity_pionNorm_FINAL.C
// Pion-normalised yields vs multiplicity with ALICE overlay
// Based on pionNorm_v3 with ALICE data integration
//
// KEY FEATURES:
// 1. Yields use WEIGHTED mean: yield = total_particles / total_events
// 2. Errors from subsampling: std dev of subsampled yields around weighted mean
// 3. Ratio errors: proper error propagation for R = A / B
// 4. ALICE pion ratio data overlaid from HEPData root file
//
// Usage in ROOT:
//   root -l
//   .L avgYields_vs_Multiplicity_pionNorm_FINAL.C
//   avgYields_vs_Multiplicity_pionNorm_FINAL("spectra_monash_100M.root",
//                                             "spectra_junctions_100M.root",
//                                             "HEPData-ins1471838-v1-root.root")

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TMath.h>
#include <TF1.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <iomanip>
#include <cmath>
#include <map>

struct ParticleDataPN {
    const char* histName;
    const char* displayName;
    int color;
    int marker;
    bool isBaryon;   // false = meson, true = baryon
    int aliceTable;  // Table number for ALICE data (-1 if not available)
};

struct SubsampledYield {
    double mean;
    double stdDev;
    std::vector<double> subsampleValues;
};

// ============================================================================
// Map <dNch/dη> to percentile
// ============================================================================
double MapMultiplicityToPercentile(TH1D* hMult, double targetMult) {
    int nbins = hMult->GetNbinsX();
    double totalEvents = hMult->Integral(1, nbins);

    int targetBin = -1;
    double minDiff = 1e9;
    for (int b = 1; b <= nbins; ++b) {
        double binCenter = hMult->GetBinCenter(b);
        double diff = std::abs(binCenter - targetMult);
        if (diff < minDiff) {
            minDiff = diff;
            targetBin = b;
        }
    }

    if (targetBin < 0) return -1;

    double eventsBelow = hMult->Integral(1, targetBin - 1);
    double percentile = 100.0 * eventsBelow / totalEvents;

    return percentile;
}

// ============================================================================
// Extract ALICE ratio data (pion-normalized)
// ============================================================================
TGraphErrors* ExtractAliceRatio(TFile* aliceFile, int tableNum, const char* graphName) {
    TDirectory* tableDir = (TDirectory*)aliceFile->Get(Form("Table %d", tableNum));
    if (!tableDir) return nullptr;
    
    auto* h = dynamic_cast<TH1*>(tableDir->Get("Hist1D_y1"));
    auto* h_e1 = dynamic_cast<TH1*>(tableDir->Get("Hist1D_y1_e1"));
    auto* h_e2 = dynamic_cast<TH1*>(tableDir->Get("Hist1D_y1_e2"));
    
    if (!h || !h_e1 || !h_e2) return nullptr;
    
    std::vector<double> vx, vy, vex, vey;
    
    int nbins = h->GetNbinsX();
    for (int b = 1; b <= nbins; ++b) {
        double x = h->GetXaxis()->GetBinCenter(b);
        double ex = h->GetXaxis()->GetBinWidth(b) / 2.0;
        double y = h->GetBinContent(b);
        double ey1 = h_e1->GetBinContent(b);
        double ey2 = h_e2->GetBinContent(b);
        double ey = std::sqrt(ey1*ey1 + ey2*ey2);
        
        if (y > 0) {
            vx.push_back(x);
            vy.push_back(y);
            vex.push_back(ex);
            vey.push_back(ey);
        }
    }
    
    if (vx.empty()) return nullptr;
    
    auto* g = new TGraphErrors((int)vx.size(), vx.data(), vy.data(), vex.data(), vey.data());
    g->SetName(graphName);
    return g;
}

// ============================================================================
// CORRECTED Subsampling with WEIGHTED mean and proper error calculation
// ============================================================================
SubsampledYield CalculateSubsampledYieldPN(
    TH2D* h2D,
    TH1D* hMult,
    int binLow,
    int binHigh,
    int nSubsamples = 20
)
{
    SubsampledYield result;
    result.subsampleValues.clear();

    TString xTitle(h2D->GetXaxis()->GetTitle());
    xTitle.ToLower();

    TH1D* fullProj = nullptr;
    if (xTitle.Contains("pt") || xTitle.Contains("p_{t}"))
        fullProj = h2D->ProjectionX("fullProj_tempPN", binLow, binHigh);
    else
        fullProj = h2D->ProjectionY("fullProj_tempPN", binLow, binHigh);

    double fullIntegral = fullProj->Integral();
    double fullEvents = hMult->Integral(binLow, binHigh);

    if (fullEvents <= 0) {
        result.mean = 0;
        result.stdDev = 0;
        delete fullProj;
        return result;
    }

    double weightedMean = fullIntegral / fullEvents;
    result.mean = weightedMean;

    int multBinsToProcess = binHigh - binLow + 1;
    int binsPerMultSubsample = std::max(1, multBinsToProcess / nSubsamples);
    int actualSubsamples =
        (multBinsToProcess + binsPerMultSubsample - 1) / binsPerMultSubsample;

    for (int s = 0; s < actualSubsamples; ++s) {
        int currentBinLow  = binLow + s * binsPerMultSubsample;
        int currentBinHigh = std::min(binLow + (s + 1) * binsPerMultSubsample - 1,
                                      binHigh);
        if (currentBinLow > currentBinHigh) break;

        double nEvSub = hMult->Integral(currentBinLow, currentBinHigh);
        if (nEvSub <= 0) continue;

        TH1D* proj = nullptr;
        if (xTitle.Contains("pt") || xTitle.Contains("p_{t}"))
            proj = h2D->ProjectionX("proj_subPN", currentBinLow, currentBinHigh);
        else
            proj = h2D->ProjectionY("proj_subPN", currentBinLow, currentBinHigh);

        double yieldSub = proj->Integral();
        double yieldPerEvent = yieldSub / nEvSub;
        result.subsampleValues.push_back(yieldPerEvent);
        delete proj;
    }

    if (!result.subsampleValues.empty()) {
        double sumSq = 0;
        for (double v : result.subsampleValues)
            sumSq += (v - weightedMean) * (v - weightedMean);

        int N = result.subsampleValues.size();
        result.stdDev = (N > 1) ? std::sqrt(sumSq / (N - 1)) : 0.0;
    } else {
        result.stdDev = 0;
    }

    delete fullProj;
    return result;
}

// ============================================================================
// Helper: Calculate error for ratio R = A / B
// ============================================================================
double CalculateRatioError(double meanA, double errA, double meanB, double errB) {
    if (meanA <= 0 || meanB <= 0) return 0;
    
    double relErrA = errA / meanA;
    double relErrB = errB / meanB;
    double relErrRatio = std::sqrt(relErrA * relErrA + relErrB * relErrB);
    
    double ratio = meanA / meanB;
    return ratio * relErrRatio;
}

// ============================================================================
// MAIN FUNCTION
// ============================================================================
void avgYields_vs_Multiplicity_pionNorm_FINAL(
    const char* monashFileName    = "spectra_monash_100M.root",
    const char* junctionsFileName = "spectra_junctions_100M.root",
    const char* aliceFileName     = "HEPData-ins1471838-v1-root.root")
{
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║ avgYields_vs_Multiplicity_pionNorm_FINAL                       ║\n";
    std::cout << "║ Pion-normalized yields with ALICE overlay                      ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n\n";

    // Open files
    TFile* fMon = TFile::Open(monashFileName, "READ");
    if (!fMon || fMon->IsZombie()) {
        std::cerr << "[ERROR] Cannot open Monash file: " << monashFileName << "\n";
        return;
    }
    TFile* fJunc = TFile::Open(junctionsFileName, "READ");
    if (!fJunc || fJunc->IsZombie()) {
        std::cerr << "[ERROR] Cannot open Junctions file: " << junctionsFileName << "\n";
        fMon->Close();
        return;
    }

    TFile* fAlice = TFile::Open(aliceFileName, "READ");
    if (!fAlice || fAlice->IsZombie()) {
        std::cerr << "[ERROR] Cannot open ALICE file\n";
    }

    auto* hMultMon  = dynamic_cast<TH1D*>(fMon->Get("fHistMultiplicity"));
    auto* hMultJunc = dynamic_cast<TH1D*>(fJunc->Get("fHistMultiplicity"));
    if (!hMultMon || !hMultJunc) {
        std::cerr << "[ERROR] Missing fHistMultiplicity in one of the files.\n";
        fMon->Close();
        fJunc->Close();
        return;
    }

    // Particle list (includes ALICE table numbers for pion-normalized yields - ABSOLUTE RATIOS)
    std::vector<ParticleDataPN> particles = {
        {"fHistPtPions",   "#pi^{#pm} (Pions)",     kBlue,      20, false, -1},   // 0 - no ratio
        {"fHistPtKaons",   "K^{#pm} (Kaons)",       kOrange+2,  22, false, -1},   // K/π - no ALICE data available
        {"fHistPtK0s",     "K^{0}_{s} (K0s)",       kGreen+2,   23, false, 36},   // K0s/π absolute ratio (Table 36)
        {"fHistPtProtons", "p/#bar{p} (Protons)",   kGray+2,    34, true,  -1},   // p/π - no ALICE data available
        {"fHistPtLambdas", "#Lambda (Lambdas)",     kRed+1,     43, true,  37},   // Λ/π absolute ratio (Table 37)
        {"fHistPtXis",     "#Xi^{#pm} (Xis)",       kViolet,    45, true,  38},   // Ξ/π absolute ratio (Table 38)
        {"fHistPtOmegas",  "#Omega^{#pm} (Omegas)", kOrange+4,  21, true,  39}    // Ω/π absolute ratio (Table 39)
    };

    const int nBins       = 10;
    const int nSubsamples = 20;
    const int nSpec       = particles.size();

    // Multiplicity class boundaries
    const int nbMult = hMultMon->GetNbinsX();
    const double total = hMultMon->Integral(1, nbMult);

    std::vector<int> binEdges;
    binEdges.push_back(nbMult + 1);
    double acc = 0.0;
    int targetIdx = 1;
    const std::array<double,10> targets{
        0.10, 0.20, 0.30, 0.40, 0.50,
        0.60, 0.70, 0.80, 0.90, 1.00
    };

    for (int b = nbMult; b >= 1 && targetIdx <= 10; --b) {
        acc += hMultMon->GetBinContent(b);
        double frac = acc / total;
        while (targetIdx <= 10 && frac >= targets[targetIdx-1]) {
            binEdges.push_back(b);
            ++targetIdx;
        }
    }
    while (binEdges.size() < 11) binEdges.push_back(1);
    std::reverse(binEdges.begin(), binEdges.end());

    std::vector<double> xValues;
    for (int i = 0; i < nBins; ++i) {
        xValues.push_back(5.0 + i * 10.0);
    }

    // Pre-load ALICE pion ratio data with PROPER percentile mapping
    std::map<std::string, TGraphErrors*> aliceRatioGraphs;
    if (fAlice) {
        for (const auto& p : particles) {
            if (p.aliceTable < 0) continue;

            TGraphErrors* gAlice_raw = ExtractAliceRatio(fAlice, p.aliceTable, Form("alice_%s", p.histName));

            if (!gAlice_raw) continue;

            // Remap ALICE x-axis to multiplicity percentages using proper mapping
            std::vector<double> x_new, y_new, ex_new, ey_new;
            int n = gAlice_raw->GetN();
            for (int i = 0; i < n; ++i) {
                double x_alice, y_alice;
                gAlice_raw->GetPoint(i, x_alice, y_alice);
                double ey_alice = gAlice_raw->GetErrorY(i);

                // Map ALICE <dNch/dη> to percentile in PYTHIA distribution
                double percentile = MapMultiplicityToPercentile(hMultMon, x_alice);

                if (percentile >= 0) {
                    x_new.push_back(percentile);
                    y_new.push_back(y_alice);
                    ex_new.push_back(0.0);  // No horizontal error bars
                    ey_new.push_back(0.0);  // No vertical error bars
                }
            }

            if (x_new.empty()) {
                delete gAlice_raw;
                continue;
            }

            auto* gAlice_mapped = new TGraphErrors((int)x_new.size(),
                                                    x_new.data(), y_new.data(),
                                                    ex_new.data(), ey_new.data());
            gAlice_mapped->SetMarkerStyle(25);
            gAlice_mapped->SetMarkerSize(0.8);
            gAlice_mapped->SetMarkerColor(p.color);
            gAlice_mapped->SetLineColor(p.color);
            gAlice_mapped->SetLineWidth(1);
            gAlice_mapped->SetLineStyle(2);  // Dashed line

            aliceRatioGraphs[p.histName] = gAlice_mapped;
            delete gAlice_raw;
        }
    }

    // Calculate yields for all species
    std::vector<std::vector<double>> yMon(nSpec), yJ(nSpec);
    std::vector<std::vector<double>> eMon(nSpec), eJ(nSpec);
    std::vector<double> pionMon(nBins), pionJ(nBins);
    std::vector<double> ePionMon(nBins), ePionJ(nBins);

    for (int ip = 0; ip < nSpec; ++ip) {
        auto* h2DMon = dynamic_cast<TH2D*>(fMon->Get(particles[ip].histName));
        auto* h2DJunc = dynamic_cast<TH2D*>(fJunc->Get(particles[ip].histName));

        if (!h2DMon || !h2DJunc) {
            std::cerr << "[WARNING] Cannot find " << particles[ip].histName << "\n";
            continue;
        }

        for (int i = 0; i < nBins; ++i) {
            int binLow = binEdges[i];
            int binHigh = binEdges[i + 1];

            auto resMon = CalculateSubsampledYieldPN(h2DMon, hMultMon, binLow, binHigh, nSubsamples);
            auto resJunc = CalculateSubsampledYieldPN(h2DJunc, hMultJunc, binLow, binHigh, nSubsamples);

            yMon[ip].push_back(resMon.mean);
            yJ[ip].push_back(resJunc.mean);
            eMon[ip].push_back(resMon.stdDev);
            eJ[ip].push_back(resJunc.stdDev);
        }
    }

    // Store pion yields for normalization
    pionMon = yMon[0];
    pionJ = yJ[0];
    ePionMon = eMon[0];
    ePionJ = eJ[0];

    // ===================== 1) All Particles =====================
    TCanvas* cAllNorm = new TCanvas("cAllNorm_pionNorm_FINAL",
        "Pion-Normalized Yields vs Multiplicity (All Particles)",
        1200, 800);
    cAllNorm->SetLogy();

    TMultiGraph* mgAllNorm = new TMultiGraph();
    TLegend* legAllNorm = new TLegend(0.65, 0.15, 0.88, 0.35);
    legAllNorm->SetBorderSize(1);
    legAllNorm->SetTextSize(0.025);
    legAllNorm->SetFillColor(kWhite);

    TGraphErrors* gRefMonAll = nullptr;
    TGraphErrors* gRefJAll = nullptr;

    for (int ip = 1; ip < nSpec; ++ip) {
        const auto& p = particles[ip];

        double rMon[10], rMonErr[10];
        double rJ[10],   rJErr[10];

        for (int i = 0; i < nBins; ++i) {
            if (pionMon[i] > 0) {
                rMon[i]    = yMon[ip][i] / pionMon[i];
                rMonErr[i] = CalculateRatioError(yMon[ip][i], eMon[ip][i], 
                                                  pionMon[i], ePionMon[i]);
            } else {
                rMon[i] = rMonErr[i] = 0.;
            }

            if (pionJ[i] > 0) {
                rJ[i]    = yJ[ip][i] / pionJ[i];
                rJErr[i] = CalculateRatioError(yJ[ip][i], eJ[ip][i], 
                                               pionJ[i], ePionJ[i]);
            } else {
                rJ[i] = rJErr[i] = 0.;
            }
        }

        TGraphErrors* gMonG = new TGraphErrors(nBins, xValues.data(), rMon, nullptr, rMonErr);
        gMonG->SetMarkerStyle(p.marker);
        gMonG->SetMarkerSize(1.2);
        gMonG->SetMarkerColor(p.color);
        gMonG->SetLineColor(p.color);
        gMonG->SetLineWidth(2);
        gMonG->SetLineStyle(1);

        TGraphErrors* gJG = new TGraphErrors(nBins, xValues.data(), rJ, nullptr, rJErr);
        gJG->SetMarkerStyle(p.marker);
        gJG->SetMarkerSize(1.2);
        gJG->SetMarkerColor(p.color);
        gJG->SetLineColor(p.color);
        gJG->SetLineWidth(2);
        gJG->SetLineStyle(2);

        // Add ALICE if available
        if (aliceRatioGraphs.count(p.histName)) {
            mgAllNorm->Add(aliceRatioGraphs[p.histName], "LPX");
        }

        mgAllNorm->Add(gMonG, "LPX");
        mgAllNorm->Add(gJG,   "LPX");
        legAllNorm->AddEntry(gMonG, p.displayName, "lp");
        if (!gRefMonAll) { gRefMonAll = gMonG; gRefJAll = gJG; }
    }

    cAllNorm->cd();
    mgAllNorm->Draw("A");
    mgAllNorm->SetTitle("Average Number of Particles per Pion (#pi^{#pm}) vs Multiplicity;"
                        "Multiplicity interval (%);"
                        "Normalized Yield N_{species}/N_{#pi^{#pm}}");
    mgAllNorm->GetXaxis()->SetLimits(0, 100);
    mgAllNorm->GetXaxis()->SetNdivisions(11);
    mgAllNorm->GetYaxis()->SetRangeUser(1e-5, 1.0);

    {
        TAxis* ax = mgAllNorm->GetXaxis();
        for (int i = 0; i <= 10; ++i)
            ax->ChangeLabel(i+1, -1, -1, -1, -1, -1,
                            Form("%d", 100 - i*10));
    }
    legAllNorm->Draw();

    if (gRefMonAll && gRefJAll) {
        TLegend* legModels = new TLegend(0.18, 0.80, 0.45, 0.90);
        legModels->SetBorderSize(0);
        legModels->SetFillStyle(0);
        legModels->SetTextSize(0.025);
        
        // Create temporary graphs with blue styling for legend
        TGraphErrors* tempMonash = new TGraphErrors();
        tempMonash->SetLineStyle(1);
        tempMonash->SetLineColor(kBlue);
        tempMonash->SetMarkerStyle(20);
        tempMonash->SetMarkerColor(kBlue);
        
        TGraphErrors* tempJunctions = new TGraphErrors();
        tempJunctions->SetLineStyle(2);
        tempJunctions->SetLineColor(kBlue);
        tempJunctions->SetMarkerStyle(20);
        tempJunctions->SetMarkerColor(kBlue);
        
        TGraphErrors* tempAlice = new TGraphErrors();
        tempAlice->SetLineStyle(2);
        tempAlice->SetLineColor(kBlue);
        tempAlice->SetMarkerStyle(25);
        tempAlice->SetMarkerColor(kBlue);
        
        legModels->AddEntry(tempMonash, "Monash", "lp");
        legModels->AddEntry(tempJunctions, "Junctions", "lp");
        legModels->AddEntry(tempAlice, "ALICE", "lp");
        
        legModels->Draw();
    }

    cAllNorm->SaveAs("avgNormYields_perPion_FINAL_all.png");

    // ===================== 2) Mesons =====================
    TCanvas* cMesons = new TCanvas("cMesons_pionNorm_FINAL",
        "Pion-Normalized Yields vs Multiplicity (Mesons)",
        1200, 800);
    cMesons->SetLogy();

    TMultiGraph* mgMesons = new TMultiGraph();
    TLegend* legMesons = new TLegend(0.65, 0.15, 0.88, 0.35);
    legMesons->SetBorderSize(1);
    legMesons->SetTextSize(0.025);
    legMesons->SetFillColor(kWhite);

    TGraphErrors* gRefMonMes = nullptr;
    TGraphErrors* gRefJMes   = nullptr;

    for (int ip = 1; ip < nSpec; ++ip) {
        const auto& p = particles[ip];
        if (p.isBaryon) continue;

        double rMon[10], rMonErr[10];
        double rJ[10],   rJErr[10];

        for (int i = 0; i < nBins; ++i) {
            if (pionMon[i] > 0) {
                rMon[i]    = yMon[ip][i] / pionMon[i];
                rMonErr[i] = CalculateRatioError(yMon[ip][i], eMon[ip][i], 
                                                  pionMon[i], ePionMon[i]);
            } else {
                rMon[i] = rMonErr[i] = 0.;
            }

            if (pionJ[i] > 0) {
                rJ[i]    = yJ[ip][i] / pionJ[i];
                rJErr[i] = CalculateRatioError(yJ[ip][i], eJ[ip][i], 
                                               pionJ[i], ePionJ[i]);
            } else {
                rJ[i] = rJErr[i] = 0.;
            }
        }

        TGraphErrors* gMonG = new TGraphErrors(nBins, xValues.data(), rMon, nullptr, rMonErr);
        gMonG->SetMarkerStyle(p.marker);
        gMonG->SetMarkerSize(1.2);
        gMonG->SetMarkerColor(p.color);
        gMonG->SetLineColor(p.color);
        gMonG->SetLineWidth(2);
        gMonG->SetLineStyle(1);

        TGraphErrors* gJG = new TGraphErrors(nBins, xValues.data(), rJ, nullptr, rJErr);
        gJG->SetMarkerStyle(p.marker);
        gJG->SetMarkerSize(1.2);
        gJG->SetMarkerColor(p.color);
        gJG->SetLineColor(p.color);
        gJG->SetLineWidth(2);
        gJG->SetLineStyle(2);

        // Add ALICE if available
        if (aliceRatioGraphs.count(p.histName)) {
            mgMesons->Add(aliceRatioGraphs[p.histName], "LPX");
        }

        mgMesons->Add(gMonG, "LPX");
        mgMesons->Add(gJG,   "LPX");
        legMesons->AddEntry(gMonG, p.displayName, "lp");
        if (!gRefMonMes) { gRefMonMes = gMonG; gRefJMes = gJG; }
    }

    cMesons->cd();
    mgMesons->Draw("A");
    mgMesons->SetTitle("Average Number of Mesons per Pion (#pi^{#pm}) vs Multiplicity;"
                       "Multiplicity interval (%);"
                       "Normalized Yield N_{species}/N_{#pi^{#pm}}");
    mgMesons->GetXaxis()->SetLimits(0, 100);
    mgMesons->GetXaxis()->SetNdivisions(11);
    mgMesons->GetYaxis()->SetRangeUser(1e-3, 1.0);

    {
        TAxis* ax = mgMesons->GetXaxis();
        for (int i = 0; i <= 10; ++i)
            ax->ChangeLabel(i+1, -1, -1, -1, -1, -1,
                            Form("%d", 100 - i*10));
    }
    legMesons->Draw();

    if (gRefMonMes && gRefJMes) {
        TLegend* legModels = new TLegend(0.18, 0.80, 0.45, 0.90);
        legModels->SetBorderSize(0);
        legModels->SetFillStyle(0);
        legModels->SetTextSize(0.025);
        
        // Create temporary graphs with blue styling for legend
        TGraphErrors* tempMonash = new TGraphErrors();
        tempMonash->SetLineStyle(1);
        tempMonash->SetLineColor(kBlue);
        tempMonash->SetMarkerStyle(20);
        tempMonash->SetMarkerColor(kBlue);
        
        TGraphErrors* tempJunctions = new TGraphErrors();
        tempJunctions->SetLineStyle(2);
        tempJunctions->SetLineColor(kBlue);
        tempJunctions->SetMarkerStyle(20);
        tempJunctions->SetMarkerColor(kBlue);
        
        TGraphErrors* tempAlice = new TGraphErrors();
        tempAlice->SetLineStyle(2);
        tempAlice->SetLineColor(kBlue);
        tempAlice->SetMarkerStyle(25);
        tempAlice->SetMarkerColor(kBlue);
        
        legModels->AddEntry(tempMonash, "Monash", "lp");
        legModels->AddEntry(tempJunctions, "Junctions", "lp");
        legModels->AddEntry(tempAlice, "ALICE", "lp");
        
        legModels->Draw();
    }

    cMesons->SaveAs("avgNormYields_perPion_FINAL_mesons.png");

    // ===================== 3) Baryons =====================
    TCanvas* cBaryons = new TCanvas("cBaryons_pionNorm_FINAL",
        "Pion-Normalized Yields vs Multiplicity (Baryons)",
        1200, 800);
    cBaryons->SetLogy();

    TMultiGraph* mgBaryons = new TMultiGraph();
    TLegend* legBaryons = new TLegend(0.65, 0.15, 0.88, 0.35);
    legBaryons->SetBorderSize(1);
    legBaryons->SetTextSize(0.025);
    legBaryons->SetFillColor(kWhite);

    TGraphErrors* gRefMonBar = nullptr;
    TGraphErrors* gRefJBar   = nullptr;

    for (int ip = 1; ip < nSpec; ++ip) {
        const auto& p = particles[ip];
        if (!p.isBaryon) continue;

        double rMon[10], rMonErr[10];
        double rJ[10],   rJErr[10];

        for (int i = 0; i < nBins; ++i) {
            if (pionMon[i] > 0) {
                rMon[i]    = yMon[ip][i] / pionMon[i];
                rMonErr[i] = CalculateRatioError(yMon[ip][i], eMon[ip][i], 
                                                  pionMon[i], ePionMon[i]);
            } else {
                rMon[i] = rMonErr[i] = 0.;
            }

            if (pionJ[i] > 0) {
                rJ[i]    = yJ[ip][i] / pionJ[i];
                rJErr[i] = CalculateRatioError(yJ[ip][i], eJ[ip][i], 
                                               pionJ[i], ePionJ[i]);
            } else {
                rJ[i] = rJErr[i] = 0.;
            }
        }

        TGraphErrors* gMonG = new TGraphErrors(nBins, xValues.data(), rMon, nullptr, rMonErr);
        gMonG->SetMarkerStyle(p.marker);
        gMonG->SetMarkerSize(1.2);
        gMonG->SetMarkerColor(p.color);
        gMonG->SetLineColor(p.color);
        gMonG->SetLineWidth(2);
        gMonG->SetLineStyle(1);

        TGraphErrors* gJG = new TGraphErrors(nBins, xValues.data(), rJ, nullptr, rJErr);
        gJG->SetMarkerStyle(p.marker);
        gJG->SetMarkerSize(1.2);
        gJG->SetMarkerColor(p.color);
        gJG->SetLineColor(p.color);
        gJG->SetLineWidth(2);
        gJG->SetLineStyle(2);

        // Add ALICE if available
        if (aliceRatioGraphs.count(p.histName)) {
            mgBaryons->Add(aliceRatioGraphs[p.histName], "LPX");
        }

        mgBaryons->Add(gMonG, "LPX");
        mgBaryons->Add(gJG,   "LPX");
        legBaryons->AddEntry(gMonG, p.displayName, "lp");
        if (!gRefMonBar) { gRefMonBar = gMonG; gRefJBar = gJG; }
    }

    cBaryons->cd();
    mgBaryons->Draw("A");
    mgBaryons->SetTitle("Average Number of Baryons per Pion (#pi^{#pm}) vs Multiplicity;"
                        "Multiplicity interval (%);"
                        "Normalized Yield N_{species}/N_{#pi^{#pm}}");
    mgBaryons->GetXaxis()->SetLimits(0, 100);
    mgBaryons->GetXaxis()->SetNdivisions(11);
    mgBaryons->GetYaxis()->SetRangeUser(1e-7, 1.0);

    {
        TAxis* ax = mgBaryons->GetXaxis();
        for (int i = 0; i <= 10; ++i)
            ax->ChangeLabel(i+1, -1, -1, -1, -1, -1,
                            Form("%d", 100 - i*10));
    }
    legBaryons->Draw();

    if (gRefMonBar && gRefJBar) {
        TLegend* legModels = new TLegend(0.18, 0.80, 0.45, 0.90);
        legModels->SetBorderSize(0);
        legModels->SetFillStyle(0);
        legModels->SetTextSize(0.025);
        
        // Create temporary graphs with blue styling for legend
        TGraphErrors* tempMonash = new TGraphErrors();
        tempMonash->SetLineStyle(1);
        tempMonash->SetLineColor(kBlue);
        tempMonash->SetMarkerStyle(20);
        tempMonash->SetMarkerColor(kBlue);
        
        TGraphErrors* tempJunctions = new TGraphErrors();
        tempJunctions->SetLineStyle(2);
        tempJunctions->SetLineColor(kBlue);
        tempJunctions->SetMarkerStyle(20);
        tempJunctions->SetMarkerColor(kBlue);
        
        TGraphErrors* tempAlice = new TGraphErrors();
        tempAlice->SetLineStyle(2);
        tempAlice->SetLineColor(kBlue);
        tempAlice->SetMarkerStyle(25);
        tempAlice->SetMarkerColor(kBlue);
        
        legModels->AddEntry(tempMonash, "Monash", "lp");
        legModels->AddEntry(tempJunctions, "Junctions", "lp");
        legModels->AddEntry(tempAlice, "ALICE", "lp");
        
        legModels->Draw();
    }

    cBaryons->SaveAs("avgNormYields_perPion_FINAL_baryons.png");

    fMon->Close();
    fJunc->Close();
    if (fAlice) fAlice->Close();

    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║ [SUCCESS] Analysis complete!                                   ║\n";
    std::cout << "║ Plots:  avgNormYields_perPion_FINAL_*.png                      ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n\n";

    delete cAllNorm;
    delete cMesons;
    delete cBaryons;
}

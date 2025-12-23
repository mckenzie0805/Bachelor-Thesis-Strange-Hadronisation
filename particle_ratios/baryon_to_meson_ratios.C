// File: pT_spectra_StrangeParticles_WithRatio.C
//
// Creates pT spectra plots for ALL particles with ratio panels
// Light: Pions, Protons
// Strange Mesons: K0s, K±
// Strange Baryons: Lambda, Xi, Omega
//
// Features:
// - 5 multiplicity classes: 0-20%, 20-40%, 40-60%, 60-80%, 80-100%
// - pT range: 0-10 GeV/c
// - Normalized to number of events: (1/Nev) × dN/dpT
// - Ratio panel below main plot (relative to 0-20% bin)
//
// Usage:
//   root -l 'pT_spectra_StrangeParticles_WithRatio.C("spectra_monash_100M.root")'

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TLine.h>
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>

static const double PT_MIN = 0.0;
static const double PT_MAX = 10.0;

struct BinRange { int lo; int hi; };

// ============================================================================
// Helper functions
// ============================================================================

static bool xIsPt(const TH2* h2) {
    auto looksLikePt = [](const char* t) {
        TString s(t ? t : "");
        s.ToLower();
        return s.Contains("p_{t}") || s.Contains("p_t") || s.Contains("pt") || s.Contains("gev");
    };
    bool xpt = looksLikePt(h2->GetXaxis()->GetTitle());
    bool ypt = looksLikePt(h2->GetYaxis()->GetTitle());
    if (xpt && !ypt) return true;
    if (!xpt && ypt) return false;
    return true;
}

// Project pT with subsampling error bars and normalize per event
static TH1D* projectPtSliceNormalized(TH2* h2, TH1D* hMult, int multLo, int multHi,
                                       const TString& name, int nSubsamples = 20) {
    if (!h2 || !hMult) return nullptr;

    // Get number of events in this multiplicity range
    double nEvents = hMult->Integral(multLo, multHi);
    if (nEvents <= 0) return nullptr;

    // Create the main projection
    TH1D* hProj = nullptr;
    if (xIsPt(h2)) {
        hProj = (TH1D*)h2->ProjectionX(name, multLo, multHi, "e");
    } else {
        hProj = (TH1D*)h2->ProjectionY(name, multLo, multHi, "e");
    }

    if (!hProj) return nullptr;

    // Subsampling for error estimation
    int multBinsToProcess = multHi - multLo + 1;
    int binsPerSubsample = std::max(1, multBinsToProcess / nSubsamples);
    int actualSubsamples = std::min(nSubsamples, multBinsToProcess);

    std::vector<TH1D*> subsamples;
    for (int s = 0; s < actualSubsamples; ++s) {
        int currentLo = multLo + s * binsPerSubsample;
        int currentHi = std::min(multLo + (s + 1) * binsPerSubsample - 1, multHi);
        if (currentLo > currentHi) break;

        TH1D* hSub = nullptr;
        if (xIsPt(h2)) {
            hSub = (TH1D*)h2->ProjectionX(Form("%s_sub%d", name.Data(), s), currentLo, currentHi, "e");
        } else {
            hSub = (TH1D*)h2->ProjectionY(Form("%s_sub%d", name.Data(), s), currentLo, currentHi, "e");
        }

        if (hSub) {
            double nEvSub = hMult->Integral(currentLo, currentHi);
            if (nEvSub > 0) {
                // Normalize each subsample
                for (int i = 1; i <= hSub->GetNbinsX(); ++i) {
                    double binWidth = hSub->GetBinWidth(i);
                    if (binWidth > 0) {
                        hSub->SetBinContent(i, hSub->GetBinContent(i) / (nEvSub * binWidth));
                    }
                }
                subsamples.push_back(hSub);
            } else {
                delete hSub;
            }
        }
    }

    // Calculate mean and standard error for each pT bin
    for (int i = 1; i <= hProj->GetNbinsX(); ++i) {
        double binWidth = hProj->GetBinWidth(i);
        if (binWidth <= 0) continue;

        std::vector<double> values;
        for (auto* hSub : subsamples) {
            double val = hSub->GetBinContent(i);
            if (val > 0) values.push_back(val);
        }

        if (!values.empty()) {
            // Calculate mean
            double mean = 0.0;
            for (double v : values) mean += v;
            mean /= values.size();

            // Calculate standard deviation
            double stddev = 0.0;
            if (values.size() > 1) {
                for (double v : values) {
                    stddev += (v - mean) * (v - mean);
                }
                stddev = std::sqrt(stddev / (values.size() - 1));
            }

            // Standard error of the mean
            double stderr = stddev / std::sqrt(values.size());

            hProj->SetBinContent(i, mean);
            hProj->SetBinError(i, stderr);
        } else {
            // No valid subsamples, use simple normalization
            double content = hProj->GetBinContent(i);
            hProj->SetBinContent(i, content / (nEvents * binWidth));
            hProj->SetBinError(i, 0);
        }
    }

    // Cleanup subsamples
    for (auto* h : subsamples) delete h;

    hProj->GetXaxis()->SetRangeUser(PT_MIN, PT_MAX);
    return hProj;
}

// Build 5 multiplicity class ranges
static std::array<BinRange,5> makeClassRanges(TH1D* hMult) {
    std::array<BinRange,5> rngs{};
    if (!hMult) return rngs;

    const int nb = hMult->GetNbinsX();
    const double total = hMult->Integral(1, nb);

    if (total <= 0) {
        int step = nb/5;
        for (int i=0; i<5; ++i) {
            int lo = i*step + 1;
            int hi = (i==4) ? nb : (i+1)*step;
            rngs[i] = {lo, hi};
        }
        return rngs;
    }

    std::vector<int> edges;
    edges.push_back(nb);

    double acc = 0.0;
    int targetIdx = 1;
    const std::array<double,5> targets{0.20, 0.40, 0.60, 0.80, 1.00};

    for (int b = nb; b >= 1 && targetIdx <= 5; --b) {
        acc += hMult->GetBinContent(b);
        const double frac = acc/total;
        while (targetIdx <= 5 && frac >= targets[targetIdx-1]) {
            edges.push_back(b);
            ++targetIdx;
        }
    }

    if (edges.size() < 6) edges.push_back(1);

    for (int i = 0; i < 5; ++i) {
        rngs[i].lo = (i+1 < (int)edges.size()) ? edges[i+1]+1 : 1;
        rngs[i].hi = edges[i];
        if (rngs[i].lo > rngs[i].hi) std::swap(rngs[i].lo, rngs[i].hi);
    }

    return rngs;
}

// ============================================================================
// Draw species with ratio panel
// ============================================================================
static void drawSpeciesWithRatio(TH2* h2,
                                 TH1D* hMult,
                                 const std::array<BinRange,5>& ranges,
                                 const char* species,
                                 const char* displayName,
                                 const char* outPng) {
    if (!h2 || !hMult) return;

    gStyle->SetOptStat(0);

    std::array<int,5> cols{ kBlue+1, kGreen+2, kOrange+1, kMagenta+1, kRed+1 };
    std::array<const char*,5> labels{
        "0-20%", "20-40%", "40-60%", "60-80%", "80-100%"
    };

    // Prepare histograms for each class
    std::array<TH1D*,5> hClass{};
    for (int i=0; i<5; ++i) {
        const BinRange br = ranges[i];
        hClass[i] = projectPtSliceNormalized(h2, hMult, br.lo, br.hi,
                                              Form("h_%s_%d", species, i));
        if (!hClass[i]) continue;

        hClass[i]->SetLineColor(cols[i]);
        hClass[i]->SetLineWidth(2);
        hClass[i]->SetMarkerStyle(1);  // No markers (invisible)
        hClass[i]->SetMarkerSize(0);
        hClass[i]->SetMarkerColor(cols[i]);
    }

    // Determine y-maximum for main plot
    double ymax = 0.0;
    for (auto* h : hClass) {
        if (h) ymax = std::max(ymax, h->GetMaximum());
    }
    if (ymax <= 0) ymax = 1.0;
    ymax *= 1.4;

    // Create canvas with two pads
    TCanvas c(Form("c_%s_all", species), Form("%s p_{T} spectra", species), 1000, 900);

    // Top pad: main spectra
    TPad* pad1 = new TPad("pad1", "Spectra", 0, 0.32, 1, 1.0);  // Changed from 0.30 to 0.32
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.15);  // Increased from 0.12 to prevent y-axis cutoff
    pad1->SetRightMargin(0.05);
    pad1->SetTopMargin(0.10);   // Increased from 0.08 to give title more space
    pad1->SetLogy();
    pad1->SetGridx();           // Add vertical gridlines
    pad1->SetGridy();           // Add horizontal gridlines
    pad1->Draw();
    pad1->cd();

    // Set light grid color
    gStyle->SetGridColor(kGray);
    gStyle->SetGridStyle(3);    // Dotted lines
    gStyle->SetGridWidth(1);

    // Draw main spectra
    TH1D* href = nullptr;
    for (int i = 0; i < 5; ++i) {
        if (!hClass[i]) continue;
        if (hClass[i]->GetMaximum() <= 0) continue;

        if (!href) {
            hClass[i]->SetTitle(Form("p_{T} spectra (%s) across multiplicity classes", displayName));
            hClass[i]->GetYaxis()->SetTitle("(1/N_{ev}) dN/dp_{T} [(GeV/c)^{-1}]");
            hClass[i]->GetYaxis()->SetTitleSize(0.05);
            hClass[i]->GetYaxis()->SetLabelSize(0.045);
            hClass[i]->GetYaxis()->SetTitleOffset(1.0);  // Bring closer to axis
            hClass[i]->SetMaximum(ymax);
            hClass[i]->SetMinimum(ymax * 1e-5);
            hClass[i]->GetXaxis()->SetLabelSize(0);  // Hide x-axis labels
            hClass[i]->Draw("HIST");
            href = hClass[i];
        } else {
            hClass[i]->Draw("HIST SAME");
        }
        hClass[i]->Draw("E SAME");  // Draw error bars
    }

    // Legend
    TLegend* leg = new TLegend(0.65, 0.55, 0.92, 0.88);
    leg->SetBorderSize(1);  // Add box around legend
    leg->SetFillStyle(1001); // Solid white background
    leg->SetFillColor(kWhite);
    leg->SetTextSize(0.04);
    for (int i=0; i<5; ++i) {
        if (hClass[i]) leg->AddEntry(hClass[i], labels[i], "l");  // Only line, no point
    }
    leg->Draw();

    c.cd();

    // Bottom pad: ratio to 0-20%
    TPad* pad2 = new TPad("pad2", "Ratio", 0, 0.0, 1, 0.32);  // Changed from 0.30 to 0.32 to match top
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);  // Match top pad
    pad2->SetRightMargin(0.05);
    pad2->SetGridx();           // Add vertical gridlines
    pad2->SetGridy();           // Add horizontal gridlines
    pad2->Draw();
    pad2->cd();

    // Create ratio histograms (relative to 0-20%)
    TH1D* hRef = hClass[0];  // 0-20% bin as reference
    if (hRef) {
        std::array<TH1D*,4> hRatio{};
        for (int i=1; i<5; ++i) {
            if (!hClass[i]) continue;
            hRatio[i-1] = (TH1D*)hClass[i]->Clone(Form("ratio_%s_%d", species, i));
            hRatio[i-1]->Divide(hRef);
            hRatio[i-1]->SetLineColor(cols[i]);
            hRatio[i-1]->SetMarkerColor(cols[i]);
            hRatio[i-1]->SetMarkerStyle(1);  // No markers
            hRatio[i-1]->SetMarkerSize(0);
            hRatio[i-1]->SetLineWidth(2);
        }

        // Set fixed y-range for ratio plot
        double ratioMin = 0.0;
        double ratioMax = 0.4;

        // Draw ratio
        TH1D* hRatioFrame = nullptr;
        for (int i=0; i<4; ++i) {
            if (!hRatio[i]) continue;
            if (!hRatioFrame) {
                hRatio[i]->SetTitle("");
                hRatio[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                hRatio[i]->GetYaxis()->SetTitle("Ratio/(0-20%)");
                hRatio[i]->GetXaxis()->SetTitleSize(0.10);
                hRatio[i]->GetXaxis()->SetLabelSize(0.10);
                hRatio[i]->GetYaxis()->SetTitleSize(0.10);
                hRatio[i]->GetYaxis()->SetLabelSize(0.09);
                hRatio[i]->GetYaxis()->SetTitleOffset(0.45);  // Move slightly away from axis
                hRatio[i]->GetXaxis()->SetTitleOffset(0.9);
                hRatio[i]->SetMinimum(ratioMin);
                hRatio[i]->SetMaximum(ratioMax);
                hRatio[i]->GetYaxis()->SetNdivisions(505);
                hRatio[i]->Draw("HIST");  // Line only, no error bars for ratio
                hRatioFrame = hRatio[i];
            } else {
                hRatio[i]->Draw("HIST SAME");
            }
        }

        // Unity line (only if it's in range)
        if (hRatioFrame && 1.0 >= ratioMin && 1.0 <= ratioMax) {
            TLine* unity = new TLine(PT_MIN, 1.0, PT_MAX, 1.0);
            unity->SetLineStyle(2);
            unity->SetLineColor(kBlack);
            unity->SetLineWidth(2);
            unity->Draw("same");
        }
    }

    c.SaveAs(outPng);
    std::cout << "[OK] Saved: " << outPng << "\n";

    // Cleanup
    for (auto* h : hClass) if (h) delete h;
}

// ============================================================================
// Main function
// ============================================================================
void pT_spectra_StrangeParticles_WithRatio(const char* filename = "spectra_monash_100M.root") {
    std::cout << "\n╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║ pT Spectra: Strange Particles with Ratio Panels              ║\n";
    std::cout << "║ Input: " << filename << "\n";
    std::cout << "║ Normalization: (1/Nev) dN/dpT                                ║\n";
    std::cout << "║ pT range: 0-10 GeV/c                                         ║\n";
    std::cout << "║ Ratio: Relative to 0-20% (highest multiplicity)             ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n\n";

    TFile* f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[ERROR] Cannot open file: " << filename << "\n";
        return;
    }

    auto* hMult = dynamic_cast<TH1D*>(f->Get("fHistMultiplicity"));
    if (!hMult) {
        std::cerr << "[ERROR] Missing fHistMultiplicity.\n";
        f->Close();
        return;
    }

    auto ranges = makeClassRanges(hMult);

    std::cout << "\nMultiplicity binning:\n";
    for (int i = 0; i < 5; ++i) {
        std::cout << "  " << (i*20) << "-" << ((i+1)*20) << "%: bins "
                  << ranges[i].lo << " to " << ranges[i].hi << "\n";
    }
    std::cout << "\n";

    // LIGHT PARTICLES (reference)
    std::cout << "\n=== Light Particles ===\n";

    auto* h2Pions = dynamic_cast<TH2D*>(f->Get("fHistPtPions"));
    if (h2Pions) {
        drawSpeciesWithRatio(h2Pions, hMult, ranges, "Pions", "#pi^{#pm}",
                            "pT_Pions_spectra_ratio.png");
    } else {
        std::cerr << "[WARNING] Missing fHistPtPions\n";
    }

    auto* h2Protons = dynamic_cast<TH2D*>(f->Get("fHistPtProtons"));
    if (h2Protons) {
        drawSpeciesWithRatio(h2Protons, hMult, ranges, "Protons", "p/#bar{p}",
                            "pT_Protons_spectra_ratio.png");
    } else {
        std::cerr << "[WARNING] Missing fHistPtProtons\n";
    }

    // STRANGE MESONS
    std::cout << "\n=== Strange Mesons ===\n";

    auto* h2K0s = dynamic_cast<TH2D*>(f->Get("fHistPtK0s"));
    if (h2K0s) {
        drawSpeciesWithRatio(h2K0s, hMult, ranges, "K0s", "K^{0}_{s}",
                            "pT_K0s_spectra_ratio.png");
    } else {
        std::cerr << "[WARNING] Missing fHistPtK0s\n";
    }

    auto* h2Kaons = dynamic_cast<TH2D*>(f->Get("fHistPtKaons"));
    if (h2Kaons) {
        drawSpeciesWithRatio(h2Kaons, hMult, ranges, "Kaons", "K^{#pm}",
                            "pT_Kaons_spectra_ratio.png");
    } else {
        std::cerr << "[WARNING] Missing fHistPtKaons\n";
    }

    // STRANGE BARYONS
    std::cout << "\n=== Strange Baryons ===\n";

    auto* h2Lambda = dynamic_cast<TH2D*>(f->Get("fHistPtLambdas"));
    if (h2Lambda) {
        drawSpeciesWithRatio(h2Lambda, hMult, ranges, "Lambda", "#Lambda/#bar{#Lambda}",
                            "pT_Lambda_spectra_ratio.png");
    } else {
        std::cerr << "[WARNING] Missing fHistPtLambdas\n";
    }

    auto* h2Xi = dynamic_cast<TH2D*>(f->Get("fHistPtXis"));
    if (h2Xi) {
        drawSpeciesWithRatio(h2Xi, hMult, ranges, "Xi", "#Xi^{-} / #bar{#Xi}^{+}",
                            "pT_Xi_spectra_ratio.png");
    } else {
        std::cerr << "[WARNING] Missing fHistPtXis\n";
    }

    auto* h2Omega = dynamic_cast<TH2D*>(f->Get("fHistPtOmegas"));
    if (h2Omega) {
        drawSpeciesWithRatio(h2Omega, hMult, ranges, "Omega", "#Omega^{-} / #bar{#Omega}^{+}",
                            "pT_Omega_spectra_ratio.png");
    } else {
        std::cerr << "[WARNING] Missing fHistPtOmegas\n";
    }

    f->Close();

    std::cout << "\n╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║ [SUCCESS] All pT spectra with ratio panels created!          ║\n";
    std::cout << "║ Output files:                                                 ║\n";
    std::cout << "║   Light particles: Pions, Protons                            ║\n";
    std::cout << "║   Strange mesons: K0s, K±                                    ║\n";
    std::cout << "║   Strange baryons: Lambda, Xi, Omega                         ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n\n";
}

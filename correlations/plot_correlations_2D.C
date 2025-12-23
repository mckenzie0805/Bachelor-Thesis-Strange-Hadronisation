// ============================================================================
//  plot_correlations_2D.C - Side-by-side 2D ΔφΔη correlation maps
// ============================================================================
// Usage: root -l -q 'plot_correlations_2D.C("corr_monash_100M.root", "corr_junctions_100M.root")'

#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <iostream>

void plot_correlations_2D(const char* monashFile, const char* junctionsFile)
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    // Open files
    TFile* fMon = TFile::Open(monashFile);
    TFile* fJunc = TFile::Open(junctionsFile);

    if (!fMon || fMon->IsZombie()) {
        std::cerr << "Cannot open Monash file: " << monashFile << std::endl;
        return;
    }
    if (!fJunc || fJunc->IsZombie()) {
        std::cerr << "Cannot open Junctions file: " << junctionsFile << std::endl;
        return;
    }

    // Get 2D histograms
    TH2D* hMon = (TH2D*)fMon->Get("hDPhiDEta");
    TH2D* hJunc = (TH2D*)fJunc->Get("hDPhiDEta");

    if (!hMon) {
        std::cerr << "Cannot find hDPhiDEta in Monash file" << std::endl;
        return;
    }
    if (!hJunc) {
        std::cerr << "Cannot find hDPhiDEta in Junctions file" << std::endl;
        return;
    }

    // Create canvas with two pads side by side
    TCanvas* c = new TCanvas("c_2D", "Correlation Maps", 1600, 700);
    c->Divide(2, 1);

    // Left pad - Monash
    c->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.12);

    hMon->SetTitle("#Xi-#bar{#Xi} Correlations (Monash);#Delta#phi (rad);#Delta#eta");
    hMon->GetXaxis()->SetTitleSize(0.05);
    hMon->GetYaxis()->SetTitleSize(0.05);
    hMon->GetXaxis()->SetLabelSize(0.04);
    hMon->GetYaxis()->SetLabelSize(0.04);
    hMon->GetZaxis()->SetLabelSize(0.035);
    hMon->Draw("COLZ");

    // Right pad - Junctions
    c->cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.12);

    hJunc->SetTitle("#Xi-#bar{#Xi} Correlations (Junctions);#Delta#phi (rad);#Delta#eta");
    hJunc->GetXaxis()->SetTitleSize(0.05);
    hJunc->GetYaxis()->SetTitleSize(0.05);
    hJunc->GetXaxis()->SetLabelSize(0.04);
    hJunc->GetYaxis()->SetLabelSize(0.04);
    hJunc->GetZaxis()->SetLabelSize(0.035);
    hJunc->Draw("COLZ");

    c->SaveAs("correlations_2D_comparison.png");

    std::cout << "Created: correlations_2D_comparison.png" << std::endl;

    fMon->Close();
    fJunc->Close();
    delete c;
}
